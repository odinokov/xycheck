#!/usr/bin/env python3
"""Fetch Umap mappability + ENCODE blacklist and build a clean chrX/chrY BED.

Run this once per genome to produce a uniquely-mappable, blacklist-subtracted
chrX/chrY track under data/<genome>/. xycheck.py then consumes that premade BED
without touching the network.
"""
import gzip
import logging
import os
import pathlib
import shutil
import sys
import tarfile
import tempfile

import argparse
import pybedtools
import pybedtools.helpers
import requests

UMAP_URL = "https://bismap.hoffmanlab.org/raw/{asm}.umap.tar.gz"
BL_URL   = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/{asm}-blacklist.v2.bed.gz"
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent


def _write_stream(response, dest: pathlib.Path) -> None:
    tmp = dest.with_name(dest.name + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with open(tmp, "wb") as fh:
        for chunk in response.iter_content(1 << 20):
            if chunk:
                fh.write(chunk)
    tmp.replace(dest)


def fetch(url: str, dest: pathlib.Path) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size:
        logging.info("cached  %s", dest)
        return dest
    logging.info("downloading %s", dest.name)
    with requests.get(url, stream=True, timeout=(10, 300)) as r:
        r.raise_for_status()
        _write_stream(r, dest)
    return dest


def _extract_member(tar, member, outdir: pathlib.Path) -> pathlib.Path:
    dst = outdir / pathlib.Path(member.name).name
    if not dst.exists():
        with tar.extractfile(member) as src, open(dst, "wb") as out:
            shutil.copyfileobj(src, out)
    return dst


def extract_umap_bed(tar_path: pathlib.Path, k: int, outdir: pathlib.Path) -> pathlib.Path:
    pats = [f"k{k}.umap.bed.gz", f"{k}bp_kmers.bed.gz", f"{k}bp_kmers.bed"]
    for pat in pats:
        dst = outdir / pat
        if dst.exists() and dst.stat().st_size:
            logging.info("cached extraction %s", dst)
            return dst
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            if any(m.name.endswith(p) for p in pats):
                return _extract_member(tar, m, outdir)
    raise RuntimeError(f"k={k} track not found in {tar_path.name}")


def gunzip(p: pathlib.Path) -> pathlib.Path:
    if p.suffix != ".gz":
        return p
    out = p.with_suffix("")
    if not out.exists():
        with gzip.open(p, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out


def _sample_bed_rows(bed_path: pathlib.Path, n: int = 10) -> list:
    opener = gzip.open if bed_path.suffix == ".gz" else open
    rows = []
    with opener(bed_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            rows.append(line.split())
            if len(rows) >= n:
                break
    return rows


def detect_score_column(bed_path: pathlib.Path) -> int:
    rows = _sample_bed_rows(bed_path)
    if not rows:
        raise RuntimeError(f"No data in {bed_path}")
    for col in range(3, len(rows[0])):
        try:
            vals = [float(row[col]) for row in rows if len(row) > col]
            if vals and all(0.0 <= v <= 1.0 for v in vals):
                logging.debug("mappability score → column %d", col)
                return col
        except ValueError:
            continue
    raise RuntimeError(
        f"No [0,1]-float column found in {bed_path}. "
        "Pass --score-col to specify it manually."
    )


def build_clean_track(
    asm: str, kmer: int, data_dir: pathlib.Path, score_col: int
) -> pathlib.Path:
    wd = data_dir / asm
    out = wd / f"clean_XY.{asm}.k{kmer}.bed"
    if out.exists() and out.stat().st_size:
        logging.info("cached BED %s", out)
        return out

    umap_tar   = fetch(UMAP_URL.format(asm=asm), wd / f"{asm}.umap.tar.gz")
    bl_gz      = fetch(BL_URL.format(asm=asm),   wd / f"{asm}.blacklist.bed.gz")
    umap_bed_gz = extract_umap_bed(umap_tar, kmer, wd)

    col = score_col if score_col is not None else detect_score_column(umap_bed_gz)

    umap_bed  = gunzip(umap_bed_gz)
    blacklist = gunzip(bl_gz)

    logging.info("filtering umap BED: chrX/Y, score[col %d] == 1.0, subtract blacklist", col)

    def _keep(x):
        try:
            return x.chrom in ("chrX", "chrY") and float(x.fields[col]) == 1.0
        except (IndexError, ValueError):
            return False

    (
        pybedtools.BedTool(umap_bed)
        .filter(_keep)
        .subtract(pybedtools.BedTool(blacklist))
        .sort()
        .merge()
        .saveas(str(out))
    )
    logging.info("clean BED written: %s", out)
    return out


def main():
    p = argparse.ArgumentParser(
        description="Fetch Umap + ENCODE blacklist and build a clean chrX/chrY BED."
    )
    p.add_argument("-g", "--genome",    required=True, choices=["hg19", "hg38"],
                   help="Reference genome")
    p.add_argument("-k", "--kmer",      default="100", choices=["100", "50", "36", "24"],
                   metavar="{100,50,36,24}",
                   help="Umap K-mer track (default: 100)")
    p.add_argument("-d", "--data-dir",  default="data", metavar="DIR",
                   help="Cache dir for downloaded/built files. Relative paths are resolved from prepare_bed.py (default: data)")
    p.add_argument("-T", "--tmp-dir",   default=tempfile.gettempdir(), metavar="DIR",
                   help="Temp dir for pybedtools scratch files (default: $TMPDIR)")
    p.add_argument("--score-col",       type=int, default=None, metavar="INT",
                   help="0-based mappability score column in umap BED (auto-detected if omitted)")
    p.add_argument("-v", "--verbose",   action="store_true",
                   help="Show download/filter progress")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    os.environ["TMPDIR"] = args.tmp_dir
    pybedtools.helpers.set_tempdir(args.tmp_dir)

    data_path = pathlib.Path(args.data_dir)
    if not data_path.is_absolute():
        data_path = SCRIPT_DIR / data_path

    try:
        bed_path = build_clean_track(args.genome, int(args.kmer), data_path, args.score_col)
        print(bed_path)
    except Exception as e:
        logging.error("%s", e)
        sys.exit(1)
    finally:
        pybedtools.helpers.cleanup(remove_all=True)


if __name__ == "__main__":
    main()
