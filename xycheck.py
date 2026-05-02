#!/usr/bin/env python3
import gzip
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import tarfile
import tempfile

import click
import pybedtools
import pybedtools.helpers
import pysam
import requests

UMAP_URL = "https://bismap.hoffmanlab.org/raw/{asm}.umap.tar.gz"
BL_URL   = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/{asm}-blacklist.v2.bed.gz"
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent


def fetch(url: str, dest: pathlib.Path) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size:
        logging.info("cached  %s", dest)
        return dest
    tmp_dest = dest.with_name(dest.name + ".tmp")
    logging.info("downloading %s", dest.name)
    if tmp_dest.exists():
        tmp_dest.unlink()
    with requests.get(url, stream=True, timeout=(10, 300)) as r:
        r.raise_for_status()
        with open(tmp_dest, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                if chunk:
                    fh.write(chunk)
    tmp_dest.replace(dest)
    return dest


def extract_umap_bed(tar_path: pathlib.Path, k: int, outdir: pathlib.Path) -> pathlib.Path:
    pats = [f"k{k}.umap.bed.gz", f"{k}bp_kmers.bed.gz", f"{k}bp_kmers.bed"]
    # Return cached extraction without opening the (potentially incomplete) tarball
    for pat in pats:
        dst = outdir / pat
        if dst.exists() and dst.stat().st_size:
            logging.info("cached extraction %s", dst)
            return dst
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            if any(m.name.endswith(p) for p in pats):
                dst = outdir / pathlib.Path(m.name).name
                if not dst.exists():
                    with tar.extractfile(m) as src, open(dst, "wb") as out:
                        shutil.copyfileobj(src, out)
                return dst
    raise RuntimeError(f"k={k} track not found in {tar_path.name}")


def gunzip(p: pathlib.Path) -> pathlib.Path:
    out = p.with_suffix("") if p.suffix == ".gz" else p
    if p.suffix == ".gz" and not out.exists():
        with gzip.open(p, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out


def detect_score_column(bed_path: pathlib.Path) -> int:
    opener = gzip.open if bed_path.suffix == ".gz" else open
    rows = []
    with opener(bed_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            rows.append(line.split())
            if len(rows) >= 10:
                break

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


def mappable_lengths(bed_path: pathlib.Path):
    lenX = lenY = 0
    with open(bed_path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            n = int(parts[2]) - int(parts[1])
            chrom = parts[0]
            if chrom in ("chrX", "X"):
                lenX += n
            elif chrom in ("chrY", "Y"):
                lenY += n
    if not lenX or not lenY:
        raise ValueError(f"Zero mappable bp: chrX={lenX} chrY={lenY}")
    return lenX, lenY


def pct_xy(bam, bed, mapq, f_inc, f_exc, threads, lenX, lenY, chromX, chromY):
    base = [
        "samtools", "view", "-c",
        "-q", str(mapq), "-f", str(f_inc), "-F", str(f_exc),
        "-@", str(threads), "-L", str(bed), bam,
    ]
    try:
        cntX = int(subprocess.run(base + [chromX], capture_output=True, text=True, check=True).stdout.strip())
        cntY = int(subprocess.run(base + [chromY], capture_output=True, text=True, check=True).stdout.strip())
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.strip() if e.stderr else "no stderr"
        logging.warning("samtools failed for %s: %s; stderr: %s", bam, e, stderr)
        return None, None
    except ValueError as e:
        logging.warning("samtools returned a non-integer count for %s: %s", bam, e)
        return None, None
    if cntX + cntY == 0:
        return None, None
    normX = cntX / lenX
    normY = cntY / lenY
    total = normX + normY
    pX = round(100 * normX / total, 2)
    return pX, round(100 - pX, 2)


def alignment_index_exists(path: pathlib.Path) -> bool:
    if path.suffix == ".bam":
        return path.with_suffix(".bai").exists() or path.with_suffix(path.suffix + ".bai").exists()
    if path.suffix == ".cram":
        return path.with_suffix(".crai").exists() or path.with_suffix(path.suffix + ".crai").exists()
    return False


def detect_sex_chromosomes(chroms):
    chrom_set = set(chroms)
    if {"chrX", "chrY"}.issubset(chrom_set):
        return True, "chrX", "chrY"
    if {"X", "Y"}.issubset(chrom_set):
        return False, "X", "Y"
    raise ValueError(
        "BAM/CRAM header does not contain chrX/chrY or X/Y contigs"
    )


@click.command()
@click.option("-b", "--bam",          required=True, type=click.Path(exists=True, dir_okay=False),
              help="Input BAM/CRAM file")
@click.option("-o", "--output",       type=click.Path(dir_okay=False))
@click.option("-g", "--genome",       required=True, type=click.Choice(["hg19", "hg38"]))
@click.option("-k", "--kmer",         default="100", show_default=True,
              type=click.Choice(["100", "50", "36", "24"]))
@click.option("-q", "--mapq",         default=60,   show_default=True)
@click.option("-f", "--include-flag", default=3,    show_default=True)
@click.option("-F", "--exclude-flag", default=3852, show_default=True)
@click.option("-t", "--threads",      default=min(8, os.cpu_count() or 1), show_default=True)
@click.option("-T", "--tmp-dir",      default=tempfile.gettempdir(), show_default=True,
              type=click.Path(file_okay=False, writable=True),
              help="Temp dir for pybedtools scratch files (set to SLURM $TMPDIR)")
@click.option("-d", "--data-dir",     default="data", show_default=True,
              type=click.Path(file_okay=False),
              help="Cache dir for umap/blacklist/clean BED files. Relative paths are resolved from xycheck.py")
@click.option("--score-col",          default=None, type=int,
              help="0-based column index of mappability score in umap BED (auto-detected if omitted)")
@click.option("--sex-threshold",      default=20.0, show_default=True, type=float,
              help="%Y cut-off: samples below are called Female, at or above are called Male")
@click.option("-v", "--verbose",      is_flag=True)
def main(bam, output, genome, kmer, mapq, include_flag, exclude_flag,
         threads, tmp_dir, data_dir, score_col, sex_threshold, verbose):

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    os.environ["TMPDIR"] = tmp_dir
    pybedtools.helpers.set_tempdir(tmp_dir)

    data_path = pathlib.Path(data_dir)
    if not data_path.is_absolute():
        data_path = SCRIPT_DIR / data_path

    try:
        bam_path = pathlib.Path(bam)

        bed_path = build_clean_track(genome, int(kmer), data_path, score_col)

        if not alignment_index_exists(bam_path):
            logging.warning("No index found for %s; samtools may fail or run slowly", bam_path)

        with pysam.AlignmentFile(str(bam_path)) as bf:
            chroms = [sq["SN"] for sq in bf.header.get("SQ", [])]
        has_chr, chromX, chromY = detect_sex_chromosomes(chroms)

        if not has_chr:
            bed_nochr = bed_path.with_suffix(".nochr.bed")
            if not bed_nochr.exists():
                logging.info("stripping chr prefix from BED to match BAM")
                with open(bed_path) as fi, open(bed_nochr, "w") as fo:
                    for line in fi:
                        fo.write(line.replace("chrX", "X").replace("chrY", "Y"))
            bed_path = bed_nochr

        lenX, lenY = mappable_lengths(bed_path)
        logging.info("mappable bp: %s=%d  %s=%d", chromX, lenX, chromY, lenY)

        header = "sample\tpct_chrX\tpct_chrY\tsex"
        x, y = pct_xy(str(bam_path), bed_path, mapq, include_flag, exclude_flag,
                      threads, lenX, lenY, chromX, chromY)
        sample = bam_path.name
        if x is None:
            line = f"{sample}\tNA\tNA\tNA"
        else:
            sex = "Female" if y < sex_threshold else "Male"
            line = f"{sample}\t{x:.2f}\t{y:.2f}\t{sex}"
        print(header)
        print(line)
        if output:
            with open(output, "w") as fh:
                fh.write(f"{header}\n{line}\n")

    except Exception as e:
        logging.error("%s", e)
        sys.exit(1)
    finally:
        pybedtools.helpers.cleanup(remove_all=True)


if __name__ == "__main__":
    main()
