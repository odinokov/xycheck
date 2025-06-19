#!/usr/bin/env python3
# ycheck.py  v1.4  (2025-06-19)
#
# MIT License
#
# Counts MAPQ-filtered, properly-paired reads on chrX and chrY for each BAM/CRAM.
# ▸ Caches Umap single-read mappability + ENCODE blacklist (hg19/hg38)
# ▸ Builds a merged, non-overlapping chrX/chrY BED  →  clean_XY.<genome>.bed
# ▸ Processes BAMs sequentially; `-t / --threads` → samtools -@ for fast BGZF
# ▸ Writes <list>.XY.tsv unless -o is given, plus a rich console table.
#
# Quick start
# -----------
#   pip install click requests rich pysam pybedtools
#   python ycheck.py -b bam_list.txt -g hg38 -t 4
# -------------------------------------------------------------------------

import os, sys, pathlib, logging, tarfile, gzip, shutil, tempfile, subprocess, shlex
import requests, click, pysam, pybedtools
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table    import Table
from rich.console  import Console

UMAP_URL = "https://bismap.hoffmanlab.org/raw/{asm}.umap.tar.gz"
BL_URL   = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/{asm}-blacklist.v2.bed.gz"
DATA_DIR = pathlib.Path("data")

# ── helpers ────────────────────────────────────────────────────────────────
def fetch(url: str, dest: pathlib.Path, prog: Progress) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    resume = tmp.stat().st_size if tmp.exists() else 0
    headers = {"Range": f"bytes={resume}-"} if resume else {}
    t = prog.add_task(f"download {dest.name}", total=None)
    with requests.get(url, stream=True, timeout=(5, 120), headers=headers) as r:
        r.raise_for_status()
        mode = "ab" if resume else "wb"
        with open(tmp, mode) as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
                prog.update(t, advance=len(chunk))
    prog.remove_task(t)
    tmp.rename(dest)
    return dest

def extract_umap_bed(tar_path: pathlib.Path, k: int, outdir: pathlib.Path) -> pathlib.Path:
    pats = [f"k{k}.umap.bed.gz", f"{k}bp_kmers.bed.gz", f"{k}bp_kmers.bed"]
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            if any(m.name.endswith(p) for p in pats):
                dst = outdir / pathlib.Path(m.name).name
                with tar.extractfile(m) as src, open(dst, "wb") as out:
                    shutil.copyfileobj(src, out)
                return dst
    raise RuntimeError(f"k={k} track not found in {tar_path.name}")

def gunzip(p: pathlib.Path) -> pathlib.Path:
    if p.suffix != ".gz":
        return p
    out = p.with_suffix('')
    if not out.exists():
        with gzip.open(p, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out

def build_clean_track(asm: str, kmer: int, prog: Progress) -> pathlib.Path:
    wd = DATA_DIR / asm
    umap_tar = fetch(UMAP_URL.format(asm=asm), wd/f"{asm}.umap.tar.gz", prog)
    bl_gz    = fetch(BL_URL .format(asm=asm), wd/f"{asm}.blacklist.bed.gz", prog)
    merged = (pybedtools
              .BedTool(gunzip(extract_umap_bed(umap_tar, kmer, wd)))
              .filter(lambda x: x.chrom in ("chrX", "chrY"))
              .subtract(pybedtools.BedTool(gunzip(bl_gz)))
              .sort()
              .merge())
    out_path = wd / f"clean_XY.{asm}.bed"
    merged.saveas(out_path)
    logging.info(f"{out_path.name} cached")
    return out_path

def pct_xy(bam, bed, mapq, f_inc, f_exc, threads):
    # counts in clean XY intervals
    xy_base = (f"samtools view -c -q {mapq} -f {f_inc} -F {f_exc} "
               f"-@ {threads} -L {shlex.quote(str(bed))} {shlex.quote(bam)}")
    cntX = int(subprocess.check_output(xy_base + " chrX", shell=True).strip())
    cntY = int(subprocess.check_output(xy_base + " chrY", shell=True).strip())
    
    if not cntX or not cntY:
        return (None, None)
    
    total = cntX + cntY

    return (round(cntX * 100 / total, 2), round(cntY * 100 / total, 2))

# ── CLI ────────────────────────────────────────────────────────────────────
@click.command()
@click.option("-b","--bam-list", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("-o","--output", type=click.Path(dir_okay=False))
@click.option("-g","--genome", type=click.Choice(["hg19","hg38"]))
@click.option("-k","--kmer", default=100, show_default=True)
@click.option("-q","--mapq", default=60, show_default=True)
@click.option("-I","--include-flag", default=3, show_default=True)
@click.option("-E","--exclude-flag", default=3852, show_default=True)
@click.option("-t","--threads", default=min(8, os.cpu_count()), show_default="all cores")
@click.option("-T","--tmp-dir", default=tempfile.gettempdir(), show_default=True,
              type=click.Path(file_okay=False, writable=True))
@click.option("-v","--verbose", is_flag=True)
def main(bam_list, output, genome, kmer, mapq,
         include_flag, exclude_flag, threads, tmp_dir, verbose):

    os.environ["TMPDIR"] = tmp_dir
    import pybedtools.helpers as _h; _h.set_tempdir(tmp_dir)

    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING,
                        format="%(levelname)s: %(message)s")
    console = Console()

    bam_paths = [p.strip() for p in pathlib.Path(bam_list).read_text().splitlines() if p.strip()]
    if not bam_paths:
        sys.exit("BAM list is empty")

    if output is None:
        output = f"{pathlib.Path(bam_list).with_suffix('')}.XY.tsv"

    console.print(f"[green]samtools -@ {threads} | TMP {tmp_dir}[/]")

    # genome autodetect (first SQ with AS tag)
    if genome is None:
        genome = "hg38"  # default fallback
        with pysam.AlignmentFile(bam_paths[0]) as bf:
            for sq in bf.header.get('SQ', []):
                asm_tag = sq.get('AS')
                if asm_tag:
                    genome = {"GRCh37": "hg19", "GRCh38": "hg38",
                              "hg19": "hg19", "hg38": "hg38"}.get(asm_tag, genome)
                    break
    console.print(f"[bold]Genome:[/] {genome}")

    # step 1: download / build clean track
    with Progress(SpinnerColumn(), BarColumn(), TextColumn("{task.description}"),
                  TimeElapsedColumn(), transient=True) as dl:
        bed_path = build_clean_track(genome, kmer, dl)

    # step 2: sequential counting
    results = []
    with Progress(SpinnerColumn(), BarColumn(), TextColumn("counting BAMs"),
                  TimeElapsedColumn()) as bar:
        task = bar.add_task("counting BAMs", total=len(bam_paths))
        for bam in bam_paths:
            x, y = pct_xy(bam, bed_path, mapq, include_flag, exclude_flag, threads)
            results.append((bam, x, y))
            bar.update(task, advance=1)

    # step 3: write tsv
    with open(output, "w") as fh:
        fh.write("bam\tpct_chrX\tpct_chrY\n")
        for bam, x, y in results:
            fh.write(f"{bam}\t{x:.2f}\t{y:.2f}\n")
    console.print(f"[cyan]Report written:[/] {output}")

    # pretty table
    tbl = Table(title=f"chrX/chrY fragment %  (MAPQ≥{mapq}, -I {include_flag}, -F {exclude_flag})")
    tbl.add_column("BAM", overflow="fold")
    tbl.add_column("%X", justify="right")
    tbl.add_column("%Y", justify="right")
    for bam, x, y in results:
        tbl.add_row(bam, f"{x:.2f}", f"{y:.2f}")
    console.print(tbl)

    # cleanup pybedtools scratch
    pybedtools.helpers.cleanup(remove_all=True)

if __name__ == "__main__":
    main()
