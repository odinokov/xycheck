#!/usr/bin/env python3
# ycheck.py  v1.6  (2025-06-20)
#
# MIT License
#
# Counts MAPQ-filtered, properly-paired reads on chrX and chrY for each BAM/CRAM.
# ▸ Caches Umap single-read mappability + ENCODE blacklist (hg19/hg38)
# ▸ Builds a merged, non-overlapping chrX/chrY BED  →  clean_XY.<genome>.k<kmer>.bed
# ▸ Normalises counts by total mappable bp per chromosome before computing %
# ▸ Processes BAMs sequentially; `-t / --threads` → samtools -@ for fast BGZF
# ▸ Writes <list>.XY.tsv unless -o is given, plus a rich console table.
#
# Quick start
# -----------
#   pip install click requests rich pysam pybedtools
#   python ycheck.py -b bam_list.txt -g hg38 -t 4
# -------------------------------------------------------------------------

import os, sys, pathlib, logging, tarfile, gzip, shutil, tempfile, subprocess
import requests, click, pysam, pybedtools
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table    import Table
from rich.console  import Console

UMAP_URL = "https://bismap.hoffmanlab.org/raw/{asm}.umap.tar.gz"
BL_URL   = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/{asm}-blacklist.v2.bed.gz"
DATA_DIR = pathlib.Path("data")

# ── helpers ────────────────────────────────────────────────────────────────

def fetch(url: str, dest: pathlib.Path, prog: Progress) -> pathlib.Path:
    """Download file with resume support and error handling"""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size:
        return dest
        
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    resume = tmp.stat().st_size if tmp.exists() else 0
    headers = {"Range": f"bytes={resume}-"} if resume else {}
    t = prog.add_task(f"download {dest.name}", total=None)
    
    try:
        with requests.get(url, stream=True, timeout=(5, 120), headers=headers) as r:
            r.raise_for_status()
            mode = "ab" if resume else "wb"
            with open(tmp, mode) as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
                    prog.update(t, advance=len(chunk))
        tmp.rename(dest)
        return dest
    except requests.RequestException as e:
        logging.error(f"Download failed: {e}")
        if tmp.exists():
            tmp.unlink()
        raise
    finally:
        prog.remove_task(t)

def extract_umap_bed(tar_path: pathlib.Path, k: int, outdir: pathlib.Path) -> pathlib.Path:
    """Extract kmer-specific BED from Umap tarball"""
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
    """Decompress .gz file if needed"""
    if p.suffix != ".gz":
        return p
    out = p.with_suffix('')
    if not out.exists():
        with gzip.open(p, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out

def build_clean_track(asm: str, kmer: int, prog: Progress) -> pathlib.Path:
    """Build cleaned chromosome X/Y BED track with kmer in filename"""
    wd = DATA_DIR / asm
    out_path = wd / f"clean_XY.{asm}.k{kmer}.bed"
    
    # Return existing file if valid
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
        
    umap_tar = fetch(UMAP_URL.format(asm=asm), wd/f"{asm}.umap.tar.gz", prog)
    bl_gz    = fetch(BL_URL .format(asm=asm), wd/f"{asm}.blacklist.bed.gz", prog)
    
    # Process with pybedtools
    umap_bed = gunzip(extract_umap_bed(umap_tar, kmer, wd))
    blacklist = gunzip(bl_gz)
    
    merged = (pybedtools
              .BedTool(umap_bed)
              .filter(lambda x: x.chrom in ("chrX", "chrY"))
              .subtract(pybedtools.BedTool(blacklist))
              .sort()
              .merge())
    merged.saveas(out_path)
    logging.info(f"{out_path.name} created")
    return out_path

def calculate_mappable_lengths(bed_path: pathlib.Path) -> tuple:
    """Calculate mappable bp for X/Y without loading entire file"""
    lenX, lenY = 0, 0
    try:
        with open(bed_path) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[:3]
                length = int(end) - int(start)
                if chrom == "chrX":
                    lenX += length
                elif chrom == "chrY":
                    lenY += length
    except Exception as e:
        logging.error(f"Error reading BED: {e}")
        raise
    
    if lenX == 0 or lenY == 0:
        raise ValueError(f"Invalid mappable lengths: chrX={lenX}, chrY={lenY}")
    
    return (lenX, lenY)

def pct_xy(bam, bed, mapq, f_inc, f_exc, threads, lenX: int, lenY: int, chromX: str, chromY: str):
    """Calculate normalized X/Y percentages with explicit chromosome names"""
    cmd_base = [
        "samtools", "view", "-c",
        "-q", str(mapq),
        "-f", str(f_inc),
        "-F", str(f_exc),
        "-@", str(threads),
        "-L", str(bed),
        bam
    ]
    
    try:
        # Count chromosome X
        resultX = subprocess.run(
            cmd_base + [chromX],
            capture_output=True, text=True, check=True
        )
        cntX = int(resultX.stdout.strip())
        
        # Count chromosome Y
        resultY = subprocess.run(
            cmd_base + [chromY],
            capture_output=True, text=True, check=True
        )
        cntY = int(resultY.stdout.strip())
    except subprocess.CalledProcessError as e:
        logging.error(f"samtools failed for {bam}: {e.stderr}")
        return (None, None)
    except ValueError:
        logging.error(f"Invalid output from samtools for {bam}")
        return (None, None)

    if (cntX + cntY) == 0:
        return (None, None)

    normX = cntX / lenX
    normY = cntY / lenY
    total = normX + normY

    pctX = round(100 * normX / total, 2)
    pctY = round(100 - pctX, 2)
    return (pctX, pctY)

# ── CLI ────────────────────────────────────────────────────────────────────
@click.command()
@click.option("-b","--bam-list", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("-o","--output", type=click.Path(dir_okay=False))
@click.option("-g","--genome", type=click.Choice(["hg19","hg38"]))
@click.option("-k","--kmer", default=100, show_default=True)
@click.option("-q","--mapq", default=60, show_default=True)
@click.option("-f","--include-flag", default=3, show_default=True)
@click.option("-F","--exclude-flag", default=3852, show_default=True)
@click.option("-t","--threads", default=min(8, os.cpu_count() or 1), show_default=True)
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

    try:
        bam_paths = [p.strip() for p in pathlib.Path(bam_list).read_text().splitlines() if p.strip()]
        if not bam_paths:
            sys.exit("BAM list is empty")

        if output is None:
            output = f"{pathlib.Path(bam_list).with_suffix('')}.XY.tsv"

        console.print(f"[green]samtools -@ {threads} | TMP {tmp_dir}[/]")

        # Genome autodetection from first BAM
        if genome is None:
            genome = "hg38"  # default fallback
            with pysam.AlignmentFile(bam_paths[0]) as bf:
                chroms = [sq['SN'] for sq in bf.header.get('SQ', [])]
                has_chr = any(c.startswith('chr') for c in chroms)
                for sq in bf.header.get('SQ', []):
                    asm_tag = sq.get('AS')
                    if asm_tag:
                        genome = {"GRCh37": "hg19", "GRCh38": "hg38",
                                  "hg19": "hg19", "hg38": "hg38"}.get(asm_tag, genome)
                        break
        console.print(f"[bold]Genome:[/] {genome}")

        # Step 1: Download/build clean track
        with Progress(SpinnerColumn(), BarColumn(), TextColumn("{task.description}"),
                      TimeElapsedColumn(), transient=True) as dl:
            bed_path = build_clean_track(genome, kmer, dl)

        # Step 1b: Detect chromosome naming convention
        with pysam.AlignmentFile(bam_paths[0]) as bf:
            chroms = [sq['SN'] for sq in bf.header.get('SQ', [])]
            has_chr = any(c.startswith('chr') for c in chroms)
            chromX = "chrX" if has_chr else "X"
            chromY = "chrY" if has_chr else "Y"

        # Step 1c: Calculate mappable lengths
        lenX, lenY = calculate_mappable_lengths(bed_path)
        console.print(f"[blue]Mappable bp → {chromX}:{lenX:,}  {chromY}:{lenY:,}[/]")

        # Step 2: Process BAMs
        results = []
        with Progress(SpinnerColumn(), BarColumn(), TextColumn("counting BAMs"),
                      TimeElapsedColumn()) as bar:
            task = bar.add_task("counting BAMs", total=len(bam_paths))
            for bam in bam_paths:
                x, y = pct_xy(
                    bam, bed_path, mapq, include_flag, exclude_flag, 
                    threads, lenX, lenY, chromX, chromY
                )
                results.append((bam, x, y))
                bar.update(task, advance=1)

        # Step 3: Write TSV output
        with open(output, "w") as fh:
            fh.write("bam\tpct_chrX\tpct_chrY\n")
            for bam, x, y in results:
                fh.write(f"{bam}\t")
                if x is None or y is None:
                    fh.write("NA\tNA\n")
                else:
                    fh.write(f"{x:.2f}\t{y:.2f}\n")
        console.print(f"[cyan]Report written:[/] {output}")

        # Step 4: Display pretty table
        tbl = Table(title=f"chrX/chrY fragment %  (MAPQ≥{mapq}, -I {include_flag}, -E {exclude_flag})")
        tbl.add_column("BAM", overflow="fold")
        tbl.add_column("%X", justify="right")
        tbl.add_column("%Y", justify="right")
        for bam, x, y in results:
            if x is None or y is None:
                tbl.add_row(bam, "NA", "NA")
            else:
                tbl.add_row(bam, f"{x:.2f}", f"{y:.2f}")
        console.print(tbl)

    except Exception as e:
        console.print(f"[bold red]Error:[/] {e}")
        sys.exit(1)
    finally:
        # Cleanup temporary files
        pybedtools.helpers.cleanup(remove_all=True)

if __name__ == "__main__":
    main()
    