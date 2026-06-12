#!/usr/bin/env python3
"""Count chrX/chrY fragments in one BAM/CRAM against a premade clean BED.

Pass the clean chrX/chrY BED built by prepare_bed.py via --bed. This script only
reads that BED and the alignment file; it does not touch the network.
"""
import argparse
import logging
import os
import pathlib
import subprocess
import sys

import pysam


def _chrom_length(parts: list) -> tuple:
    if len(parts) < 3:
        return 0, 0
    n = int(parts[2]) - int(parts[1])
    if parts[0] in ("chrX", "X"):
        return n, 0
    if parts[0] in ("chrY", "Y"):
        return 0, n
    return 0, 0


def mappable_lengths(bed_path: pathlib.Path):
    lenX = lenY = 0
    with open(bed_path) as fh:
        for line in fh:
            dx, dy = _chrom_length(line.split())
            lenX += dx
            lenY += dy
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


def _strip_chr_from_bed(bed_path: pathlib.Path) -> pathlib.Path:
    bed_nochr = bed_path.with_suffix(".nochr.bed")
    if not bed_nochr.exists():
        logging.info("stripping chr prefix from BED to match BAM")
        with open(bed_path) as fi, open(bed_nochr, "w") as fo:
            for line in fi:
                fo.write(line.replace("chrX", "X").replace("chrY", "Y"))
    return bed_nochr


def main():
    p = argparse.ArgumentParser(
        description="Count chrX/chrY fragments in a BAM/CRAM against a premade clean BED."
    )
    p.add_argument("-b", "--bam",          required=True, metavar="FILE",
                   help="Input BAM/CRAM file")
    p.add_argument("-B", "--bed",          required=True, metavar="FILE",
                   help="Premade clean chrX/chrY BED (built by prepare_bed.py)")
    p.add_argument("-o", "--output",       metavar="FILE",
                   help="Optional TSV output path (always printed to stdout)")
    p.add_argument("-q", "--mapq",         type=int, default=30, metavar="INT",
                   help="Minimum MAPQ (default: 30)")
    p.add_argument("-f", "--include-flag", type=int, default=3, metavar="INT",
                   help="samtools -f flag (default: 3)")
    p.add_argument("-F", "--exclude-flag", type=int, default=3852, metavar="INT",
                   help="samtools -F flag (default: 3852)")
    p.add_argument("-t", "--threads",      type=int,
                   default=min(8, os.cpu_count() or 1), metavar="INT",
                   help="samtools threads (default: min(8, nproc))")
    p.add_argument("--sex-threshold",      type=float, default=20.0, metavar="FLOAT",
                   help="%%Y cut-off: below → Female, at or above → Male (default: 20.0)")
    p.add_argument("-v", "--verbose",      action="store_true",
                   help="Show progress messages")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    try:
        bam_path = pathlib.Path(args.bam)
        bed_path = pathlib.Path(args.bed)

        if not bam_path.exists():
            raise FileNotFoundError(f"BAM/CRAM not found: {bam_path}")
        if not bed_path.exists():
            raise FileNotFoundError(f"BED not found: {bed_path}")

        if not alignment_index_exists(bam_path):
            logging.warning("No index found for %s; samtools may fail or run slowly", bam_path)

        with pysam.AlignmentFile(str(bam_path)) as bf:
            chroms = [sq["SN"] for sq in bf.header.get("SQ", [])]
        has_chr, chromX, chromY = detect_sex_chromosomes(chroms)

        if not has_chr:
            bed_path = _strip_chr_from_bed(bed_path)
        lenX, lenY = mappable_lengths(bed_path)
        logging.info("mappable bp: %s=%d  %s=%d", chromX, lenX, chromY, lenY)

        header = "sample\tpct_chrX\tpct_chrY\tsex"
        x, y = pct_xy(str(bam_path), bed_path, args.mapq, args.include_flag, args.exclude_flag,
                      args.threads, lenX, lenY, chromX, chromY)
        sample = bam_path.name
        if x is None:
            line = f"{sample}\tNA\tNA\tNA"
        else:
            sex = "Female" if y < args.sex_threshold else "Male"
            line = f"{sample}\t{x:.2f}\t{y:.2f}\t{sex}"
        print(header)
        print(line)
        if args.output:
            with open(args.output, "w") as fh:
                fh.write(f"{header}\n{line}\n")

    except Exception as e:
        logging.error("%s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
