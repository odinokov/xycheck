# xycheck

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.7+](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/)
[![samtools](https://img.shields.io/badge/requires-samtools-green.svg)](https://www.htslib.org/)

Counts chrX / chrY reads in one BAM/CRAM file, normalises by the uniquely-mappable X/Y territory, and reports the relative X/Y percentage with an inferred sex call. Useful as a quick sample-identity / sex-mismatch check for WGS, WES and cfDNA libraries.

## Quick start

```bash
pip install pysam                 # plus samtools on PATH
python xycheck.py -b sample.bam -B data/hg38/clean_XY.hg38.k100.bed
```

```
sample	pct_chrX	pct_chrY	sex
sample.bam	99.91	0.09	Female
```

Premade BEDs for hg19 and hg38 ship in `data/`, so no download is needed for the common case.

## Install

**`xycheck.py`** (run step):
- Python 3 with `pysam`
- `samtools`

**`prepare_bed.py`** (one-time BED build):
- additionally `requests`, `pybedtools`
- `bedtools` available to `pybedtools`

```bash
pip install -r requirements.txt          # pysam requests pybedtools
conda install -c bioconda samtools bedtools
```

Python 3.7+ (uses `subprocess.run(capture_output=...)`).

## Usage

The workflow is two steps: build the clean BED once per genome with
`prepare_bed.py`, then run `xycheck.py` per BAM/CRAM. The repo ships premade
clean BEDs for **hg19** and **hg38** under `data/`, so step 1 is only needed for
a different genome, k-mer, or to rebuild.

### 1. Prepare the clean BED (network step, once per genome)

```
python prepare_bed.py -g hg38
```

This downloads and caches Umap mappability (several hundred MB per genome) +
the ENCODE blacklist under `data/<genome>/`, then builds a uniquely-mappable
chrX/chrY BED (score == 1.0, blacklist subtracted) named
`clean_XY.<genome>.k<kmer>.bed`.

| flag              | default       | note                                                       |
| ----------------- | ------------- | ---------------------------------------------------------- |
| `-g/--genome`     | **required**  | hg19 / hg38                                                |
| `-k/--kmer`       | 100           | umap K-mer track (100 / 50 / 36 / 24)                      |
| `-d/--data-dir`   | data          | cache dir, relative to `prepare_bed.py` unless absolute    |
| `-T/--tmp-dir`    | $TMPDIR       | temp dir for pybedtools scratch                            |
| `--score-col`     | auto          | 0-based mappability score column in the umap BED           |
| `-v/--verbose`    |               | show download/filter progress                              |

### 2. Run the check (offline, per sample)

```
python xycheck.py -b sample.bam -B data/hg38/clean_XY.hg38.k100.bed -t 4
```

| flag              | default       | note                                           |
| ----------------- | ------------- | ---------------------------------------------- |
| `-b/--bam`        | **required**  | input BAM/CRAM file                            |
| `-B/--bed`        | **required**  | premade clean chrX/chrY BED from `prepare_bed.py` |
| `-o/--output`     | stdout        | optional TSV report path                       |
| `-t/--threads`    | min(8, nproc) | passed to `samtools -@`                        |
| `-q` `-f` `-F`    | 30 3 3852     | MAPQ, FLAG filters                             |
| `--sex-threshold` | 20.0          | %Y cut-off: below → Female, at or above → Male |
| `-v/--verbose`    |               | show progress                                  |

`xycheck.py` counts properly-paired MAPQ-filtered reads from one BAM/CRAM
overlapping the given clean BED on X and Y, normalises each count by the
mappable bp of its chromosome, and reports the share of X vs Y:

```
pct_chrX = 100 * (nX / lenX) / (nX / lenX + nY / lenY)
pct_chrY = 100 - pct_chrX
```

It does not touch the network. If the alignment uses `X`/`Y` contigs (no
`chr` prefix) the BED is stripped to match, writing a `.nochr.bed` alongside.

Requirements and caveats:

- The BAM/CRAM must be indexed (`.bai` / `.crai`); `samtools view -L` needs it.
- The default `-f 3` keeps only properly-paired reads, so single-end data
  yields zero counts. Use `-f 0` for single-end libraries.
- If samtools fails or no reads pass the filters, the row is
  `sample	NA	NA	NA` and the exit code is 0; on a missing file or unknown
  contigs the script exits 1.

Run one `xycheck.py` process per BAM/CRAM for large cohorts, e.g.

```bash
for b in *.bam; do python xycheck.py -b "$b" -B data/hg38/clean_XY.hg38.k100.bed -o "${b%.bam}.XY.tsv"; done
```

**Output** — TSV is printed to stdout. Use `-o/--output` to also write it to a file.
Sex is called `Female` when `%Y < --sex-threshold` (default 20), `Male` otherwise. Tune `--sex-threshold` if your assay, reference, or filtering strategy shifts the expected XX/XY clusters.

## License

MIT — see [LICENSE](LICENSE).

## References

 - Jeong S, Kim J, Park W, Jeon H, Kim N. SEXCMD: Development and validation of sex marker sequences for whole-exome/genome and RNA sequencing. PLOS ONE 12(9): e0184087 (2017). https://doi.org/10.1371/journal.pone.0184087
 - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
 - Mehran Karimzadeh, Carl Ernst, Anshul Kundaje, Michael M Hoffman, Umap and Bismap: quantifying genome and methylome mappability, Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120, https://doi.org/10.1093/nar/gky677
