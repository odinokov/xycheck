# xycheck

Counts chrX / chrY fragments in one BAM/CRAM file, adjusts counts by the callable X/Y territory, and reports the relative X/Y percentage with an inferred sex call.

## Install

Dependencies:

- Python 3 with `click`, `requests`, `pysam`, `pybedtools`
- `samtools`
- `bedtools` available to `pybedtools`

```bash
pip install click requests pysam pybedtools
```

## Usage

The workflow is two steps: build the clean BED once per genome with
`prepare_bed.py`, then run `xycheck.py` per BAM/CRAM. The repo ships premade
clean BEDs for **hg19** and **hg38** under `data/`, so step 1 is only needed for
a different genome, k-mer, or to rebuild.

### 1. Prepare the clean BED (network step, once per genome)

```
python prepare_bed.py -g hg38
```

This downloads and caches Umap mappability + the ENCODE blacklist under
`data/<genome>/`, then builds a uniquely-mappable chrX/chrY BED (score == 1.0,
blacklist subtracted) named `clean_XY.<genome>.k<kmer>.bed`.

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
| `-t/--threads`    | all cores     | passed to `samtools -@`                        |
| `-q` `-f` `-F`    | 30 3 3852     | MAPQ, FLAG filters                             |
| `--sex-threshold` | 20.0          | %Y cut-off: below → Female, at or above → Male |
| `-v/--verbose`    |               | show progress                                  |

`xycheck.py` counts properly-paired MAPQ-filtered reads from one BAM/CRAM
overlapping the given clean BED on X and Y, adjusts by callable bp, and calls
sex. It does not touch the network. If the alignment uses `X`/`Y` contigs
(no `chr` prefix) the BED is stripped to match, writing a `.nochr.bed` alongside.

Run one `xycheck.py` process per BAM/CRAM for large cohorts.

**Output** — TSV is printed to stdout. Use `-o/--output` to also write it to a file:

```
sample	pct_chrX	pct_chrY	sex
sample.sorted.bam	99.91	0.09	Female
```

Sex is called `Female` when `%Y < --sex-threshold` (default 20), `Male` otherwise. Tune `--sex-threshold` if your assay, reference, or filtering strategy shifts the expected XX/XY clusters.

## Citation

 - Jeong S, Kim J, Park W, Jeon H, Kim N. SEXCMD: Development and validation of sex marker sequences for whole-exome/genome and RNA sequencing. PLOS ONE 12(9): e0184087 (2017). https://doi.org/10.1371/journal.pone.0184087
 - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
 - Mehran Karimzadeh, Carl Ernst, Anshul Kundaje, Michael M Hoffman, Umap and Bismap: quantifying genome and methylome mappability, Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120, https://doi.org/10.1093/nar/gky677
