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

```
python xycheck.py -b sample.bam -g hg38 -t 4
```

| flag                | default       | note                                           |
| ------------------- | ------------- | ---------------------------------------------- |
| `-b/--bam`          | **required**  | input BAM/CRAM file                            |
| `-g/--genome`       | **required**  | hg19 / hg38                                    |
| `-o/--output`       | stdout        | optional TSV report path                       |
| `-t/--threads`      | all cores     | passed to `samtools -@`                        |
| `-k` `-q` `-f` `-F` | 100 60 3 3852 | umap K-mers, MAPQ, FLAG filters                |
| `-T`                | $TMPDIR       | temp dir for pybedtools scratch                |
| `-d/--data-dir`     | data          | cache dir, relative to `xycheck.py` unless absolute |
| `--sex-threshold`   | 20.0          | %Y cut-off: below → Female, at or above → Male |
| `-v/--verbose`      |               | show download/filter progress                  |

**What it does**

1. Downloads and caches Umap mappability + ENCODE blacklist under `data/`.
2. Builds a uniquely-mappable chrX/chrY BED once (score == 1.0, blacklist subtracted).
3. Counts properly-paired MAPQ-filtered reads from one BAM/CRAM overlapping that BED on X and Y, adjusts by callable bp, and calls sex.

Run one process per BAM/CRAM for large cohorts. By default, the shared `data/` cache lives next to `xycheck.py`, so it is reused even when the script is run from another directory.

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
