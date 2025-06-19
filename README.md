# ycheck
Counts chrX / chrY fragments in BAM/CRAM files and reports the length‑normalised percentage for each chromosome.

```
python ycheck.py -b bam_list.txt -g hg38 -t 4 -T .
```

| flag                | default       | note                               |
| ------------------- | ------------- | ---------------------------------- |
| `-b/--bam-list`     | **required**  | text file, one BAM/CRAM per line   |
| `-g/--genome`       | **required**  | hg19 / hg38                        |
| `-o/--output`       | <list>.XY.tsv | TSV report path                    |
| `-t/--threads`      | all cores     | passed to `samtools -@`            |
| `-k` `-q` `-I` `-E` | 100 60 3 3852 | umap K-mers, MAPQ, FLAG filters    |
| `-T`                | \$TMPDIR      | temp dir                           |

**What it does**

1. Caches Umap mappability + ENCODE blacklist (`data/`).
2. Builds clean chrX/chrY BED (once).
3. For each BAM outputs %X, %Y.


**If you use it, please cite:**

 - Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z
 - Mehran Karimzadeh, Carl Ernst, Anshul Kundaje, Michael M Hoffman, Umap and Bismap: quantifying genome and methylome mappability, Nucleic Acids Research, Volume 46, Issue 20, 16 November 2018, Page e120, https://doi.org/10.1093/nar/gky677
   
