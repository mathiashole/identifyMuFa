# Blast-and-HMM-searches

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=R&logoColor=white&labelColor=101010)](https://www.r-project.org/about.html)
[![bash](https://img.shields.io/badge/bash-4EAA25?style=for-the-badge&logo=gnubash&logoColor=white&labelColor=101010)](https://www.gnu.org/software/bash/)
[![Perl](https://img.shields.io/badge/Perl-blue?style=for-the-badge&logo=perl&logoColor=white&labelColor=101010)](https://www.perl.org)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mathiashole/identifyMuFa?style=for-the-badge&labelColor=101010&color=white)
![GitHub](https://img.shields.io/github/license/mathiashole/identifyMuFa?color=%23179287&style=for-the-badge&logoColor=white&labelColor=101010&cacheSeconds=60)

`Blast-and-HMM-searches` 

A modular pipeline for homology-based gene family identification, re-annotation, pseudogene detection, and profile HMM searches across genomes and proteomes. The workflow combines BLAST, ORF prediction, HMMER, and custom R scripts to identify complete genes, pseudogenes, and distant homologs from genomic or proteomic datasets.

## :book: Features

- Identification of gene family members using BLAST.
- Reconstruction of fragmented loci.
- Automatic ORF prediction.
- Gene and pseudogene classification based on customizable thresholds.
- Profile HMM construction from multiple sequence alignments.
- Detection of distant homologs using HMMER.
- Support for genomes with or without annotated proteomes.
- Fully configurable for any gene family.

## Modules

### OE1 — Gene Identification and Re-annotation

Identifies homologous loci using BLAST searches, reconstructs fragmented genomic regions, predicts open reading frames (ORFs), and classifies candidates as complete genes or pseudogenes according to user-defined similarity, coverage, and ORF-length thresholds.

#### Main tasks

- BLASTn/BLASTp homology searches
- Merge nearby BLAST hits
- ORF prediction
- Gene boundary reconstruction
- Gene/pseudogene classification
- Export annotated coordinates and sequences

---

### OE2 — Profile HMM Search

Builds profile Hidden Markov Models (HMMs) from multiple sequence alignments and searches genomes or proteomes for distant homologs. For genomes lacking protein annotations, ORFs are automatically extracted prior to HMM screening.

#### Main tasks

- Multiple sequence alignment
- HMM construction
- HMMER searches
- Domain detection
- Candidate protein reconstruction

---

## Dependencies

The pipeline requires the following software:

| Software | Purpose |
|----------|---------|
| BLAST+ | Sequence similarity searches |
| EMBOSS (`transeq`, `getorf`) | ORF prediction and translation |
| HMMER | Profile HMM construction and searches |
| MAFFT | Multiple sequence alignment |
| R (≥ 4.2 recommended) | Data processing and visualization |
| Bash | Workflow execution |


## :wrench: Usage

### Usage OE1

#### 1. Using the Terminal
-  Navigate to the folder `OE1` containing `main.R` and execute the script with the necessary arguments:

```{bash, eval = FALSE}
Rscript main.R --no_gff <config.tsv>
```

#### Example of `config.tsv`

| /Path/to/genome | /Path/to/gene-searches | 900 | 60 |
| /Path/to/other-genome | /Path/to/gene-searches (same or other) | 500 | 80 |
| /Path/to/other-genome | /Path/to/gene-searches (same or other) | 1900 | 50 |

## :hammer: in progress ...


## :sparkling_heart: Contributing

- :octocat: [Pull requests](https://github.com/mathiashole/JustMAQT/pulls) and :star2: stars are always welcome.
- For major changes, please open an [issue](https://github.com/mathiashole/JustMAQT/issues) first to discuss what you would like to change.
- Please make sure to update tests as appropriate.

## :mega: Contact

