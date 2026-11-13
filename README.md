# Population Genetics in R

R scripts from my second-semester coursework in population and conservation genetics. They cover genetic diversity, Hardy–Weinberg tests, population structure (PCA), genomic relatedness, F<sub>ST</sub> and runs of homozygosity (ROH) in livestock and wild species.

## Scripts

| Script | Topic | Data / species |
|---|---|---|
| `01_genepop_HWE_diversity_FST_PCA.R` | HWE tests, gene diversity and F<sub>IS</sub> (`genepop`), pairwise Weir & Cockerham F<sub>ST</sub> (`hierfstat`), PCA (`ade4`) | Hedgehog microsatellites (SSR) |
| `02_PCA_plink.R` | IBS distance matrix and PCA with PLINK, PCA plot | SNP data after QC |
| `03_genomic_relatedness_IBD.R` | LD pruning, IBD / PI_HAT estimation with PLINK | Goats (ADAPTmap) |
| `04_multipopulation_FST_goats.R` | Pairwise F<sub>ST</sub> between breeds (`StAMPP`), neighbour-joining trees (`ape`) | Saanen, Toggenburg, Rangeland, Teddi goats |
| `05_FST_dairy_vs_beef_candidate_genes.R` | Breed-wise QC, merging, genome-wide F<sub>ST</sub> Dairy vs Beef, Manhattan plots, candidate genes from Ensembl (`biomaRt`) | Holstein, Jersey, Brown Swiss vs Angus, Brahman, Wagyu cattle |
| `06_ROH_inbreeding_brown_swiss.R` | ROH detection at 1–16 Mb (cgaTOH), F<sub>ROH</sub> inbreeding, ROH-island Manhattan plots, common-ROH analysis | Brown Swiss cattle |

## Tools and packages

- **External:** [PLINK 1.9](https://www.cog-genomics.org/plink/), cgaTOH (`TOH_ClusteringSuite`)
- **R:** `tidyverse`, `ggplot2`, `qqman`, `genepop`, `adegenet`, `hierfstat`, `ade4`, `StAMPP`, `ape`, `biomaRt`

## Data

The genotype data are **not included**. They came from course material and public datasets (e.g. ADAPTmap). To run a script, put the matching PLINK/Genepop files in your working directory and update `setwd()`.

## Credits

Scripts 02–04 started from lecture and *Genomics Boot Camp* tutorial code, and I adapted them for the course exercises.
