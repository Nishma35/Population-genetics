# genomic relatedness


# Clear workspace
rm(list = ls())
# Set working directory
setwd("C:/analysis/2026_dataAnalysisCourse")
# 
library(tidyverse)

## Run specifics
breed1 <- "ANG"
## Create "breed1.txt" and "breed2.txt" files based on R input
writeLines(breed1, "breed1.txt")

##################### 
# LD pruning
#####################
system("plink --bfile data/ADAPTmap_genotypeTOP_20160222_full --chr-set 29 --keep-fam breed1.txt --indep-pairwise 50 5 0.3 --out breed1")

##################### 
# Create subsets
#####################
# Breed 1
system("plink --bfile data/ADAPTmap_genotypeTOP_20160222_full --chr-set 29 --keep-fam breed1.txt --exclude breed1.prune.out --mind 0.1 --geno 0.1 --genome --out breed1")


# visualize IBD values
genome <- read_table("breed1.genome")
hist(genome$PI_HAT)

# Genomic relationship = IBD = PI_HAT

##################### 
# Explanations
#####################
# 
#
#   
#   
#   Probabilities summed up to 1
#   
# Biological interpretation
# 
# For a given SNP:
#   
# IBD = 0 → neither allele came from the same ancestral chromosome.
# IBD = 1 → one allele is inherited from a common ancestor.
# IBD = 2 → both alleles are inherited from the same ancestral source.
# 
# 
#   | Relationship           | Z0   | Z1   | Z2   | PI_HAT |
#   | ---------------------- | ---- | ---- | ---- | ------ |
#   | Parent–offspring       | 0    | 1    | 0    | 0.50   |
#   | Full siblings          | 0.25 | 0.50 | 0.25 | 0.50   |
#   | Half siblings          | 0.50 | 0.50 | 0    | 0.25   |
#   | Grandparent–grandchild | 0.50 | 0.50 | 0    | 0.25   |
#   | First cousins          | 0.75 | 0.25 | 0    | 0.125  |
#   | Unrelated              | ~1   | ~0   | ~0   | ~0     |
#   
#
#
#   | Parameter | Meaning                                            |
#   | --------- | -------------------------------------------------- |
#   | **Z0**    | Probability that the pair shares **0 IBD alleles** |
#   | **Z1**    | Probability that the pair shares **1 IBD allele**  |
#   | **Z2**    | Probability that the pair shares **2 IBD alleles** |
