# Clear workspace
rm(list = ls())

# set up the working directory
setwd("c:/analysis/2020_GenomicsBootCamp_Demo/")

library(tidyverse)

#install.packages("ape", dependencies = TRUE)
library(ape)
#install.packages("StAMPP", dependencies = TRUE)
library(StAMPP)
#install.packages("adegenet", dependencies = TRUE)
library(adegenet)


# select breeds to analyze
# Breeds: Saanen, Toggenburg, Rangeland, Teddi
breeds = c("SAA", "TOG", "RAN", "TED")

write.table(breeds, 'selectedBreeds.txt', row.names = F, col.names = F, quote =F)


# extract breeds and quality control
system(str_c("plink --bfile ADAPTmap_genotypeTOP_20160222_full --chr-set 29 --autosome ",
             "--keep-fam selectedBreeds.txt --nonfounders ",
             "--geno 0.1 --mind 0.1 ",
             "--recode A --out afterQC"))

# prepares a "afterQC.raw" file as output >>> genotype coding

#convert to genlight format
PopToGenlight=read.PLINK(file="afterQC.raw",parallel=FALSE)

#import to StAMPP
goatBreed=stamppConvert(PopToGenlight, "genlight")

# Calculate pairwise Fst values between each population
# takes a LONG time
goatBreed.fst <- stamppFst(goatBreed, 10, 95, 1)

goatBreed.fst

# Plot Phylogenies
arbol <- nj(as.dist(goatBreed.fst))

plot(arbol)

plot.phylo(arbol, "phylogram", edge.width = 2, font = 3, tip.color="blue", cex=1.7)

plot.phylo(arbol, "cladogram", edge.width = 2, font = 3, tip.color="blue", cex=1.7)

plot.phylo(arbol, "fan", edge.width = 2, font = 3, tip.color="blue")

plot.phylo(arbol, "unrooted", edge.width = 2, font = 3, tip.color="blue")

plot.phylo(arbol, "radial", edge.width = 2, font = 3, tip.color="blue", cex=1.7)

plot.phylo(arbol, "tidy", edge.width = 2, font = 3, tip.color="blue", cex=1.7)


##################

# rerun with bootstrapping - 10 iterations to show, too low for real analysis - 100 is default!!!


