install.packages("genepop")
install.packages("adegenet")
install.packages("hierfstat")
install.packages("ade4")
update.packages(ask = FALSE, checkBuilt = TRUE)

library(adegenet)
library(genepop)
library(hierfstat)
library(ade4)

#Define name of input Genepop file
file <- "Hedgehog_SSR_data.genepop.txt"
list.files()

# Estimate deviations from HWE into the file "NileTilapiaMatrix.genepop.gen.P"
test_HW("Hedgehog_SSR_data.genepop.txt")

# Calculate genetic diversity into the file "NileTilapiaMatrix.genepop.gen.txt.DIV"
genedivFis("Hedgehog_SSR_data.genepop.txt")




# Calculate pairwise FST
library(adegenet)
library(hierfstat)

# Load file using adegenet. Hierfstat requires the file to have a .gen ending
## ncode corresponds to the way genotypes are coded. If an allele corresponds to
## the number of digits used to code an allele (e.g. 12 -> ncode = 2; 243 -> ncode = 3)
gen <- read.genepop("Hedgehog_SSR_data.genepop.gen", ncode = 2)  # or ncode = 3

# Transform genepop to hierfstat format using adegenet
hf <- genind2hierfstat(gen)

#Calculate pairwise FST using hierfstat and save it into a matrix
pairwise_fst <- pairwise.WCfst(hf)
pairwise_fst

#Save pairwise FST into the file "pairwise_FST.csv"
write.csv(pairwise_fst, "pairwise_FST.csv")


#Make PCA
library(ade4)

#Convert the genotype data into a matrix
mat <- tab(gen, NA.method = "mean")

#Make a PCA using ade4
pca <- dudi.pca(
          mat,
          center = TRUE,
          scale = FALSE,
          scannf = FALSE,
          nf = 10
)
nLoc(gen)
nInd(gen)
#Plot PCA
library(ggplot2)

#Make dataframe out of the Eigenvalues
df <- data.frame(
          PC1 = pca$li[,1],
          PC2 = pca$li[,2],
          Population = pop(gen)
)
#Plot scatter plot colored by population
ggplot(df, aes(PC1, PC2, color = Population)) +
          geom_point(size = 3) +
          stat_ellipse() +
          theme_classic()