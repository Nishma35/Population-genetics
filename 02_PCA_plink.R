# Simple PCA analysis with PLINK

# Clear workspace
rm(list = ls())
# Set working directory
setwd("XXXXXXXXXXXXX")
library(tidyverse)


# Compute a similarity matrix
system("plink --bfile afterQC --cow --distance square ibs --out distanceMatrix_IBS")

# do the PCA with PLINK
system("plink --bfile afterQC --chr-set 29 --pca --out plinkPCA")

###
# Visualize PCA results
###

# read in result files
eigenValues <- read_delim("plinkPCA.eigenval", delim = " ", col_names = F)
eigenVectors <- read_delim("plinkPCA.eigenvec", delim = " ", col_names = F)

## Proportion of variation captured by each vector
eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

# PCA plot
ggplot(data = eigenVectors) +
  geom_point(mapping = aes(x = X3, y = X4, color = X1, shape = X1), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA plot",
      x = paste0("Principal component 1 (",eigen_percent[1,1]," %)"),
      y = paste0("Principal component 2 (",eigen_percent[2,1]," %)"),
      colour = "Breeds", shape = "Breeds") +
  theme_minimal()

