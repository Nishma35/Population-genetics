# FST Analysis with Individual Breed QC (No MAF/HWE filtering)
# Then merge as Production Types (Dairy vs Beef)

rm(list = ls())
library(tidyverse)
library(qqman)

# ============================================================================
# STEP 1: Prepare breed-specific files
# ============================================================================
print("Step 1: Preparing breed-specific QC files...")

cattle <- read.table("allcattle.fam", header = FALSE)

# Create separate files for each breed
breeds <- c("HO", "JER", "BSW", "AN", "BR", "WAGY")

for (breed in breeds) {
          breed_samples <- cattle %>%
                    filter(V1 == breed) %>%
                    select(V1, V2)
          
          write_delim(breed_samples,
                      paste0(breed, "_samples.txt"),
                      col_names = FALSE, delim = "\t")
          
          cat("✓ Created", breed, "sample file\n")
}

# ============================================================================
# STEP 2: Individual QC per breed (ONLY missing data filters)
# NO MAF, NO HWE
# ============================================================================
print("\nStep 2: Running breed-specific QC (no MAF/HWE)...")

for (breed in breeds) {
          cmd <- str_c(
                    "plink --bfile allcattle --chr-set 29 --autosome ",
                    "--keep ", breed, "_samples.txt --nonfounders ",
                    "--geno 0.1 --mind 0.1 ",
                    "--make-bed --out QC_", breed
          )
          system(cmd)
          cat("✓ QC completed for", breed, "\n")
}

# ============================================================================
# STEP 3: Merge all breed-QC files
# ============================================================================
print("\nStep 3: Merge breeds within group (Dairy & Beef)...")

# ===== MERGE DAIRY BREEDS =====
cat("\n--- Merging Dairy Breeds ---\n")

# Merge HO + BSW
system("plink --bfile QC_HO --cow --autosome --bmerge QC_BSW --nonfounders --allow-no-sex --make-bed --out Dairy_HO_BSW")

# Merge HO+BSW + JER
system("plink --bfile Dairy_HO_BSW --cow --autosome --bmerge QC_JER --nonfounders --allow-no-sex --make-bed --out mergedDairy")

cat("✓ Dairy breeds merged (HO + BSW + JER)\n")

# ===== MERGE BEEF BREEDS =====
cat("\n--- Merging Beef Breeds ---\n")

# Merge AN + BR
system("plink --bfile QC_AN --cow --autosome --bmerge QC_BR --nonfounders --allow-no-sex --make-bed --out Beef_AN_BR")

# Merge AN+BR + WAGY
system("plink --bfile Beef_AN_BR --cow --autosome --bmerge QC_WAGY --nonfounders --allow-no-sex --make-bed --out mergedBeef")

cat("✓ Beef breeds merged (AN + BR + WAGY)\n")

# ===== MERGE DAIRY + BEEF (ROBUST APPROACH) =====
cat("\n--- Merging Dairy + Beef (Robust) ---\n")

# Step 1: Initial merge
system("plink --bfile mergedDairy --cow --autosome --bmerge mergedBeef --nonfounders --allow-no-sex --make-bed --out DairyBeef_Combined")

# Step 2: Identify duplicates
system("plink --bfile mergedDairy --cow --autosome --bmerge mergedBeef --list-duplicate-vars ids-only suppress-first --nonfounders --allow-no-sex --make-bed --out DairyBeef_Combined")

# Step 3: Remove duplicates
system("plink --bfile DairyBeef_Combined --cow --autosome --exclude DairyBeef_Combined.dupvar --nonfounders --allow-no-sex --make-bed --out DairyBeef_Combined_Clean")

cat("✓ Dairy + Beef merged (duplicates removed)\n")
cat("✓ Ready for FST: use DairyBeef_Combined_Clean\n")




# compute Fst
dairy_fam <- read.table("mergedDairy.fam") %>% 
          select(V1, V2) %>%
          mutate(V2 = as.character(V2))
dairy_fam$population <- 1

beef_fam <- read.table("mergedBeef.fam") %>% 
          select(V1, V2) %>%
          mutate(V2 = as.character(V2))
beef_fam$population <- 2

# Combine
pop_file <- bind_rows(dairy_fam, beef_fam)

# Save
write.table(pop_file, "prodType_forFST.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cat("✓ Created prodType_forFST.txt\n")
head(pop_file)
#####################
system("plink --cow --bfile DairyBeef_Combined_Clean --fst --within prodType_forFST.txt --out FstResults")

# Load results
fstResults <- read.delim("FstResults.fst") %>% drop_na()

cat("FST Summary:\n")
cat("SNPs:", nrow(fstResults), "\n")
cat("Mean FST:", round(mean(fstResults$FST), 4), "\n")
cat("Max FST:", round(max(fstResults$FST), 4), "\n\n")

# Manhattan plot
png("FST_Manhattan.png", width = 1400, height = 700)
manhattan(fstResults, chr = "CHR", bp = "POS", p = "FST", logp = FALSE,
          main = "Genome-wide FST: Dairy vs Beef Cattle",
          ylab = "FST")
dev.off()

cat("✓ Plot saved\n")


# Top 0.0307% (FST >= 0.6)
threshold <- 0.6
top_snps <- fstResults %>% filter(FST >= threshold) %>% arrange(desc(FST))

cat("Top 0.0307% SNPs (FST ≥ 0.6):\n")
cat("Count:", nrow(top_snps), "\n\n")
print(top_snps)

# Add threshold line to plot
png("FST_Manhattan_WithThreshold.png", width = 1400, height = 700)
manhattan(fstResults, chr = "CHR", bp = "POS", p = "FST", logp = FALSE,
          genomewideline = threshold,
          main = "Genome-wide FST: Dairy vs Beef Cattle (Top 0.0307%)",
          ylab = "FST")
dev.off()

cat("\n✓ Plot with threshold saved\n")

############################Huntfor genes#####################################
######################################
# Explore genomics regions           #
######################################

# Script file for the video "Hunt for genes! | How to extract gene content of genomic regions" on the Genomics Boot Camp YouTube channel

# Clear workspace
rm(list = ls())
#########################
# Regions of interest (RoI)
#########################
# Create 500kb windows around each FST SNP (±250kb)

segmentsBiomart <- fstResults %>%
          filter(FST >= 0.6) %>%
          mutate(RoI = paste(CHR, POS-250000, POS+250000, sep = ":")) %>%
          dplyr::select(RoI) %>%
          t()

write_delim(as.data.frame(t(segmentsBiomart)), "RoI4BioMart.txt", col_names = FALSE, delim = ",")

cat("✓ Created RoI file with", length(segmentsBiomart), "regions\n\n")

################################
# Extract genomic regions with biomaRt
################################
if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

BiocManager::install("biomaRt", quietly = TRUE)
library("biomaRt")

# Connect to Ensembl
cat("Connecting to Ensembl...\n")
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl")

# Query genes in regions of interest
cat("Extracting genes from candidate regions...\n")
biomaRt_results <- getBM(
          attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                         "chromosome_name", "start_position", "end_position", "description"),
          filters = c("chromosomal_region", "biotype"),
          values = list(chromosomal_region = segmentsBiomart, biotype = "protein_coding"),
          mart = ensembl
)

cat("✓ Found", nrow(biomaRt_results), "protein-coding genes\n\n")

# Write full results
write_delim(biomaRt_results, "resultsFromBioMart_FST.txt", delim = ";")
cat("✓ Saved: resultsFromBioMart_FST.txt\n")

################################
# Post-processing: Extract gene names
################################
geneNames <- biomaRt_results %>%
          drop_na(external_gene_name) %>%
          dplyr::select(external_gene_name, chromosome_name, start_position, end_position, description) %>%
          unique() %>%
          arrange(chromosome_name)

write_delim(geneNames, "listOfGenes_FST.txt", delim = ";")

cat("\n=== CANDIDATE GENES (FST ≥ 0.6) ===\n")
print(geneNames)

########################print genes#############################

# Summary by chromosome
gene_summary <- geneNames %>%
          group_by(chromosome_name) %>%
          summarise(
                    count = n(),
                    genes = paste(external_gene_name, collapse = ", ")
          ) %>%
          arrange(desc(count))

cat("\n=== CANDIDATE GENES BY CHROMOSOME ===\n")
print(gene_summary)

# Highlight key genes (known to affect production traits)
key_genes <- c("FGF10", "NF2", "LIF", "OSM", "TFEB", "PKM", "CCND3", "PLIN3", "PLIN4", "PLIN5")

interesting <- geneNames %>%
          filter(external_gene_name %in% key_genes)

cat("\n=== KEY CANDIDATE GENES (Production/Metabolic Traits) ===\n")
print(interesting)

# Save all results
write.csv(geneNames, "Candidate_Genes_FST_Full.csv", row.names = FALSE)
write.csv(gene_summary, "Candidate_Genes_Summary.csv", row.names = FALSE)

cat("\n✓ Saved: Candidate_Genes_FST_Full.csv\n")
cat("✓ Saved: Candidate_Genes_Summary.csv\n")



# Top 5 SNPs with highest FST
top5_snps <- fstResults %>%
          filter(FST >= 0.6) %>%
          arrange(desc(FST)) %>%
          head(5)

cat("=== TOP 5 SNPS BY FST ===\n")
print(top5_snps)

# Genes from TOP 5 SNPs only
top5_snp_chrs <- top5_snps$CHR

genes_top5 <- geneNames %>%
          filter(chromosome_name %in% top5_snp_chrs) %>%
          arrange(chromosome_name)

cat("\n=== GENES FROM TOP 5 SNPS ===\n")
cat("Chromosomes:", paste(unique(top5_snp_chrs), collapse = ", "), "\n")
cat("Total genes:", nrow(genes_top5), "\n\n")
print(genes_top5)

# Save top5 genes
write.csv(genes_top5, "Candidate_Genes_Top5_SNPs.csv", row.names = FALSE)
write.csv(top5_snps, "Top5_SNPs_FST.csv", row.names = FALSE)

cat("\n✓ Saved: Candidate_Genes_Top5_SNPs.csv\n")
cat("✓ Saved: Top5_SNPs_FST.csv\n")












