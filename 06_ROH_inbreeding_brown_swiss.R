
# ============================================
# COMPLETE ROH ANALYSIS - BROWN SWISS (BSW) BREED
# ============================================
# Includes: ROH calculation, inbreeding stats, boxplots,
# Manhattan plots, common ROH analysis, and visualizations

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(qqman)

# Set working directory
# setwd("YOUR_PATH_HERE")

inputData <- "QC_BSW"

# ============================================
# SECTION 1: QUALITY CONTROL & ROH CALCULATION
# ============================================

cat("=== QUALITY CONTROL ===\n")
system(paste0("plink --bfile ",inputData ," --cow --autosome --nonfounders --allow-no-sex --recode --out arrForROH"))

cat("\n=== LAUNCHING cgaTOH ===\n")
system(paste("TOH_ClusteringSuite_v1_0.exe"))

cat("\n=== ROH DETECTION (5 THRESHOLDS) ===\n")
system(paste("TOH_ClusteringSuite_v1_0.exe -map arrForROH -p arrForROH -o inbreeding/", inputData, "roh_1 -l 15 -min_length 1000000 ",
             "-max_gap 1000000 -max_missing 0 -max_hetero 0 -skip_clustering -force_proceed",sep=""))
system(paste("TOH_ClusteringSuite_v1_0.exe -map arrForROH -p arrForROH -o inbreeding/", inputData, "roh_2 -l 15 -min_length 2000000 ",
             "-max_gap 1000000 -max_missing 0 -max_hetero 0 -skip_clustering -force_proceed",sep=""))
system(paste("TOH_ClusteringSuite_v1_0.exe -map arrForROH -p arrForROH -o inbreeding/", inputData, "roh_4 -l 15 -min_length 4000000 ",
             "-max_gap 1000000 -max_missing 1 -max_hetero 0 -skip_clustering -force_proceed",sep=""))
system(paste("TOH_ClusteringSuite_v1_0.exe -map arrForROH -p arrForROH -o inbreeding/", inputData, "roh_8 -l 15 -min_length 8000000 ",
             "-max_gap 1000000 -max_missing 2 -max_hetero 0 -skip_clustering -force_proceed",sep=""))
system(paste("TOH_ClusteringSuite_v1_0.exe -map arrForROH -p arrForROH -o inbreeding/", inputData, "roh_16 -l 15 -min_length 16000000 ",
             "-max_gap 1000000 -max_missing 4 -max_hetero 1 -skip_clustering -force_proceed",sep=""))

# ============================================
# SECTION 2: INBREEDING ANALYSIS
# ============================================

cat("\n=== INBREEDING COEFFICIENT ANALYSIS ===\n")

tohLengths <- c(1, 2, 4, 8, 16)

# Read BSW ROH results
tmp_BSW <- read.csv(paste0("inbreeding/", inputData, "roh_1.homozygousruns"))
tohResults_BSW <- data.frame(nr=1:length(unique(tmp_BSW$Label)), cow=unique(tmp_BSW$Label))

for (i in 1:length(tohLengths)) {
          tohRuns <- read.csv(paste0("inbreeding/", inputData, "roh_", tohLengths[i], ".homozygousruns"))
          tohSum <- aggregate(tohRuns$PositionLength ~ tohRuns$Label, tohRuns, sum)
          tohSum$inbreedingLevel <- round(tohSum[,2]/2700000000, digits=4)
          colnames(tohSum) <- c("cow", paste("tohTotal_", tohLengths[i], sep=""),
                                paste("inbreedingLevel_", tohLengths[i], sep=""))
          tohResults_BSW <- merge(tohResults_BSW, tohSum, by="cow", sort=F, all.x=T)
}

tohResultsBSW=tohResults_BSW[order(tohResults_BSW$nr),]
write.table(tohResultsBSW,"allTohResultsBSW.txt",quote = F,row.names = F,col.names = T)
cat("✓ ROH analysis complete for BSW (all 5 thresholds)\n")
head(tohResultsBSW)

# Summary statistics
cat("\n=== BSW (BROWN SWISS) INBREEDING SUMMARY ===\n")
for(j in 1:5) {
          col_name <- paste("inbreedingLevel_", tohLengths[j], sep="")
          cat("\nThreshold", tohLengths[j], "Mb:\n")
          cat("  Min:", round(min(tohResults_BSW[[col_name]], na.rm=T), 4), "\n")
          cat("  Q1:", round(quantile(tohResults_BSW[[col_name]], 0.25, na.rm=T), 4), "\n")
          cat("  Median:", round(median(tohResults_BSW[[col_name]], na.rm=T), 4), "\n")
          cat("  Mean:", round(mean(tohResults_BSW[[col_name]], na.rm=T), 4), "\n")
          cat("  Q3:", round(quantile(tohResults_BSW[[col_name]], 0.75, na.rm=T), 4), "\n")
          cat("  Max:", round(max(tohResults_BSW[[col_name]], na.rm=T), 4), "\n")
          cat("  SD:", round(sd(tohResults_BSW[[col_name]], na.rm=T), 4), "\n")
}

# ============================================
# SECTION 3: BOXPLOT
# ============================================

inbreeding_data_BSW <- data.frame(
          threshold = rep(c("1Mb", "2Mb", "4Mb", "8Mb", "16Mb"), each=nrow(tohResults_BSW)),
          inbreeding = c(tohResults_BSW$inbreedingLevel_1,
                         tohResults_BSW$inbreedingLevel_2,
                         tohResults_BSW$inbreedingLevel_4,
                         tohResults_BSW$inbreedingLevel_8,
                         tohResults_BSW$inbreedingLevel_16)
)
inbreeding_data_BSW <- na.omit(inbreeding_data_BSW)
inbreeding_data_BSW$threshold <- factor(inbreeding_data_BSW$threshold,
                                        levels=c("1Mb", "2Mb", "4Mb", "8Mb", "16Mb"))
p_BSW <- ggplot(inbreeding_data_BSW, aes(x=threshold, y=inbreeding, fill=threshold)) +
          geom_boxplot(alpha=0.7, width=0.6) +
          geom_point(position="jitter", size=2, alpha=0.4, color="darkblue") +
          theme_minimal() +
          labs(title="Inbreeding Trend: ROH Threshold 1-16Mb (BROWN SWISS)",
               x="ROH Length Threshold",
               y="Inbreeding Coefficient") +
          theme(legend.position="none")
ggsave("ROH_Inbreeding_BSW_1to16Mb.png", p_BSW, width=12, height=6, dpi=300)
cat("\n✓ Boxplot saved: ROH_Inbreeding_BSW_1to16Mb.png\n")

# ============================================
# SECTION 4: SELECTION SIGNATURES (MANHATTAN)
# ============================================

roi_combined_BSW <- data.frame()
for(thresh in tohLengths) {
          roi <- read.table(paste0("inbreeding/", inputData, "roh_", thresh, ".incidenceofcommonrunsperSNP"),
                            header=TRUE, stringsAsFactors=FALSE)
          
          roi_parsed <- roi %>%
                    separate(SNP.Chr.Position.Column.Runs,
                             c("SNP", "CHR", "POS", "COL", "Incidence"),
                             sep=",") %>%
                    mutate(CHR = as.numeric(CHR),
                           POS = as.numeric(POS),
                           Incidence = as.numeric(Incidence),
                           Threshold = thresh) %>%
                    select(SNP, CHR, POS, Incidence, Threshold)
          
          roi_combined_BSW <- bind_rows(roi_combined_BSW, roi_parsed)
}

# Manhattan plots
for(thresh in tohLengths) {
          roi_data <- roi_combined_BSW %>% filter(Threshold == thresh)
          
          png(paste0("ROH_Manhattan_BSW_", thresh, "Mb.png"), width=1400, height=700)
          manhattan(roi_data, chr="CHR", bp="POS", p="Incidence", logp=FALSE,
                    main=paste("ROH Incidence - BSW Selection Signatures (", thresh, "Mb)"),
                    ylab="ROH Incidence")
          dev.off()
          
          cat("✓ Manhattan plot: ROH_Manhattan_BSW_", thresh, "Mb.png\n")
}
cat("\n✓ BSW (Brown Swiss) ROH analysis complete!\n")

# ============================================
# SECTION 5: SAVE SUMMARY STATISTICS
# ============================================

summary_stats_BSW <- data.frame()
for(j in 1:length(tohLengths)) {
          col_name <- paste("inbreedingLevel_", tohLengths[j], sep="")
          
          stats <- data.frame(
                    Threshold_Mb = tohLengths[j],
                    Min = round(min(tohResultsBSW[[col_name]], na.rm=T), 4),
                    Q1 = round(quantile(tohResultsBSW[[col_name]], 0.25, na.rm=T), 4),
                    Median = round(median(tohResultsBSW[[col_name]], na.rm=T), 4),
                    Mean = round(mean(tohResultsBSW[[col_name]], na.rm=T), 4),
                    Q3 = round(quantile(tohResultsBSW[[col_name]], 0.75, na.rm=T), 4),
                    Max = round(max(tohResultsBSW[[col_name]], na.rm=T), 4),
                    SD = round(sd(tohResultsBSW[[col_name]], na.rm=T), 4),
                    N = sum(!is.na(tohResultsBSW[[col_name]]))
          )
          summary_stats_BSW <- rbind(summary_stats_BSW, stats)
}

write.csv(summary_stats_BSW, "ROH_Summary_Statistics_BSW.csv", row.names=FALSE)
cat("\n✓ Saved: ROH_Summary_Statistics_BSW.csv\n")
print(summary_stats_BSW)

# ============================================
# SECTION 6: COMMON ROH ANALYSIS (≥50%)
# ============================================

cat("\n=== COMMON ROH ANALYSIS - BROWN SWISS ===\n")

n_animals <- nrow(tohResultsBSW)
threshold_pct <- 0.50
threshold_count <- ceiling(n_animals * threshold_pct)

cat("Total animals: ", n_animals, "\n")
cat("50% threshold: ≥", threshold_count, " animals\n\n")

common_roh_summary <- data.frame()

for (thresh in tohLengths) {
          cat("Processing threshold:", thresh, "Mb...\n")
          
          roi_file <- paste0("inbreeding/", inputData, "roh_", thresh, ".incidenceofcommonrunsperSNP")
          roi <- read.table(roi_file, header=TRUE, stringsAsFactors=FALSE)
          
          roi_parsed <- roi %>%
                    separate(SNP.Chr.Position.Column.Runs,
                             c("SNP", "CHR", "POS", "COL", "Incidence"),
                             sep=",") %>%
                    mutate(CHR = as.numeric(CHR),
                           POS = as.numeric(POS),
                           Incidence = as.numeric(Incidence),
                           Incidence_Pct = (Incidence / n_animals) * 100)
          
          common_snps <- roi_parsed %>% filter(Incidence >= threshold_count)
          chromosomes_with_common <- unique(common_snps$CHR)
          
          n_snps_common <- nrow(common_snps)
          n_chromosomes <- length(chromosomes_with_common)
          
          cat("  SNPs with ≥50% incidence:", n_snps_common, "\n")
          cat("  Chromosomes affected:", n_chromosomes, "\n")
          cat("  Chr list:", paste(sort(chromosomes_with_common), collapse=", "), "\n\n")
          
          top_snps <- common_snps %>%
                    arrange(desc(Incidence)) %>%
                    head(10)
          
          cat("  Top 10 SNPs by incidence:\n")
          print(top_snps %>% select(SNP, CHR, POS, Incidence, Incidence_Pct))
          cat("\n")
          
          summary_row <- data.frame(
                    Threshold_Mb = thresh,
                    SNPs_50pct = n_snps_common,
                    Chromosomes = n_chromosomes,
                    Mean_Incidence = round(mean(common_snps$Incidence), 2),
                    Max_Incidence = max(common_snps$Incidence),
                    Chr_List = paste(sort(chromosomes_with_common), collapse=",")
          )
          
          common_roh_summary <- rbind(common_roh_summary, summary_row)
}

cat("\n=== SUMMARY TABLE: SNPs with ≥50% Animal Incidence ===\n")
print(common_roh_summary)

write.csv(common_roh_summary, "Common_ROH_50pct_Summary_BSW.csv", row.names=FALSE)
cat("\n✓ Saved: Common_ROH_50pct_Summary_BSW.csv\n")

# Chromosome-level analysis
cat("\n=== CHROMOSOME-LEVEL ANALYSIS ===\n")
for (thresh in tohLengths) {
          roi_file <- paste0("inbreeding/", inputData, "roh_", thresh, ".incidenceofcommonrunsperSNP")
          roi <- read.table(roi_file, header=TRUE, stringsAsFactors=FALSE)
          
          roi_parsed <- roi %>%
                    separate(SNP.Chr.Position.Column.Runs,
                             c("SNP", "CHR", "POS", "COL", "Incidence"),
                             sep=",") %>%
                    mutate(CHR = as.numeric(CHR),
                           POS = as.numeric(POS),
                           Incidence = as.numeric(Incidence),
                           Incidence_Pct = (Incidence / n_animals) * 100)
          
          common_snps <- roi_parsed %>% filter(Incidence >= threshold_count)
          
          if (nrow(common_snps) > 0) {
                    cat("\n", thresh, "Mb threshold - SNPs by Chromosome:\n")
                    chr_summary <- common_snps %>%
                              group_by(CHR) %>%
                              summarise(
                                        N_SNPs = n(),
                                        Mean_Incidence = round(mean(Incidence), 2),
                                        Max_Incidence = max(Incidence),
                                        .groups = 'drop'
                              ) %>%
                              arrange(desc(N_SNPs))
                    
                    print(chr_summary)
          } else {
                    cat("\n", thresh, "Mb threshold - No SNPs with ≥50% incidence\n")
          }
}

# ============================================
# SECTION 7: VISUALIZATION PLOTS
# ============================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# Load summary data
common_roh_summary <- read.csv("Common_ROH_50pct_Summary_BSW.csv")

# PLOT 1: SNPs and Chromosomes by Threshold
summary_long <- common_roh_summary %>%
          select(Threshold_Mb, SNPs_50pct, Chromosomes) %>%
          pivot_longer(cols = c(SNPs_50pct, Chromosomes),
                       names_to = "Metric", values_to = "Count")

p1 <- ggplot(summary_long, aes(x = as.factor(Threshold_Mb), y = Count, fill = Metric)) +
          geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
          theme_minimal() +
          labs(title = "Common ROH: SNPs and Chromosomes Affected (≥50% Incidence)",
               x = "ROH Length Threshold (Mb)",
               y = "Count",
               fill = "Metric") +
          scale_fill_manual(values = c("SNPs_50pct" = "#1f77b4", "Chromosomes" = "#ff7f0e"),
                            labels = c("SNPs_50pct" = "SNPs", "Chromosomes" = "Chromosomes")) +
          theme(text = element_text(size = 12, face = "bold"))

ggsave("Common_ROH_SNPs_Chromosomes_BSW.png", p1, width = 10, height = 6, dpi = 300)
cat("✓ Plot 1 saved: Common_ROH_SNPs_Chromosomes_BSW.png\n")

# PLOT 2: Threshold Decline
p2 <- ggplot(common_roh_summary, aes(x = as.factor(Threshold_Mb), y = SNPs_50pct)) +
          geom_line(group = 1, linewidth = 1.5, color = "#d62728") +
          geom_point(size = 5, color = "#d62728") +
          geom_text(aes(label = SNPs_50pct), vjust = -1, size = 4) +
          theme_minimal() +
          labs(title = "Decline in Common ROH SNPs with Threshold",
               x = "ROH Length Threshold (Mb)",
               y = "Number of SNPs with ≥50% Incidence",
               subtitle = "Fewer animals share very old inbreeding (16Mb vs 1Mb)") +
          ylim(0, max(common_roh_summary$SNPs_50pct) * 1.2) +
          theme(text = element_text(size = 12, face = "bold"))

ggsave("Common_ROH_Decline_by_Threshold_BSW.png", p2, width = 10, height = 6, dpi = 300)
cat("✓ Plot 2 saved: Common_ROH_Decline_by_Threshold_BSW.png\n")

# PLOT 3: Chromosome Heatmap
chr_data_all <- data.frame()

for (thresh in tohLengths) {
          roi_file <- paste0("inbreeding/", inputData, "roh_", thresh, ".incidenceofcommonrunsperSNP")
          roi <- read.table(roi_file, header = TRUE, stringsAsFactors = FALSE)
          
          roi_parsed <- roi %>%
                    separate(SNP.Chr.Position.Column.Runs,
                             c("SNP", "CHR", "POS", "COL", "Incidence"),
                             sep = ",") %>%
                    mutate(CHR = as.numeric(CHR),
                           POS = as.numeric(POS),
                           Incidence = as.numeric(Incidence))
          
          common_snps <- roi_parsed %>% filter(Incidence >= threshold_count)
          
          chr_summary <- common_snps %>%
                    group_by(CHR) %>%
                    summarise(N_SNPs = n(), .groups = 'drop') %>%
                    mutate(Threshold_Mb = thresh)
          
          chr_data_all <- rbind(chr_data_all, chr_summary)
}

p3 <- ggplot(chr_data_all, aes(x = as.factor(Threshold_Mb), y = as.factor(CHR), fill = N_SNPs)) +
          geom_tile(color = "white", size = 0.5) +
          scale_fill_gradient(low = "white", high = "#2ca02c", name = "SNPs") +
          theme_minimal() +
          labs(title = "Chromosome × Threshold Heatmap: Common ROH SNPs",
               x = "ROH Length Threshold (Mb)",
               y = "Chromosome",
               subtitle = "Darker = more SNPs with ≥50% animal incidence") +
          theme(text = element_text(size = 11, face = "bold"),
                axis.text.y = element_text(size = 9))

ggsave("Common_ROH_Heatmap_BSW.png", p3, width = 10, height = 12, dpi = 300)
cat("✓ Plot 3 saved: Common_ROH_Heatmap_BSW.png\n")

# PLOT 4: Chromosome Ranking
chr_total <- chr_data_all %>%
          group_by(CHR) %>%
          summarise(Total_SNPs = sum(N_SNPs), .groups = 'drop') %>%
          arrange(desc(Total_SNPs)) %>%
          head(10) %>%
          mutate(CHR = as.factor(CHR))

p4 <- ggplot(chr_total, aes(x = reorder(CHR, Total_SNPs), y = Total_SNPs, fill = Total_SNPs)) +
          geom_bar(stat = "identity", alpha = 0.8) +
          scale_fill_gradient(low = "#1f77b4", high = "#d62728", name = "Total SNPs") +
          theme_minimal() +
          labs(title = "Top 10 Chromosomes by Common ROH SNPs (All Thresholds)",
               x = "Chromosome",
               y = "Total SNPs with ≥50% Incidence",
               subtitle = "Brown Swiss Dairy Cattle") +
          coord_flip() +
          theme(text = element_text(size = 12, face = "bold"))

ggsave("Common_ROH_Top_Chromosomes_BSW.png", p4, width = 10, height = 6, dpi = 300)
cat("✓ Plot 4 saved: Common_ROH_Top_Chromosomes_BSW.png\n")

cat("\n=== TOP CHROMOSOMES BY TOTAL COMMON ROH SNPs ===\n")
print(chr_total)

cat("\n✓ ALL BROWN SWISS ROH ANALYSIS COMPLETE!\n")





