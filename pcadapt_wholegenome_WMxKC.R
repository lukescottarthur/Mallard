# this is a script made to be submitted in bash, intended to analyze
# whole genome data and output files along the way

# GFMxWM

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)
library(CMplot)

# read data
WMxKC <- read.pcadapt("/home/las80898/mallard_wholegenome_data/WMxKC.bed", type = "bed")

# initial analysis
x1 <- pcadapt(WMxKC, K = 2, LD.clumping = list(size = 5000, thr = 0.1))


# plotting with qqman
#make dataframe with values from pcadapt
# Read bim
bim <- read.table("/home/las80898/mallard_wholegenome_data/WMxKC.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))


# Step 3: recode and filter scaffolds
bim$CHR <- as.character(bim$CHR)
bim$CHR[bim$CHR == "chrZ" | bim$CHR == "Z"] <- "30"
bim$CHR <- suppressWarnings(as.numeric(bim$CHR))
scaffold_keep <- !is.na(bim$CHR)

# Step 4: filter bim AND pvalues with the same index
bim_filtered   <- bim[scaffold_keep, ]
pvals_filtered <- x1$pvalues[scaffold_keep]

# Step 5: filter NAs and replace zero p-values
pval_keep <- !is.na(pvals_filtered)
SNP <- bim_filtered$SNP[pval_keep]
CHR <- bim_filtered$CHR[pval_keep]
BP  <- bim_filtered$BP[pval_keep]
P   <- pvals_filtered[pval_keep]
P[P == 0] <- .Machine$double.xmin

# Step 6: build qqdf — outlier_snps already set above, no renaming needed
qqdf_WMxKC <- data.frame(SNP, CHR, BP, P)


# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_4/WMxKC_qqman_T2.png", width = 1800, height = 850, units = "px", pointsize = 14)
manhattan(qqdf_WMxKC, 
          cex.axis = 0.8, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          genomewideline = FALSE,
          xlab = "Chromosome number", 
          cex = 0.6, 
          ylim = c(0, 150))
dev.off()

# list of sig snps
significant_snps <- qqdf_WMxKC %>%
  filter(-log10(P) > 7.301)

# ── New Output 1: SNPs grouped and listed by chromosome ──────────────────────
sink("/scratch/las80898/pcadapt_output_4/WMxKC_outliers_by_chr.txt")

cat("Significant SNPs grouped by chromosome\n")
cat("Threshold: -log10(P) > 7.301\n")
cat("Total significant SNPs:", nrow(significant_snps), "\n")
cat(strrep("=", 60), "\n\n")

chr_list <- sort(unique(significant_snps$CHR))

for (chr in chr_list) {
  chr_snps <- significant_snps %>%
    filter(CHR == chr) %>%
    arrange(BP)

  chr_label <- ifelse(chr == 30, "Z", as.character(chr))

  cat(sprintf("Chromosome %s  (%d SNP%s)\n",
              chr_label, nrow(chr_snps), ifelse(nrow(chr_snps) == 1, "", "s")))
  cat(strrep("-", 40), "\n")
  cat(sprintf("  %-25s %12s %10s\n", "SNP", "BP", "-log10(P)"))

  for (i in seq_len(nrow(chr_snps))) {
    cat(sprintf("  %-25s %12d %10.3f\n",
                chr_snps$SNP[i],
                chr_snps$BP[i],
                -log10(chr_snps$P[i])))
  }
  cat("\n")
}

sink()

# ── New Output 2: inferred SNP clusters per chromosome ───────────────────────
#
# Clustering logic:
#   SNPs are sorted by BP within each chromosome. A new cluster begins
#   whenever the gap between consecutive SNPs exceeds `cluster_gap`.
#   Each cluster is summarised by its start BP, end BP, span, and SNP count.
#
cluster_gap <- 50000   # 50 kb — adjust to taste

sink("/scratch/las80898/pcadapt_output_4/WMxKC_outlier_clusters.txt")

cat("Inferred SNP clusters by chromosome\n")
cat("Threshold : -log10(P) > 7.301\n")
cat(sprintf("Cluster gap: %s bp (%.0f kb)\n", format(cluster_gap, big.mark=","), cluster_gap/1e3))
cat("Total significant SNPs:", nrow(significant_snps), "\n")
cat(strrep("=", 60), "\n\n")

all_clusters <- list()   # collect for a summary table at the end

for (chr in sort(unique(significant_snps$CHR))) {
  chr_snps <- significant_snps %>%
    filter(CHR == chr) %>%
    arrange(BP)

  chr_label <- ifelse(chr == 30, "Z", as.character(chr))

  # assign cluster IDs
  gaps      <- c(Inf, diff(chr_snps$BP))   # Inf ensures first SNP starts cluster 1
  cluster_id <- cumsum(gaps > cluster_gap)

  chr_snps$cluster <- cluster_id

  # summarise clusters
  cluster_summary <- chr_snps %>%
    group_by(cluster) %>%
    summarise(
      n_snps   = n(),
      start_bp = min(BP),
      end_bp   = max(BP),
      span_bp  = max(BP) - min(BP),
      top_snp  = SNP[which.min(P)],
      top_neglogP = -log10(min(P)),
      .groups  = "drop"
    ) %>%
    mutate(CHR = chr_label)

  all_clusters[[chr_label]] <- cluster_summary

  cat(sprintf("Chromosome %s — %d cluster%s\n",
              chr_label,
              nrow(cluster_summary),
              ifelse(nrow(cluster_summary) == 1, "", "s")))
  cat(strrep("-", 60), "\n")

  for (i in seq_len(nrow(cluster_summary))) {
    cl <- cluster_summary[i, ]
    cat(sprintf("  Cluster %d\n", i))
    cat(sprintf("    SNPs      : %d\n", cl$n_snps))
    cat(sprintf("    Start BP  : %s\n", format(cl$start_bp, big.mark=",")))
    cat(sprintf("    End BP    : %s\n", format(cl$end_bp,   big.mark=",")))
    cat(sprintf("    Span      : %s bp (%.2f kb)\n",
                format(cl$span_bp, big.mark=","), cl$span_bp / 1e3))
    cat(sprintf("    Top SNP   : %s  (-log10P = %.3f)\n",
                cl$top_snp, cl$top_neglogP))
    cat("\n")
  }
}

# ── Summary table across all chromosomes ─────────────────────────────────────
cat(strrep("=", 60), "\n")
cat("Summary table\n")
cat(strrep("=", 60), "\n")

summary_df <- bind_rows(all_clusters) %>%
  select(CHR, cluster, n_snps, start_bp, end_bp, span_bp, top_snp, top_neglogP)

cat(sprintf("%-6s %8s %8s %14s %14s %14s  %-25s %10s\n",
            "CHR", "Cluster", "N_SNPs", "Start_BP", "End_BP", "Span_BP",
            "Top_SNP", "-log10P"))
cat(strrep("-", 100), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-6s %8d %8d %14s %14s %14s  %-25s %10.3f\n",
              r$CHR, r$cluster, r$n_snps,
              format(r$start_bp, big.mark=","),
              format(r$end_bp,   big.mark=","),
              format(r$span_bp,  big.mark=","),
              r$top_snp,
              r$top_neglogP))
}

sink()
