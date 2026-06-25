# WMxKC - pcadapt whole genome analysis

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)

# read data
WMxKC <- read.pcadapt("/home/las80898/mallard_wholegenome_data/WMxKC.bed", type = "bed")

x1 <- pcadapt(WMxKC, K = 2, LD.clumping = list(size = 5000, thr = 0.1))

bim <- read.table("/home/las80898/mallard_wholegenome_data/WMxKC.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))

bim$CHR <- as.character(bim$CHR)
bim$CHR[bim$CHR == "chrZ" | bim$CHR == "Z"] <- "30"
bim$CHR <- suppressWarnings(as.numeric(bim$CHR))
scaffold_keep <- !is.na(bim$CHR)

bim_filtered   <- bim[scaffold_keep, ]
pvals_filtered <- x1$pvalues[scaffold_keep]

pval_keep <- !is.na(pvals_filtered)
SNP <- bim_filtered$SNP[pval_keep]
CHR <- bim_filtered$CHR[pval_keep]
BP  <- bim_filtered$BP[pval_keep]
P   <- pvals_filtered[pval_keep]
P[P == 0] <- .Machine$double.xmin

qqdf_WMxKC <- data.frame(SNP, CHR, BP, P)

# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_4/WMxKC_qqman_T2.png", 
    width = 1800, height = 850, units = "px", pointsize = 14)
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

# helper to relabel chr 30 back to Z in output
chr_label <- function(chr) ifelse(chr == 30, "Z", as.character(chr))

# ── Output 1: SNPs listed by chromosome ──────────────────────────────────────
sink("/scratch/las80898/pcadapt_output_4/WMxKC_outliers_by_chr.txt")

cat("Significant SNPs grouped by chromosome\n")
cat("Threshold: -log10(P) > 7.301\n")
cat("Total significant SNPs:", nrow(significant_snps), "\n")
cat(strrep("=", 60), "\n\n")

for (chr in sort(unique(significant_snps$CHR))) {
  chr_snps <- significant_snps %>%
    filter(CHR == chr) %>%
    arrange(BP)

  cat(sprintf("Chromosome %s  (%d SNP%s)\n",
              chr_label(chr), nrow(chr_snps), ifelse(nrow(chr_snps) == 1, "", "s")))
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

# ── Output 2: SNP clusters per chromosome (filtered, no top SNP/logP) ────────
#
# A new cluster begins whenever the gap between consecutive SNPs exceeds
# cluster_gap bp. Clusters where start == end (span = 0) are excluded.
#
cluster_gap <- 50000   # 50 kb — adjust to taste

all_clusters <- list()

for (chr in sort(unique(significant_snps$CHR))) {
  chr_snps <- significant_snps %>%
    filter(CHR == chr) %>%
    arrange(BP)

  gaps       <- c(Inf, diff(chr_snps$BP))
  chr_snps$cluster <- cumsum(gaps > cluster_gap)

  cluster_summary <- chr_snps %>%
    group_by(cluster) %>%
    summarise(
      n_snps   = n(),
      start_bp = min(BP),
      end_bp   = max(BP),
      span_bp  = max(BP) - min(BP),
      .groups  = "drop"
    ) %>%
    mutate(CHR = chr_label(chr)) %>%
    filter(span_bp > 0)   # drop singleton clusters

  all_clusters[[chr_label(chr)]] <- cluster_summary
}

summary_df <- bind_rows(all_clusters) %>%
  select(CHR, cluster, n_snps, start_bp, end_bp, span_bp)

# ── Write cluster output ──────────────────────────────────────────────────────
sink("/scratch/las80898/pcadapt_output_4/WMxKC_outlier_clusters.txt")

cat("Inferred SNP clusters by chromosome\n")
cat("Threshold  : -log10(P) > 7.301\n")
cat(sprintf("Cluster gap: %s bp (%.0f kb)\n", format(cluster_gap, big.mark=","), cluster_gap/1e3))
cat("Filter     : Span_BP > 0 (singletons excluded)\n")
cat("Total significant SNPs:", nrow(significant_snps), "\n")
cat(strrep("=", 60), "\n\n")

# per-chromosome detail blocks
for (chr in sort(unique(significant_snps$CHR))) {
  label <- chr_label(chr)
  cl    <- all_clusters[[label]]
  if (is.null(cl) || nrow(cl) == 0) next

  cat(sprintf("Chromosome %s — %d cluster%s\n",
              label, nrow(cl), ifelse(nrow(cl) == 1, "", "s")))
  cat(strrep("-", 60), "\n")

  for (i in seq_len(nrow(cl))) {
    cat(sprintf("  Cluster %d\n",       cl$cluster[i]))
    cat(sprintf("    SNPs    : %d\n",   cl$n_snps[i]))
    cat(sprintf("    Start BP: %s\n",   format(cl$start_bp[i], big.mark=",")))
    cat(sprintf("    End BP  : %s\n",   format(cl$end_bp[i],   big.mark=",")))
    cat(sprintf("    Span    : %s bp (%.2f kb)\n",
                format(cl$span_bp[i], big.mark=","), cl$span_bp[i] / 1e3))
    cat("\n")
  }
}

# summary table
cat(strrep("=", 70), "\n")
cat("Summary table\n")
cat(strrep("=", 70), "\n")
cat(sprintf("%-6s %8s %8s %14s %14s %14s\n",
            "CHR", "Cluster", "N_SNPs", "Start_BP", "End_BP", "Span_BP"))
cat(strrep("-", 70), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-6s %8d %8d %14s %14s %14s\n",
              r$CHR, r$cluster, r$n_snps,
              format(r$start_bp, big.mark=","),
              format(r$end_bp,   big.mark=","),
              format(r$span_bp,  big.mark=",")))
}
# ── Output 3: filtered cluster summary table only ─────────────────────────────
sink("/scratch/las80898/pcadapt_output_4/WMxKC_outlier_clusters_filtered.txt")

cat("Filtered SNP cluster summary\n")
cat("Threshold  : -log10(P) > 7.301\n")
cat(sprintf("Cluster gap: %s bp (%.0f kb)\n", format(cluster_gap, big.mark=","), cluster_gap/1e3))
cat("Filter     : Span_BP > 0 (singletons excluded)\n")
cat(sprintf("Total clusters: %d\n", nrow(summary_df)))
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-6s %8s %8s %14s %14s %14s\n",
            "CHR", "Cluster", "N_SNPs", "Start_BP", "End_BP", "Span_BP"))
cat(strrep("-", 70), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-6s %8d %8d %14s %14s %14s\n",
              r$CHR, r$cluster, r$n_snps,
              format(r$start_bp, big.mark=","),
              format(r$end_bp,   big.mark=","),
              format(r$span_bp,  big.mark=",")))
}

sink()