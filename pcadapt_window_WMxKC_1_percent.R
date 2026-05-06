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

# ── Sliding window analysis on pcadapt results ────────────────────────────────
#
# For each chromosome, a 50 kb window slides in 5 kb steps across all SNPs.
# Windows with mean -log10(P) > 7.301 are retained and written to CSV.
#
window_size <- 50000   # 50 kb
step_size   <-   5000   #   5 kb

# helper to relabel chr 30 back to Z in output
chr_label <- function(chr) ifelse(chr == 30, "Z", as.character(chr))

window_results <- list()

for (chr in sort(unique(qqdf_WMxKC$CHR))) {

  chr_snps <- qqdf_WMxKC %>%
    filter(CHR == chr) %>%
    arrange(BP)

  if (nrow(chr_snps) == 0) next

  bp_min <- min(chr_snps$BP)
  bp_max <- max(chr_snps$BP)

  starts <- seq(bp_min, bp_max, by = step_size)

  chr_windows <- data.frame(
    CHR          = character(),
    window_start = integer(),
    window_end   = integer(),
    mean_neglogP = numeric()
  )

  for (w_start in starts) {
    w_end <- w_start + window_size - 1

    snps_in_window <- chr_snps %>%
      filter(BP >= w_start & BP <= w_end)

    if (nrow(snps_in_window) == 0) next

    mean_neglogP <- mean(-log10(snps_in_window$P))

    chr_windows <- rbind(chr_windows, data.frame(
      CHR          = chr_label(chr),
      window_start = w_start,
      window_end   = w_end,
      mean_neglogP = round(mean_neglogP, 4)
    ))
  }

  if (nrow(chr_windows) > 0) {
    window_results[[chr_label(chr)]] <- chr_windows
  }
}

# combine all windows across chromosomes
window_df <- bind_rows(window_results)

# calculate top 1% threshold across all non-empty windows and filter
top1_threshold <- quantile(window_df$mean_neglogP, probs = 0.99)

window_df_filtered <- window_df %>%
  filter(mean_neglogP >= top1_threshold)

write.csv(window_df_filtered,
          file = "/scratch/las80898/pcadapt_output_4/WMxKC_sliding_windows_1_percent.csv",
          row.names = FALSE,
          quote = FALSE)

cat(sprintf("Sliding window analysis complete\n"))
cat(sprintf("  Total windows scored  : %d\n", nrow(window_df)))
cat(sprintf("  Top 1%% threshold      : %.4f\n", top1_threshold))
cat(sprintf("  Windows retained      : %d\n", nrow(window_df_filtered)))