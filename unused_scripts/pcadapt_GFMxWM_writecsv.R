# GFMxKC - pcadapt whole genome analysis

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)

# read data
GFMxWM <- read.pcadapt("/home/las80898/mallard_wholegenome_data/GFMxWM.bed", type = "bed")

x1 <- pcadapt(GFMxWM, K = 2, LD.clumping = list(size = 5000, thr = 0.1))

bim <- read.table("/home/las80898/mallard_wholegenome_data/GFMxWM.bim",
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

qqdf_GFMxWM <- data.frame(SNP, CHR, BP, P)

write.csv(qqdf_GFMxWM,
          file = "/scratch/las80898/pcadapt_output_4/GFMxWM_dataframe.csv",
          row.names = FALSE,
          quote = FALSE)