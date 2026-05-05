# this is a script made to be submitted in bash, intended to analyze
# whole genome data and output files along the way

#####GFMxKC

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)
library(CMplot)

# read data
GFMxWM <- read.pcadapt("/home/las80898/mallard_wholegenome_data/GFMxWM.bed", type = "bed")

# initial analysis
x1 <- pcadapt(GFMxKC, K = 2, LD.clumping = list(size = 50000, thr = 0.1)

# outlier adjustment
padjbonf <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.00000005
outliersbonf <- which(padjbonf < alpha)


# plotting with qqman
#make dataframe with values from pcadapt
# Read bim
bim <- read.table("/home/las80898/mallard_wholegenome_data/GFMxKC.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))

# Step 2: get outlier SNP IDs from ORIGINAL unfiltered bim
outlier_snps <- bim$SNP[outliersbonf]

# Start redirecting output to a file
sink("GFMxKC_outliers.txt")
# Your print statement
print(outlier_snps)
# Close the file connection
sink()  

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
qqdf_GFMxWM <- data.frame(SNP, CHR, BP, P)


# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_4/GFMxKC_qqman_T1.png", width = 1024, height = 768, units = "px", pointsize = 14)
manhattan(qqdf_GFMxWM, 
          cex.axis = 0.8, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          genomewideline = -log10(5e-08),
          xlab = "Chromosome number", 
          cex = 0.5, 
          ylim = c(0, 120))
dev.off()









