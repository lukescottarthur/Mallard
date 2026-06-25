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

GFMxKC <- read.pcadapt("/home/las80898/mallard_wholegenome_data/GFMxKC.bed", type = "bed")

GFMxWM <- read.pcadapt("/home/las80898/mallard_wholegenome_data/GFMxWM.bed", type = "bed")

# initial analysis
x1 <- pcadapt(WMxKC, K = 2, LD.clumping = list(size = 5000, thr = 0.1))
x2 <- pcadapt(GFMxKC, K = 2, LD.clumping = list(size = 5000, thr = 0.1))
x3 <- pcadapt(GFMxWM, K = 2, LD.clumping = list(size = 5000, thr = 0.1))

# plotting with qqman
#make dataframe with values from pcadapt
# Read bim
bim <- read.table("/home/las80898/mallard_wholegenome_data/WMxKC.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))
bim2 <- read.table("/home/las80898/mallard_wholegenome_data/GFMxKC.bim",
                  header = FALSE, col.names = c("CHR2","SNP2","CM2","BP2","A12","A22"))
bim3 <- read.table("/home/las80898/mallard_wholegenome_data/GFMxWM.bim",
                  header = FALSE, col.names = c("CHR3","SNP3","CM3","BP3","A13","A23"))

# Step 3: recode and filter scaffolds
bim$CHR <- as.character(bim$CHR)
bim$CHR[bim$CHR == "chrZ" | bim$CHR == "Z"] <- "30"
bim$CHR <- suppressWarnings(as.numeric(bim$CHR))
scaffold_keep <- !is.na(bim$CHR)

bim2$CHR2 <- as.character(bim2$CHR2)
bim2$CHR2[bim2$CHR2 == "chrZ" | bim2$CHR2 == "Z"] <- "30"
bim2$CHR2 <- suppressWarnings(as.numeric(bim2$CHR2))
scaffold_keep2 <- !is.na(bim2$CHR2)

bim3$CHR3 <- as.character(bim3$CHR3)
bim3$CHR3[bim3$CHR == "chrZ" | bim3$CHR3 == "Z"] <- "30"
bim3$CHR3 <- suppressWarnings(as.numeric(bim3$CHR3))
scaffold_keep3 <- !is.na(bim3$CHR3)

# Step 4: filter bim AND pvalues with the same index
bim_filtered   <- bim[scaffold_keep, ]
pvals_filtered <- x1$pvalues[scaffold_keep]

bim_filtered2   <- bim2[scaffold_keep2, ]
pvals_filtered2 <- x2$pvalues[scaffold_keep2]

bim_filtered3   <- bim3[scaffold_keep3, ]
pvals_filtered3 <- x3$pvalues[scaffold_keep3]

# Step 5: filter NAs and replace zero p-values
pval_keep <- !is.na(pvals_filtered)
SNP <- bim_filtered$SNP[pval_keep]
CHR <- bim_filtered$CHR[pval_keep]
BP  <- bim_filtered$BP[pval_keep]
P   <- pvals_filtered[pval_keep]
P[P == 0] <- .Machine$double.xmin

pval_keep2 <- !is.na(pvals_filtered2)
SNP2 <- bim_filtered2$SNP2[pval_keep2]
CHR2 <- bim_filtered2$CHR2[pval_keep2]
BP2  <- bim_filtered2$BP2[pval_keep2]
P2   <- pvals_filtered2[pval_keep2]
P2[P2 == 0] <- .Machine$double.xmin

pval_keep3 <- !is.na(pvals_filtered3)
SNP3 <- bim_filtered3$SNP3[pval_keep3]
CHR3 <- bim_filtered3$CHR3[pval_keep3]
BP3  <- bim_filtered3$BP3[pval_keep3]
P3   <- pvals_filtered3[pval_keep3]
P3[P3 == 0] <- .Machine$double.xmin

# Step 6: build qqdf — outlier_snps already set above, no renaming needed
WMxKC <- data.frame(SNP, CHR, BP, P)

GFMxKC <- data.frame(SNP2, CHR2, BP2, P2)

GFMxWM <- data.frame(SNP3, CHR3, BP3, P3)

# subset data
filtered_data <- subset(WMxKC, CHR %in% c(1, 2, 3, 5, 12, 30))

filtered_data2 <- subset(GFMxKC, CHR2 %in% c(1, 2, 3, 5, 12, 30))

filtered_data3 <- subset(GFMxWM, CHR3 %in% c(1, 2, 3, 5, 12, 30))


# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_4/WMxKC_subset.png", width = 1800, height = 850, units = "px", pointsize = 14)
manhattan(filtered_data, 
          cex.axis = 0.8, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          genomewideline = FALSE,
          xlab = "WMxKC - Chromosome Number", 
          cex = 0.6, 
          ylim = c(0, 150))
dev.off()

png(filename = "/scratch/las80898/pcadapt_output_4/GFMxKC_subset.png", width = 1800, height = 850, units = "px", pointsize = 14)
manhattan(filtered_data2, 
          cex.axis = 0.8, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          genomewideline = FALSE,
          xlab = "GFMxKC - Chromosome Number", 
          cex = 0.6, 
          ylim = c(0, 250),
          chr = "CHR2",
          bp ="BP2",
          snp = "SNP2",
          p = "P2")
dev.off()

png(filename = "/scratch/las80898/pcadapt_output_4/GFMxWM_subset.png", width = 1800, height = 850, units = "px", pointsize = 14)
manhattan(filtered_data3, 
          cex.axis = 0.8, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          genomewideline = FALSE,
          xlab = "GFMxWM - Chromosome Number", 
          cex = 0.6, 
          ylim = c(0, 50),
          chr = "CHR3",
          bp ="BP3",
          snp = "SNP3",
          p = "P3")
dev.off()