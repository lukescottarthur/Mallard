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
GFMxWM <- read.pcadapt("/home/las80898/mallard_wholegenome_data/GFMxWM.bed", type = "bed")

# initial analysis
x1 <- pcadapt(GFMxWM, K = 2)

png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_initial_manhattan.png")
plot(x1 , option = "manhattan")
dev.off()

png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_qqplot.png")
plot(x1, option = "qqplot")
dev.off()

png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_pvalue_histogram.png")
hist(x1$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()

png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_stats_distribution.png")
plot(x1, option = "stat.distribution")
dev.off()

# outlier adjustment
padjbonf <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.0000001
outliersbonf <- which(padjbonf < alpha)
sink("/scratch/las80898/pcadapt_output_2/GFMxWM_outliers.txt")
print(outliersbonf)
sink()

# LD
png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_LD.png",
    width = 800, height = 400 * 2)
par(mfrow = c(2, 1))
for (i in 1:2)
  plot(x1$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
dev.off()

#association between pc and outliers
snp_pc <- get.pc(x1, outliersbonf)
sink("/scratch/las80898/pcadapt_output_2/GFMxWM_snp_pc_associations.txt")
print(snp_pc)
sink()

# plotting with qqman
#make dataframe with values from pcadapt
bim <- read.table("/home/las80898/mallard_wholegenome_data/GFMxWM.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))

message("Unique CHR values in bim: ")
print(unique(bim$CHR))
message("Total SNPs in bim: ", nrow(bim))
message("Total pvalues: ", length(x1$pvalues))

# Get SNP IDs of outliers BEFORE any filtering
outlier_snps <- bim$SNP[outliersbonf]

# Convert CHR and recode
bim$CHR <- as.character(bim$CHR)
bim$CHR[bim$CHR == "chrZ" | bim$CHR == "Z"] <- "30"
bim$CHR <- suppressWarnings(as.numeric(bim$CHR))

# Track which rows to keep (non-scaffold rows)
scaffold_keep <- !is.na(bim$CHR)

# Filter bim AND pvalues together using the same index
bim_filtered  <- bim[scaffold_keep, ]
pvals_filtered <- x1$pvalues[scaffold_keep]

# Now filter NAs from pvalues, again keeping bim aligned
pval_keep <- !is.na(pvals_filtered)
SNP <- bim_filtered$SNP[pval_keep]
CHR <- bim_filtered$CHR[pval_keep]
BP  <- bim_filtered$BP[pval_keep]
P   <- pvals_filtered[pval_keep]

# Sanity check - all should be equal and > 0
message("Lengths after filtering - SNP: ", length(SNP), 
        " CHR: ", length(CHR), 
        " BP: ", length(BP), 
        " P: ", length(P))
message("Any NA in P: ", any(is.na(P)))
message("Any infinite in P: ", any(!is.finite(P)))

qqdf_GFMxWM <- data.frame(SNP, CHR, BP, P)

# verify outliers
sink("/scratch/las80898/pcadapt_output_2/GFMxWM_snp_pc_associations.txt")
message("Number of outliers: ", length(outliersbonf))
print("Number of outliers: ", length(outliersbonf))
message("Number of valid outlier SNP IDs: ", sum(!is.na(outlier_snps)))
print("Number of valid outlier SNP IDs: ", sum(!is.na(outlier_snps)))
sink()

# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_2/GFMxWM_qqman.png")
manhattan(qqdf_GFMxWM, main = " pcadapt SNP Outliers", cex.axis = 0.8, cex.main = .8, annotatePval = 0.0000001, suggestiveline = F , annotateTop = FALSE, xlab = "Chromosome number", cex = 0.3, highlight = outlier_snps)
dev.off()

# plotting with ggplot
# Prepare the dataset
don <- qqdf_GFMxWM %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(qqdf_GFMxWM, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% outlier_snps, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(P)>8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
GFMxWM_ggplot_manhattan <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black"))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), 
                    size=2, max.overlaps = 15) +
  #horizontal line
  geom_hline(yintercept = 8, col="red", alpha = 0.5) +
  
  labs(title = "pcadapt SNP outliers", 
       x = "Chromosome number") +
  
  # Customize theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# saving ggplot
ggsave("/scratch/las80898/pcadapt_output_2/GFMxWM_ggplot_manhattan.png", GFMxWM_ggplot_manhattan, width = 8, height = 6, dpi = 600)

###circular plot - need to amend chr.labels i think
setwd("/scratch/las80898/pcadapt_output_2")
CMplot(qqdf_GFMxWM, plot.type="c", r=1.6,
       outward=TRUE, cir.chr.h=.1, chr.den.col="orange",
       file.name="GFMxWM_circular_manhattan",
       file="jpg", dpi=600, chr.labels=seq(1, 31))
