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

# outlier adjustment
padjbonf <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.0000001
outliersbonf <- which(padjbonf < alpha)


# plotting with qqman
#make dataframe with values from pcadapt
# Read bim
bim <- read.table("/home/las80898/mallard_wholegenome_data/GFMxWM.bim",
                  header = FALSE, col.names = c("CHR","SNP","CM","BP","A1","A2"))

bim$SNP <- paste0("snp", seq_len(nrow(bim)))

message("Unique CHR values in bim: ")
print(unique(bim$CHR))
message("Total SNPs in bim: ", nrow(bim))
message("Total pvalues: ", length(x1$pvalues))

# Step 2: get outlier SNP IDs from ORIGINAL unfiltered bim
outlier_snps <- bim$SNP[outliersbonf]

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


# verify outliers
message("Number of outliers: ", length(outliersbonf))
message("Number of valid outlier SNP IDs: ", sum(!is.na(outlier_snps)))

message("nrow of qqdf: ", nrow(qqdf_GFMxWM))
message("CHR range: ", min(CHR), " to ", max(CHR))
message("BP range: ", min(BP), " to ", max(BP))
message("P range: ", min(P), " to ", max(P))
message("CHR unique values: ")
print(sort(unique(CHR)))
message("Sample SNP IDs: ")
print(head(SNP, 10))

# build manhattan qqman
png(filename = "/scratch/las80898/pcadapt_output_3/GFMxWM_qqman_2.png", width = 2400, height = 1600, res = 600,)
manhattan(qqdf_GFMxWM, 
          main = "pcadapt SNP Outliers", 
          cex.axis = 0.8, cex.main = .8,
          annotatePval = 0.0000001, 
          suggestiveline = FALSE,
          annotateTop = FALSE, 
          xlab = "Chromosome number", 
          cex = 0.3, 
          highlight = outlier_snps,
          ylim = c(0, 320))   # 320 safely covers -log10(.Machine$double.xmin)
dev.off()


# plotting with ggplot
# Prepare the dataset
don <- qqdf_GFMxWM %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(qqdf_GFMxWM, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(logP = pmin(-log10(P), 320)) %>%        # move up here
  mutate(is_highlight=ifelse(SNP %in% outlier_snps, "yes", "no")) %>%
  mutate(is_annotate=ifelse(logP > 8, "yes", "no"))  # use logP, not -log10(P)

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
GFMxWM_ggplot_manhattan <- ggplot(don, aes(x=BPcum, y=logP)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey60"), length.out = 30)) +
  
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
ggsave("/scratch/las80898/pcadapt_output_3/GFMxWM_ggplot_manhattan_2.png", GFMxWM_ggplot_manhattan, width = 16, height = 12, dpi = 600)

###circular plot - need to amend chr.labels i think
setwd("/scratch/las80898/pcadapt_output_3")
CMplot(qqdf_GFMxWM, plot.type="c", r=1.6,
       outward=TRUE, cir.chr.h=.1, chr.den.col="orange",
       file.name="GFMxWM_circular_manhattan_2",
       file="jpg", dpi=600, chr.labels=seq(1, 30))
