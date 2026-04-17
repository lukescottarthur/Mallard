# this is a script made to be submitted in bash, intended to analyze
# whole genome data and output files along the way

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)
library(CMplot)

# read data
GFMxWM <- read.pcadapt("GFMxWM.bed", type = "bed")

~~~~use the following code to save what I want.
# save ggplot using ggsave
# make a plot first with a name, then save it
ggsave(
    "/scratch/las80898/pcadapt_output/*.png",
    combo_plot
)

# or base r plot using png, pdf, etc.
png("/scratch/las80898/pcadapt_output/*.png")
heatmap.2(mtrx, scale = "row", Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none")
dev.off()

# Example: Saving a base plot to a PNG file
png(filename = "output_plot.png")
plot(1:10, 1:10, col = "blue", pch = 18, main = "Base Plot")
dev.off()   
~~~~

# initial analysis
x1 <- pcadapt(GFMxWM, K = 2)
summary(x1)
plot(x1 , option = "manhattan")
plot(x1, option = "qqplot")
hist(x1$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x1, option = "stat.distribution")

# outlier adjustment
padjbonf <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.0001
outliersbonf <- which(padjbonf < alpha)
length(outliersbonf)
outliersbonf

# LD
for (i in 1:2)
  plot(x1$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

#association between pc and outliers
snp_pc <- get.pc(x1, outliersbonf)
snp_pc

# plotting with qqman

#make dataframe with values from pcadapt
SNP <- x1$pass
CHR <- rep(17, times = 9093) 
BP <- x1$pass
P <- x1$pvalues[!is.na(x1$pvalues)]

qqdf_GFMxWM <- data.frame(SNP, CHR, BP, P)

manhattan(qqdf_GFMxWM, main = " pcadapt SNP Outliers", cex.axis = 0.8, cex.main = .8, 
          annotatePval = 0.0000001, suggestiveline = F , annotateTop = FALSE, 
          xlab = "Chromosome number", cex = 0.3, highlight = outliersbonf)

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
  mutate( is_highlight=ifelse(SNP %in% outliersbonf, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(P)>8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black"))) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
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

# saving plots
ggsave(".....png", chr17gg, width = 8, height = 6, dpi = 600)
ggsave("......png", plot = my_plot, width = 8, height = 6, dpi = 300)   


###circular plot
CMplot(qqdf_GFMxWM, plot.type="c", r=1.6,
       outward=TRUE, cir.chr.h=.1 ,chr.den.col="orange", file="jpg",
       dpi=300, chr.labels=seq(1,22))
