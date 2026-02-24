#install.packages("pcadapt")
#install.packages("CMplot")
#install_github("jdstorey/qvalue")
#install.packages("qqman")
library(ggrepel)
library(pcadapt)
library(ggplot2)
library(dplyr)
library(CMplot)
library(devtools)
library(qvalue)
library(qqman)
#library(grid)
#library(gridGraphics)



test17 <- read.pcadapt("chr17.bed", type = "bed")
test17
class(test17)


x <- pcadapt(input = test17, K = 5)
class(x)
summary(x)

plot(x, option = "screeplot")

#K=2

#computing test statistic
x1 <- pcadapt(test17, K = 2)
summary(x1)
plot(x1 , option = "manhattan")
plot(x1, option = "qqplot")
hist(x1$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x1, option = "stat.distribution")


#choosing method for cutoff of outlier detection

#least conservative
qval <- qvalue(x1$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)

#moderately conservative
padjBH <- p.adjust(x1$pvalues,method="BH")
alpha <- 0.1
outliersBH <- which(padjBH < alpha)
length(outliersBH)
outliersBH

#most conservative
padjbonf <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.0001
outliersbonf <- which(padjbonf < alpha)
length(outliersbonf)
outliersbonf

#LD thinning

#Recommendation from authors:display the loadings (contributions of each SNP to the PC) and to evaluate 
#if the loadings are clustered in a single or several genomic regions
for (i in 1:2)
  plot(x1$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))


#doesn't look like LD is affecting PCs

#association between pc and outliers
snp_pc <- get.pc(x1, outliersbonf)
snp_pc


#converting pcadapt object contents into qqman format
#manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
#need: chromosome number (1), bp location number (9093), snp identifier, p value, snps of interest
summary(x1)

SNP <- x1$pass
CHR <- rep(17, times = 9093)   
BP <- x1$pass
P <- x1$pvalues[!is.na(x1$pvalues)]


qqdf17 <- data.frame(SNP, CHR, BP, P)

####this one
qq17plot <- manhattan(qqdf17, main = " Pcadapt SNP Outliers", cex.axis = 0.8, cex.main = .8, annotatePval = 0.0000001, suggestiveline = F , annotateTop = FALSE, xlab = "Chromosome 17", cex = 0.3, highlight = outliersbonf)

######convert qq data frame to ggplot

don <- qqdf17 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(qqdf17, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=1.3) +
  scale_color_manual(values = rep(c("black") )) +
  
  #horizontal line
  geom_hline(yintercept = 8, col="red", alpha = .7) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

##highlight snps in ggplot
# List of SNPs to highlight are in the outliersbonf object
# We will use ggrepel for the annotation


# Prepare the dataset
don <- qqdf17 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(qqdf17, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% outliersbonf, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(P)>8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
chr17gg <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=1.3) +
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
       x = "Chromosome 17 SNPs") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#ggsave("chr17gg.png", chr17gg, width = 8, height = 6, dpi = 600)
#ggsave("my_plot.png", plot = my_plot, width = 8, height = 6, dpi = 300)   


###circular plot
CMplot(qqdf17, plot.type="c", r=1.6,
       outward=TRUE, cir.chr.h=.1 ,chr.den.col="orange", file="jpg",
       dpi=300, chr.labels=seq(1,22))

