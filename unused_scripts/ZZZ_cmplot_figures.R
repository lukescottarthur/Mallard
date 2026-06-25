# GFMxWM - sliding window figure

# load libraries
library(pcadapt)
library(ggplot2)
library(qqman)
library(ggrepel)
library(dplyr)
library(devtools)
library(CMplot)

# load data
WMAxGFM <- read.csv("/scratch/las80898/mallard_cmplot/cmplot_WMAxGFM.csv")
GFMxKC <- read.csv("/scratch/las80898/mallard_cmplot/cmplot_GFMxKC.csv")
WMAxKC <- read.csv("/scratch/las80898/mallard_cmplot/cmplot_WMAxKC.csv")

# make plots
CMplot(WMAxGFM, plot.type="m", multracks=TRUE, file="png", file.name=cmplot_WMAxGFM_1,dpi=500, file.output=TRUE, verbose=FALSE, ylab=c("Fst","-log10(p) pcadapt","SS-X12"))