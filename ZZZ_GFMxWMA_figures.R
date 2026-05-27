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
Fst <- read.csv("/scratch/las80898/pcadapt_output_5/Fst_windows.csv")
neg_log_p <- read.csv("/scratch/las80898/pcadapt_output_5/GFMxWM_revised_windows.csv")

# make plots
aphgpoashgso