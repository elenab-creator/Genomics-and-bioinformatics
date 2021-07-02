# RStudio version 1.4.1103
# R R version 4.0.4

# Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pathview")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

install.packages("ggplot2")
install.packages("pheatmap")
install.packages("readxl")
install.packages("DESeq2")
install.packages("gridExtra")
install.packages("DOSE")
install.packages("AnnotationHub")
install.packages("org.Mm.eg.db")

library("ggplot2")
library("pheatmap")
library("readxl")
library("DESeq2") 
library(gridExtra)
library(DOSE)
library(pathview)
library(clusterProfiler) 
library(AnnotationHub) 
library (org.Mm.eg.db)
library(data.table)

