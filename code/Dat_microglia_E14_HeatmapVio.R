setwd("./Raw_data/Dat_microglia/Dat_microglia_raw")
rm(list = ls())
library(dplyr)
library(Seurat)
library(ggplot2)

E14 <- readRDS("./Dat_E14/processed_dat/E14_SeuratObject.rds")

# Heatmap plots
gene_ls1 <- read.csv("./Dat_E14/processed_dat/heatmap1.csv")
# For heatmap1 gene list, use gene name F11r, instead of F11r.1 or F11r.2
#     use gene name Spi1, instead of Spi1a/b
genes1 <- gene_ls1$Gene_names
sum(genes1 %in% rownames(E14)) # 22/24 are in the data: Gpr34b and Selenop are not in the gene names

gene_ls2 <- read.csv("./Dat_E14/processed_dat/heatmap2.csv")
genes2 <- gene_ls2$Gene_names
sum(genes2 %in% rownames(E14)) # 94/94 are in the data

pdf("./Dat_E14/res/Feb_0219_2024/heatmap1.pdf", width = 10, height = 7)
DoHeatmap(E14, features = genes1, angle = 60, raster = FALSE, size = 5) + 
  theme(axis.text.y = element_text(size = 10))
dev.off()

pdf("./Dat_E14/res/Feb_0219_2024/heatmap2.pdf", width = 10, height = 20)
DoHeatmap(E14, features = genes2, angle = 60, raster = FALSE, size = 5) + 
  theme(axis.text.y = element_text(size = 10))
dev.off()

# Violin plots
gene_ls <- c("Cmtm7", "Cndp2", "F11r", "Gas6", "Irf8", "Mpp1") # Replace F11r.1 with F11r
all(gene_ls %in% rownames(E14)) # TRUE
pdf("./Dat_E14/res/Feb_0219_2024/violin.pdf", height = 8, width = 12)
VlnPlot(E14, features = gene_ls, slot = "counts", log = TRUE)
dev.off()

# UMAP plots
gene_ls <- c("Cmtm7", "Cndp2", "F11r", "Gas6", "Irf8", "Mpp1") # Replace F11r.1 with F11r
pdf("./Dat_E14/res/Feb_0219_2024/FeaturePlot.pdf", height = 12, width = 10)
FeaturePlot(E14, features = gene_ls)
dev.off()