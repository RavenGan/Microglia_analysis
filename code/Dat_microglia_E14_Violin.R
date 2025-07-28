setwd("/Users/david/Desktop/scRNA_seq_Neuron/Raw_data/Dat_microglia/Dat_microglia_raw")
rm(list = ls())
library(dplyr)
library(Seurat)

E14 <- readRDS("./Dat_E14/E14_SeuratObject.rds")

gene_ls <- c("Tmem119", "P2ry12", "Mrc1", "F13a1",
             "Ccr", "Ccr2", "Cldn5", "Vtn", "Pecam1",
             "Neurod6", "Nfib", "Elavl3", "Folr2",
             "F13a1", "Cd163", "Ccl2", "Ms4a7")

x11()
pdf("./Dat_E14/E14_violin.pdf", height = 20, width = 20)
VlnPlot(E14, features = gene_ls, slot = "counts", log = TRUE)
dev.off()


gene_ls2 <- c("Lamp1", "Cd63", "Casp8", "Ccr2")

x11()
pdf("./Dat_E14/Ana_Mar_09/E14_violin_new.pdf", height = 8, width = 20)
VlnPlot(E14, features = gene_ls2, slot = "counts", log = TRUE, ncol = 4)
dev.off()

gene_ls3 <- c("Sox17")

x11()
pdf("./Dat_E14/Ana_Mar_09/E14_violin_sox17.pdf", height = 5, width = 5)
VlnPlot(E14, features = gene_ls3, slot = "counts", log = TRUE)
dev.off()