setwd("/Users/david/Desktop/scRNA_seq_Neuron/Raw_data/Dat_microglia/Dat_microglia_raw")
rm(list = ls())
library(dplyr)
library(Seurat)

E14_raw <- readRDS("./Dat_E14/E14_count.rds")
meta_data <- readRDS("./Dat_E14/E14_meta.rds")

# Initialize the Seurat object with the raw (non-normalized data).
E14 <- CreateSeuratObject(counts = E14_raw, project = "E14Cells", min.cells = 3, min.features = 200)
E14

E14[["percent.mt"]] <- PercentageFeatureSet(E14, pattern = "^MT-")

# Visualize QC metrics as a violin plot
x11()
VlnPlot(E14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalizing the data
E14 <- NormalizeData(E14, normalization.method = "LogNormalize", scale.factor = 10000)

# highly variable features
E14 <- FindVariableFeatures(E14, selection.method = "vst", nfeatures = 2000)

# scaling
all.genes <- rownames(E14)
E14 <- ScaleData(E14, features = all.genes)

# Linear dimension reduction
E14 <- RunPCA(E14, features = VariableFeatures(object = E14))

# Decide the number of PCs
x11()
ElbowPlot(E14) # keep all 20 PCs

# cluster cells
E14 <- FindNeighbors(E14, dims = 1:20)
E14 <- FindClusters(E14, resolution = 0.5)

# UMAP
E14 <- RunUMAP(E14, dims = 1:20)

x11()
pdf("./Dat_E14/E14_umap.pdf", height = 8, width = 8)
DimPlot(E14, reduction = "umap")
dev.off()

# Change ident of E14
# Idents(E14) <- meta_data$cell_batch
# x11()
# pdf("./Dat_E14/E14_umap_batch_effects.pdf", height = 8, width = 8)
# DimPlot(E14, reduction = "umap")
# dev.off()

saveRDS(E14, file = "./Dat_E14/E14_SeuratObject.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
E14.markers <- FindAllMarkers(E14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Marker_gene_lst <- E14.markers %>%
  group_by(cluster) %>%
  slice_max(n = Inf, order_by = avg_log2FC)

n_cluster <- length(unique(Marker_gene_lst$cluster))

# Save marker genes for each cluster
for (i in 1:n_cluster) {
  chosen_mg_lst <- Marker_gene_lst[which(Marker_gene_lst$cluster == (i-1)), ]
  save_path <- paste0("./Dat_E14/Marker_gene_ls/Cluster_", (i-1), "_mg_lst.csv")
  write.csv(chosen_mg_lst, save_path)
}
