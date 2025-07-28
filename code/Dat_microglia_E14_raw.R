setwd("/Users/david/Desktop/scRNA_seq_Neuron/Raw_data/Dat_microglia/Dat_microglia_raw")
rm(list = ls())
library(dplyr)
library(Seurat)
# We only consider E14 files
# Consider male first
# Add suffix for each cell due to replicated cells
E14_M_B8 <- read.table("./unzipped_raw/GSM3442013_E14_M_B8.dge.txt", header = TRUE)
rownames(E14_M_B8) <- E14_M_B8$GENE
E14_M_B8 <- E14_M_B8[, -1] # dim 14356*2594
colnames(E14_M_B8) <- paste0(colnames(E14_M_B8), "M_B8") 


E14_M_B7 <- read.table("./unzipped_raw/GSM3442012_E14_M_B7.dge.txt", header = TRUE)
rownames(E14_M_B7) <- E14_M_B7$GENE
E14_M_B7 <- E14_M_B7[, -1] # dim 15849*3180
colnames(E14_M_B7) <- paste0(colnames(E14_M_B7), "M_B7")


E14_M_B9 <- read.table("./unzipped_raw/GSM3442010_E14_M_B9.dge.txt", header = TRUE)
View(E14_M_B9[1:5, 1:5])
rownames(E14_M_B9) <- E14_M_B9$GENE
E14_M_B9 <- E14_M_B9[, -1] # dim 16000*3906
dim(E14_M_B9)
colnames(E14_M_B9) <- paste0(colnames(E14_M_B9), "M_B9")

E14_M_B11 <- read.table("./unzipped_raw/GSM3442007_E14_M_B11.dge.txt", header = TRUE)
View(E14_M_B11[1:5, 1:5])
rownames(E14_M_B11) <- E14_M_B11$GENE
E14_M_B11 <- E14_M_B11[, -1] # dim 14816*2242
dim(E14_M_B11)
colnames(E14_M_B11) <- paste0(colnames(E14_M_B11), "M_B11")


# Consider female
E14_F_B6 <- read.table("./unzipped_raw/GSM3442011_E14_F_B6.dge.txt", header = TRUE)
View(E14_F_B6[1:5, 1:5])
rownames(E14_F_B6) <- E14_F_B6$GENE
E14_F_B6 <- E14_F_B6[, -1] # dim 14664*2553
dim(E14_F_B6)
colnames(E14_F_B6) <- paste0(colnames(E14_F_B6), "F_B6")

E14_F_C1 <- read.table("./unzipped_raw/GSM3442009_E14_F_C1.dge.txt", header = TRUE)
View(E14_F_C1[1:5, 1:5])
rownames(E14_F_C1) <- E14_F_C1$GENE
E14_F_C1 <- E14_F_C1[, -1] # dim 14598*2657
dim(E14_F_C1)
colnames(E14_F_C1) <- paste0(colnames(E14_F_C1), "F_C1")

E14_F_B12 <- read.table("./unzipped_raw/GSM3442008_E14_F_B12.dge.txt", header = TRUE)
View(E14_F_B12[1:5, 1:5])
rownames(E14_F_B12) <- E14_F_B12$GENE
E14_F_B12 <- E14_F_B12[, -1] # dim 16442*4248
dim(E14_F_B12)
colnames(E14_F_B12) <- paste0(colnames(E14_F_B12), "F_B12")

E14_F_B10 <- read.table("./unzipped_raw/GSM3442006_E14_F_B10.dge.txt", header = TRUE)
View(E14_F_B10[1:5, 1:5])
rownames(E14_F_B10) <- E14_F_B10$GENE
E14_F_B10 <- E14_F_B10[, -1] # dim 15588*2732
dim(E14_F_B10)
colnames(E14_F_B10) <- paste0(colnames(E14_F_B10), "F_B10")

# Union genes
union_genes <- union(rownames(E14_M_B8), rownames(E14_M_B7)) %>%
  union(rownames(E14_M_B9)) %>%
  union(rownames(E14_M_B11)) %>%
  union(rownames(E14_F_B6)) %>%
  union(rownames(E14_F_C1)) %>%
  union(rownames(E14_F_B12)) %>%
  union(rownames(E14_F_B10))
length(union_genes) # 18541 genes

# Create new data using union genes
E14_M_B8_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_M_B8))
rownames(E14_M_B8_new) <- union_genes
colnames(E14_M_B8_new) <- colnames(E14_M_B8)
E14_M_B8_new[rownames(E14_M_B8), ] <- as.matrix(E14_M_B8)


E14_M_B7_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_M_B7))
rownames(E14_M_B7_new) <- union_genes
colnames(E14_M_B7_new) <- colnames(E14_M_B7)
E14_M_B7_new[rownames(E14_M_B7), ] <- as.matrix(E14_M_B7)


E14_M_B9_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_M_B9))
rownames(E14_M_B9_new) <- union_genes
colnames(E14_M_B9_new) <- colnames(E14_M_B9)
E14_M_B9_new[rownames(E14_M_B9), ] <- as.matrix(E14_M_B9)

E14_M_B11_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_M_B11))
rownames(E14_M_B11_new) <- union_genes
colnames(E14_M_B11_new) <- colnames(E14_M_B11)
E14_M_B11_new[rownames(E14_M_B11), ] <- as.matrix(E14_M_B11)



E14_F_B6_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_F_B6))
rownames(E14_F_B6_new) <- union_genes
colnames(E14_F_B6_new) <- colnames(E14_F_B6)
E14_F_B6_new[rownames(E14_F_B6), ] <- as.matrix(E14_F_B6)

E14_F_C1_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_F_C1))
rownames(E14_F_C1_new) <- union_genes
colnames(E14_F_C1_new) <- colnames(E14_F_C1)
E14_F_C1_new[rownames(E14_F_C1), ] <- as.matrix(E14_F_C1)

E14_F_B12_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_F_B12))
rownames(E14_F_B12_new) <- union_genes
colnames(E14_F_B12_new) <- colnames(E14_F_B12)
E14_F_B12_new[rownames(E14_F_B12), ] <- as.matrix(E14_F_B12)

E14_F_B10_new <- matrix(0, nrow = length(union_genes), ncol = ncol(E14_F_B10))
rownames(E14_F_B10_new) <- union_genes
colnames(E14_F_B10_new) <- colnames(E14_F_B10)
E14_F_B10_new[rownames(E14_F_B10), ] <- as.matrix(E14_F_B10)

E14 <- cbind(E14_M_B8_new, E14_M_B7_new)%>%
  cbind(E14_M_B9_new) %>%
  cbind(E14_M_B11_new) %>%
  cbind(E14_F_B6_new) %>%
  cbind(E14_F_C1_new) %>%
  cbind(E14_F_B12_new) %>%
  cbind(E14_F_B10_new)
dim(E14)

# Create meta for 8 batches
batch_idx <- c("M_B8", "M_B7", "M_B9", "M_B11",
               "F_B6", "F_C1", "F_B12", "F_B10")
rep_num <- c(ncol(E14_M_B8_new), ncol(E14_M_B7_new), ncol(E14_M_B9_new), ncol(E14_M_B11_new),
             ncol(E14_F_B6_new), ncol(E14_F_C1_new), ncol(E14_F_B12_new), ncol(E14_F_B10_new))
meta_data <- data.frame(cell_id = colnames(E14),
                        cell_batch = rep(batch_idx, rep_num))

saveRDS(E14, file = "./Dat_E14/E14_count.rds")
saveRDS(meta_data, file = "./Dat_E14/E14_meta.rds")
