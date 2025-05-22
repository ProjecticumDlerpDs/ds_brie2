# script voor clusterlabels uit seurat halen 
# werken in Project-BRIE2/seurat_pipeline/bewerkte_data
# conda envirment project_brie2 openen 
# R opstarten in terminal 

# openen van packages 
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(Seurat)
library(Matrix)

# seurat object inlezen 
seurat_obj <- readRDS("umap_seurat_embryoE85.rds")

# meta data bekijken 
head(seurat_obj@meta.data)

# design matrix maken 
cluster_labels <- seurat_obj@meta.data[, c("orig.ident", "seurat_clusters")]
write.table(cluster_labels, file = "design_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

# kolom namen veranderen 
meta <- seurat_obj@meta.data
meta$cell_id <- rownames(meta)
design_df <- meta[, c("cell_id", "seurat_clusters")]
colnames(design_df) <- c("cell_id", "cluster")

# design matrix maken voor cluster 5 en 9 
write.table(design_df, "design_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
subset_cells <- subset(seurat_obj, idents = c("5", "9"))
design_df <- data.frame(
  cell_id = rownames(subset_cells@meta.data),
  cluster = as.character(subset_cells@meta.data$seurat_clusters)
)
write.table(design_df, "design_matrix_5vs9.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# design_matrix.tsv een design_matrix_5vs9.tsv opgeslagen in seurat_pipeline/bewerkte_data