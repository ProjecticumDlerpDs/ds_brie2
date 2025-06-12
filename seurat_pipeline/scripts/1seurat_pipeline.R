#1seurat_pipeline.R
#script voor het uitvoeren van de seurat filtering 

#packages laden 
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(Seurat)
library(Matrix)
library(patchwork)
library(here)


#data laden 
features <- read.csv(here("seurat_pipeline", "raw_data", "e85_feature_metadata.csv.gz"))
samples <- read.csv(here("seurat_pipeline", "raw_data", "e85_sample_metadata.csv"))


#MTX file laden 
counts <- ReadMtx(here("seurat_pipeline", "raw_data", "e85_count_matrix.mtx.gz"),
                  here("seurat_pipeline", "raw_data", "e85_sample_metadata.csv"),
                  here("seurat_pipeline", "raw_data", "e85_feature_metadata.csv.gz"),
                  feature.sep = ",",
                  cell.sep = ",",
                  cell.column = 3,
                  feature.column = 1,
                  skip.cell = 1,
                  skip.feature = 1
)

#seurat object maken 
seurat <- CreateSeuratObject(counts = counts, 
                             project = "mouse_embryo", 
                             min.cells = 3, 
                             min.features = 200)
seurat

#percent.MT wordt toegevoegd
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

#QC aantal features, countes en percentage MT
png("plot.png", width = 800, height = 600)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()  

#shatterplot van QC
png("scatter_plots.png", width = 1200, height = 800)
cnt_ftr <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(cnt_ftr)  

cnt_mt <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
print(cnt_mt)  # Print de plot naar het bestand

dev.off()  

#subset van de data
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normaliseren van de data
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

#10 meest variabele genen uitzoeken 
top10 <- head(VariableFeatures(seurat), 10)

#data plotten zonder label 
png("variable_feature_plot.png", width = 1200, height = 800)
plot1 <- VariableFeaturePlot(seurat)
print(plot1)  

vrbl <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
print(vrbl)  

dev.off()  

#iets nog 
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

#PCA analyse
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

#plotten 
png("elbowplot.png", width = 1200, height = 800)
ElbowPlot(seurat)
dev.off()  

png("dimplot.png", width = 1200, height = 800)
dimplot <- DimPlot(seurat, reduction = "pca") + NoLegend()
print(dimplot)  
dev.off()

#data clusteren
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)

#UMAP maken 
png("umap_plot.png", width = 1200, height = 800)
seurat <- RunUMAP(seurat, dims = 1:10)

umap_plot <- DimPlot(seurat, reduction = "umap")
print(umap_plot)  

dev.off()  

#seurat_object opslaan
saveRDS(seurat, file = "Project-BRIE2/seurat_pipeline/bewerkte_data/umap_seurat_embryoE85.rds")

#markers vinden en allen de positieve rapporteren
markers <- FindAllMarkers(seurat, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#heatmap maken van markers
png("heatmap_top10_features.png", width = 1200, height = 800)

heatmap_plot <- DoHeatmap(seurat, features = top10$gene) + NoLegend()
print(heatmap_plot)  

dev.off()  



