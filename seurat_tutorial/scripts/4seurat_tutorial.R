# 4seurat_tutorial.R
# Script om Seurat te laden in een actieve conda omgeving 
# en Seurat tutorial te volgen 

# laden van libraries
library(Seurat)
library(dplyr)
library(patchwork)


#QC en cellen selecteren voor verdere analyse 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#QC metrics eerste 5 cellen 
head(pbmc@meta.data, 5)

#Visualisatie QC metrics in violin plot 
pdf("combined_vlnplot.pdf", width = 12, height = 4)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


#Visualisatie FeatureScatter 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

png("scatterplots.png", width = 1000, height = 500)
print(combined)
dev.off()

#QC filter >200 gene expression, <2500 cellen met veel genen, <5 mitochondriale expression
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalisering data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#features selection 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#identificering 10 meest variable genen 
top10 <- head(VariableFeatures(pbmc), 10)

#variable features plotten met en zonder labels 
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png("variable_features_labeled.png", width = 800, height = 600)
print(plot1 + plot2)
dev.off()

#scaling the data 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#linear dimensionale reduction 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#visualisatie PCA results 
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

png("pca_loadings.png", width = 1000, height = 500)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

png("pca_dimplot.png", width = 1000, height = 500)
DimPlot(pbmc, reduction = "pca") + NoLegend()
dev.off()

pdf("dimheatmap.pdf", width = 10, height = 8)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#dimentie dataset bepalen
png("elbowplot.png")
ElbowPlot(pbmc)
dev.off()

#cellen clusteren
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#eerste 5 IDs clusters 
head(Idents(pbmc), 5)

#UMAP maken 
pbmc <- RunUMAP(pbmc, dims = 1:10)

png("umap.png")
DimPlot(pbmc, reduction = "umap")
dev.off()

#object opslaan
saveRDS(pbmc, file = "~/Project-BRIE2/seurat_tutorial/bewerkte_data/pbmc_tutorial.rds")

