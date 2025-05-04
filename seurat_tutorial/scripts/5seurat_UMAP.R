#5seurat_UMAP.R
#script om de effecten van instellingen te zien op UMAP

# laden van libraries
library(Seurat)
library(dplyr)
library(patchwork)

# Verschillende QC-varianten maken
pbmc_strict <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc_mild <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
pbmc_nofilter <- pbmc  # geen filter toegepast

# functie seurat 
process_seurat <- function(qc_filter_obj, label) {
  qc_filter_obj <- NormalizeData(qc_filter_obj)
  qc_filter_obj <- FindVariableFeatures(qc_filter_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(qc_filter_obj)
  qc_filter_obj <- ScaleData(qc_filter_obj, features = all.genes)
  qc_filter_obj <- RunPCA(qc_filter_obj, features = VariableFeatures(qc_filter_obj))
  qc_filter_obj <- FindNeighbors(qc_filter_obj, dims = 1:10)
  qc_filter_obj <- FindClusters(qc_filter_obj, resolution = 0.5)
  qc_filter_obj <- RunUMAP(qc_filter_obj, dims = 1:10)
  qc_filter_obj$filter_label <- label
  return(qc_filter_obj)
}

#QC uitvoeren voor varianten 
pbmc_strict <- process_seurat(pbmc_strict, "Strikt")
pbmc_mild <- process_seurat(pbmc_mild, "Mild")
pbmc_nofilter <- process_seurat(pbmc_nofilter, "Geen filter")

#individule UMAP plots maken 
umap_strict <- DimPlot(pbmc_strict, reduction = "umap")
umap_mild <- DimPlot(pbmc_mild, reduction = "umap")
umap_nofilter <- DimPlot(pbmc_nofilter, reduction = "umap")

#UMAP titel geven 
combined_umap <- (umap_strict + umap_mild + umap_nofilter) +
  plot_annotation(title = "Vergelijking van filtering effect op UMAP", 
                  subtitle = "Strikt, Mild, Geen filtering")

#UMAP Strikt, Mild, Geen filtering opslaan
png("umap_qc_variaties.png", width = 1800, height = 600)
print(combined_umap)
dev.off()
