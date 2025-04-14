# 3seurat_object.R
# script voor het maken van het seurat object

#Seurat object maken met raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

##An object of class Seurat 
##13714 features across 2700 samples within 1 assay 
##Active assay: RNA (13714 features, 0 variable features)
##2 layers present: counts, data

#examine eerste 3 genen in de eerste 30 cellen 
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#0 values wegfilteren en een compacter versie object opslaan 
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size



