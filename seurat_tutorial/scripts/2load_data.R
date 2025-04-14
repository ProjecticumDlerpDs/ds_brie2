# 2load_data.R
# script om data in te laden 

#downloaden van de data 
#de data is gedownload van 10x Genomics. Het bestand is een `tar.gz`bestand dat UMI count matrices bevat. 
#bestand lokaal gedownload en geupload op locatie Project-BRIE2/seurat_tutorial/raw_data

#bestand uitpakken in terminal 
(project_brie2) 1843893@hu-p-dataanalysis01:~/Project-BRIE2/seurat_tutorial/raw_data$ tar -xzvf pbmc3k_filtered_gene_bc_matrices.tar.gz


#data inlezen
pbmc.data <- Read10X(data.dir = "~/Project-BRIE2/seurat_tutorial/raw_data/filtered_gene_bc_matrices/hg19/")

