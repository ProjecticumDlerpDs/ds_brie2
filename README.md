# Data science project BRIE2
Een data science project voor de analyse van alternatieve splicing in VASA-sequencing van muisembryo's met BRIE2. 

#Beschrijving
Het doel van dit project is om een pipeline te schrijven voor het analyseren van VASA-sequencing met BRIE2 op alternative splicing. 
In deze Githup repository zullen alle scripts en data worden opgeslagen. Niet alle data kan op githup worden geplaats vanwege de grote van de bestanden. Deze bestanden zijn terug te vinden in het gitgnore bestand, en hieronder zal worden beschreven waar de data op de HU Rstudio server is opgeslagen. 

#Data
Het orginele dataset is opgeslagen op de Rstudio Server "/home/data/projecticum/splicing/data"

#Stap 1 Seurat 
Eerst vind inspectie en clustering van het data set E8.5 plaats met de package seurat van [@satijalab/seurat] (https://github.com/satijalab/seurat)
Hier is een pipeline voor geschreven:[seurat_pipeline](https://github.com/ProjecticumDlerpDs/ds_brie2/tree/main/seurat_pipeline) 

#Stap 2 BRIE2
voor verdere analyse met is gebruik gemaakt van package BRIE2 [@huangyh09/brie](https://github.com/huangyh09/brie)

- stap 2 is gestart 

#Orginele data 
- Artikel: High-throughput total RNA sequencing in single cells using VASA-seq (https://doi.org/10.1101/2021.09.15.460240)
- Github:  [@hemberg-lab/VASAseq_2022](https://github.com/hemberg-lab/VASAseq_2022)

#Versies

R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

# Auteur
Anne Brussaard [@annebrussaard](https://github.com/annebrussaard)