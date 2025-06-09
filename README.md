# Data science project BRIE2
Een data science project voor de analyse van alternatieve splicing in VASA-sequencing van muisembryo's met BRIE2. 

## Beschrijving
Het doel van dit project is om een pipeline te schrijven voor het analyseren van VASA-sequencing met BRIE2 op alternative splicing. 
In deze Githup repository zullen alle scripts en data worden opgeslagen. Niet alle data kan op githup worden geplaats vanwege de grote van de bestanden. Deze bestanden zijn terug te vinden in het gitgnore bestand, en hieronder zal worden beschreven waar de data op de HU Rstudio server is opgeslagen.

## Reproduceerbaarheid 
Alle benodigde informatie over de scripts en het uitvoeren hiervan is te vinden in de map scripts_info. Lees deze eerst door voor het uitvoeren van scritps.

## Workflow 
Deelvraag 1: Hoe werkt Seurat en wat leer ik via de tutorial?
- Seurat installeren en tutorial volgen 
- Script: 1setup_conda_env.sh en Tutorial_seurat.pdf

Deelvraag 2: Hoe kan de ruwe data worden gefilterd?
- Filtering testen 
- Script: Seurat_filtering.pdf

Deelvraag 3: Hoe kan eigen data worden geclusterd?
- Kwaliteit van data bepalen
- Bepalen welke PC worden meegenomen 
- Script: seurat_pipeline.pdf

Deelvraag 4: Welke clusters gaan we analyseren met BRIE2?
- Selecteer relevante clusters voor splicing-analyse.
- Output: seurat_pipeline.pfd

Deelvraag 5: Hoe werkt BRIE2 en wat leer ik via de tutorial?
- Installeren BRIE2 
- Doorloop de BRIE2 tutorial met voorbeelddata.
- Script: 1setup_conda_envs.sh en brie2_tutorial.pdf

Deelvraag 6: Hoe wordt data geschikt voor input van BRIE2?
- data overzetten van seurat naar brie2
- Script: seurat_pipeline/bewerkte_data/designmatrix_5vs9.tsv en 2download_fastq.sh

Deelvraag 7: Hoe kan de data geanalyseerd worden in BRIE2

Deelvraag 8: Kan de output van BRIE2 count gekwalificeerd worden?

Deelvraag 9: Hoe kan de data van BRIE2 gevisualiseerd worden?

## Data
Het orginele dataset is opgeslagen op de Rstudio Server "/home/data/projecticum/splicing/data"

## Stap 1 Seurat 
Eerst is er een seurat tutorial uitgevoerd. (https://github.com/satijalab/seurat/blob/HEAD/vignettes/pbmc3k_tutorial.Rmd). 

Daarna vind inspectie en clustering van het data set E8.5 plaats met de package seurat van [@satijalab/seurat] (https://github.com/satijalab/seurat)

De filtering met seurat is getest en de resultaten hiervan zijn beschreven in een Rmarkdown (https://github.com/ProjecticumDlerpDs/ds_brie2/tree/main/seurat_pipeline/scripts)

Hier is een analyse pipeline voor geschreven:[seurat_pipeline](https://github.com/ProjecticumDlerpDs/ds_brie2/tree/main/seurat_pipeline) 

## Stap 2 BRIE2
voor verdere analyse met is gebruik gemaakt van package BRIE2 [@huangyh09/brie](https://github.com/huangyh09/brie)

De tutorial van brie2 is uitgevoerd. (https://github.com/ProjecticumDlerpDs/ds_brie2/tree/main/brie2_tutorial/scripts)

Nu wordt gewerkt aan het overzetten van de data vanuit Seurat naar een format voor BRIE2

## Orginele data 
- Artikel: High-throughput total RNA sequencing in single cells using VASA-seq (https://doi.org/10.1101/2021.09.15.460240)
- Github:  [@hemberg-lab/VASAseq_2022](https://github.com/hemberg-lab/VASAseq_2022)

## Versies
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

### Structuur repository 
Project-BRIE2/
├── seurat_tutorial/ # R-project met Seurat tutorial
├── seurat_pipeline/ # R-project met de volledige Seurat pipeline
├── brie2_tutorial/ # Python tutorial voor BRIE2
├── brie2_pipeline/ # Volledige BRIE2 analysepipeline
├── scripts_info/ # Instructies voor runnen van de scripts

### Benodigheden
- Conda (voor Python-omgevingen: TFProb, briekit en voor R-omgeving project_brie2)
- R en RStudio (voor Seurat-analyses)
- Git (voor versiebeheer)

### Environments  
|  Map               | Taal    | Environment        |  Versies                                
| `seurat_tutorial/` | R       | `project_brie2`    | seurat-4.4.0                            
| `seurat_pipeline/` | R       | `project_brie2`    | seurat-4.4.0                            
| `brie2_tutorial/`  | Python  | `TFProb`           | python-3.11.11                           
| `brie2_pipeline/`  | Python  | `briekit`, `TFProb`| briekit = python-2.7 en TFProb = python-3.11.11


## Auteur
Anne Brussaard [@annebrussaard](https://github.com/annebrussaard)