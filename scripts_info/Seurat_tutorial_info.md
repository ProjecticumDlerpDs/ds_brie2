# Seurat Tutorial - Script info 

## Benodigde omgeving 
- Conda environment: `project_brie2`
- R versie: 4.2.3
- Belangrijkste packages:
  - Seurat v4.4.0
  
## Bestanden 
in Project-BRIE2/seurat_tutorial
-`scripts/1setup_conda_env.sh`
-`scripts/2load_data.R`
-`scripts/3seurat_object.R`
-`scripts/4seurat_tutorial.R`
-`scripts/5seurat_UMAP.R`
-`Tutorial_seurat.Rmd`

## Scrips
Script `Tutorial_seurat.Rmd` kan worden gerund om alle onderstaande stappen uit te voeren. 

Of
1. Voer `scripts/1setup_conda_env.sh`uit om de conda omgeving te maken
2. Gebruik `scripts/2load_data.R`voor het downloaden en laden van de data
3. Ga verder met `scripts/3seurat_object.R`voor het maken van het seurat object
4. Volg `scripts/4seurat_tutorial.R`voor het verder uitvoeren van de tutorial 
5. Sluit af met `scripts/5seurat_UMAP.R`voor het visualiseren van de data

## Runnen 
1. conda activate project_brie2
2. activeer R in de terminal 
3. Voor `Tutorial_seurat.Rmd` voer uit in de terminal rmarkdown::render("seurat_tutorial/tutorial_analysis.Rmd")
4. voor losse scripts voer de commands uit in de terminal

### Output 
- `Tutorial_seurat.pdf`