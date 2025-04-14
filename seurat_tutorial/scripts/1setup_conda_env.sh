#!/bin/bash

# 1setup_conda_env.sh
# script om conda omgeving aan te maken en Seurat te installeren

# Stap 1: nieuwe conda omgeving aanmaken
conda create -n project_brie2
   j
# Stap 2: omgeving activeren
# commanod handmatig runnen via terminal:
conda activate project_brie2

# Stap 3: versie 4 van Seurat installeren via conda
conda install -c bioconda r-seurat 