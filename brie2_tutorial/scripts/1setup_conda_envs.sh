#!/bin/bash

# 1setup_conda_envs.sh
# script om conda omgevingen aan te maken en BRIE2 en briekit te installeren

# Stap 1: nieuwe conda omgeving aanmaken
# voor BRIE2
conda create -n TFProb python=3.11

#voor briekit
conda create -n briekit python=2.7 numpy=1.13.0

# Stap 2: omgeving activeren
# commanod handmatig runnen via terminal:
# voor BRIE2
conda activate TFProb

#v voor BRIEkit
conda activate briekit

# Stap 3: BRIE2 installeren via pip
# voor BRIE2
pip install -U git+https://github.com/huangyh09/brie  

# voor briekit
pip install briekit

# Stap 4: Instalatie testen 
# voor BRIE2
brie-quant 


