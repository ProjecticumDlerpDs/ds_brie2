#!/bin/bash

# 1setup_conda_envs.sh
# script om conda omgevingen aan te maken en BRIE2 en briekit te installeren

# Stap 1: nieuwe conda omgeving aanmaken
# voor BRIE2
conda create -n TFProb python=3.11 -y

#voor briekit
conda create -n briekit python=2.7 numpy=1.13.0 -y

# Stap 2: installeren via pip
# packages installeren voor BRIE2 in omgeving TFProb
conda run -n TFProb pip install -U git+https://github.com/huangyh09/brie  

# packages installeren voor briekit in omgeving briekit
conda run -n briekit pip install briekit

# Stap 4: Instalatie testen 
# voor BRIE2 in omgeving TFProb
conda run -n TFProb brie-quant 


