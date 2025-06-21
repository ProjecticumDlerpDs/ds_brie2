#!/bin/bash

# SRA-accession nummers met bijbehorende sample namen
# Format: SampleName SRRnummer

samples=(
  "E8.5-19_i36 SRR14783084"
  "E8.5-1_i7 SRR14783085"
  "E8.5-20_i37 SRR14783086"
  "E8.5-21_i38 SRR14783087"
  "E8.5-22_i39 SRR14783088"
  "E8.5-23_i40 SRR14783089"
  "E8.5-4_i8 SRR14783092"
  "E8.5-5_i18 SRR14783093"
  "E8.5-8_i21 SRR14783096"
  "E8.5-9_i33 SRR14783097"
 )

mkdir -p ../raw_data/fastq

for sample in "${samples[@]}"
do
name=$(echo $sample | cut -d' ' -f1)
srr=$(echo $sample | cut -d' ' -f2)

echo "Processing sample: $name with SRR: $srr"

# Download de .sra file
prefetch $srr

# Zet om naar fastq en sla op in fastq/ map, naam prefix met sample naam
fastq-dump $srr --split-files --gzip -O ../raw_data/fastq/  2>> ../raw_data/fastq_errors.log
  
# Hernoem bestanden naar herkenbare namen
if [[ -f ../raw_data/fastq/${srr}_1.fastq.gz && -f ../raw_data/fastq/${srr}_2.fastq.gz ]]; then
    mv ../raw_data/fastq/${srr}_1.fastq.gz ../raw_data/fastq/${name}_1.fastq.gz
    mv ../raw_data/fastq/${srr}_2.fastq.gz ../raw_data/fastq/${name}_2.fastq.gz
    echo "$srr succesvol geconverteerd en hernoemd."
    
# Zoek en verwijder het .sra bestand
    sra_file=$(find ~/.ncbi/public/sra -name "${srr}.sra" 2>/dev/null | head -n 1)
    if [[ -f "$sra_file" ]]; then
      rm "$sra_file"
      echo "🗑️  .sra bestand verwijderd: $sra_file"
    fi
  else
    echo "FASTQ-bestanden ontbreken voor $srr. Conversie mislukt?"
  fi
done

echo "Klaar met verwerken van alle samples. Controleer 'fastq_errors.log' indien nodig."


#uitvoerbaar maken chmod +x 2download_fastq.sh (eenmalig)
#runnen met ./2download_fastq.sh

#E8.5-19_i36 SRR14783084 .sra bestand is 3.3 GB fastq files zijn ~ 6 GB voor 10 sampels is dat ~ 30 GB aan sra files en 60 GB aan fastq files. 
#Er is dus ongeveer 100 GB nodig voor downloaden. Op dit moment 21JUN25 is er 287 GB open op de server. 
#Dit script is getest op het eerste sample en werkt, inverband met opslag nog niet alles gedownload. 
