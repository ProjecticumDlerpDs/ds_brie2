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

mkdir -p fastq

for sample in "${samples[@]}"
do
name=$(echo $sample | cut -d' ' -f1)
srr=$(echo $sample | cut -d' ' -f2)

echo "Processing sample: $name with SRR: $srr"

# Download de .sra file
prefetch $srr

# Zet om naar fastq en sla op in fastq/ map, naam prefix met sample naam
fasterq-dump $srr --split-files --gzip -O fastq/
  
# Hernoem bestanden naar herkenbare namen
mv fastq/${srr}_1.fastq.gz fastq/${name}_1.fastq.gz
mv fastq/${srr}_2.fastq.gz fastq/${name}_2.fastq.gz
done


#uitvoerbaar maken chmod +x download_fastq.sh (eenmalig)
#runnen met /download_fastq.sh

