#!/bin/bash

# SRA-accession nummers met bijbehorende sample namen
# Format: SampleName SRRnummer

samples=(
  "E8.5-19_i36 SRR14783084"
  "E8.5-1_i7 SRR14783085")

mkdir -p fastq

for sample in "${samples[@]}"
do
name=$(echo $sample | cut -d' ' -f1)
srr=$(echo $sample | cut -d' ' -f2)

echo "Processing sample: $name with SRR: $srr"

# Download de .sra file
prefetch $srr

# Zet om naar fastq en sla op in fastq/ map, naam prefix met sample naam
fasterq-dump $srr --split-files --gzip -O fastq/ -t ./tmp -e 4 2>> fasterq_errors.log
  
# Hernoem bestanden naar herkenbare namen
if [[ -f fastq/${srr}_1.fastq.gz && -f fastq/${srr}_2.fastq.gz ]]; then
    mv fastq/${srr}_1.fastq.gz fastq/${name}_1.fastq.gz
    mv fastq/${srr}_2.fastq.gz fastq/${name}_2.fastq.gz
    echo "$srr succesvol geconverteerd en hernoemd."
    
# Zoek en verwijder het .sra bestand
    sra_file=$(find ~/.ncbi -name "${srr}.sra" 2>/dev/null | head -n 1)
    if [[ -f "$sra_file" ]]; then
      rm "$sra_file"
      echo "🗑️  .sra bestand verwijderd: $sra_file"
    fi
  else
    echo "FASTQ-bestanden ontbreken voor $srr. Conversie mislukt?"
  fi
done

echo "Klaar met verwerken van alle samples. Controleer 'fasterq_errors.log' indien nodig."

.
#uitvoerbaar maken chmod +x download_fastq.sh (eenmalig)
#runnen met ./2download_fastq.sh

