#!/bin/bash

## example to show for one dataset
## text file input contains all the SRA run IDs
HSCFile='GSE48968_accessions.txt'
echo Starting Cortex Fastq
while read -r p; do
    fastq-dump -I --split-files -O GSE48968/ $p
done < $HSCFile

