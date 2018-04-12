#!/bin/bash

rsem-prepare-reference --gtf ./ref/homo/gencode.v27.annotation.gtf --bowtie2 ./ref/homo/GRCh38.primary_assembly.genome.fa ./ref/homo/rsem-ref

rsem-prepare-reference --gtf ./ref/mus/gencode.vM16.annotation.gtf --bowtie2 ./ref/mus/GRCm38.primary_assembly.genome.fa ./ref/mus/rsem-ref

mkdir OUTPUT
mkdir OUTPUT/mus_cortex
mkdir OUTPUT/homo_lymphoid
mkdir OUTPUT/mus_embryo
cortexFile='fastq/cortex_sample.txt'
echo Starting Cortex Alignment
while read -r p; do
    echo "Aligning " $p
    rsem-calculate-expression -p 4 --bowtie2 --estimate-rspd --append-names --single-cell-prior --time --output-genome-bam fastq/cortex/${p}* ./ref/mus/rsem-ref OUTPUT/mus_cortex/${p}
    echo ${p} " done!"
done < $cortexFile

embryoFile='fastq/embryo_ids.txt'
echo Starting Embyro Alignment
while read -r p; do
    echo "Aligning " $p
    rsem-calculate-expression -p 4 --bowtie2 --estimate-rspd --append-names --single-cell-prior --time --output-genome-bam fastq/embryo/${p}* ./ref/mus/rsem-ref OUTPUT/mus_embryo/${p}
    echo ${p} " done!"
done < $embryoFile

echo "Aligning SRR2088075"
rsem-calculate-expression -p 4 --bowtie2 --estimate-rspd --append-names --single-cell-prior --time --output-genome-bam fastq/GSE70580/SRR2088075* ./ref/homo/rsem-ref OUTPUT/homo_lymphoid/SRR2088075