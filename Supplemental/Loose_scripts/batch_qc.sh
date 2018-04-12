#!/bin/bash

echo Starting GSE70580
COUNTER=1
TOTAL=$(ls ~/FYP/fastq/GSE70580 | grep .fastq | wc -l)
for file in ~/FYP/fastq/GSE70580/*.fastq
do
fastqc -f fastq -o ~/FYP/qc_reports/GSE70580 ${file}
echo $COUNTER of $TOTAL done
COUNTER=$((COUNTER+1))
done

