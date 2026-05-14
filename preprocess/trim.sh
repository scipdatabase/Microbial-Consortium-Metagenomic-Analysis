#!/bin/bash

mkdir fastqc_trimmed/

# For trimming adaptor

for f1 in raw_reads/*_R1.fastq.gz;
do

f2=${f1/_R1.fastq.gz/_R2.fastq.gz}

base=$(basename $f1 _R1.fastq.gz)

echo "===Trimming sample using fastp: $base"

fastp -w 6 \
-i "$f1" -I "$f2" \
-o "trimmed_reads/${base}_R1_paired.fq.gz" -O "trimmed_reads/${base}_R2_paired.fq.gz" \
--unpaired1 "${base}_R1_unpaired.fq.gz" --unpaired2 "${base}_R2_unpaired.fq.gz" \
-h "trimmed_reads/${base}fastp.html" -j "trimmed_reads/${base}fastp.json" 
done

#for fastqc
echo "===Running fastqc...."
fastqc -t 6 trimmed_reads/*_paired.fq.gz -o fastqc_trimmed/

#Multiqc
echo "===Generating multiqc report...."
multiqc fastqc_trimmed/ -o trimmed_multiqc_report/



