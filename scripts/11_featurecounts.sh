#!/bin/bash

featureCounts \
  -a prokka_global/global.gff \
  -o featurecounts/gene_counts.txt \
  -t CDS -g ID \
  -p --countReadPairs \
  -T 60 \
  -M --fraction \
  bowtie2_mapping_global/*.sorted.bam
