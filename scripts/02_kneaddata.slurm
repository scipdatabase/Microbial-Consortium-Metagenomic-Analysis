#!/bin/bash
#SBATCH --job-name=knead_host
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1-3
#SBATCH --time=24:00:00
#SBATCH --array=1-27
#SBATCH --output=log/knead_%A_%a.out
#SBATCH --error=log/knead_%A_%a.err

# Define sample ID using SLURM array task ID
SAMPLES=(
  L1-1 L1-2 L1-3 L1-4 L1-5 L1-6 L1-7 L1-8 L1-9
  L3-10 L3-11 L3-12 L3-13 L3-14 L3-15 L3-16 L3-17 L3-18
  L5-19 L5-20 L5-21 L5-22 L5-23 L5-24 L5-25 L5-26 L5-27
)

SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}

# Input files
R1="reads/${SAMPLE}_R1_paired.fq.gz"
R2="reads/${SAMPLE}_R2_paired.fq.gz"

# Output folder
OUTDIR="host_rem/host_rem_${SAMPLE}"

# Path to chickpea bowtie2 index
REF_DB="ref/ICC4958_index"

# Create output directory
mkdir -p ${OUTDIR}

echo "Processing sample: ${SAMPLE}"

kneaddata \
  --input1 ${R1} \
  --input2 ${R2} \
  --reference-db ${REF_DB} \
  --output ${OUTDIR} \
  --threads 16 \
  --bypass-trim \
  --remove-intermediate-output

echo "Finished sample: ${SAMPLE}"
