#!/bin/bash
#SBATCH --job-name=kraken_vriksha
#SBATCH --output=log/kraken_%A_%a.out
#SBATCH --error=log/kraken_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --partition=vriksha
#SBATCH --array=1-27

source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken_env

# --- PATHS ---
BASE_DIR="/lustre/scratch/skmlab107/amina/trim_files"
DB_PATH="${BASE_DIR}/db"
OUT_DIR="${BASE_DIR}/output"
mkdir -p "$OUT_DIR"

# --- SAMPLE MAPPING ---
if [ $SLURM_ARRAY_TASK_ID -le 9 ]; then
    LOC="L1"
elif [ $SLURM_ARRAY_TASK_ID -le 18 ]; then
    LOC="L3"
else
    LOC="L5"
fi

SAMPLE="${LOC}-${SLURM_ARRAY_TASK_ID}"
SAMPLE_DIR="host_rem_${SAMPLE}"

# Define input paths (Relative to your host_rem directory)
R1="${BASE_DIR}/host_rem/${SAMPLE_DIR}/${SAMPLE}_R1_paired_kneaddata_paired_1.fastq"
R2="${BASE_DIR}/host_rem/${SAMPLE_DIR}/${SAMPLE}_R1_paired_kneaddata_paired_2.fastq"

# --- ANALYSIS ---
# Use 0.05 as your primary confidence for rhizosphere samples
CONFIDENCES=(0.1)

if [[ -f "$R1" && -f "$R2" ]]; then
    for CONF in "${CONFIDENCES[@]}"; do
        
        CONF_OUT_DIR="${OUT_DIR}/conf_${CONF}"
        mkdir -p "$CONF_OUT_DIR"

        REPORT="${CONF_OUT_DIR}/${SAMPLE}_report.txt"
        OUTFILE="${CONF_OUT_DIR}/${SAMPLE}.kraken"
        BRACKEN="${CONF_OUT_DIR}/${SAMPLE}_bracken.txt"

        echo "Processing ${SAMPLE} | Confidence: ${CONF} | Node: $SLURMD_NODENAME"

        # Kraken2 Run
        kraken2 --db "$DB_PATH" --threads "$SLURM_CPUS_PER_TASK" --paired "$R1" "$R2" \
                --report "$REPORT" --output "$OUTFILE" --confidence "$CONF"

        # Bracken Run (Ensure -r 150 matches your Illumina run length)
        if [ -s "$REPORT" ]; then
            bracken -d "$DB_PATH" -i "$REPORT" -o "$BRACKEN" -r 150 -l S
        fi
        
        # Space Saving: The .kraken output is huge (~20GB+). 
        # Delete it after the report is generated to stay under your 2TB quota.
        rm "$OUTFILE"
    done
else
    echo "ERROR: Files for ${SAMPLE} not found at $R1"
    exit 1
fi
