# Rhizosphere Microbiome of Chickpea Under Combined Stress and Microbial Consortium Inoculation

> Whole-genome shotgun metagenomics of chickpea (*Cicer arietinum*) roots across three field locations, comparing healthy controls, combined-stress plants, and combined-stress plants inoculated with a microbial consortium.

---

## Table of Contents

- [Study Overview](#study-overview)
- [Data Availability](#data-availability)
- [Repository Structure](#repository-structure)
- [Pipeline Overview](#pipeline-overview)
- [Requirements](#requirements)
- [Step-by-Step Usage](#step-by-step-usage)
  - [1. Quality Control](#1-quality-control)
  - [2. Host Read Removal](#2-host-read-removal)
  - [3. Taxonomic Profiling](#3-taxonomic-profiling)
  - [4. Co-assembly](#4-co-assembly)
  - [5. Read Mapping](#5-read-mapping)
  - [6. Contig Taxonomy (CAT)](#6-contig-taxonomy-cat)
  - [7. Gene Prediction](#7-gene-prediction)
  - [8. Functional Annotation](#8-functional-annotation)
  - [9. CAZyme Annotation](#9-cazyme-annotation)
  - [10. Read Counting](#10-read-counting)
  - [11. Downstream R Analysis](#11-downstream-r-analysis)
- [Sample Design](#sample-design)
- [Software Versions](#software-versions)
- [Citation](#citation)

---

## Study Overview

Root samples were collected during the 2024–25 rabi season at 50% flowering stage from three field locations (L1, L3, L5). Three treatment groups were profiled:

| Treatment | Code | Description |
|---|---|---|
| Control | CL | Healthy plants, no stress or inoculation |
| Combined Stress | CS | Disease & drought stress, no inoculation |
| Consortium + Stress | MC | Combined stress with microbial consortium |

Nine samples per location (3 treatments × 3 replicates) = **27 samples total**.  
Sequencing: Illumina NovaSeq 6000, S4 flow cell, 2×150 bp paired-end, ~5–10 GB per sample.

---

## Data Availability

Raw sequencing data are deposited at NCBI under BioProject **[PRJNA1464083](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1464083)**.

---

## Repository Structure

```
.
├── scripts/
│   ├── 01_trim.sh                   
│   ├── 02_kneaddata.slurm           
│   ├── 03_kraken_loc.slurm          
│   ├── 04_global_assembly.pbs       
│   ├── 05_global_bowtie_ind.pbs     
│   ├── 06_global_bowtie_map.pbs     
│   ├── 07_global_cat.pbs            
│   ├── 08_global_prokka.pbs         
│   ├── 09_global_eggnog.pbs         
│   └── 10_global_dbcan.pbs          
│
├── analysis/
│   ├── alpha_beta_analysis.Rmd      
│   ├── maaslin_eggnog.R             # KEGG pathway differential abundance (MaAsLin2)
│   ├── network_treatment.R          # Co-occurrence network per treatment (NetCoMi)
│   ├── ternary_plot.R               
│   └── distribution_script.R        # Abundance distribution visualisation
│
└── README.md
```

---
## Sample Design

| Sample | Location | Treatment |
|---|---|---|
| L1-1, L1-2, L1-3 | Location 1 | Control (CL) |
| L1-4, L1-5, L1-6 | Location 1 | Combined Stress (CS) |
| L1-7, L1-8, L1-9 | Location 1 | CS + Consortium (MC) |
| L3-10 … L3-12 | Location 3 | Control (CL) |
| L3-13 … L3-15 | Location 3 | Combined Stress (CS) |
| L3-16 … L3-18 | Location 3 | CS + Consortium (MC) |
| L5-19 … L5-21 | Location 5 | Control (CL) |
| L5-22 … L5-24 | Location 5 | Combined Stress (CS) |
| L5-25 … L5-27 | Location 5 | CS + Consortium (MC) |

---


## Pipeline Overview

```
Raw reads (27 samples, 2×150 bp)
        │
        ▼
[1] fastp — adapter trimming, quality filtering
        │
        ▼
[2] KneadData — host (chickpea ICC4958) read removal
        │
        ├──────────────────────────────────────────────────────┐
        ▼                                                      ▼
[3] Kraken2 + Bracken                              [4] MEGAHIT co-assembly
    taxonomic profiling                                (all 27 samples pooled)
    (Kraken2 core_nt DB,                               --presets meta-sensitive
     confidence = 0.1)                                 contigs ≥ 1000 bp
        │                                                      │
        ▼                                             ┌────────┴────────┐
  phyloseq object                                     ▼                 ▼
  alpha/beta diversity                      [5] Bowtie2 mapping    [7] Prokka
  PERMANOVA                                    per sample              ORF prediction
  MaAsLin2 KEGG analysis                       (--very-sensitive)      (--metagenome)
  co-occurrence networks                       → sorted BAM               │
                                                    │              ┌───────┴──────┐
                                               [6] CAT_pack        ▼              ▼
                                               contig taxonomy  [8] eggNOG    [9] dbCAN
                                                               KEGG/COG/GO    CAZymes
                                                                    │
                                                               [10] featureCounts
                                                               per-gene counts
                                                               (27 samples × genes)
                                                                    │
                                                               R downstream analysis
```

---

## Requirements


## Software Versions

```
fastp           0.20.1
KneadData       (host removal)
Kraken2         2.1.2
Bracken         (species-level re-estimation, r=150)
MEGAHIT         1.2.9
Bowtie2         2.5.5
SAMtools        1.18
CAT_pack        (NR database: 20241212)
Prokka          1.14.6
eggNOG-mapper   2.1.13  (eggNOG 5.0 database)
DIAMOND         (bundled with eggNOG-mapper)
dbCAN           5.2.8   (CAZy database)
featureCounts   2.1.1   (Subread package)

R               4.5.2
  phyloseq      1.54.2
  vegan         2.7-2
  MaAsLin2      1.24.1
  NetCoMi       1.2.0
  edgeR         (TMM normalisation)
  FAPROTAX      1.2.12
  FUNGuild      1.1
```


---

## Method

### 1. Quality Control

```bash
# Trims adapters, filters low-quality reads, runs FastQC + MultiQC
# Place raw reads in raw_reads/ as *_R1.fastq.gz / *_R2.fastq.gz

bash scripts/01_trim.sh
```

**Output:** `trimmed_reads/`, `fastqc_trimmed/`, `trimmed_multiqc_report/`

---

### 2. Host Read Removal

```bash
# Removes chickpea (ICC4958) host reads using KneadData
# Requires Bowtie2 index of chickpea genome at ref/ICC4958_index

# SLURM array job — one job per sample (27 total)
sbatch scripts/02_kneaddata.slurm
```

---

### 3. Taxonomic Profiling

```bash
# Kraken2 classification against core_nt database
# Bracken re-estimation at species level (read length = 150 bp)
# Confidence threshold = 0.1

sbatch scripts/03_kraken_loc.slurm
```

> Downstream taxonomic analysis (alpha/beta diversity, PERMANOVA, network analysis) is performed in R — see `analysis/alpha_beta_analysis.Rmd` and `analysis/network_treatment.R`.

---

### 4. Co-assembly

```bash
# Pools all 27 samples for co-assembly with MEGAHIT
# --presets meta-sensitive, --min-contig-len 1000 bp, 20 threads

qsub scripts/04_global_assembly.pbs
```


---

### 5. Read Mapping

```bash
# Build Bowtie2 index from co-assembly
qsub scripts/05_global_bowtie_ind.pbs

# Map each sample back to co-assembly (array job, 5 concurrent)
# --very-sensitive, --no-unal, sorted BAM + flagstat per sample
qsub scripts/06_global_bowtie_map.pbs
```

---

### 6. Contig Taxonomy (CAT)

```bash
# Assigns taxonomy to contigs using CAT_pack against NR database

qsub scripts/07_global_cat.pbs
```

---

### 7. Gene Prediction

```bash
# Predicts ORFs from co-assembly contigs using Prokka
# --metagenome mode, 20 threads

qsub scripts/08_global_prokka.pbs
```

---

### 8. Functional Annotation

```bash
# Annotates predicted proteins with eggNOG-mapper
# DIAMOND search, sensitive mode, eggNOG 5.0 database
# E-value < 1e-3, bit score >= 60, GO terms: non-electronic only

qsub scripts/09_global_eggnog.pbs
```

> KEGG pathway differential abundance analysis is in `analysis/maaslin_eggnog.R`

---

### 9. CAZyme Annotation

```bash
# Annotates CAZymes using dbCAN v5.2.8
# Both DIAMOND and HMM methods against CAZy database

qsub scripts/10_global_dbcan.pbs
```

---

### 10. Read Counting

```bash
# featureCounts counts reads per CDS feature per sample
# Paired-end mode, multi-mappers fractionally distributed
bash scripts/11_featurecounts.sh

```

---

### 11. Downstream R Analysis


```r
# Alpha/beta diversity, PERMANOVA, compositional analysis
rmarkdown::render("analysis/alpha_beta_analysis.Rmd")

# KEGG pathway differential abundance (MaAsLin2)
Rscript analysis/maaslin_eggnog.R

# Co-occurrence networks per treatment (NetCoMi, CLR + location regression)
# Set TREATMENT variable to "Control", "CS", or "CS with consortium"
Rscript analysis/network_treatment.R

# Abundance distribution
Rscript analysis/distribution_script.R
```

#### Network Analysis Details

Co-occurrence networks were constructed per treatment using `network_treatment.R`:

- **Normalisation:** CLR transformation with multiplicative zero replacement (zCompositions CZM)
- **Location correction:** Location effect regressed out from CLR values before correlation
- **Correlation:** Pairwise Pearson on CLR residuals
- **Significance:** Bootstrap permutation test (999 permutations)
- **Threshold:** |r| > 0.82, p < 0.01
- **Topology:** Computed via NetCoMi v1.2.0 (degree, betweenness, clustering coefficient, modularity, natural connectivity)
- **Module detection:** Fast-greedy algorithm
- **Node roles:** Zi–Pi framework (Network Hubs, Module Hubs, Connectors, Peripherals)
- **Visualisation:** Exported as GraphML for Cytoscape / Gephi v0.10.1

---

## Sample Design

| Sample | Location | Treatment |
|---|---|---|
| L1-1, L1-2, L1-3 | Location 1 | Control (CL) |
| L1-4, L1-5, L1-6 | Location 1 | Combined Stress (CS) |
| L1-7, L1-8, L1-9 | Location 1 | CS + Consortium (MC) |
| L3-10 … L3-12 | Location 3 | Control (CL) |
| L3-13 … L3-15 | Location 3 | Combined Stress (CS) |
| L3-16 … L3-18 | Location 3 | CS + Consortium (MC) |
| L5-19 … L5-21 | Location 5 | Control (CL) |
| L5-22 … L5-24 | Location 5 | Combined Stress (CS) |
| L5-25 … L5-27 | Location 5 | CS + Consortium (MC) |

---

