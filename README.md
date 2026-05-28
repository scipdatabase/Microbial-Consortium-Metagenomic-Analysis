# \# Rhizosphere Microbiome of Chickpea Under Combined Stress and Microbial Consortium Inoculation

# 

# > Whole-genome shotgun metagenomics of chickpea (\*Cicer arietinum\*) roots across three field locations, comparing healthy controls, combined-stress plants, and combined-stress plants inoculated with a microbial consortium.

# 

# \---

# 

# \## Table of Contents

# 

# \- \[Study Overview](#study-overview)

# \- \[Data Availability](#data-availability)

# \- \[Repository Structure](#repository-structure)

# \- \[Pipeline Overview](#pipeline-overview)

# \- \[Requirements](#requirements)

# \- \[Step-by-Step Usage](#step-by-step-usage)

# &#x20; - \[1. Quality Control](#1-quality-control)

# &#x20; - \[2. Host Read Removal](#2-host-read-removal)

# &#x20; - \[3. Taxonomic Profiling](#3-taxonomic-profiling)

# &#x20; - \[4. Co-assembly](#4-co-assembly)

# &#x20; - \[5. Read Mapping](#5-read-mapping)

# &#x20; - \[6. Contig Taxonomy (CAT)](#6-contig-taxonomy-cat)

# &#x20; - \[7. Gene Prediction](#7-gene-prediction)

# &#x20; - \[8. Functional Annotation](#8-functional-annotation)

# &#x20; - \[9. CAZyme Annotation](#9-cazyme-annotation)

# &#x20; - \[10. Read Counting](#10-read-counting)

# &#x20; - \[11. Downstream R Analysis](#11-downstream-r-analysis)

# \- \[Sample Design](#sample-design)

# \- \[Software Versions](#software-versions)

# \- \[Citation](#citation)

# 

# \---

# 

# \## Study Overview

# 

# Root samples were collected during the 2024–25 rabi season at 50% flowering stage from three field locations (L1, L3, L5). Three treatment groups were profiled:

# 

# | Treatment | Code | Description |

# |---|---|---|

# | Control | CL | Healthy plants, no stress or inoculation |

# | Combined Stress | CS | Disease \& drought stress, no inoculation |

# | Consortium + Stress | MC | Combined stress with microbial consortium |

# 

# Nine samples per location (3 treatments × 3 replicates) = \*\*27 samples total\*\*.  

# Sequencing: Illumina NovaSeq 6000, S4 flow cell, 2×150 bp paired-end, \~5–10 GB per sample.

# 

# \---

# 

# \## Data Availability

# 

# Raw sequencing data are deposited at NCBI under BioProject \*\*\[PRJNA1464083](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1464083)\*\*.

# 

# \---

# 

# \## Repository Structure

# 

# ```

# .

# ├── scripts/

# │   ├── 01\_trim.sh                   

# │   ├── 02\_kneaddata.slurm           

# │   ├── 03\_kraken\_loc.slurm          

# │   ├── 04\_global\_assembly.pbs       

# │   ├── 05\_global\_bowtie\_ind.pbs     

# │   ├── 06\_global\_bowtie\_map.pbs     

# │   ├── 07\_global\_cat.pbs            

# │   ├── 08\_global\_prokka.pbs         

# │   ├── 09\_global\_eggnog.pbs         

# │   └── 10\_global\_dbcan.pbs          

# │

# ├── analysis/

# │   ├── alpha\_beta\_analysis.Rmd      

# │   ├── maaslin\_eggnog.R             # KEGG pathway differential abundance (MaAsLin2)

# │   ├── network\_treatment.R          # Co-occurrence network per treatment (NetCoMi)

# │   ├── ternary\_plot.R               

# │   └── distribution\_script.R        # Abundance distribution visualisation

# │

# └── README.md

# ```

# 

# \---

# \## Sample Design

# 

# | Sample | Location | Treatment |

# |---|---|---|

# | L1-1, L1-2, L1-3 | Location 1 | Control (CL) |

# | L1-4, L1-5, L1-6 | Location 1 | Combined Stress (CS) |

# | L1-7, L1-8, L1-9 | Location 1 | CS + Consortium (MC) |

# | L3-10 … L3-12 | Location 3 | Control (CL) |

# | L3-13 … L3-15 | Location 3 | Combined Stress (CS) |

# | L3-16 … L3-18 | Location 3 | CS + Consortium (MC) |

# | L5-19 … L5-21 | Location 5 | Control (CL) |

# | L5-22 … L5-24 | Location 5 | Combined Stress (CS) |

# | L5-25 … L5-27 | Location 5 | CS + Consortium (MC) |

# 

# \---



# 

# \## Pipeline Overview

# 

# ```

# Raw reads (27 samples, 2×150 bp)

# &#x20;       │

# &#x20;       ▼

# \[1] fastp — adapter trimming, quality filtering

# &#x20;       │

# &#x20;       ▼

# \[2] KneadData — host (chickpea ICC4958) read removal

# &#x20;       │

# &#x20;       ├──────────────────────────────────────────────────────┐

# &#x20;       ▼                                                      ▼

# \[3] Kraken2 + Bracken                              \[4] MEGAHIT co-assembly

# &#x20;   taxonomic profiling                                (all 27 samples pooled)

# &#x20;   (Kraken2 core\_nt DB,                               --presets meta-sensitive

# &#x20;    confidence = 0.1)                                 contigs ≥ 1000 bp

# &#x20;       │                                                      │

# &#x20;       ▼                                             ┌────────┴────────┐

# &#x20; phyloseq object                                     ▼                 ▼

# &#x20; alpha/beta diversity                      \[5] Bowtie2 mapping    \[7] Prokka

# &#x20; PERMANOVA                                    per sample              ORF prediction

# &#x20; MaAsLin2 KEGG analysis                       (--very-sensitive)      (--metagenome)

# &#x20; co-occurrence networks                       → sorted BAM               │

# &#x20;                                                   │              ┌───────┴──────┐

# &#x20;                                              \[6] CAT\_pack        ▼              ▼

# &#x20;                                              contig taxonomy  \[8] eggNOG    \[9] dbCAN

# &#x20;                                                              KEGG/COG/GO    CAZymes

# &#x20;                                                                   │

# &#x20;                                                              \[10] featureCounts

# &#x20;                                                              per-gene counts

# &#x20;                                                              (27 samples × genes)

# &#x20;                                                                   │

# &#x20;                                                              R downstream analysis

# ```

# 

# \---

# 

# \## Requirements

# 

# 

# \## Software Versions

# 

# ```

# fastp           0.20.1

# KneadData       (host removal)

# Kraken2         2.1.2

# Bracken         (species-level re-estimation, r=150)

# MEGAHIT         1.2.9

# Bowtie2         2.5.5

# SAMtools        1.18

# CAT\_pack        (NR database: 20241212)

# Prokka          1.14.6

# eggNOG-mapper   2.1.13  (eggNOG 5.0 database)

# DIAMOND         (bundled with eggNOG-mapper)

# dbCAN           5.2.8   (CAZy database)

# featureCounts   2.1.1   (Subread package)

# 

# R               4.5.2

# &#x20; phyloseq      1.54.2

# &#x20; vegan         2.7-2

# &#x20; MaAsLin2      1.24.1

# &#x20; NetCoMi       1.2.0

# &#x20; edgeR         (TMM normalisation)

# &#x20; FAPROTAX      1.2.12

# &#x20; FUNGuild      1.1

# ```

# 



# \---

# 

# \## Method

# 

# \### 1. Quality Control

# 

# ```bash

# \# Trims adapters, filters low-quality reads, runs FastQC + MultiQC

# \# Place raw reads in raw\_reads/ as \*\_R1.fastq.gz / \*\_R2.fastq.gz

# 

# bash scripts/01\_trim.sh

# ```

# 

# \*\*Output:\*\* `trimmed\_reads/`, `fastqc\_trimmed/`, `trimmed\_multiqc\_report/`

# 

# \---

# 

# \### 2. Host Read Removal

# 

# ```bash

# \# Removes chickpea (ICC4958) host reads using KneadData

# \# Requires Bowtie2 index of chickpea genome at ref/ICC4958\_index

# 

# \# SLURM array job — one job per sample (27 total)

# sbatch scripts/02\_kneaddata.slurm

# ```

# 

# \---

# 

# \### 3. Taxonomic Profiling

# 

# ```bash

# \# Kraken2 classification against core\_nt database

# \# Bracken re-estimation at species level (read length = 150 bp)

# \# Confidence threshold = 0.1

# 

# sbatch scripts/03\_kraken\_loc.slurm

# ```

# 

# > Downstream taxonomic analysis (alpha/beta diversity, PERMANOVA, network analysis) is performed in R — see `analysis/alpha\_beta\_analysis.Rmd` and `analysis/network\_treatment.R`.

# 

# \---

# 

# \### 4. Co-assembly

# 

# ```bash

# \# Pools all 27 samples for co-assembly with MEGAHIT

# \# --presets meta-sensitive, --min-contig-len 1000 bp, 20 threads

# 

# qsub scripts/04\_global\_assembly.pbs

# ```

# 

# 

# \---

# 

# \### 5. Read Mapping

# 

# ```bash

# \# Build Bowtie2 index from co-assembly

# qsub scripts/05\_global\_bowtie\_ind.pbs

# 

# \# Map each sample back to co-assembly (array job, 5 concurrent)

# \# --very-sensitive, --no-unal, sorted BAM + flagstat per sample

# qsub scripts/06\_global\_bowtie\_map.pbs

# ```

# 

# \---

# 

# \### 6. Contig Taxonomy (CAT)

# 

# ```bash

# \# Assigns taxonomy to contigs using CAT\_pack against NR database

# 

# qsub scripts/07\_global\_cat.pbs

# ```

# 

# \---

# 

# \### 7. Gene Prediction

# 

# ```bash

# \# Predicts ORFs from co-assembly contigs using Prokka

# \# --metagenome mode, 20 threads

# 

# qsub scripts/08\_global\_prokka.pbs

# ```

# 

# \---

# 

# \### 8. Functional Annotation

# 

# ```bash

# \# Annotates predicted proteins with eggNOG-mapper

# \# DIAMOND search, sensitive mode, eggNOG 5.0 database

# \# E-value < 1e-3, bit score >= 60, GO terms: non-electronic only

# 

# qsub scripts/09\_global\_eggnog.pbs

# ```

# 

# > KEGG pathway differential abundance analysis is in `analysis/maaslin\_eggnog.R`

# 

# \---

# 

# \### 9. CAZyme Annotation

# 

# ```bash

# \# Annotates CAZymes using dbCAN v5.2.8

# \# Both DIAMOND and HMM methods against CAZy database

# 

# qsub scripts/10\_global\_dbcan.pbs

# ```

# 

# \---

# 

# \### 10. Read Counting

# 

# ```bash

# \# featureCounts counts reads per CDS feature per sample

# \# Paired-end mode, multi-mappers fractionally distributed

# bash scripts/11\_featurecounts.sh

# 

# ```

# 

# \---

# 

# \### 11. Downstream R Analysis

# 

# 

# ```r

# \# Alpha/beta diversity, PERMANOVA, compositional analysis

# rmarkdown::render("analysis/alpha\_beta\_analysis.Rmd")

# 

# \# KEGG pathway differential abundance (MaAsLin2)

# Rscript analysis/maaslin\_eggnog.R

# 

# \# Co-occurrence networks per treatment (NetCoMi, CLR + location regression)

# \# Set TREATMENT variable to "Control", "CS", or "CS with consortium"

# Rscript analysis/network\_treatment.R

# 

# \# Abundance distribution

# Rscript analysis/distribution\_script.R

# ```

# 

# \#### Network Analysis Details

# 

# Co-occurrence networks were constructed per treatment using `network\_treatment.R`:

# 

# \- \*\*Normalisation:\*\* CLR transformation with multiplicative zero replacement (zCompositions CZM)

# \- \*\*Location correction:\*\* Location effect regressed out from CLR values before correlation

# \- \*\*Correlation:\*\* Pairwise Pearson on CLR residuals

# \- \*\*Significance:\*\* Bootstrap permutation test (999 permutations)

# \- \*\*Threshold:\*\* |r| > 0.82, p < 0.01

# \- \*\*Topology:\*\* Computed via NetCoMi v1.2.0 (degree, betweenness, clustering coefficient, modularity, natural connectivity)

# \- \*\*Module detection:\*\* Fast-greedy algorithm

# \- \*\*Node roles:\*\* Zi–Pi framework (Network Hubs, Module Hubs, Connectors, Peripherals)

# \- \*\*Visualisation:\*\* Exported as GraphML for Cytoscape / Gephi v0.10.1

# 

# \---

# 

# \## Sample Design

# 

# | Sample | Location | Treatment |

# |---|---|---|

# | L1-1, L1-2, L1-3 | Location 1 | Control (CL) |

# | L1-4, L1-5, L1-6 | Location 1 | Combined Stress (CS) |

# | L1-7, L1-8, L1-9 | Location 1 | CS + Consortium (MC) |

# | L3-10 … L3-12 | Location 3 | Control (CL) |

# | L3-13 … L3-15 | Location 3 | Combined Stress (CS) |

# | L3-16 … L3-18 | Location 3 | CS + Consortium (MC) |

# | L5-19 … L5-21 | Location 5 | Control (CL) |

# | L5-22 … L5-24 | Location 5 | Combined Stress (CS) |

# | L5-25 … L5-27 | Location 5 | CS + Consortium (MC) |

# 

# \---

# 

# 

