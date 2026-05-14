library(readxl)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(vegan)
library(patchwork)
library(pals)
library(RColorBrewer)
library(scales)


path = "D:/SKM/Metagenome/dwn_analysis"
setwd(path)
set.seed(124)

#====================================
#Loading data and metadata
#====================================

raw_data <- read.csv("D:/SKM/metagenome/analysis/kraken/Newest_table_No_Outliers.csv") %>%
  mutate(Kingdom = str_replace(Kingdom,                  
                               "Unclassified", "Others"))

#Metadata
metadata <- data.frame(
  Sample_ID = c("L1.1","L1.2","L1.3","L1.4","L1.5","L1.6","L1.7","L1.8","L1.9",
                "L3.10","L3.11","L3.12","L3.13","L3.14","L3.15","L3.16","L3.17","L3.18",
                "L5.19","L5.20","L5.21","L5.22","L5.23","L5.24","L5.25","L5.26","L5.27"),
  Location = c(rep("Location 1", 9), rep("Location 3", 9), rep("Location 5", 9)),
  Treatment = c(rep("CL", 3), rep("DD", 3), rep("MC", 3),
                rep("CL", 3), rep("DD", 3), rep("MC", 3),
                rep("CL", 3), rep("DD", 3), rep("MC", 3))
)

#====================================
#Phyloseq object
#====================================
otu_mat <- as.data.frame(raw_data[, 9:ncol(raw_data)])
rownames(otu_mat) <- raw_data[[1]] 
OTU <- otu_table(as.matrix(otu_mat), taxa_are_rows = TRUE)

tax_mat <- as.data.frame(raw_data[, 1:8])
rownames(tax_mat) <- tax_mat[[1]]
tax_mat_final <- as.matrix(tax_mat[, -1])
TAX <- tax_table(tax_mat_final)

metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$Sample_ID
META <- sample_data(metadata_df)

# Check the names in your OTU table vs Metadata
print(colnames(OTU))
print(rownames(META))

# See which ones are missing from the metadata
setdiff(colnames(OTU), rownames(META))
ps <- phyloseq(OTU, TAX, META)
saveRDS(ps, "ps_raw.rds")

#====================================
#Prevalence filtering
#====================================
PREV_THRESHOLD <- 0.10

ps_filt <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= PREV_THRESHOLD * nsamples(ps),
  prune = TRUE
)

ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x))

get_ra_long <- function(ps_obj, rank = "Phylum", kingdom_filter = NULL) {
  
  df <- psmelt(ps_obj)                                   # phyloseq -> long data frame
  
  if (!is.null(kingdom_filter)) {
    df <- df %>% filter(Kingdom == kingdom_filter)
  }
  
  # Re-normalise within kingdom after subsetting so bars sum to 100%
  df <- df %>%
    group_by(Sample, Location, Treatment) %>%
    mutate(RA_within_kingdom = Abundance / sum(Abundance) * 100) %>%
    ungroup()
  
  # Rename the target rank column to "Taxon" for uniform downstream use
  df <- df %>% rename(Taxon = all_of(rank))
  
  return(df)
}


prev_df <- data.frame(
  Taxon      = taxa_names(ps),
  Prevalence = apply(otu_table(ps), 1, function(x) sum(x > 0)),
  TotalReads = taxa_sums(ps)
) %>%
  mutate(PrevPct = Prevalence / nsamples(ps) * 100)

cat("--- Prevalence summary (all taxa, raw) ---\n")
cat("Total taxa:                  ", nrow(prev_df), "\n")
cat("Taxa in >=5% samples:        ", sum(prev_df$PrevPct >= 5), "\n")
cat("Taxa in >=10% samples:       ", sum(prev_df$PrevPct >= 10), "\n")
cat("Taxa in >=20% samples:       ", sum(prev_df$PrevPct >= 20), "\n")
cat("Reads retained at 20% cutoff:", 
    round(sum(prev_df$TotalReads[prev_df$PrevPct >= 20]) / sum(prev_df$TotalReads) * 100, 1), "%\n")

colour_pal <- rev(c("#7FB07D", "#74171B", "#F42B1F", "#FFCF1D", "#755819", "#003B7B", "#FA972F", "#287B1E",
                "#D272DC", "#8D099B","#d3b3b0", "#6D48B8","#67C8F5","#7436BD","#2F5E5F", "#FF2573" ))


plot_top_n_composition <- function(ps_obj,
                                   kingdom_filter,
                                   rank       = "Phylum",
                                   top_n      = 15,
                                   out_file   = NULL) {
  
  df <- get_ra_long(ps_obj, rank = rank, kingdom_filter = kingdom_filter)
  
 
  top_taxa <- df %>%
    group_by(Taxon) %>%
    summarise(grand_total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(!Taxon %in% c("Others", "Unclassified", "NA", NA)) %>%
    # slice_max(order_by = grand_total, n = top_n, with_ties = FALSE) %>%
    pull(Taxon)
  

  df_plot <- df %>%
    mutate(Taxon = if_else(Taxon %in% top_taxa, as.character(Taxon), "Others")) %>%
    # First sum per sample (handling multi-OTU assignments to same rank)
    group_by(Location, Treatment, Sample, Taxon) %>%
    summarise(RA_per_sample = sum(RA_within_kingdom), .groups = "drop") %>%
    # Then average the replicates (n=3) for the final bar
    group_by(Location, Treatment, Taxon) %>%
    summarise(RA = mean(RA_per_sample, na.rm = TRUE), .groups = 'drop')
  

  taxon_order <- df_plot %>%
    group_by(Taxon) %>%
    summarise(total_ra = sum(RA), .groups = "drop") %>%
    arrange(desc(total_ra)) %>%
    mutate(is_others = Taxon == "Others") %>%
    arrange(is_others, desc(total_ra)) %>%
    pull(Taxon)
  
  df_plot$Taxon <- factor(df_plot$Taxon, levels = taxon_order)
  
  # 4. COLOR MAPPING: Match ordered levels to palette
  taxa_names <- levels(df_plot$Taxon)
  core_taxa  <- taxa_names[taxa_names != "Others"]
  
  # Assign colors from your wheel to core taxa
  my_colors <- setNames(colour_pal[1:length(core_taxa)], core_taxa)
  
  if ("Others" %in% taxa_names) {
    my_colors["Others"] <- "grey70"
  }
  
  # PLOTTING
  p <- ggplot(df_plot, aes(x = Treatment, y = RA, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack",
             color = "white", linewidth = 0.15, width = 0.8) +
    facet_grid(~ Location, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = my_colors) +
    # Use expand = c(0,0) to remove the empty space gap
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) + 
    theme_minimal() +
    labs(y = "Relative Abundance (%)", x = NULL, fill = rank) +
    theme(
      axis.text.x        = element_text(face = "bold", size = 9, colour = "black"),
      axis.text.y        = element_text(face = "bold", size = 9, colour = "black"),
      axis.title.y       = element_text(face = "bold", size = 9, colour = "black"),
      strip.text         = element_text(face = "bold", size = 9),
      panel.border       = element_rect(color = "gray80", fill = NA, linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      legend.key.size    = unit(0.7, "lines"),
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8)
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = 5, height = 5, dpi = 300, bg = "white")
  }
  
  return(invisible(p))
}



# Phylum level
plot_top_n_composition(ps_rel,
                       kingdom_filter = "Bacteria",
                       rank           = "Phylum",
                       top_n          = 15,
                       out_file       = "new_top15_bact_phy.png")

#Class level
plot_top_n_composition(ps_rel,
                       kingdom_filter = "Bacteria",
                       rank           = "Class",
                       top_n          = 15,
                       out_file       = "final_top15_bact_class.png")


plot_top_n_composition(ps_rel,
                       kingdom_filter = "Fungi",
                       rank           = "Phylum",
                       top_n          = 10,               
                       out_file       = "new_plot/final_top15_Phy_genus.png")

plot_top_n_composition(ps_rel,
                       kingdom_filter = "Fungi",
                       rank           = "Class",
                       top_n          = 10,
                       out_file       = "new_plot/final_top15_fungi_class.png")

