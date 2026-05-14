# install.packages("ggtern")

library(readxl)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(vegan)
library(patchwork)
library(pals)
library(RColorBrewer)
library(scales)


path = "."
setwd(path)
set.seed(124)

#====================================
#Loading data and metadata
#====================================

raw_data <- read.csv("fcr_abundance_table_2404.csv")

#Metadata
metadata <- data.frame(
  Sample_ID = c("CS.1", "CS.2" ,"CS.3", "P.1", "P.2",  "P.3",  "C.1" ,"C.2","C.3"),
  Treatment = c(rep("Disease and Drought", 3), rep("Disease", 3), rep("Control", 3))
)

family <- c("Botryosphaeriaceae",
            "Nectriaceae",
            "Nitrobacteraceae",
            "Burkholderiaceae",
            "Devosiaceae",
            "Phyllobacteriaceae",
            "Sphingomonadaceae",
            "Rhizobiaceae",
            "Glomeraceae")

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
PREV_THRESHOLD <- 0.05

ps_filt <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= PREV_THRESHOLD * nsamples(ps),
  prune = TRUE
)
ps_family <- subset_taxa(ps_filt, Family %in% family)

ps_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))

get_ra_long <- function(ps_obj, rank = "Family", kingdom_filter = NULL) {
  
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



library(ggtern)

tern_df <- ps_rel %>%
  psmelt() %>%
  filter(Kingdom %in% c("Fungi", "Bacteria")) %>%
  filter(Family != "" & !is.na(Family) & Family != " ") %>% 
  group_by(Treatment, Kingdom, Family, Genus) %>% 
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = Abundance, values_fill = 0) %>%
  mutate(Total_RA = Control + Disease + `Disease and Drought`)

colour_pal <- rev(c("#7FB07D", "#74171B", "#F42B1F", "#8D099B", "#755819", "#003B7B", 
                 "#287B1E", "#D272DC", "#FFCF1D", "#d3b3b0", "#FA972F", 
                "#67C8F5", "#7436BD", "#2F5E5F", "#FF2573"))

p <- ggtern(data = tern_df, aes(x = Control, y = Disease, z = `Disease and Drought`)) +
  geom_point(aes(color = Family, size = Total_RA, shape = Kingdom), alpha = 0.8) +
  theme_rgbw() +
  theme_hidetitles() + 
  scale_color_manual(values = colour_pal) + 
  labs(
    size  = "Relative Abundance (%)",
    x     = "Control", 
    y     = "Disease", 
    z     = "Disease and Drought"
  ) +
  theme_showarrows() +
  theme(
    legend.position = "bottom",
    legend.text     = element_blank(),
    legend.title    = element_blank(),
    tern.axis.arrow.text = element_text(size = 12)
  )

print(p)

ggsave("ternary_plot.png", plot = p, width = 7, height = 5, dpi = 300, bg = "white")
