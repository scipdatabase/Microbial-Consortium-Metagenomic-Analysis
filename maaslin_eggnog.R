# ============================================================================
# Full pipeline: eggNOG → KEGG matrix → MaAsLin2 → Volcano plots
# ============================================================================

library(tidyverse)
library(KEGGREST)
library(Maaslin2)
library(ggrepel)
library(patchwork)

path <- "./eggnog_global/"
setwd(path)
egg <- read.csv("eggnog_count.csv")

#Single-letter COG categories
egg <- egg %>% filter(nchar(COG_category) == 1)

sample_cols <- c("L1.1","L1.2","L1.3","L1.4","L1.5","L1.6","L1.7","L1.8","L1.9",
                 "L3.10","L3.11","L3.12","L3.13","L3.14","L3.15","L3.16","L3.17","L3.18",
                 "L5.19","L5.20","L5.21","L5.22","L5.23","L5.24","L5.25","L5.26","L5.27")

# Metadata 
metadata <- data.frame(
  Sample_ID = sample_cols,
  Location  = c(rep("L1",9), rep("L3",9), rep("L5",9)),
  Treatment = c(rep("CL",3), rep("DD",3), rep("MC",3),
                rep("CL",3), rep("DD",3), rep("MC",3),
                rep("CL",3), rep("DD",3), rep("MC",3)),
  stringsAsFactors = FALSE
) %>%
  mutate(Treatment = factor(Treatment, levels = c("CL","DD","MC")))

# Sequencing depth 
seq_depth <- egg %>%
  select(all_of(sample_cols)) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(),
               names_to  = "Sample_ID",
               values_to = "depth") %>%
  mutate(log_depth = log(depth))

# Depth distribution
cat("Depth summary:\n")
print(seq_depth)

# Plot depth
p_depth <- ggplot(seq_depth,
                  aes(x    = reorder(Sample_ID, log_depth),
                      y    = log_depth,
                      fill = sub("\\..*", "", Sample_ID))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(seq_depth$log_depth),
             linetype = "dashed", color = "red") +
  labs(x     = "Sample",
       y     = "log(sequencing depth)",
       fill  = "Location",
       title = "Sequencing depth per sample") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_depth)
ggsave("sequencing_depth.png", p_depth, width = 10, height = 6, dpi = 300)

# Depth summary by treatment
depth_by_treatment <- seq_depth %>%
  left_join(metadata, by = "Sample_ID") %>%
  group_by(Treatment) %>%
  summarise(
    mean_depth = mean(depth),
    sd_depth   = sd(depth),
    min_depth  = min(depth),
    max_depth  = max(depth),
    .groups    = "drop"
  )
cat("\nDepth by treatment:\n")
print(depth_by_treatment)

# Identify outliers 
outlier_threshold <- mean(seq_depth$log_depth) + 3*sd(seq_depth$log_depth)
outlier_samples   <- seq_depth %>%
  filter(log_depth > outlier_threshold) %>%
  pull(Sample_ID)
cat("\nOutlier samples (>3 SD):", outlier_samples, "\n")

meta_maaslin <- metadata %>%
  left_join(seq_depth %>% select(Sample_ID, log_depth),
            by = "Sample_ID") %>%
  column_to_rownames("Sample_ID")

cat("\nMetadata preview:\n")
print(head(meta_maaslin))

# KEGG pathway 
kegg_names <- keggList("pathway") %>%
  as.data.frame() %>%
  rownames_to_column("pathway_id") %>%
  rename(pathway_name = ".") %>%
  mutate(
    pathway_id   = gsub("path:", "", pathway_id),
    pathway_name = gsub(" - .*", "", pathway_name)   # remove organism suffix
  )

cat("KEGG pathways fetched:", nrow(kegg_names), "\n")

#KEGG abundance matrix 
global_maps <- c("map01100","map01110","map01120","map01130")

# Keep only maps starting with 000 (metabolism) and 020 (environmental info processing)
kegg_matrix <- egg %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(KEGG_Pathway = trimws(KEGG_Pathway)) %>%
  filter(grepl("^map", KEGG_Pathway)) %>%
  filter(!KEGG_Pathway %in% global_maps) %>%
  filter(grepl("^map0[0-2]", KEGG_Pathway)) %>%   
  filter(KEGG_Pathway %in% kegg_names$pathway_id)  %>%
  pivot_longer(cols      = all_of(sample_cols),
               names_to  = "SampleID",
               values_to = "Abundance") %>%
  group_by(KEGG_Pathway, SampleID) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from  = SampleID,
              values_from = Abundance,
              values_fill = 0)

cat("KEGG matrix dimensions (pathways x samples):",
    nrow(kegg_matrix), "x", ncol(kegg_matrix)-1, "\n")


kegg_mat <- kegg_matrix %>%
  column_to_rownames("KEGG_Pathway") %>%
  t() %>%
  as.data.frame()

cat("kegg_mat dimensions (samples x pathways):",
    nrow(kegg_mat), "x", ncol(kegg_mat), "\n")

# MaAsLin2 
run_maaslin <- function(mat, meta, output_dir, fixed_fx) {
  Maaslin2(
    input_data     = mat,
    input_metadata = meta,
    output         = output_dir,
    fixed_effects  = fixed_fx,
    reference      = c("Treatment,CL", "Location,L1"),  
    normalization  = "TSS",
    transform      = "LOG",
    min_prevalence = 0.2,
    min_abundance  = 0.0001,
    cores          = 4
  )
}
#Treatment comparison
# CL vs DD
meta_dd <- meta_maaslin %>%
  filter(Treatment %in% c("DD","CL")) %>%
  mutate(Treatment = relevel(factor(Treatment), ref = "CL"))
mat_dd <- kegg_mat[rownames(meta_dd), ]

# CL vs MC
meta_mc <- meta_maaslin %>%
  filter(Treatment %in% c("CL","MC")) %>%
  mutate(Treatment = relevel(factor(Treatment), ref = "CL"))
mat_mc <- kegg_mat[rownames(meta_mc), ]

# Without outliers
meta_dd_noout <- meta_dd %>% filter(!rownames(.) %in% outlier_samples)
mat_dd_noout  <- kegg_mat[rownames(meta_dd_noout), ]

meta_mc_noout <- meta_mc %>% filter(!rownames(.) %in% outlier_samples)
mat_mc_noout  <- kegg_mat[rownames(meta_mc_noout), ]

cat("\nRunning MaAsLin2...\n")


# CL vs DD 
fit_dd_v2 <- run_maaslin(mat_dd, meta_dd,
                         "maaslin_dd_v2_withdepth",
                         c("Treatment","Location","log_depth"))

# CL vs DD — without outliers
fit_dd_v3 <- run_maaslin(mat_dd_noout, meta_dd_noout,
                         "maaslin_dd_v3_nooutliers",
                         c("Treatment","Location","log_depth"))

# CL vs MC 
fit_mc_v2 <- run_maaslin(mat_mc, meta_mc,
                         "maaslin_mc_v2_withdepth",
                         c("Treatment","Location","log_depth"))

# CL vs MC — without outliers
fit_mc_v3 <- run_maaslin(mat_mc_noout, meta_mc_noout,
                         "maaslin_mc_v3_nooutliers",
                         c("Treatment","Location","log_depth"))

# ── STEP 11: Compare models with and without depth ───────────────────────────
compare_models <- function(fit_v1, fit_v2, label) {
  v1 <- fit_v1$results %>%
    filter(metadata == "Treatment") %>%
    select(feature, coef_v1 = coef, qval_v1 = qval) %>%
    mutate(sig_v1 = qval_v1 < 0.05)
  
  v2 <- fit_v2$results %>%
    filter(metadata == "Treatment") %>%
    select(feature, coef_v2 = coef, qval_v2 = qval) %>%
    mutate(sig_v2 = qval_v2 < 0.05)
  
  combined <- left_join(v1, v2, by = "feature") %>%
    mutate(status = case_when(
      sig_v1 &  sig_v2 ~ "Robust",
      sig_v1 & !sig_v2 ~ "Lost with depth correction",
      !sig_v1 & sig_v2 ~ "Gained with depth correction",
      TRUE             ~ "NS in both"
    ))
  
  cat("\n===", label, "===\n")
  print(table(combined$status))
  return(combined)
}

comp_dd <- compare_models(fit_dd_v1, fit_dd_v2, "DD vs CL")
comp_mc <- compare_models(fit_mc_v1, fit_mc_v2, "MC vs CL")

# Save comparison tables
write.csv(comp_dd, "model_comparison_dd.csv", row.names = FALSE)
write.csv(comp_mc, "model_comparison_mc.csv", row.names = FALSE)

# ── STEP 12: Volcano plot function ───────────────────────────────────────────
make_volcano <- function(maaslin_fit, title) {
  
  res <- maaslin_fit$results %>%
    filter(metadata == "Treatment") %>%
    left_join(kegg_names, by = c("feature" = "pathway_id")) %>%
    mutate(
      pathway_name = ifelse(is.na(pathway_name), feature, pathway_name),
      Significant  = qval < 0.05,
      Direction    = case_when(
        qval < 0.05 & coef > 0 ~ "Enriched",
        qval < 0.05 & coef < 0 ~ "Depleted",
        TRUE                   ~ "NS"
      )
    )
  
  # Count significant features
  n_enriched <- sum(res$Direction == "Enriched")
  n_depleted <- sum(res$Direction == "Depleted")
  cat(title, "— Enriched:", n_enriched, "| Depleted:", n_depleted, "\n")
  
  # Symmetric x-axis
  x_lim <- max(abs(res$coef), na.rm = TRUE)
  y_max <- max(-log10(res$qval), na.rm = TRUE)
  
  ggplot(res, aes(x = coef,
                  y = -log10(qval),
                  color = Direction)) +
    
    geom_point(size = 5, alpha = 1) +
    
    # Label only significant (same as your logic, just cleaner repel)
    # geom_text_repel(
    #   data = res %>% filter(Significant),
    #   aes(label = pathway_name),
    #   size = 3,
    #   box.padding = 0.4,
    #   max.overlaps = 20,
    #   segment.color = "grey70",
    #   show.legend = FALSE
    # ) +
    # 
    # Threshold lines (unchanged meaning)
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "grey40",
               linewidth = 0.5) +
    
    geom_vline(xintercept = 0,
               color = "grey75",
               linewidth = 0.4) +
    
    # Balanced axes
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    
    # Cleaner annotation placement
    # annotate("text",
    #          x = x_lim * 0.95,
    #          y = y_max * 0.95,
    #          label = paste0("Enriched: ", n_enriched),
    #          hjust = 1,
    #          color = "#d6604d",
    #          size = 3.5) +
    # 
    # annotate("text",
    #          x = -x_lim * 0.95,
    #          y = y_max * 0.95,
    #          label = paste0("Depleted: ", n_depleted),
    #          hjust = 0,
    #          color = "#2166ac",
    #          size = 3.5) +
    
    scale_color_manual(
      values = c("Enriched" = "#d6604d",
                 "Depleted" = "#2166ac",
                 "NS"       = "#bdbdbd"),
      name = NULL
    ) +
    
    labs(
      x = "Coefficient (effect size)",
      y = "-Log10 Adjusted P-value"
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none" ,
     
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
      axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
      axis.text.y = element_text(size = 18, colour = "black", face = "bold")
    )
}

# Primary plots (with depth correction — V2)
p_dd <- make_volcano(fit_dd_v2, "DD vs CL")
p_mc <- make_volcano(fit_mc_v2, "MC vs CL")

# Sensitivity plots (without outliers — V3)
p_dd_sens <- make_volcano(fit_dd_v3, "DD vs CL (outliers removed)")
p_mc_sens <- make_volcano(fit_mc_v3, "MC vs CL (outliers removed)")

# Save primary
ggsave("volcano_dd_vs_cl2.png",   p_dd,   width = 7,  height = 7, dpi = 300)
ggsave("volcano_mc_vs_cl2.png",   p_mc,   width = 7,  height = 7, dpi = 300)

# Save sensitivity
ggsave("volcano_dd_vs_cl_sensitivity2.png", p_dd_sens, width = 7,  height = 7, dpi = 300)
ggsave("volcano_mc_vs_cl_sensitivity2.png", p_mc_sens, width = 7,  height = 7, dpi = 300)


# MaAsLin2 
run_maaslin <- function(mat, meta, output_dir, fixed_fx) {
  Maaslin2(
    input_data     = mat,
    input_metadata = meta,
    output         = output_dir,
    fixed_effects  = fixed_fx,
    reference      = c("Treatment,DD", "Location,L1"),  # added Location reference
    normalization  = "TSS",
    transform      = "LOG",
    min_prevalence = 0.2,
    min_abundance  = 0.0001,
    cores          = 4
  )
}
#Treatment comparison

# CL vs MC
meta_mc <- meta_maaslin %>%
  filter(Treatment %in% c("DD","MC")) %>%
  mutate(Treatment = relevel(factor(Treatment), ref = "DD"))
mat_mc <- kegg_mat[rownames(meta_mc), ]

# Without outliers

meta_mc_noout <- meta_mc %>% filter(!rownames(.) %in% outlier_samples)
mat_mc_noout  <- kegg_mat[rownames(meta_mc_noout), ]

## DD vs MC 
fit_mc_v2 <- run_maaslin(mat_mc, meta_mc,
                         "maaslin_ddcl_v2_withdepth",
                         c("Treatment","Location","log_depth"))

# DD vs MC  without outliers
fit_mc_v3 <- run_maaslin(mat_mc_noout, meta_mc_noout,
                         "maaslin_ddcl_v3_nooutliers",
                         c("Treatment","Location","log_depth"))

compare_models <- function(fit_v1, fit_v2, label) {
  v1 <- fit_v1$results %>%
    filter(metadata == "Treatment") %>%
    select(feature, coef_v1 = coef, qval_v1 = qval) %>%
    mutate(sig_v1 = qval_v1 < 0.05)
  
  v2 <- fit_v2$results %>%
    filter(metadata == "Treatment") %>%
    select(feature, coef_v2 = coef, qval_v2 = qval) %>%
    mutate(sig_v2 = qval_v2 < 0.05)
  
  combined <- left_join(v1, v2, by = "feature") %>%
    mutate(status = case_when(
      sig_v1 &  sig_v2 ~ "Robust",
      sig_v1 & !sig_v2 ~ "Lost with depth correction",
      !sig_v1 & sig_v2 ~ "Gained with depth correction",
      TRUE             ~ "NS in both"
    ))
  
  cat("\n===", label, "===\n")
  print(table(combined$status))
  return(combined)
}

comp_mc <- compare_models(fit_mc_v1, fit_mc_v2, "DD vs MC")


write.csv(comp_mc, "model_comparison_ddmc.csv", row.names = FALSE)

make_volcano <- function(maaslin_fit, title) {
  
  res <- maaslin_fit$results %>%
    filter(metadata == "Treatment") %>%
    left_join(kegg_names, by = c("feature" = "pathway_id")) %>%
    mutate(
      pathway_name = ifelse(is.na(pathway_name), feature, pathway_name),
      Significant  = qval < 0.05,
      Direction    = case_when(
        qval < 0.05 & coef > 0 ~ "Enriched",
        qval < 0.05 & coef < 0 ~ "Depleted",
        TRUE                   ~ "NS"
      )
    )
  
  
  n_enriched <- sum(res$Direction == "Enriched")
  n_depleted <- sum(res$Direction == "Depleted")
  cat(title, "— Enriched:", n_enriched, "| Depleted:", n_depleted, "\n")
  
  x_lim <- max(abs(res$coef), na.rm = TRUE)
  y_max <- max(-log10(res$qval), na.rm = TRUE)
  
  ggplot(res, aes(x = coef,
                  y = -log10(qval),
                  color = Direction)) +
    
    geom_point(size = 5, alpha = 1) +
    
    # geom_text_repel(
    #   data = res %>% filter(Significant),
    #   aes(label = pathway_name),
    #   size = 3,
    #   box.padding = 0.4,
    #   max.overlaps = 20,
    #   segment.color = "grey70",
    #   show.legend = FALSE
    # ) +
    # 
 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "grey40",
               linewidth = 0.5) +
    
    geom_vline(xintercept = 0,
               color = "grey75",
               linewidth = 0.4) +
    
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    
    # annotate("text",
    #          x = x_lim * 0.95,
    #          y = y_max * 0.95,
    #          label = paste0("Enriched: ", n_enriched),
    #          hjust = 1,
    #          color = "#d6604d",
    #          size = 3.5) +
    # 
    # annotate("text",
    #          x = -x_lim * 0.95,
    #          y = y_max * 0.95,
    #          label = paste0("Depleted: ", n_depleted),
    #          hjust = 0,
    #          color = "#2166ac",
    #          size = 3.5) +
    
    scale_color_manual(
      values = c("Enriched" = "#d6604d",
                 "Depleted" = "#2166ac",
                 "NS"       = "#bdbdbd"),
      name = NULL
    ) +
    
    labs(
      x = "Coefficient (effect size)",
      y = "-Log10 Adjusted P-value"
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none" ,
      
      axis.line = element_line(color = "black"),
      axis.title.x = element_text(size = 18, colour = "black", face = "bold"),
      axis.title.y = element_text(size = 18, colour = "black", face = "bold"),
      axis.text.x = element_text(size = 18, colour = "black", face = "bold"),
      axis.text.y = element_text(size = 18, colour = "black", face = "bold")
    )
}

p_mc <- make_volcano(fit_mc_v2, "DD vs MC")


p_mc_sens <- make_volcano(fit_mc_v3, "DD vs MC (outliers removed)")

ggsave("volcano_dd_vs_mc2.png",   p_mc,   width = 7,  height = 7, dpi = 300)

ggsave("volcano_dd_vs_mc_sensitivity2.png", p_mc_sens, width = 7,  height = 7, dpi = 300)


cat("\nDone. Check your working directory for output files.\n")