#================================================

local_lib <- Sys.getenv("R_LIBS_USER")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(NetCoMi)
library(phyloseq)
library(readxl)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

TREATMENT <- "Control" 


BASE_DIR <- Sys.getenv("ANALYSIS_DIR",
                       unset = "D:/SKM/metagenome/dwn_analysis/")
setwd(BASE_DIR)

OUT <- paste0("networks/", TREATMENT, "/")
dir.create(OUT,            showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(OUT, "Cytoscape"), showWarnings = FALSE)
dir.create(paste0(OUT, "plots"),     showWarnings = FALSE)
dir.create(paste0(OUT, "rds"),       showWarnings = FALSE)

cat("\n======================================\n")
cat(" Treatment:", TREATMENT, "\n")
cat(" Output dir:", OUT, "\n")
cat("======================================\n\n")

#=========================
#LOADDATA
#=========================
raw_data <- read_excel("Abundance_table.csv")
raw_data$`L5-28` <- NULL
raw_data <- as.data.frame(raw_data)

# Taxonomy
tax_df <- raw_data[, 2:8]
rownames(tax_df) <- raw_data[[1]]
colnames(tax_df) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# OTU matrix (ASVs x samples)
otu_mat <- as.matrix(raw_data[, 9:ncol(raw_data)])
rownames(otu_mat) <- raw_data[[1]]

# Metadata
metadata <- data.frame(
  Sample_ID = c("L1-1","L1-2","L1-3","L1-4","L1-5","L1-6","L1-7","L1-8","L1-9",
                "L3-10","L3-11","L3-12","L3-13","L3-14","L3-15","L3-16","L3-17","L3-18",
                "L5-19","L5-20","L5-21","L5-22","L5-23","L5-24","L5-25","L5-26","L5-27"),
  Location  = c(rep("L1",9), rep("L3",9), rep("L5",9)),
  Treatment = rep(c(rep("Control",3), rep("CS",3), rep("CS with consortium",3)), 3),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$Sample_ID

#=========================
# PHYLOSEQ OBJECT
#=========================
ps_full <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data(metadata),
  tax_table(as.matrix(tax_df))
)

#=========================
# FILTERING
#=========================
total_reads    <- sum(otu_table(ps_full))
keep_abund     <- taxa_sums(ps_full) / total_reads > 0.0001
keep_prev      <- rowSums(otu_table(ps_full) > 0) / nsamples(ps_full) >= 0.20
ps_filt        <- prune_taxa(keep_abund & keep_prev, ps_full)
cat("After filtering:", ntaxa(ps_filt), "taxa retained\n")

#=========================
# TREATMENT
#=========================
ps_trt <- subset_samples(ps_filt, Treatment == TREATMENT)
ps_trt <- prune_taxa(taxa_sums(ps_trt) > 0, ps_trt)   # drop zero-sum taxa
cat("Treatment subset:", ntaxa(ps_trt), "taxa |", nsamples(ps_trt), "samples\n")

# Sample order check
cat("Samples:", sample_names(ps_trt), "\n")
cat("Locations:", sample_data(ps_trt)$Location, "\n")

#=========================
# 7. CLR TRANSFORMATION
#=========================
count_mat <- t(as.matrix(otu_table(ps_trt)))   # samples x taxa

# Multiplicative zero replacement 
if (!requireNamespace("zCompositions", quietly = TRUE)) {
  install.packages("zCompositions",
                   lib   = local_lib,
                   repos = "https://cloud.r-project.org")
}
library(zCompositions)

count_mat_nz <- cmultRepl(count_mat,
                          label  = 0,
                          method = "CZM",    # Bayesian multiplicative replacement
                          output = "p-counts",
                          suppress.print = TRUE)

clr_mat <- t(apply(count_mat_nz, 1, function(x) {
  log(x) - mean(log(x))
  
  #=========================
  # REGRESS OUT LOCATION EFFECT
  #=========================
  
  location_vec <- factor(sample_data(ps_trt)$Location)
  
  clr_resid <- apply(clr_mat, 2, function(y) {
    fit <- lm(y ~ location_vec)
    resid(fit)
  })
  rownames(clr_resid) <- rownames(clr_mat)
  
  cat("Residuals matrix dim:", dim(clr_resid), "\n")
  
  write.csv(as.data.frame(clr_mat),
            paste0(OUT, "CLR_matrix_", TREATMENT, ".csv"))
  write.csv(as.data.frame(clr_resid),
            paste0(OUT, "CLR_residuals_location_regressed_", TREATMENT, ".csv"))
  
  #=========================
  # NETWORK CONSTRUCTION 
  #=========================
  t_start <- proc.time()
  
  n_cpus <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
  
  cor_obs <- cor(clr_resid, method = "pearson")
  
  cat(sprintf("  Running 999 permutations on %d core(s)...\n", n_cpus))
  set.seed(123)
  n_boot       <- 999
  exceed_count <- matrix(0L,
                         nrow = ncol(clr_resid), ncol = ncol(clr_resid),
                         dimnames = list(colnames(clr_resid), colnames(clr_resid)))
  
  if (n_cpus > 1 && requireNamespace("parallel", quietly = TRUE)) {
    library(parallel)
    cl <- makeCluster(n_cpus)
    clusterExport(cl, c("clr_resid", "cor_obs"), envir = environment())
    counts_list <- parLapply(cl, seq_len(n_boot), function(b) {
      perm_mat <- apply(clr_resid, 2, sample)
      cor_boot <- cor(perm_mat, method = "pearson")
      (abs(cor_boot) >= abs(cor_obs)) * 1L
    })
    stopCluster(cl)
    for (ct in counts_list) exceed_count <- exceed_count + ct
  } else {
    for (b in seq_len(n_boot)) {
      perm_mat     <- apply(clr_resid, 2, sample)
      cor_boot     <- cor(perm_mat, method = "pearson")
      exceed_count <- exceed_count + (abs(cor_boot) >= abs(cor_obs)) * 1L
      if (b %% 200 == 0) cat("  Bootstrap", b, "/ 999\n")
    }
  }
  
  p_matrix       <- (exceed_count + 1) / (n_boot + 1)
  diag(p_matrix) <- 1
  
  write.csv(cor_obs,  paste0(OUT, "Pearson_CLR_correlation_", TREATMENT, ".csv"))
  write.csv(p_matrix, paste0(OUT, "Bootstrap_pvalues_",       TREATMENT, ".csv"))
  
  cor_thresh <- cor_obs
  cor_thresh[abs(cor_obs) <= 0.82 | p_matrix >= 0.01] <- 0
  diag(cor_thresh) <- 0
  cat(sprintf("  Edges retained after thresholding: %d\n",
              sum(cor_thresh[upper.tri(cor_thresh)] != 0)))
  
  net_obj <- netConstruct(
    data        = cor_thresh,      
    dataType    = "correlation",   
    sparsMethod = "none",          
    normMethod  = "none",
    seed        = 123,
    verbose     = TRUE
  )
  
  
  #=========================
  #NETWORK ANALYSIS
  #=========================
  ana_obj <- netAnalyze(
    net_obj,
    clustMethod  = "cluster_fast_greedy",
    hubPar       = "degree",           
    hubQuant     = 0.95,
    lnormFit     = TRUE,
    normDeg      = TRUE,
    normBetw     = TRUE,
    normClose    = TRUE,
    normEigen    = TRUE
  )
  saveRDS(ana_obj, paste0(OUT, "rds/ana_", TREATMENT, ".rds"))
  
  #=========================
  # TOPOLOGY TABLE
  #=========================
  lcc  <- ana_obj$globalPropsLCC
  adj  <- ana_obj$input$adjaMat1
  asso <- ana_obj$input$assoEst1
  
  adj_upper  <- adj[upper.tri(adj)]
  asso_upper <- asso[upper.tri(asso)]
  
  topology_df <- data.frame(
    Treatment              = TREATMENT,
    Nodes                  = sum(rowSums(adj) > 0),
    Edges                  = sum(adj_upper != 0),
    Positive_edges         = sum(adj_upper != 0 & asso_upper > 0),
    Negative_edges         = sum(adj_upper != 0 & asso_upper < 0),
    LCC_size               = lcc$lccSize1,
    LCC_size_relative      = round(lcc$lccSizeRel1, 4),
    Average_degree         = round(mean(rowSums(adj != 0)), 4),
    Average_path_length    = round(lcc$avPath1, 4),
    Clustering_coefficient = round(lcc$clustCoef1, 4),
    Modularity             = round(lcc$modularity1, 4),
    Graph_density          = round(lcc$density1, 6),
    Natural_connectivity   = round(lcc$natConnect1, 6),
    Positive_edge_percent  = round(lcc$pep1, 2)
  )
  print(topology_df)
  write.csv(topology_df, paste0(OUT, "Topology_", TREATMENT, ".csv"), row.names = FALSE)
  
  #=========================
  # 12. Zi-Pi & HUBS BY DEGREE
  #=========================
  cent     <- ana_obj$centralities
  degree_v <- cent$degree1
  betw_v   <- cent$between1
  close_v  <- cent$close1
  eigenv_v <- cent$eigenv1
  
  clust    <- ana_obj$clusteringLCC$clust1
  all_asvs <- names(degree_v)
  raw_deg  <- rowSums(adj != 0)
  
  mod_all  <- rep(0L, length(all_asvs))
  names(mod_all) <- all_asvs
  mod_all[names(clust)] <- clust
  
  
  zi <- sapply(all_taxid, function(taxid) {
    m <- mod_all[taxid]
    if (m == 0) return(NA_real_)
    peers  <- names(mod_all)[mod_all == m]
    k_i    <- raw_deg[taxid]
    k_mean <- mean(raw_deg[peers], na.rm = TRUE)
    k_sd   <- sd(raw_deg[peers],   na.rm = TRUE)
    if (is.na(k_sd) || k_sd == 0) return(0)
    (k_i - k_mean) / k_sd
  })
  
  # Pi (participation coefficient)
  pi_val <- sapply(all_taxid, function(taxid) {
    k_i <- raw_deg[asv]
    if (k_i == 0) return(NA_real_)
    neighbours <- names(which(adj[taxid, ] != 0))
    if (length(neighbours) == 0) return(NA_real_)
    k_is_vec <- table(mod_all[neighbours])
    1 - sum((k_is_vec / k_i)^2)
  })
  
  # Node roles
  node_role <- case_when(
    is.na(zi) | is.na(pi_val)  ~ "Peripheral",
    zi >= 2.5 & pi_val >= 0.62 ~ "Network Hub",
    zi >= 2.5 & pi_val <  0.62 ~ "Module Hub",
    zi <  2.5 & pi_val >= 0.62 ~ "Connector",
    TRUE                        ~ "Peripheral"
  )
  
  # Full node table with taxonomy
  tax_micro <- as.data.frame(tax_table(ps_trt))
  node_df <- data.frame(
    Treatment   = TREATMENT,
    taxid         = all_taxid,
    Zi          = round(zi,     4),
    Pi          = round(pi_val, 4),
    Node_role   = node_role,
    Degree      = round(degree_v[all_asvs], 4),
    Betweenness = round(betw_v[all_asvs],   6),
    Closeness   = round(close_v[all_asvs],  6),
    Eigenvector = round(eigenv_v[all_asvs], 6),
    Module      = mod_all[all_asvs],
    row.names   = NULL
  ) %>%
    left_join(tax_micro %>% mutate(taxid = rownames(tax_micro)), by = "taxid") %>%
    dplyr::select(Treatment, taxid,
                  Kingdom, Phylum, Class, Order, Family, Genus, Species,
                  Degree, Betweenness, Closeness, Eigenvector,
                  Module, Zi, Pi, Node_role)
  
  write.csv(node_df,
            paste0(OUT, "Node_Taxonomy_ZiPi_", TREATMENT, ".csv"),
            row.names = FALSE)
  
  # TOP 20 HUBS BY DEGREE 
  top20_hubs <- node_df %>%
    filter(!is.na(Degree)) %>%
    arrange(desc(Degree), desc(Betweenness)) %>%
    slice_head(n = 20) %>%
    mutate(Hub_rank = row_number())
  
  cat("\n=== Top 20 Hub Microbes (", TREATMENT, ") ===\n")
  print(top20_hubs %>% dplyr::select(Hub_rank, taxid, Phylum, Genus, Degree, Betweenness, Node_role))
  write.csv(top20_hubs,
            paste0(OUT, "Top20_Hubs_", TREATMENT, ".csv"),
            row.names = FALSE)
  
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  
  node_names <- V(g)$name
  nd_sub     <- node_df[match(node_names, node_df$ASV), ]
  
  V(g)$Phylum    <- nd_sub$Phylum
  V(g)$Genus     <- nd_sub$Genus
  V(g)$Degree    <- nd_sub$Degree
  V(g)$Zi        <- nd_sub$Zi
  V(g)$Pi        <- nd_sub$Pi
  V(g)$Node_role <- nd_sub$Node_role
  V(g)$Module    <- nd_sub$Module
  V(g)$Top20_hub <- node_names %in% top20_hubs$ASV
  E(g)$sign      <- ifelse(E(g)$weight > 0, "positive", "negative")
  
  write_graph(g,
              paste0(OUT, "Cytoscape/", TREATMENT, "_network.graphml"),
              format = "graphml")
  
  el <- as_edgelist(g)
  write.csv(data.frame(source = el[,1], target = el[,2],
                       weight = E(g)$weight, sign = E(g)$sign),
            paste0(OUT, "Cytoscape/", TREATMENT, "_edges.csv"),
            row.names = FALSE)
  write.csv(nd_sub,
            paste0(OUT, "Cytoscape/", TREATMENT, "_nodes.csv"),
            row.names = FALSE)
  
  t_total <- proc.time() - t_start
  cat(sprintf("\n====== %s COMPLETE in %.1f minutes ======\n",
              TREATMENT, t_total["elapsed"] / 60))
}
  
  
  