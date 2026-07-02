# =============================================================================
# Script 4: Heatmap & Phylogenetic Tree Visualization
# =============================================================================
# Description:
#   Builds log10-transformed SRA hit heatmaps for RHO and HAR1A, overlaid on
#   Jaccard-distance-based phylogenetic trees. Three heatmap backends are
#   provided for HAR1A pending a final choice.
#
# Dependencies (must be run before this script):
#   - Script 2 (metadata): provides rho_filtered, har_filtered_clean
#   - Script 3 (clustering): provides clean_rho_hc, build_jaccard_matrix()
#
# Inputs:  rho_filtered, har_filtered_clean, clean_rho_hc, build_jaccard_matrix()
# Outputs: Interactive plots only (no files written)
# =============================================================================
# =============================================================================
# SECTION 1: RHO — Heatmap + Phylogenetic Tree
# =============================================================================

# 1. Convert empty organism strings to NA and remove them
rho_filtered <- rho_filtered %>%
  dplyr::mutate(organism = dplyr::na_if(organism, ""))

# 2. Get the top N most frequent NCBI organisms
top_ncbi_rho <- rho_filtered %>%
  dplyr::filter(!is.na(organism)) %>%
  dplyr::count(organism) %>%
  dplyr::slice_max(n, n = top_n_rho) %>%
  dplyr::pull(organism)

# 3. Build the heatmap matrix (rows = query species, cols = NCBI organisms)
rho_heatmap_data <- rho_filtered %>%
  dplyr::filter(organism %in% top_ncbi_rho) %>%
  dplyr::group_by(Query_Species, organism) %>%
  dplyr::summarise(hits = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = organism, values_from = hits, values_fill = 0)

# 4. Convert to matrix and apply log10 transform (+1 prevents log(0) errors)
rho_heatmap_matrix <- as.data.frame(rho_heatmap_data)
rownames(rho_heatmap_matrix) <- rho_heatmap_matrix$Query_Species
rho_heatmap_matrix$Query_Species <- NULL
rho_heatmap_log <- log10(rho_heatmap_matrix + 1)

# 5. Convert Jaccard clustering object to phylogenetic tree
#    Requires: clean_rho_hc from Script 3
rho_phylo <- as.phylo(clean_rho_hc)

# 6. Draw the base phylogenetic tree
p_rho <- ggtree(rho_phylo) +
  geom_tiplab(size = 4, align = TRUE, linesize = 0.5) +
  hexpand(0.3)

# 7. Overlay the heatmap on the tree
gheatmap(p_rho, rho_heatmap_log,
         offset           = 0.2,
         width            = 1.5,
         colnames         = TRUE,
         colnames_angle   = 45,
         colnames_offset_y = -0.5,
         font.size        = 3) +
  scale_fill_viridis_c(name = "log10(SRA Hits)", option = "plasma") +
  theme(legend.position = "right") +
  labs(
    title    = "RHO: Evolutionary Sequence Clustering vs. SRA Metadata",
    subtitle = "Warmer colors indicate higher abundance in that NCBI dataset"
  )


# =============================================================================
# SECTION 2: HAR1A — Heatmap + Phylogenetic Tree
# =============================================================================
# Requires: har_filtered_clean from Script 2 (additional cleaning step)

# 1. Convert empty organism strings to NA and remove them
har_filtered_clean <- har_filtered_clean %>%
  dplyr::mutate(organism = dplyr::na_if(organism, ""))

# 2. Get the top N most frequent NCBI organisms
top_ncbi_har <- har_filtered_clean %>%
  dplyr::filter(!is.na(organism)) %>%
  dplyr::count(organism, sort = TRUE) %>%
  dplyr::slice_head(n = top_n_har) %>%
  dplyr::pull(organism)

# 3. Build the heatmap matrix (rows = variations, cols = NCBI organisms)
har_heatmap_data <- har_filtered_clean %>%
  dplyr::filter(organism %in% top_ncbi_har) %>%
  dplyr::group_by(Query_Species, organism) %>%
  dplyr::summarise(hits = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = organism, values_from = hits, values_fill = 0)

# 4. Convert to matrix and apply log10 transform (+1 prevents log(0) errors)
har_heatmap_matrix <- as.data.frame(har_heatmap_data)
rownames(har_heatmap_matrix) <- har_heatmap_matrix$Query_Species
har_heatmap_matrix$Query_Species <- NULL
har_heatmap_log <- log10(har_heatmap_matrix + 1)

# 5. Build Jaccard distance matrix using biosample as the unit of comparison
#    Requires: build_jaccard_matrix() from Script 3
har_split <- split(har_filtered_clean$biosample, har_filtered_clean$Query_Species)
har_ready_for_tree <- lapply(har_split, function(biosample_vector) {
  setNames(rep(1, length(biosample_vector)), biosample_vector)
})

clean_har_dist <- build_jaccard_matrix(har_ready_for_tree)
clean_har_hc   <- hclust(clean_har_dist, method = "average")

# 6. Convert to phylogenetic tree object
har_phylo <- as.phylo(clean_har_hc)

# 7. Draw the base phylogenetic tree
p_har <- ggtree(har_phylo) +
  geom_tiplab(size = 3.5, align = TRUE, linesize = 0.5) +
  hexpand(0.4)

# 8. Overlay the heatmap on the tree
gheatmap(p_har, har_heatmap_log,
         offset            = 0.05,
         width             = 2,
         colnames          = TRUE,
         colnames_angle    = 90,
         colnames_offset_y = -0.5,
         font.size         = 2.5) +
  scale_fill_viridis_c(name = "log10(SRA Hits)", option = "plasma") +
  theme(legend.position = "right") +
  labs(title = "HAR1A: Evolutionary Sequence Clustering vs. SRA Metadata")


# =============================================================================
# SECTION 3: HAR1A — Alternative Heatmap Backends (CHOOSE ONE)
# =============================================================================
# TODO: Select one of the three implementations below and remove the others.
#       All three use the same Jaccard clustering (clean_har_hc) for rows and
#       Euclidean clustering for columns.
# =============================================================================

# --- Option A: ggtree + gheatmap (already shown above in Section 2) ----------
# (Retained above — best for combining with phylogenetic tree layout)

# --- Option B: ComplexHeatmap -------------------------------------------------
Heatmap(
  as.matrix(har_heatmap_log),
  name             = "log10(Hits)",
  cluster_rows     = clean_har_hc,   # Jaccard-based clustering
  cluster_columns  = TRUE,           # Euclidean clustering for NCBI organisms
  col              = viridis(100, option = "plasma"),
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 10),
  column_title     = "HAR1A: Biosample Co-occurrence (Jaccard) vs. SRA Metadata"
)

# --- Option C: pheatmap -------------------------------------------------------
pheatmap(
  mat          = har_heatmap_log,
  cluster_rows = clean_har_hc,       # Jaccard-based clustering
  cluster_cols = TRUE,               # Set FALSE to preserve original column order
  color        = viridis(100, option = "plasma"),
  fontsize_row = 8,
  fontsize_col = 8,
  main         = "HAR1A: Biosample Co-occurrence (Jaccard) vs. SRA Metadata"
)

