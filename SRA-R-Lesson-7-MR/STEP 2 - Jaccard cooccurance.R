# =============================================================================
# sra_jaccard_clustering.R
#
# Purpose : Compute pairwise Jaccard distances between groups based on shared
#           biosample accessions, cluster via UPGMA, and plot dendrograms for
#           RHO (by species) and HAR1A (by variation).
#
# Depends : sra_mapping.R and sra_metadata_filter.R must be sourced first so
#           that `rho_json_clean`, `har1a_json_clean`, `rho_filtered`, and
#           `har1a_filtered` exist in the session.
#
# Inputs  : rho_json_clean    (list, from sra_mapping.R)
#           har1a_json_clean  (list, from sra_mapping.R)
#           rho_filtered      (data frame, from sra_metadata_filter.R)
#           har1a_filtered    (data frame, from sra_metadata_filter.R)
#
# Outputs : plot_rho_dendrogram.png
#           plot_har1a_dendrogram.png
#
# Requires: dplyr, pheatmap, viridis
# =============================================================================

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

main <- function() {
  
  # Guard: confirm upstream objects exist before proceeding
  required <- c("rho_json_clean", "har1a_json_clean",
                "rho_filtered",   "har1a_filtered")
  missing_objs <- setdiff(required, ls(envir = .GlobalEnv))
  if (length(missing_objs) > 0L) {
    stop(
      "Required object(s) not found in global environment: ",
      paste(missing_objs, collapse = ", "),
      ".\nPlease source sra_mapping.R and sra_metadata_filter.R first."
    )
  }
  
  
  # ---- RHO: Jaccard clustering by species ----------------------------------
  
  rho_input    <- prepare_jaccard_input(rho_filtered, group_col = "Query_Species")
  rho_dist     <- build_jaccard_matrix(rho_input)
  rho_hc       <- hclust(rho_dist, method = "average")
  
  save_dendrogram_png(
    hc       = rho_hc,
    file     = "plot_rho_dendrogram.png",
    title    = "RHO Biosample Co-occurrence (Metagenomes Removed)",
    use_dend = TRUE,
    mar      = c(8, 4, 4, 2)
  )
  cat("Plot saved: plot_rho_dendrogram.png\n")
  
  
  # ---- HAR1A: filter to genomic assays and known vertebrates ---------------
  
  har1a_filtered_clean <- har1a_filtered |>
    filter(assay_type %in% GENOMIC_ASSAYS) |>
    filter(!organism  %in% EXCLUDE_ORGANISMS)
  
  cat(sprintf(
    "HAR1A : %d rows retained after genomic-assay filter and organism exclusions\n",
    nrow(har1a_filtered_clean)
  ))
  
  
  # ---- HAR1A: Jaccard clustering by variation ------------------------------
  
  har1a_input  <- prepare_jaccard_input(har1a_filtered_clean, group_col = "Query_Species")
  har1a_dist   <- build_jaccard_matrix(har1a_input)
  har1a_hc     <- hclust(har1a_dist, method = "average")
  
  save_dendrogram_png(
    hc       = har1a_hc,
    file     = "plot_har1a_dendrogram.png",
    title    = "HAR1A Biosample Co-occurrence Clustering (by Variation)",
    use_dend = FALSE,
    mar      = c(10, 4, 4, 2),
    cex      = 0.7
  )
  cat("Plot saved: plot_har1a_dendrogram.png\n")
  
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

main()
