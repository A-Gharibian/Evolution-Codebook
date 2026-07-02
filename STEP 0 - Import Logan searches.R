# =============================================================================
# sra_mapping.R
#
# Purpose : Load per-gene JSON files, extract unique SRA accessions, build a
#           master mapping table, export results, and produce QC plots.
#
# Inputs  : RHO.json, HAR1A.json  (in the working directory)
# Outputs : rho_target_sras.txt, har_target_sras.txt, target_sras.txt,
#           SRA_Mapping_Table.csv, RHO_cleaned_names.json,
#           HAR1A_cleaned_names.json, plot_rho_hits.png, plot_har1a_hits.png
#
# Requires: jsonlite, ggplot2, scales
# =============================================================================

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

main <- function() {
  
  # ---- RHO: load, clean, extract ------------------------------------------
  
  rho_json <- tryCatch(
    read_json(RHO_INPUT_FILE),
    error = function(e) stop("Could not read '", RHO_INPUT_FILE, "': ", e$message)
  )
  
  rho_json_clean <- rho_json[!(names(rho_json) %in% RHO_SPECIES_DROP)]
  
  # Rename non-standard keys to human-readable species names
  names(rho_json_clean)[names(rho_json_clean) == "X."]      <- "Frog"
  names(rho_json_clean)[names(rho_json_clean) == "Painted"] <- "Turtle"
  
  rho_table       <- flatten_sra_data(rho_json_clean, "RHO")
  unique_rho_sras <- unique(unlist(lapply(rho_json_clean, names)))
  
  writeLines(unique_rho_sras, "rho_target_sras.txt")
  cat(sprintf("RHO   : %d unique SRA accessions written to rho_target_sras.txt\n",
              length(unique_rho_sras)))
  
  
  # ---- HAR1A: load, clean, extract -----------------------------------------
  
  har1a_json <- tryCatch(
    read_json(HAR1A_INPUT_FILE),
    error = function(e) stop("Could not read '", HAR1A_INPUT_FILE, "': ", e$message)
  )
  
  har1a_json_clean <- har1a_json[!(names(har1a_json) %in% HAR1A_VARS_DROP)]
  
  har1a_table       <- flatten_sra_data(har1a_json_clean, "HAR1A")
  unique_har1a_sras  <- unique(har1a_table$SRA_Accession)
  
  writeLines(as.character(unique_har1a_sras), "har_target_sras.txt")
  cat(sprintf("HAR1A : %d unique SRA accessions written to har_target_sras.txt\n",
              length(unique_har1a_sras)))
  
  
  # ---- Master mapping: combine and export -----------------------------------
  
  master_mapping <- rbind(rho_table, har1a_table)
  write.csv(master_mapping, "SRA_Mapping_Table.csv", row.names = FALSE)
  
  unique_sras <- unique(master_mapping$SRA_Accession)
  writeLines(as.character(unique_sras), "target_sras.txt")
  cat(sprintf("TOTAL : %d unique SRAs written to target_sras.txt\n",
              length(unique_sras)))
  
  # Export cleaned JSON for downstream use
  write_json(rho_json_clean,   "RHO_cleaned_names.json",   pretty = TRUE)
  write_json(har1a_json_clean, "HAR1A_cleaned_names.json",  pretty = TRUE)
  
  
  # ---- Plot 1: RHO hits by species (log10 scale) ---------------------------
  
  rho_df <- data.frame(
    Species = names(rho_json_clean),
    Hits    = as.numeric(lengths(rho_json_clean))
  )
  
  p_rho <- ggplot(rho_df, aes(x = reorder(Species, -Hits), y = Hits)) +
    geom_col(fill = "steelblue", color = "black") +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(
      title = "RHO: SRA Hits by Species",
      x     = "Species",
      y     = "Number of SRA Hits (log\u2081\u2080)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("plot_rho_hits.png", plot = p_rho, width = 8, height = 5, dpi = 300)
  cat("Plot saved: plot_rho_hits.png\n")
  
  
  # ---- Plot 2: HAR1A hits by variation (log10 scale) -----------------------
  
  har1a_df <- data.frame(
    Variation = names(har1a_json_clean),
    Hits      = as.numeric(lengths(har1a_json_clean))
  )
  
  p_har1a <- ggplot(har1a_df, aes(x = reorder(Variation, -Hits), y = Hits)) +
    geom_col(fill = "darkred", color = "black") +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(
      title = "HAR1A: SRA Hits by Variation",
      x     = "Variation",
      y     = "Number of SRA Hits (log\u2081\u2080)"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave("plot_har1a_hits.png", plot = p_har1a, width = 10, height = 5, dpi = 300)
  cat("Plot saved: plot_har1a_hits.png\n")
  
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

main()

