# =============================================================================
# sra_metadata_filter.R
#
# Purpose : Merge SRA accessions from cleaned JSON objects with downloaded
#           metadata CSVs, remove metagenome hits, and normalise organism
#           names for downstream analysis.
#
# Depends : sra_mapping.R must be sourced first so that `rho_json_clean` and
#           `har1a_json_clean` exist in the session.
#
# Inputs  : rho_json_clean    (list, from sra_mapping.R)
#           har1a_json_clean  (list, from sra_mapping.R)
#           final_rho_table_clean.csv
#           final_har_table_clean.csv
#
# Outputs : rho_filtered  (data frame, in session)
#           har1a_filtered (data frame, in session)
#           rho_filtered.csv
#           har1a_filtered.csv
#
# Requires: dplyr
# =============================================================================

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

main <- function() {
  
  # Guard: confirm upstream objects exist before proceeding
  missing_objs <- setdiff(
    c("rho_json_clean", "har1a_json_clean"),
    ls(envir = .GlobalEnv)
  )
  if (length(missing_objs) > 0L) {
    stop(
      "Required object(s) not found in global environment: ",
      paste(missing_objs, collapse = ", "),
      ".\nPlease source sra_mapping.R before running this script."
    )
  }
  
  # -----------------------------------------------------------------------------
  
  rho_json_clean <-  RHO_INPUT_FILE[!(names(RHO_INPUT_FILE) %in% RHO_SPECIES_DROP)]
  har1a_json_clean <- HAR1A_INPUT_FILE[!(names(HAR1A_INPUT_FILE) %in% HAR1A_VARS_DROP)]
  
  
  # ---- RHO: merge metadata and filter metagenomes -------------------------
  
  rho_filtered <- load_and_filter_meta(rho_json_clean, RHO_META_CSV )
  
  cat(sprintf("RHO   : %d rows after merging and dropping metagenomes\n",
              nrow(rho_filtered)))
  
  
  # ---- HAR1A: merge metadata and filter metagenomes -----------------------
  
  har1a_filtered <- load_and_filter_meta(har1a_json_clean, HAR1A_META_CSV)
  
  cat(sprintf("HAR1A : %d rows after merging and dropping metagenomes\n",
              nrow(har1a_filtered)))
  

  
  # ---- HAR1A: normalise organism names ------------------------------------
  har1a_filtered <- har1a_filtered %>%
    filter(
      !is.na(organism),
      trimws(organism) != "",
      #Total_SRAs_in_DB >= MIN_TOTAL_SRA,
      !grepl(paste(exclude_patterns, collapse = "|"), organism, ignore.case = TRUE)
    ) %>%
    # Add new rules in the case_when block below.
    # Exact-match rules (==) take precedence; partial-match rules (grepl)
    # follow. The final TRUE clause preserves any name not otherwise matched.
    mutate(organism = apply_name_mapping(organism, exact_map, pattern_map))
  
  cat(sprintf("HAR1A : %d rows after organism name normalisation\n",
              nrow(har1a_filtered)))
  
  
  # ---- Export filtered tables ---------------------------------------------
  
  write.csv(rho_filtered,   "rho_filtered.csv",   row.names = FALSE)
  write.csv(har1a_filtered, "har1a_filtered.csv",  row.names = FALSE)
  cat("Exported rho_filtered.csv and har1a_filtered.csv\n")
  
  
  # Assign results to the global environment so downstream scripts can use them
  assign("rho_filtered",   rho_filtered,   envir = .GlobalEnv)
  assign("har1a_filtered", har1a_filtered, envir = .GlobalEnv)
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

main()

