# =============================================================================
# Script 8: Tajima's D Analysis — HAR1A
# =============================================================================
# Description:
#   Computes population genetic statistics (S, pi, Tajima's D) per organism
#   from aligned HAR1A sequences, merges results with vertebrate life history
#   traits, and plots Tajima's D against effective population size (Ne).
#
# Dependencies (must be run before this script):
#   - Script X (not yet finalized): provides har_mapped_data_clean
#     (contains aligned sequence columns pos_1 to pos_31 and organism metadata)
#
# Inputs:
#   - har_mapped_data_clean  : data frame with aligned HAR1A sequences
#   - Vertebrate mutation rates.csv : life history traits including Ne
#
# Outputs:
#   - Output/har1a_combined_dataset.csv  : Tajima's D merged with life history
#   - Output/har1a_stats_with_D.csv      : Tajima's D stats only
#   - Interactive plot: Tajima's D vs. Ne
# =============================================================================

# =============================================================================
# SECTION 2: Compute Per-Organism Statistics
# =============================================================================
# RHO


# HAR1A

har1a_stats_with_D <- har_mapped_data_clean %>%
  dplyr::filter(organism != "", !is.na(organism)) %>%
  dplyr::group_by(organism) %>%
  dplyr::summarise(
    # Sample size
    n = dplyr::n(),
    
    # Number of segregating sites (S): positions with more than one unique base
    S = sum(sapply(pick(all_of(seq_cols)), function(col) {
      length(unique(col[!is.na(col)])) > 1
    })),
    
    # Average pairwise diversity (pi): sum of per-site heterozygosity
    pi_val = sum(sapply(pick(all_of(seq_cols)), function(col) {
      valid <- col[!is.na(col)]
      if (length(valid) > 1) {
        freqs         <- table(valid) / length(valid)
        heterozygosity <- 1 - sum(freqs^2)
        return(heterozygosity * (length(valid) / (length(valid) - 1)))
      } else {
        return(0)
      }
    })),
    
    .groups = "drop"
  ) %>%
  
  # Apply minimum sample size and segregating site filters
  dplyr::filter(n >= min_n, S > 0) %>%
  
  # Compute Watterson's theta and Tajima's D per row
  dplyr::rowwise() %>%
  dplyr::mutate(
    a_n       = sum(1 / (1:(n - 1))),
    theta_W   = S / a_n,
    Tajimas_D = calc_tajimas_d(n, S, pi_val)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(organism)


# =============================================================================
# SECTION 3: Merge with Vertebrate Life History Traits
# =============================================================================

Vert_mut_rates <- read.csv("Vertebrate mutation rates.csv")

# Join on organism name — check CSV 'Species.name' column matches SRA organism
har1a_combined_dataset <- har1a_stats_with_D %>%
  dplyr::inner_join(Vert_mut_rates, by = c("organism" = "Species.name"))

print(paste("Successfully merged species:", nrow(har1a_combined_dataset)))


# =============================================================================
# SECTION 4: Plot Tajima's D vs. Effective Population Size
# =============================================================================

ggplot(har1a_combined_dataset, aes(x = Ne, y = Tajimas_D, label = Common.name)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(color = "blue", size = 3) +
  geom_text_repel(box.padding = 0.5, point.padding = 0.5) +
  theme_minimal() +
  labs(
    title = "Tajima's D vs. Effective Population Size at HAR1A Locus",
    x     = "Effective Population Size (Ne)",
    y     = "Tajima's D"
  )


# =============================================================================
# SECTION 5: Export Results
# =============================================================================
# Note: write.csv2 uses semicolon separators (European locale standard)

write.csv2(har1a_combined_dataset,
           file      = "Output/har1a_combined_dataset.csv",
           row.names = FALSE)

write.csv(har1a_stats_with_D,
          file      = "Output/har1a_stats_with_D.csv",
          row.names = FALSE)