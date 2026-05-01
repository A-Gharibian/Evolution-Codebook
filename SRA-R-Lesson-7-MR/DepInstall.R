# -----------------------------------------------------------------------------
# Package Management
# -----------------------------------------------------------------------------
getOption("repos")
options(repos = c(CRAN = "https://cran.rstudio.com"))
# CRAN bootstrap
if (!require("pacman", quietly = TRUE)) install.packages("pacman")

# Bioconductor bootstrap — required for ComplexHeatmap
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("ComplexHeatmap")
bioc_missing  <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(bioc_missing) > 0) {
  message("Installing missing Bioconductor packages: ", paste(bioc_missing, collapse = ", "))
  BiocManager::install(bioc_missing, ask = FALSE)
}

# Core packages required for the full pipeline
pacman::p_load(
  # Data I/O & wrangling
  jsonlite, readr, dplyr, tidyr, forcats, scales,
  
  # Biological & phylogenetic analysis
  seqinr, pegas, ape,
  
  # Plotting & tree visualization
  ggplot2, ggtree, ggrepel, ggalluvial,
  
  # Heatmaps & color scales
  pheatmap, circlize, viridis
)

# Verify all packages loaded successfully and report any failures
not_loaded <- pacman::p_loaded()[!sapply(pacman::p_loaded(), requireNamespace, quietly = TRUE)]
if (length(not_loaded) > 0) {
  warning("The following packages failed to load: ", paste(not_loaded, collapse = ", "))
} else {
  cat("All packages loaded successfully.\n")
}

# -----------------------------------------------------------------------------
# readr helpers — consistent wrappers for pipeline I/O
# -----------------------------------------------------------------------------

read_csv_clean <- function(path, col_names_override = NULL) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!is.null(col_names_override)) {
    stopifnot(
      "col_names_override length must match number of columns" =
        length(col_names_override) == ncol(df)
    )
    colnames(df) <- col_names_override
  }
  df
}

