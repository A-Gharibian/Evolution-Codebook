# =============================================================================
# Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# Flatten a named JSON list into a three-column data frame.
# -----------------------------------------------------------------------------
#' @param json_list  Named list where each element is itself a named list of
#'                   SRA accessions.
#' @param gene_name  Character scalar; gene label added to every row as the
#'                   Gene column.
#'
#' @return A data frame with columns Gene, Variation_ID, and SRA_Accession,
#'         or an empty data frame if every element of json_list is empty.
flatten_sra_data <- function(json_list, gene_name) {
  stopifnot(
    is.list(json_list),
    is.character(gene_name),
    length(gene_name) == 1L
  )
  
  df_list <- lapply(names(json_list), function(var_name) {
    sra_names <- names(json_list[[var_name]])
    
    if (length(sra_names) == 0L) return(NULL)
    
    data.frame(
      Gene          = gene_name,
      Variation_ID  = var_name,
      SRA_Accession = sra_names,
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(df_list)
}


# -----------------------------------------------------------------------------
# Apply exact and pattern-based name mappings to a character vector.
# -----------------------------------------------------------------------------
#' @param x           Character vector of organism names.
#' @param exact_map   Named character vector; keys are the original names,
#'                    values are the replacements.
#' @param pattern_map Named character vector; keys are regex patterns,
#'                    values are replacements (applied in order).
#' @return Character vector of normalised names, same length as x.
apply_name_mapping <- function(x, exact_map, pattern_map) {
  x <- ifelse(x %in% names(exact_map), exact_map[x], x)
  
  for (pattern in names(pattern_map)) {
    x[grepl(pattern, x, ignore.case = TRUE)] <- pattern_map[[pattern]]
  }
  
  x
}


# -----------------------------------------------------------------------------
# Load metadata CSV, merge with flattened JSON accessions, and filter.
# -----------------------------------------------------------------------------
#' @param json_list Named list (loaded from a *_cleaned_names.json file).
#' @param csv_path  Path to the metadata CSV produced by the bash step.
#' @param gene_name Character scalar passed to flatten_sra_data().
#' @param col_names Character vector of column names to impose on the CSV.
#'                  Defaults to META_COL_NAMES defined in config.R.
#' @return A joined and metagenome-filtered data frame.
load_and_filter_meta <- function(json_list, csv_path, gene_name,
                                 col_names = META_COL_NAMES) {
  stopifnot(
    is.list(json_list),
    is.character(csv_path),    length(csv_path)  == 1L,
    is.character(gene_name),   length(gene_name) == 1L,
    is.character(col_names)
  )
  
  meta <- tryCatch(
    readr::read_csv(csv_path, show_col_types = FALSE),
    error = function(e) stop("Could not read '", csv_path, "': ", e$message)
  )
  
  stopifnot(
    "col_names length must match number of columns in CSV" =
      length(col_names) == ncol(meta)
  )
  colnames(meta) <- col_names
  
  query_df <- flatten_sra_data(json_list, gene_name)
  if (nrow(query_df) == 0L) stop("Flattened JSON is empty for: ", gene_name)
  
  dplyr::inner_join(query_df, meta, by = c("SRA_Accession" = "acc")) |>
    dplyr::filter(!grepl("metagenome", organism, ignore.case = TRUE))
}


# -----------------------------------------------------------------------------
# Build a pairwise Jaccard distance matrix from a named list of sets.
# -----------------------------------------------------------------------------
#' Each element of json_list is treated as a set whose members are the
#' element's names (SRA / biosample accessions). The Jaccard distance between
#' two sets A and B is defined as 1 - |A ∩ B| / |A ∪ B|, with the
#' convention that the distance is 0 when both sets are empty.
#'
#' @param json_list Named list; names become row/column labels. Each element
#'                  must itself be a named vector or list (names are the set
#'                  members).
#'
#' @return An object of class "dist".
build_jaccard_matrix <- function(json_list) {
  stopifnot(
    is.list(json_list),
    length(json_list) >= 2L,
    !is.null(names(json_list))
  )
  
  n      <- length(json_list)
  groups <- names(json_list)
  mat    <- matrix(0, nrow = n, ncol = n, dimnames = list(groups, groups))
  
  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      set_i         <- names(json_list[[i]])
      set_j         <- names(json_list[[j]])
      intersect_len <- length(intersect(set_i, set_j))
      union_len     <- length(union(set_i, set_j))
      
      j_dist <- if (union_len == 0L) 0 else 1 - (intersect_len / union_len)
      
      mat[i, j] <- j_dist
      mat[j, i] <- j_dist
    }
  }
  
  as.dist(mat)
}


# -----------------------------------------------------------------------------
# Convert a filtered data frame into the named-list format expected by
# build_jaccard_matrix().
# -----------------------------------------------------------------------------
#' @param df            Data frame containing at least one grouping column and
#'                      one accession column.
#' @param group_col     Name of the grouping column (character scalar).
#' @param accession_col Name of the accession column (character scalar).
#'                      Defaults to "biosample".
#'
#' @return Named list; each element is a named integer vector of 1s whose
#'         names are biosample accessions.
prepare_jaccard_input <- function(df, group_col, accession_col = "biosample") {
  stopifnot(
    is.data.frame(df),
    group_col     %in% colnames(df),
    accession_col %in% colnames(df)
  )
  
  split_list <- split(df[[accession_col]], df[[group_col]])
  
  lapply(split_list, function(accessions) {
    setNames(rep(1L, length(accessions)), accessions)
  })
}


# -----------------------------------------------------------------------------
# Render a dendrogram to a PNG file.
# -----------------------------------------------------------------------------
#' @param hc        An hclust object.
#' @param file      Output file path (character scalar).
#' @param title     Plot title.
#' @param use_dend  Logical; if TRUE coerce to dendrogram for angled leaf
#'                  labels, otherwise plot the hclust object directly.
#' @param mar       Numeric vector of length 4 passed to par(mar = ...).
#' @param cex       Character expansion factor for leaf labels.
#' @param width,height,res  PNG dimensions and resolution.
save_dendrogram_png <- function(hc,
                                file,
                                title,
                                use_dend = TRUE,
                                mar      = c(10, 4, 4, 2),
                                cex      = 0.85,
                                width    = PNG_WIDTH,
                                height   = PNG_HEIGHT,
                                res      = PNG_RES) {
  grDevices::png(file, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off())
  
  graphics::par(mar = mar, cex = cex)
  
  if (use_dend) {
    plot(
      as.dendrogram(hc),
      main = title,
      ylab = "Jaccard Distance",
      lwd  = 2
    )
  } else {
    plot(
      hc,
      main = title,
      ylab = "Jaccard Distance",
      xlab = "",
      sub  = ""
    )
  }
}


# -----------------------------------------------------------------------------
# Count SRA hits along each unique combination of grouping columns.
# -----------------------------------------------------------------------------
#' @param df       A filtered SRA metadata data frame.
#' @param ...      Unquoted column names to group by (passed to
#'                 dplyr::group_by()).
#' @param min_freq Integer; rows with fewer hits than this value are dropped.
#'                 Default 1L (keep all).
#'
#' @return A data frame with the grouping columns plus Freq.
aggregate_flow <- function(df, ..., min_freq = 1L) {
  stopifnot(is.data.frame(df), is.numeric(min_freq), min_freq >= 1L)
  
  df |>
    dplyr::group_by(...) |>
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop") |>
    dplyr::filter(Freq >= min_freq)
}


# -----------------------------------------------------------------------------
# Extract and normalise variation labels from seqinr FASTA annotations.
# -----------------------------------------------------------------------------
#' @note  Expects annotations of the form "> Variation N ...". For sequences
#'        that do not match this pattern the full cleaned annotation string is
#'        returned as-is, so the function degrades gracefully for other FASTA
#'        formats.
#'
#' @param fasta_list  Named list returned by seqinr::read.fasta().
#'
#' @return Character vector of clean variation names, one per sequence.
extract_variation_names <- function(fasta_list) {
  stopifnot(is.list(fasta_list), length(fasta_list) > 0L)
  
  raw_annotations <- sapply(fasta_list, function(seq) attr(seq, "Annot"))
  cleaned         <- gsub("^>\\s*|\\s*$", "", raw_annotations)
  
  has_variation <- grepl("^Variation\\s*\\d+", cleaned, ignore.case = TRUE)
  
  cleaned[has_variation] <- gsub(
    "\\s+", "",
    regmatches(
      cleaned[has_variation],
      regexpr("^Variation\\s*\\d+", cleaned[has_variation], ignore.case = TRUE)
    )
  )
  
  cleaned
}


# -----------------------------------------------------------------------------
# Convert a seqinr FASTA list into a data frame of per-position bases.
# -----------------------------------------------------------------------------
#' @param fasta_list   Named list returned by seqinr::read.fasta().
#' @param clean_names  Character vector of row labels (length must equal
#'                     length(fasta_list)).
#' @param key_col      Name of the column added to hold row labels.
#'                     Default "Variation_ID".
#'
#' @return Data frame with one row per sequence, columns pos_1 … pos_N,
#'         and the key column appended as the final column.
fasta_to_dataframe <- function(fasta_list, clean_names, key_col = "Variation_ID") {
  stopifnot(
    is.list(fasta_list),
    is.character(clean_names),
    length(clean_names) == length(fasta_list)
  )
  
  seq_matrix            <- do.call(rbind, fasta_list)
  seq_matrix            <- toupper(seq_matrix)
  rownames(seq_matrix)  <- clean_names
  
  seq_df            <- as.data.frame(seq_matrix, stringsAsFactors = FALSE)
  colnames(seq_df)  <- paste0("pos_", seq_len(ncol(seq_matrix)))
  seq_df[[key_col]] <- rownames(seq_df)
  
  seq_df
}


# =============================================================================
# Population Genetics Helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Per-site nucleotide diversity (pi).
# -----------------------------------------------------------------------------
#' Uses the unbiased estimator π = k/(k-1) * (1 - Σ p_i²), where k is the
#' number of valid (non-gap, non-ambiguous) bases and p_i are their
#' frequencies. Sites with fewer than two valid bases return 0.
#'
#' @param base_vec  Character vector of raw base calls (upper or lower case;
#'                  gaps, Ns, and NAs are silently removed).
#'
#' @return Numeric scalar in [0, 1].
calc_pi_site <- function(base_vec) {
  base_clean <- toupper(base_vec)
  base_clean <- base_clean[base_clean %in% c("A", "C", "T", "G")]
  
  k <- length(base_clean)
  if (k <= 1L) return(0)
  
  allele_freqs <- table(base_clean) / k
  (k / (k - 1L)) * (1 - sum(allele_freqs ^ 2))
}


# -----------------------------------------------------------------------------
# Harmonic number a(n) used in Watterson's estimator.
# -----------------------------------------------------------------------------
#' Computes a_n = Σ(i=1 to n-1) 1/i.
#' Returns 0 for n ≤ 1 (undefined for a sample of one sequence).
#'
#' @param n  Positive integer sample size.
#'
#' @return Numeric scalar.
calc_harmonic <- function(n) {
  if (n <= 1L) return(0)
  sum(1 / seq_len(n - 1L))
}


# =============================================================================
# Tajima's D
# =============================================================================

#' Calculate Tajima's D for a single organism.
#'
#' @param n       Integer. Number of sequences sampled.
#' @param S       Integer. Number of segregating (variable) sites.
#' @param pi_val  Numeric. Average pairwise nucleotide diversity.
#' @return        Numeric. Tajima's D statistic, or NA if inputs are invalid.
calc_tajimas_d <- function(n, S, pi_val) {
  if (n < 2L || S == 0L || is.na(S)) return(NA_real_)
  
  a1 <- sum(1 / seq_len(n - 1L))
  a2 <- sum(1 / seq_len(n - 1L) ^ 2)
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- (2 * (n ^ 2 + n + 3)) / (9 * n * (n - 1))
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / a1 ^ 2)
  e1 <- c1 / a1
  e2 <- c2 / (a1 ^ 2 + a2)
  
  variance <- (e1 * S) + (e2 * S * (S - 1))
  if (variance <= 0) return(NA_real_)
  
  theta_W <- S / a1
  (pi_val - theta_W) / sqrt(variance)
}


# =============================================================================
# Chord Diagram Helpers
# =============================================================================

#' Compute the pairwise overlap matrix for a named JSON list.
#'
#' @param json_list Named list; element names are group labels, sub-element
#'                  names are SRA accessions.
#' @return A symmetric numeric matrix with zero diagonal.
build_overlap_matrix <- function(json_list) {
  stopifnot(
    is.list(json_list),
    length(json_list) >= 2L,
    !is.null(names(json_list))
  )
  
  groups <- names(json_list)
  n      <- length(groups)
  mat    <- matrix(0L, nrow = n, ncol = n, dimnames = list(groups, groups))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      mat[i, j] <- length(intersect(
        names(json_list[[i]]),
        names(json_list[[j]])
      ))
    }
  }
  
  mat
}


#' Render a chord diagram to a PNG file.
#'
#' @param overlap_mat  Symmetric numeric matrix produced by build_overlap_matrix().
#' @param colors       Named character vector of sector colours.
#' @param file         Output file path.
#' @param title        Plot title.
#' @param gap_after    Gap between sectors in degrees.
#' @param label_cex    Character expansion for sector labels.
#' @param track_height Height of the label track.
#' @param width,height,res  PNG dimensions and resolution.
save_chord_png <- function(overlap_mat,
                           colors,
                           file,
                           title        = "",
                           gap_after    = 4,
                           label_cex    = 0.9,
                           track_height = 0.15,
                           width        = PNG_WIDTH,
                           height       = PNG_HEIGHT,
                           res          = PNG_RES) {
  grDevices::png(file, width = width, height = height, units = "in", res = res)
  on.exit(grDevices::dev.off())
  
  circlize::circos.clear()
  circlize::circos.par(gap.after = gap_after, start.degree = 90)
  
  circlize::chordDiagram(
    overlap_mat,
    grid.col          = colors,
    transparency      = 0.4,
    annotationTrack   = "grid",
    preAllocateTracks = list(track.height = track_height)
  )
  
  circlize::circos.track(
    track.index = 1,
    panel.fun   = function(x, y) {
      circlize::circos.text(
        circlize::CELL_META$xcenter,
        circlize::CELL_META$ylim[1],
        circlize::CELL_META$sector.index,
        facing     = "clockwise",
        niceFacing = TRUE,
        adj        = c(0, 0.5),
        cex        = label_cex,
        font       = 2
      )
    },
    bg.border = NA
  )
}


# =============================================================================
# SFS Helpers
# =============================================================================

#' Compute SFS shape descriptors and minor allele counts.
#'
#' @param seq_df A data frame containing only the sequence (position) columns.
#' @return A 1-row data frame with columns Segregating_Sites, Mean_MAC,
#'         Singleton_Prop, and SFS_Raw.
calc_sfs_metrics <- function(seq_df) {
  stopifnot(is.data.frame(seq_df))
  
  macs <- vapply(seq_df, function(col) {
    valid_bases <- toupper(col)
    valid_bases <- valid_bases[valid_bases %in% c("A", "C", "T", "G")]
    
    if (length(unique(valid_bases)) < 2L) return(NA_integer_)
    
    freqs <- table(valid_bases)
    as.integer(sum(freqs) - max(freqs))
  }, NA_integer_)
  
  macs <- macs[!is.na(macs)]
  S    <- length(macs)
  
  if (S == 0L) {
    return(data.frame(
      Segregating_Sites = 0L,
      Mean_MAC          = NA_real_,
      Singleton_Prop    = NA_real_,
      SFS_Raw           = NA_character_
    ))
  }
  
  sfs_table <- table(macs)
  
  data.frame(
    Segregating_Sites = S,
    Mean_MAC          = mean(macs),
    Singleton_Prop    = sum(macs == 1L) / S,
    SFS_Raw           = paste(names(sfs_table), sfs_table, sep = ":", collapse = " | ")
  )
}
