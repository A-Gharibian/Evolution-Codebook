# =============================================================================
# config.R
# Purpose: Define global paths, constants, and drop rules.
# =============================================================================

# File Paths
RHO_INPUT_PATH     <- "Input/RHO.json"
HAR1A_INPUT_PATH   <- "Input/HAR1A.json"
RHO_META_CSV       <- "Input/final_rho_table_clean.csv"
HAR1A_META_CSV     <- "Input/final_har_table_clean.csv"
VERT_MUT_RATES_CSV <- "Input/Vertebrate mutation rates.csv"
HAR1A_FASTA_FILE   <- "Input/trimmed_43_HAR1A_31bp.fasta"
RHO_FASTA_FILE     <- "Input/trimmed_12_RHO_31bp.fasta"

# -----------------------------------------------------------------------------
# Columns & Drop Rules
# -----------------------------------------------------------------------------
META_COL_NAMES <- c("acc", "organism", "assay_type", "biosample", "bioproject", "mbases")

RHO_SPECIES_DROP <- c("Cat", "Green")
HAR1A_VARS_DROP  <- c(
  "Variation3",  "Variation7",  "Variation8",  "Variation9",  "Variation10",
  "Variation16", "Variation19", "Variation23", "Variation25", "Variation26",
  "Variation31", "Variation32", "Variation38", "Variation40", "Variation41",
  "Variation42"
)

# -----------------------------------------------------------------------------
# Constants & Mappings
# -----------------------------------------------------------------------------
GENOMIC_ASSAYS <- c("WGS", "WGA", "WCS")

EXCLUDE_ORGANISMS <- c(
  "SARS-CoV-2",
  "Arabidopsis thaliana",
  "Drosophila melanogaster",
  "Plasmodium falciparum",
  "Anopheles gambiae"
)

exclude_patterns <- c(
  "metagenome",
  "synthetic construct",
  "bacterium",
  "mixed sample",
  "uncultured eukaryote",
  "uncultured bacterium",
  "unidentified"
)

pattern_map <- c(
  "Mus musculus" = "Mus musculus",
  "Severe acute respiratory syndrome coronavirus 2" = "SARS-CoV-2"
)

exact_map <- c(
  "Gallus"                  = "Gallus gallus",
  "Mus musculus castaneus"  = "Mus musculus",
  "Mus musculus domesticus" = "Mus musculus"
)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
top_n_rho  <- 16L
top_n_har  <- 16L

MIN_TOTAL_SRA <- 8192L

RHO_MIN_FREQ   <- 128L
HAR1A_MIN_FREQ <- 128L

MIN_VARIATION_FREQ <- 0.001  # Relative MAF cutoff (0.1%)

# Column names holding per-site base calls; must match sra_sequence_mapping.R
N_POSITIONS <- 31L
SEQ_COLS    <- paste0("pos_", seq_len(N_POSITIONS))

MIN_N <- 512L  # Minimum valid SRA hits required for statistical evaluation
