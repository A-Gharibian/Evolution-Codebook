# A pipeline for exploring the evolutionary history and population genetics of the **HAR1A** and **RHO** genes using public SRA sequencing data, Jaccard-based phylogenetic clustering, and Tajima's D analysis.

## Table of Contents

* [Overview](#overview)
* [Pipeline Architecture](#pipeline-architecture)
* [Repository Structure](#repository-structure)
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Usage](#usage)
* [Inputs](#inputs)
* [Outputs](#outputs)
* [Citation](#citation)
* [Contributing](#contributing)
* [License](#license)

## Overview

This project investigates the comparative genomics and population genetics of two genes:

* **HAR1A (Human Accelerated Region 1A)** — a non-coding RNA locus with unusually rapid evolution in the human lineage.
* **RHO (Rhodopsin)** — a well-characterized visual system gene used as a comparative reference across vertebrates.

The pipeline queries the NCBI Sequence Read Archive (SRA) for publicly available sequencing data, clusters species/variations by SRA co-occurrence using Jaccard distance, visualizes results as phylogenetic heatmaps, and computes Tajima's D to test for signatures of selection at the HAR1A locus.

## Pipeline Architecture

```text
JSON Files (per gene)
        │
        ▼
┌─────────────────────┐
│  Script 1           │  Load & clean JSON, filter variations,
│  json_processing.R  │  export SRA accession lists and mapping table
└─────────────────────┘
        │
        │  target_sras.txt / rho_target_sras.txt / har_target_sras.txt
        ▼
┌─────────────────────┐
│  Bash Pipeline      │  Query SRA index, extract metadata
│  (external)         │  → final_rho_table_clean.csv
└─────────────────────┘    final_har_table_clean.csv
        │
        ▼
┌─────────────────────┐
│  Script 2           │  Load metadata CSVs, merge with SRA mapping,
│  metadata.R         │  filter metagenomes, clean organism names
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Script 3           │  Compute pairwise Jaccard distances,
│  clustering.R       │  hierarchical clustering (average linkage)
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Script 4           │  Phylogenetic tree + SRA heatmap overlays
│  heatmap_viz.R      │  (ggtree, ComplexHeatmap, or pheatmap)
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Script 5           │  Tajima's D per organism, merge with
│  tajimas_d.R        │  vertebrate life history traits, plot vs Ne
└─────────────────────┘
```

## Repository Structure

```text
.
├── data/
│   ├── RHO.json                          # Raw RHO variation → SRA mapping
│   ├── HAR1A_trim.json                   # Preprocessed HAR1A variation → SRA mapping
│   └── Vertebrate mutation rates.csv     # Life history traits including Ne
│
├── bash/
│   └── extract_metadata.sh               # SRA index query and metadata extraction
│
├── R/
│   ├── script1_json_processing.R
│   ├── script2_metadata.R
│   ├── script3_clustering.R
│   ├── script4_heatmap_visualization.R
│   └── script5_tajimas_d.R
│
├── Output/                               # Generated output files (not tracked by git)
│   ├── SRA_Mapping_Table.csv
│   ├── RHO_cleaned_names.json
│   ├── HAR_cleaned_names.json
│   ├── har1a_combined_dataset.csv
│   └── har1a_stats_with_D.csv
│
└── README.md
```

## Dependencies

### R (≥ 4.1.0)

| Package        | Source       | Purpose                         |
| -------------- | ------------ | ------------------------------- |
| jsonlite       | CRAN         | JSON parsing                    |
| dplyr          | CRAN         | Data manipulation               |
| tidyr          | CRAN         | Data reshaping                  |
| ggplot2        | CRAN         | Plotting                        |
| ggrepel        | CRAN         | Non-overlapping plot labels     |
| ape            | CRAN         | Phylogenetic tree conversion    |
| ggtree         | Bioconductor | Phylogenetic tree visualization |
| ComplexHeatmap | Bioconductor | Heatmap rendering (Option B)    |
| pheatmap       | CRAN         | Heatmap rendering (Option C)    |
| viridis        | CRAN         | Color palettes                  |

### System

| Tool | Purpose                                  |
| ---- | ---------------------------------------- |
| bash | SRA metadata extraction pipeline         |
| awk  | Column extraction from SRA metadata CSVs |

## Installation

### R packages

```r
# CRAN packages
install.packages(c(
  "jsonlite", "dplyr", "tidyr", "ggplot2",
  "ggrepel", "ape", "pheatmap", "viridis"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ComplexHeatmap"))
```

## Usage

Run the scripts in order. Each script depends on objects produced by the previous one.

### 1. Process JSON files

```r
source("R/script1_json_processing.R")
```

### 2. Query the SRA index (Bash)

```bash
bash bash/extract_metadata.sh
```

### 3. Process metadata

```r
source("R/script2_metadata.R")
```

### 4. Cluster by SRA co-occurrence

```r
source("R/script3_clustering.R")
```

### 5. Visualize heatmaps

```r
source("R/script4_heatmap_visualization.R")
```

### 6. Compute Tajima's D

```r
source("R/script5_tajimas_d.R")
```

## Inputs

| File                            | Description                                                              |
| ------------------------------- | ------------------------------------------------------------------------ |
| `RHO.json`                      | One JSON file per gene mapping variation/species names to SRA accessions |
| `HAR1A_trim.json`               | Preprocessed HAR1A JSON with biologically redundant variations removed   |
| `final_rho_table_clean.csv`     | SRA metadata for RHO                                                     |
| `final_har_table_clean.csv`     | SRA metadata for HAR1A                                                   |
| `Vertebrate mutation rates.csv` | Life history traits per species                                          |

## Outputs

| File                                | Description                                       |
| ----------------------------------- | ------------------------------------------------- |
| `SRA_Mapping_Table.csv`             | Combined gene × variation × SRA accession mapping |
| `target_sras.txt`                   | Combined accession list                           |
| `rho_target_sras.txt`               | RHO-only accession list                           |
| `har_target_sras.txt`               | HAR1A-only accession list                         |
| `RHO_cleaned_names.json`            | Cleaned RHO JSON                                  |
| `HAR_cleaned_names.json`            | Cleaned HAR1A JSON                                |
| `Output/har1a_stats_with_D.csv`     | Tajima's D results                                |
| `Output/har1a_combined_dataset.csv` | Merged dataset                                    |

## Citation

[Placeholder] — Manuscript in preparation.

## Contributing

* Follow existing script structure
* Use `dplyr::` prefixes
* Document functions with roxygen comments
* Test on both datasets

## License

GNU Affero General Public License v3.0 (AGPL-3.0).
[https://www.gnu.org/licenses/agpl-3.0](https://www.gnu.org/licenses/agpl-3.0)
