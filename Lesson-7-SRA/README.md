Here is the complete text from start to finish, fully formatted as Markdown in a single, continuous code block. You can copy the content inside the block and paste it directly into your `README.md` file.

```markdown
# HAR1A & RHO — Comparative Genomics and Population Genetics Pipeline

A pipeline for exploring the evolutionary history and population genetics of the HAR1A and RHO genes using public SRA sequencing data, Jaccard-based phylogenetic clustering, and Tajima's D analysis.

## Table of Contents
- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

This project investigates the comparative genomics and population genetics of two genes:
* **HAR1A (Human Accelerated Region 1A)** — a non-coding RNA locus with unusually rapid evolution in the human lineage
* **RHO (Rhodopsin)** — a well-characterized visual system gene used as a comparative reference across vertebrates

The pipeline queries the NCBI Sequence Read Archive (SRA) for publicly available sequencing data, clusters species/variations by SRA co-occurrence using Jaccard distance, visualizes results as phylogenetic heatmaps, and computes Tajima's D to test for signatures of selection at the HAR1A locus.

---

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

---

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

---

## Dependencies

### R (≥ 4.1.0)

| Package | Source | Purpose |
| :--- | :--- | :--- |
| `jsonlite` | CRAN | JSON parsing |
| `dplyr` | CRAN | Data manipulation |
| `tidyr` | CRAN | Data reshaping |
| `ggplot2` | CRAN | Plotting |
| `ggrepel` | CRAN | Non-overlapping plot labels |
| `ape` | CRAN | Phylogenetic tree conversion |
| `ggtree` | Bioconductor | Phylogenetic tree visualization |
| `ComplexHeatmap` | Bioconductor | Heatmap rendering (Option B) |
| `pheatmap` | CRAN | Heatmap rendering (Option C) |
| `viridis` | CRAN | Color palettes |

### System

| Tool | Purpose |
| :--- | :--- |
| `bash` | SRA metadata extraction pipeline |
| `awk` | Column extraction from SRA metadata CSVs |

---

## Installation

### R packages

```R
# CRAN packages
install.packages(c(
  "jsonlite", "dplyr", "tidyr", "ggplot2",
  "ggrepel", "ape", "pheatmap", "viridis"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ComplexHeatmap"))
```

---

## Usage

Run the scripts in order. Each script depends on objects produced by the previous one.

### 1. Process JSON files
```R
source("R/script1_json_processing.R")
```
> Produces `target_sras.txt`, `rho_target_sras.txt`, `har_target_sras.txt`, `SRA_Mapping_Table.csv`, and cleaned JSON files.

### 2. Query the SRA index (Bash)
```bash
bash bash/extract_metadata.sh
```
> Uses the `.txt` accession lists to query a local SRA index and produces `final_rho_table_clean.csv` and `final_har_table_clean.csv`. See the script header for index path configuration.

### 3. Process metadata
```R
source("R/script2_metadata.R")
```
> Merges SRA metadata with gene mappings, filters metagenomes, and standardizes organism names.

### 4. Cluster by SRA co-occurrence
```R
source("R/script3_clustering.R")
```
> Computes pairwise Jaccard distances using biosample co-occurrence and performs hierarchical clustering.

### 5. Visualize heatmaps
```R
source("R/script4_heatmap_visualization.R")
```
> Renders phylogenetic trees with SRA metadata heatmap overlays for both genes.  
> **Note:** Three heatmap backends are provided (`ggtree`, `ComplexHeatmap`, `pheatmap`). Select one and remove the others before final use — see the TODO block in the script.

### 6. Compute Tajima's D
```R
source("R/script5_tajimas_d.R")
```
> Computes S, π, and Tajima's D per organism from aligned HAR1A sequences and plots results against vertebrate effective population sizes.

---

## Inputs

| File | Description |
| :--- | :--- |
| `RHO.json` | One JSON file per gene mapping variation/species names to SRA accessions |
| `HAR1A_trim.json` | Preprocessed HAR1A JSON with biologically redundant variations removed |
| `final_rho_table_clean.csv` | SRA metadata for RHO (columns: `acc`, `organism`, `assay_type`, `biosample`, `sra_study`, `bioproject`, `mbases`) |
| `final_har_table_clean.csv` | SRA metadata for HAR1A (same schema) |
| `Vertebrate mutation rates.csv` | Life history traits per species, including effective population size (Ne) |

---

## Outputs

| File | Description |
| :--- | :--- |
| `SRA_Mapping_Table.csv` | Combined gene × variation × SRA accession mapping |
| `target_sras.txt` | Combined RHO + HAR1A SRA accession list for full index query |
| `rho_target_sras.txt` | RHO-only SRA accession list for gene-specific index query |
| `har_target_sras.txt` | HAR1A-only SRA accession list for gene-specific index query |
| `RHO_cleaned_names.json` | RHO JSON with corrected variation names |
| `HAR_cleaned_names.json` | HAR1A JSON with biologically redundant variations removed |
| `Output/har1a_stats_with_D.csv` | Per-organism S, π, θ_W, and Tajima's D for HAR1A |
| `Output/har1a_combined_dataset.csv` | Tajima's D merged with vertebrate life history traits |

---

## Citation

**[Placeholder]** — A manuscript describing this pipeline is in preparation.  
If you use this pipeline in your work, please check back here for an updated citation or contact the authors directly.

---

## Contributing

Contributions, bug reports, and suggestions are welcome. Please open an issue or submit a pull request.

When contributing:
* Follow the existing script structure and header format.
* Add `dplyr::` namespace prefixes for all `dplyr` calls.
* Document any new functions with `roxygen`-style comments.
* Test changes against both RHO and HAR1A datasets before submitting.

---

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).  
See `LICENSE` for the full license text or visit [https://www.gnu.org/licenses/agpl-3.0](https://www.gnu.org/licenses/agpl-3.0).
```
