# HAR1A and RHO Evolutionary History Pipeline

A pipeline for exploring the evolutionary history and population genetics of the **HAR1A** and **RHO** genes using public SRA sequencing data, Jaccard-based phylogenetic clustering, and Tajima's D analysis.

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

## Overview

This project investigates the comparative genomics and population genetics of two genes:

- **HAR1A (Human Accelerated Region 1A)** — a non-coding RNA locus with unusually rapid evolution in the human lineage.
- **RHO (Rhodopsin)** — a well-characterized visual system gene used as a comparative reference across vertebrates.

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
