# SalmoSideroScan

**A pipeline to assess conservation of siderophore-related genes across _Salmonella_ genomes using BLAST.**

## Workflow Overview

1. **Genomes QC**
   - Filter genomes based on CheckM completeness (>90%) and contamination (<5%)
   - Remove genomes with anomalous CDS counts via linear regression filtering

2. **Filtering by Size and MLST**
   - Keep only genomes ≥ 4 Mbp (`remove4000kbp.py`)
   - Retain one strain per MLST type to maximize diversity (`MLST_filter.py`)

3. **Gene Detection**
   - Run ABRicate using a custom siderophore-related gene database
   - Genes detected at ≥75% identity and ≥60% coverage are considered present

4. **Phylogenetic Tree Construction**
   - Use PhyloPhlAn to reconstruct maximum likelihood tree
   - Annotate and visualize with iTOL
