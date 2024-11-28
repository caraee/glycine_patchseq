scRNAseq Analysis with pySCENIC and mixOmics

This repository contains R and Python scripts for analyzing single-cell RNA sequencing (scRNAseq) data. The workflow includes data preprocessing,
running pySCENIC for gene regulatory network (GRN) inference, processing SCENIC output, and performing integrative analysis with mixOmics. 
This project specifically investigates glycine receptor signaling and associated transcription factors in pancreatic islet cells.

Workflow Overview

    Initial Data Processing and Cleaning (R):
        Performs decontamination using celda.
        Annotates cell type.
        Filters cells and genes.
        Normalizes, and scales data for downstream analyses.
        Generates visualizations like PCA, UMAP, and tSNE for exploratory analysis.

    pySCENIC GRN Inference (Python):
        This step should be performed on a high performance cluster.
        Converts processed data to Loom format.
        Executes pySCENIC to identify active regulons and their targets.
        Aggregates results from multiple pySCENIC runs for consensus GRN inference.

    SCENIC Output Processing (R):
        Loads pySCENIC output into an R environment.
        Integrates regulon activity scores with transcriptomic data for downstream analysis.

    Integrative Analysis with mixOmics (R):
        Combines transcriptomic, SCENIC regulon activity, and electrophysiology datasets.
        Explores correlations and relationships between gene regulatory networks and cellular electrophysiological properties.

Acknowledgments

This project was developed using publicly available tools and libraries. Special thanks to the developers of pySCENIC, Seurat, and mixOmics for their contributions to bioinformatics. 
The counts and metadata for these analyses are available at GEO under the accession number GSE280267
