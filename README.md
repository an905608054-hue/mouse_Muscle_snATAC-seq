# mouse_Muscle_snATAC-seq
This repository contains key R scripts for the analysis of  snATAC-seq and snRNA-seq data from mouse muscle, presented in our study. These scripts cover quality control, clustering, pseudotime analysis, differential analysis, and multi - omics integration.

## Script Overview

| Script Name          | Purpose                                                                 |
|----------------------|-------------------------------------------------------------------------|
| `QC+clustering.r`    | Performs quality control (QC) on snATAC-seq data and conducts cell clustering. |
| `pseudotime.r`       | Infers pseudotemporal ordering of Type IIb myonuclei to capture developmental or dynamic processes. |
| `DAR.txt`              | Identifies Differentially Accessible Regions (DARs) from ATAC-seq data. |
| `Integrate_RNA_ATAC.r` | Integrates snATAC-seq and snRNA-seq datasets to uncover multi-omics relationships. |
