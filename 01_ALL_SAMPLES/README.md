# Characterization of NK Cells Infiltrating Lung Tumor Bed in Multi-omic scATACseq

## Overview

This project aims to describe gene regulation networks (GRN) and differentiation trajectories of NK cells within the lung tumor niche. Using integrated single-cell ATAC-seq and RNA-seq data, we analyze gene expression, chromatin accessibility, and regulatory networks, focusing on the NK cell population.

---

## Repository Structure

The analysis scripts are organized into distinct modules, each corresponding to specific analyses performed during the study. Below is a branching view of the `03_Script` directory:

```plaintext
03_Script/
├── 01_Insight_From_Tang_Data            # Initial observations using public datasets
├── 02_Preparing_DataTables_For_Scoring  # Extract gene signatures for scoring analyses
├── 03_First_Look_at_Data                # Data exploration: clustering, signatures, QC
├── 04_SCENICplus_Analysis               # GRN inference with SCENIC+ pipeline
├── 05_MultiVelocity                     # scRNA/ATAC-based trajectory analysis
├── 06_Cytosig                           # Cytokine response inference
├── 07_MultiVelocity_NK3                 # Trajectory analysis focused on NK3 population
├── 08_SCENICplus_Analysis_NK3           # GRN inference focused on NK3 population
└── 09_Label_Transfer                    # Transfer of NK cell labels using label transfer methods

```
---

## Major Results

1. **Mapping NK Cell Clusters** 
   Classification of NK cells into well-known subsets using both scATAC and scRNA data.

2. **NK Cell Trajectories** 
   Identification of differentiation trajectories of NK cells in the tumor bed.

3. **Characterization of GRNs** 
   Identification of major gene regulatory networks in NK cell subsets within the tumor niche.

4. **Modeling Perturbation Effects** 
   Simulation of GRN perturbations and their impact on NK cell trajectories.

---

## How to Use This Repository

1. **Data Exploration** 
   Begin with `03_First_Look_at_data` for an overview of the dataset and preliminary insights.

2. **GRN Analysis** 
   Use `04_SCENICplus_Analysis` or `08_SCENICplus_Analysis_NK3` for GRN inference.

3. **Trajectory Analysis** 
   Dive into `05_MultiVelocity` or `07_MultiVelocity_NK3` for velocity-based trajectory analysis.

4. **Cytokine Response** 
   Refer to `06_Cytosig` for cytokine-specific responses inference.

5. **Label Transfer** 
   Use `09_Label_Transfer` to explore label transfer methodologies for NK subsets.

---

This repository provides a comprehensive pipeline for the analysis of NK cells in the lung tumor microenvironment, using tools to explore data, infer regulatory networks, and model cellular dynamics.


