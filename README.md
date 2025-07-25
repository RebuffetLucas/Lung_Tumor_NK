# Tumor-associated NK Cells with Cytotoxic Potential in Lung Cancer

This repository accompanies the study:

**"Integrated Single Cell Analysis Identifies CD39+ Tumor-Associated NK Cells with Cytotoxic Potential in Lung Cancer"**  
*Clara Serger<sup>1</sup>, Lucas Rebuffet<sup>2,11</sup>, Michael T. Sandholzer<sup>1,11</sup>, Irene Fusi<sup>1</sup>, Nicole Oelgarth<sup>1</sup>, Sofia Tundo<sup>1,3</sup>, Aljaz Hojski<sup>5</sup>, Didier Lardinois<sup>5</sup>, Marcel Trefny<sup>1</sup>, Nicole Kirchhammer<sup>1</sup>, Marina Natoli<sup>1</sup>, Matthias Matter<sup>4</sup>, Karin Schaeuble<sup>1</sup>, Eric Vivier<sup>2,6-9,12</sup>, Andrea Romagnani<sup>3,12</sup>, Alfred Zippelius<sup>1,10,12,13</sup>*

---

## 🧬 Summary

Despite growing interest in NK cell-targeting immunotherapies, the functional states and transcriptional programs of tumor-infiltrating NK cells remain incompletely characterized. Here, we used matched single-nucleus RNA and ATAC sequencing (snRNA-seq, snATAC-seq) to deconvolve the heterogeneity of intratumoral NK cells in non-small cell lung cancer (NSCLC).

We identified two tumor-associated NK (taNK) cell subsets (CD103⁺ and CD49a⁺) with features of circulating NK3 cells, but also signatures of tissue residency, dysfunction, and adaptive reprogramming. Trajectory and regulon analyses showed that inflammatory signals drive their differentiation from early GZMK⁺ NK3 precursors towards a CD39⁺ (ENTPD1⁺) effector phenotype, marked by enhanced interferon signaling and cytotoxic capacity. Functional assays further revealed that CD39⁺ taNK cells retain and even boost cytotoxicity following cytokine stimulation and NKG2A blockade, highlighting them as promising targets for immunotherapy.

---

## 📁 Repository structure

This repository contains the full pipeline used in the study, structured into modular subfolders corresponding to different stages of the analysis.

### Core folders

- `01_QC/`: Quality control for RNA and ATAC data, object creation
- `02_GlobalHeterogeneity/`: NK cell extraction, integration, clustering, signature scoring, functional annotation, and regulon analysis
- `07_MultiVelocity_NK3/`: Velocity-based inference of dynamic NK3 trajectories on a patient-by-patient basis
- `08_SCENICplus_Analysis_NK3/`: Regulatory network reconstruction using SCENIC+
- `10_DiffusionMap/`: Diffusion map pseudotime reconstruction on RNA and ATAC modalities
- `11_MultiVelocity_NK3_Full_Reprocessed/`: Reprocessed multi-velocity pipeline for joint patient analysis
- `12_Trajectory_Characterization/`: Mapping of NK3 trajectories based on gene activity, gene expression, and regulons

### Configuration files

- `globalParams.R` / `globalParams.py`: Global constants and paths
- `sampleParams.R` / `sampleParams.py`: Sample-specific metadata

---

## 🛠️ Requirements

Most scripts are written in R and Python, with a few notebooks (`.ipynb`) for exploratory and multi-step analyses. See each folder’s `analysisParams.*` for required packages.

- R ≥ 4.2
- Python ≥ 3.8
- Dependencies: Seurat, Signac, SCENIC+, Scanpy, pycistopic, destiny, chromVAR, BayesPrism...

We recommend using `renv` (for R), Docker or Singularity, and `conda` or `virtualenv` (for Python) to isolate environments.

---

## 📌 Key insights

- Discovery of **two major taNK cell states** (CD103⁺ and CD49a⁺) with tissue-residency and adaptive features.
- CD39⁺ NK3 cells emerge as **cytotoxic and targetable effectors** in NSCLC tumors.
- Joint **snRNA-seq and snATAC-seq** resolves NK3 differentiation continuum.
- **Regulon and velocity analyses** reveal drivers of dysfunction and effector adaptation.
- **Therapeutic implications**: cytokine priming and NKG2A blockade may boost NK responses in lung tumors.

---

## 👥 Authors and affiliations

<sup>1</sup>Department of Biomedicine, University Hospital and University of Basel, Switzerland  
<sup>2</sup>CIML, Aix Marseille Université, CNRS, INSERM, France  
<sup>3</sup>Cancer Cell Dependencies, Roche Pharma Research and Early Development, Switzerland  
<sup>4</sup>Department of Pathology, University Hospital Basel, Switzerland  
<sup>5</sup>Department of Thoracic Surgery, University Hospital Basel, Switzerland  
<sup>6–9</sup>Innate Pharma, APHM, Paris-Saclay Cancer Cluster, INSERM U1015  
<sup>10</sup>Medical Oncology, University Hospital Basel, Switzerland  
<sup>11</sup>These authors contributed equally  
<sup>12</sup>Senior authors, contributed equally  
<sup>13</sup>Lead contact


---

## 📬 Contact

For questions or data access requests, please reach out to the lead contact:

**Prof. Alfred Zippelius**  
University Hospital Basel, Switzerland  


---

## 📖 Citation

> Please cite our paper once published (citation pending). A preprint or DOI link will be added here once available.

---


