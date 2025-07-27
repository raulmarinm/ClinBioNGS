# 🧬 Overview

**ClinBioNGS** is an open-source, panel-agnostic bioinformatics pipeline for the analysis of somatic NGS cancer panels using tumor-only DNA and RNA data. Designed for clinical and translational research, ClinBioNGS enables accurate, reproducible, and interpretable detection of a wide range of genomic alterations. Built with Nextflow and containerized environments, the pipeline ensures full portability, modularity, and flexibility across computing platforms.

ClinBioNGS integrates consensus small variant calling, panel-specific CNA and MSI reference models, automated annotation and prioritization modules, internal QC systems, and a variant database for longitudinal tracking. Final results are compiled into interactive, visual HTML reports, facilitating multidisciplinary interpretation. The pipeline has been validated on SEQC2 reference datasets and thousands of real-world clinical samples, showing performance comparable to commercial tools while offering greater transparency and extended biomarker support.

## 🚀 Key Features

- 🔬 **Comprehensive variant detection**:
  - Small variants (SNVs/indels)
  - Copy-number alterations (CNAs)
  - Gene fusions and splice variants (RNA-seq)
  - Complex biomarkers: TMB and MSI

- 🧩 **Panel-agnostic and modular**:
  - Supports multiple targeted NGS panels
    - Full analysis: Illumina TruSight™ Oncology 500 (TSO500), Thermo Fisher Oncomine™ Precision Assay (OPA), and Thermo Fisher Oncomine™ Comprehensive Assay
  - Adaptable to other hybrid-capture or amplicon panels

- ⚙️ **Reproducible and portable**:
  - Built with Nextflow and containerized using Singularity images
  - Compatible with SLURM, SGE, and local environments

- 🧠 **Smart integration**:
  - Consensus small variant calling with multiple callers (e.g. Mutect2, Pisces, VarDict, Octopus, and TVC)
  - Panel-specific reference models for CNAs (e.g, TSO500, OPA, OCA) and MSI (e.g., TSO500)
  - Longitudinal variant tracking with internal database

- 🩺 **Clinically oriented outputs**:
  - Automated variant annotation and prioritization
  - Interactive, user-friendly HTML reports
  - Quality control summaries and flagging systems

- 🧪 **Validated and benchmarked**:
  - High analytical accuracy in small variant detection with SEQC2 public datasets
  - Real-world performance evaluated in over 2000 clinical tumor samples

- 🆓 **Open and accessible**:
  - Freely available for research and academic use
  - Fully documented and customizable
