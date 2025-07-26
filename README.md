# ClinBioNGS

![License: Academic Use Only](https://img.shields.io/badge/license-Academic%20Use%20Only-blue)
![Nextflow](https://img.shields.io/badge/nextflow-24.10.1-23aa62?logo=nextflow)

## 🧬 Pipeline Summary

**ClinBioNGS** is an open-source, panel-agnostic bioinformatics pipeline for the analysis of somatic NGS cancer panels using tumor-only DNA and RNA data. Designed for clinical and translational research, ClinBioNGS enables accurate, reproducible, and interpretable detection of a wide range of genomic alterations. Built with Nextflow and containerized environments, the pipeline ensures full portability, modularity, and flexibility across computing platforms.

ClinBioNGS integrates consensus small variant calling, panel-specific CNA and MSI reference models, automated annotation and prioritization modules, internal QC systems, and a variant database for longitudinal tracking. Final results are compiled into interactive, visual HTML reports, facilitating multidisciplinary interpretation. The pipeline has been validated on SEQC2 reference datasets and thousands of real-world clinical samples, showing performance comparable to commercial tools while offering greater transparency and extended biomarker support.

![ClinBioNGS Workflow](docs/images/workflow_metromap.png)

## 📂 Documentation

See the [docs/](docs) folder for full documentation:

- Overview
- Installation
- Configuration
- Input formats
- Output structure
- Test run

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

## 📦 Requirements

- Nextflow ≥ 24.10.1
- Apptainer ≥ 1.4.1
- Linux environment (HPC with SLURM/SGE or local)

## 🔧 Installation

See [nextflow.io](https://www.nextflow.io/docs/latest/install.html) and [apptainer.org](https://apptainer.org/docs/user/latest/quick_start.html#installation) for instructions.

Clone the repository:

```bash
nextflow clone raulmarinm/ClinBioNGS # or using git
cd ClinBioNGS
chmod +x bin/*
```

The Singularity Image Format (SIF) images are expected to be found on the local computer before running the analysis. This will download the SIF images (**8.7Gb**) on their specific folder (`./resources/singularity`).

```bash
nextflow run main.nf --prepareImages --runName setup # Profiles (e.g., sge, slurm) can be included
```

Resources files can be previously generated to avoid delays during the first analysis. This will perform only the preparation of general resources (**~200Gb**). Panel-specific resources (e.g., manifest files) can be generated with the appropiated configuration.

```bash
nextflow run main.nf --resourcesOnly --runName setup

rm -r work # optional: delete intermediate files to save space (resources files and logs are copied to the ClinBioNGS folder)
```

This step can also be performed during the anaylsis. ClinBioNGS checks internally if each resource exists in the specific folder (`./resources/<...>`) and generates it in case of not finding it.

Now the pipeline is ready for launching the analysis.

## 🧪 Example Run

### TSO500 analysis on a SLURM cluster

```bash
nextflow run main.nf -profile slurm,tso500 \
  --runName TSO500_RUN \
  --projectDir /mnt/projects/ClinBioNGS/output \
  --dataDir /mnt/projects/ClinBioNGS/data \
  --startingDataDir /mnt/illumina_runs/TSO500_Run/BclDirectory
```

### Custom panel analysis on an SGE cluster

```bash
nextflow run main.nf -profile sge,custom \
  --runName customPanel_RUN \
  --projectDir /mnt/projects/ClinBioNGS/output \
  --dataDir /mnt/projects/ClinBioNGS/data \
  --startingDataDir /mnt/data/custom_samples \
  --sampleSheet ./resources/sampleSheets/SampleSheet_customPanel_RUN.csv
```

## 📄 License

Licensed for research and academic use only. Commercial use requires prior approval. See [LICENSE.md](LICENSE.md).

## 📣 Citation

