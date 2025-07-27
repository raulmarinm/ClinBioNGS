# 🛠️ Installation

## 1. Prerequisites

- Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ v.24.10.1)
- Install [Apptainer](https://apptainer.org/docs/user/latest/quick_start.html#installation) (≥ v.1.4.1)

## 2. Clone ClinBioNGS

```bash
nextflow clone raulmarinm/ClinBioNGS # or using git
cd ClinBioNGS
chmod +x bin/*
```

## 3. Download SIF container images

The Singularity Image Format (SIF) images are expected to be found on the local computer before running the analysis.

```bash
nextflow run main.nf --prepareImages --runName setup
```

This will download the SIF images (8.7 GB) on their specific folder (`resources/singularity`). SIF images can also be downloaded manually and placed in their folder. 

See the [nextflow.config](../nextflow.config) file for consulting or modifiying the tool versions and links. More details and the list of software tools are available in the [Software](software.md) section.

## 4. Prepare pipeline resources

To avoid delays on first run, resources files can be previously generated:

```bash
nextflow run main.nf --resourcesOnly --runName setup
```

This will download general resources (~200 GB). Panel-specific resources (e.g., manifest files) can be generated automatically on the first panel run or also pre-generated with the appropriate configuration.

This step can also be performed during the anaylsis. ClinBioNGS checks internally if each resource exists in the specific folder (`resources/<...>`) and generates it in case of not finding it.

See the [nextflow.config](../nextflow.config) file for consulting or modifiying the resources versions and links. More details and the list of resources are available in the [Resources](resources.md) section.
