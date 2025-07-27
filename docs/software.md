# 💻 Software

All software tools required by **ClinBioNGS** are encapsulated within Singularity container images. 

This design ensures:

* ✅ Reproducibility across systems
* ⚙️ Minimal manual intervention
* 📦 Standardized analysis across computing environments

The selection of software prioritizes:

* Open-source availability
* Broad accessibility
* Active maintenance
* Community-wide adoption

See the [nextflow.config](../nextflow.config) file to review or customize the versions and download links of each image.

## 🔧 Containerization

Software tools are distributed via container images that are:

* Downloaded directly or built from publicly available sources, including [BioContainers](https://biocontainers.pro/registry), [Galaxy Project](https://depot.galaxyproject.org/singularity/), and [Docker Hub](https://hub.docker.com/)
* Managed internally by the pipeline and used transparently via Apptainer (or Singularity)

### 🐳 Custom Docker Images

To address compatibility and dependency issues in certain environments, the following custom Docker images were created and are publicly available on [Docker Hub](https://hub.docker.com/u/raulmarinm):

* [Pisces](https://hub.docker.com/r/raulmarinm/pisces-5.3.0.0) variant caller
* [Octopus](https://hub.docker.com/r/raulmarinm/octopus-0.7.4) variant caller
* [R environment](https://hub.docker.com/r/raulmarinm/r-4.0.4) with all required packages

## 📦 Software Table

| Tool             | Version | Role in ClinBioNGS                                                                             |
| ---------------- | ------- | ---------------------------------------------------------------------------------------------- |
| Abra2            | 2.24    | Realignment of DNA reads around target regions                                                 |
| ASCETS           | 1.1.2   | Estimate arm-level CNAs from inferred segment copy ratio                                       |
| AWS CLI          | 2.15.32 | Download of reference genomes from the AWS iGenomes repository                                 |
| Bcftools         | 1.16    | VCF processing and filtering                                                                   |
| BCL Convert      | 4.0.3   | Conversion of BCL to FASTQ                                                                     |
| Bedtools         | 2.31.0  | Processing and manipulation of BED files                                                       |
| Bioawk           | 1.0     | Text manipulation for biological data                                                          |
| BWA-MEM2         | 2.2.1   | Alignment of DNA reads to the reference genome                                                 |
| CNVkit           | 0.9.9   | Copy-number variation analysis from DNA reads                                                  |
| CTAT-splicing    | 0.0.2   | Detection and annotation of splicing variants from RNA data                                    |
| Ensembl VEP      | 113     | Annotation of small variants                                                                   |
| FastP            | 0.23.4  | Pre-processing and quality filtering of FASTQ files                                            |
| FastQC           | 0.12.1  | QC of raw and processed FASTQ files                                                            |
| GATK4            | 4.3.0.0 | File manipulation, deduplication (MarkDuplicates), small variant calling (Mutect2), and BAM QC |
| Gencore          | 0.17.2  | UMI-aware deduplication for paired-end reads                                                   |
| Mosdepth         | 0.3.3   | Per-base and region-level read coverage calculation                                            |
| MSIsensor-pro    | 1.2.0   | MSI detection from DNA reads                                                                   |
| MultiQC          | 1.22.3  | Aggregation and visualization of QC metrics                                                    |
| Octopus          | 0.7.4   | Small variant calling from DNA reads                                                           |
| Pisces           | 5.3.0.0 | Small variant calling from DNA reads                                                           |
| R                | 4.0.4   | Data manipulation, statistical analysis, and generation of reports                             |
| Samtools         | 1.18    | Processing and manipulation of BAM and FASTQ files                                             |
| STAR-Fusion      | 1.13.0  | Detection of fusion transcripts from RNA data                                                  |
| TMAP / TVC       | 5.12.1  | Alignment and small variant calling for Ion Torrent DNA reads                                  |
| UCSC bigBedToBed | 377     | Conversion of bigBed files to BED format                                                       |
| UCSC liftOver    | 377     | Conversion of genomic coordinates between genome builds (e.g., hg19 to hg38)                   |
| UMI-tools        | 1.1.2   | UMI-aware deduplication for single-end reads                                                   |
| UMI-transfer     | 1.0.0   | Transfer of UMI info from separate FASTQ files into read headers                               |
| VarDict          | 1.8.3   | Small variant calling from DNA reads                                                           |
| Vt               | 0.57721 | VCF decomposition and normalization of InDels                                                  |
| Xengsort         | 1.1.0   | Filtering of mouse-derived reads in xenograft sequencing data                                  |
