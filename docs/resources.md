# 📚 Resources

**ClinBioNGS** relies on a comprehensive set of genomic and annotation resources to perform accurate and reproducible somatic variant analysis. All required files are either:

* Already included in the repository under the `resources/` folder.
* Automatically downloaded and processed internally by the pipeline from trusted external sources.

This design ensures:

* ✅ Reproducibility across systems
* ⚙️ Minimal manual intervention
* 📦 Standardized analysis across computing environments

The selection of resources prioritizes:

* Open-source availability
* Broad accessibility
* Active maintenance
* Community-wide adoption

See the [nextflow.config](../nextflow.config) file to review or customize the versions and download links of each resource.

## 🧾 Summary of Resource Types

* User-defined metadata files
* Reference genomes and genome resources
* MANE annotation files
* Target region files
* VCF headers
* Gene role and oncogenicity resources
* VEP-related resources
* Cancer hotspot resources
* Problematic and high-confidence regions
* GENIE cancer registry
* Clinical evidence files (CIViC)
* RNA-specific resources
* Panel-recurrent small variants (TSO500, OPA, OCA)
* Panel-specific CNA baselines (TSO500, OPA, OCA)
* Panel-specific MSI baselines (TSO500)

## 📦 Resources Table

| Category                                | Resource                                 | Version     | Role in ClinBioNGS                                            |
| --------------------------------------- | ---------------------------------------- | ----------- | ------------------------------------------------------------- |
| User-defined metadata files             | Sample sheet                             | -           | Specify DNA and RNA samples                                   |
|                                         | *SampleInfo.csv*                         | -           | Provide sample metadata (e.g., sex, tumor type, purity, DOID) |
|                                         | *TumorNames.csv*                         | -           | Define specific tumor names and associated DOIDs              |
|                                         | *WhitelistGenes.csv*                     | -           | Define tumor-specific gene lists for prioritization           |
| Genome resources                        | Human reference genome                   | GRCh38      | Reference genome for DNA analysis                             |
|                                         | BWA-MEM2 indexed genome                  | -           | DNA alignment with BWA-MEM2                                   |
|                                         | TMAP indexed genome                      | -           | DNA alignment with TMAP                                       |
|                                         | Mouse reference genome                   | GRCm38      | Mouse read filtering (e.g., PDX)                              |
|                                         | Xengsort indexed genome                  | -           | Support for Xengsort-based mouse filtering                    |
|                                         | UCSC cytoband                            | hg38        | Annotate gene cytobands                                       |
|                                         | UCSC arm coordinates                     | hg38        | Arm-level CNA annotations                                     |
|                                         | UCSC liftover chain files                | -           | Convert coordinates (hg19/hg38)                               |
| MANE annotation files                   | MANE GTF                                 | 1.4         | Transcript structure definition                               |
|                                         | MANE summary and regions                 | -           | Region annotation and transcript linking                      |
| Target region files                     | Panel manifests                          | hg19        | Define panel capture regions                                  |
|                                         | Target BED files (4-column)              | hg38        | Define padded analysis regions                                |
|                                         | Target genes (HGNC)                      | -           | Annotate gene content per panel                               |
|                                         | Off-target BED                           | hg38        | Support BAM clipping in amplicon panels                       |
|                                         | Target chromosomes                       | -           | Enable per-chromosome variant calling                         |
| VCF headers                             | All header templates                     | 4.2         | Standardize output VCF formatting                             |
| Gene role and oncogenicity              | Network of Cancer Genes (NCG)            | 1.7         | Annotate oncogenes and TSGs                                   |
|                                         | Catalog of Validated Oncogenic Mutations | 20180130    | Validated oncogenic variants                                  |
|                                         | CIViC oncogenic evidence                 | 01-Nov-2024 | Support for oncogenic classification                          |
|                                         | GENIE oncogenic mutations                | 16.1        | Previously classified oncogenic variants                      |
|                                         | ClinGen/CGC/VICC SOP                     | -           | Previously classified oncogenic variants                      |
| VEP-related resources                   | VEP cache                                | 113         | Variant annotation                                            |
|                                         | gnomAD                                   | 4.1         | Population frequencies                                        |
|                                         | CADD                                     | 1.7         | Pathogenicity prediction                                      |
|                                         | REVEL                                    | 1.3         | Pathogenicity prediction                                      |
|                                         | AlphaMissense                            | hg38        | Pathogenicity prediction                                      |
|                                         | ClinVar                                  | 20241103    | Clinical annotations                                          |
|                                         | CIViC VCF                                | 01-Nov-2024 | Clinical evidence annotations                                 |
| Cancer hotspots                         | Panel-specific hotspot BED               | hg38        | Define user-specific hotspot regions                          |
|                                         | AACR GENIE whitelist BED                 | hg19        | Known somatic hotspots                                        |
|                                         | Cancer Hotspots                          | V2          | Statistically enriched mutations                              |
| Problematic / high-confidence regions   | Panel-specific blacklist BED             | hg38        | Define user-specific problematic regions                      |
|                                         | UCSC problematic regions                 | 20240606    | Difficult regions for variant calling                         |
|                                         | GIAB stratification BED                  | 3.5         | Difficult regions for variant calling                         |
|                                         | CTR regions                              | hg38        | High-confidence callable regions                              |
| GENIE cancer registry                   | GENIE variant data (mutation/CNA/fusion) | 16.1        | Annotate variant recurrence in cancer                         |
| Clinical evidence (CIViC)               | CIViC variant summaries (raw)            | 01-Nov-2024 | List of CIViC variants                                        |
|                                         | CIViC molecular profiles (raw)           | 01-Nov-2024 | Map variants to evidence profiles                             |
|                                         | CIViC clinical evidence (raw)            | 01-Nov-2024 | Map molecular profiles to clinical evidences                  |
|                                         | CIViC processed evidence                 | 01-Nov-2024 | Tumor-specific curated clinical evidence                      |
| RNA resources                           | CTAT library                             | Oct292023   | Fusion/splicing reference bundles                             |
|                                         | CTAT splicing DB                         | Jun232020   | Annotate cancer-enriched splice events                        |
|                                         | Mitelman database                        | 20241105    | Annotate fusion recurrence in cancer                          |
|                                         | Fusion/splicing whitelists               | -           | Curated list of known fusion/splice variants                  |
| Panel-specific recurrence and baselines | Recurrent variants (TSO500/OPA/OCA)      | 2024XX      | Flag panel-specific recurrent small variants                  |
|                                         | CNA baselines (TSO500/OPA/OCA)           | 2024XX      | CNA analysis reference profiles                               |
|                                         | MSI baseline (TSO500)                    | 20230124    | Panel-specific MSI baseline                                   |
| Other resources                         | MSigDB MMR gene sets                     | 2024.1      | MSI-related annotation support                                |
|                                         | TVC parameters file                      | -           | Ion Torrent caller configuration                              |
|                                         | Callable genome (CNAs)                   | hg38        | Define CNA-accessible regions                                 |
|                                         | CNA problematic regions (GIAB)           | hg38        | Regions with low CNA reliability                              |
|                                         | ASCETS resources files                   | 1.1.2       | Arm-level CNA detection                                       |
|                                         | Microsatellite loci                      | hg38        | MSI locus identification (10–20bp)                            |
