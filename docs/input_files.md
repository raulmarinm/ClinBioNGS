# 📄 Input Files

## Sample sheet

A CSV file with sample identifiers. There are some templates in [`resources/sampleSheets/`](../resources/sampleSheets/).

Apart from Illumina BCL or Ion Torrent uBAM analyses that have specific sample sheet formats, this is the common format:

|           |                  |
|:---------:|:----------------:|
|  [Header] |                  |
|  RunName  | ClinBioNGS-RUN01 |
| ...       | ...              |
| [Data]    |                  |
| Sample_ID | Sample_Type      |
| SAMPLE01  | DNA              |
| SAMPLE01  | RNA              |

## Metadata files

### [`SampleInfo.csv`](../resources/SampleInfo.csv)

A CSV file that provides sample-level metadata such as sex, age, tumor type, and estimated tumor purity. Example:

|        Run       |  Sample  | Sex | Age | Purity | Tumor | DOID | Whitelist | LibraryDna | LibraryRna |
|:----------------:|:--------:|:---:|-----|--------|:-----:|:----:|:---------:|------------|------------|
| ClinBioNGS-RUN01 | SAMPLE01 |  ND | ND  | ND     | NSCLC | 3908 |    LUNG   | yes        | no         |

### [`WhitelistGenes.csv`](../resources/WhitelistGenes.csv)

A CSV file that defines tumor-specific or general whitelist genes for prioritization. Example:

|  LUNG | TUMOR |
|:-----:|:-----:|
| GENE1 |       |
|       | GENE2 |

### [`TumorNames.csv`](../resources/TumorNames.csv)

A CSV file that maps user-defined tumor names to DOIDs, top-level ontology nodes, and OncoTree tumor codes, ensuring compatibility with clinical evidence annotations from CIViC. Example:

|           TUMOR_NAME          | TUMOR_DOID | TUMOR_CODE | TOPNODE_NAME | TOPNODE_DOID |
|:-----------------------------:|:----------:|------------|--------------|--------------|
| lung non-small cell carcinoma |    3908    | NSCLC      | lung cancer  | 1324         |

## Panel-specific files

## Data files
