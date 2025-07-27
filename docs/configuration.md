# ⚙️ Configuration

The configuration system in **ClinBioNGS** is modular, flexible, and designed to support both standard and custom setups across a variety of computing environments and sequencing panels.

Parameters can be passed using both the command-line options and config files. The following priority to override default parameters is defined from highest to lowest:

- Parameters specified on the command-line (`--something value`).
- Parameters provided in a custom config file (called as `-c <FILE>`).
- Parameters defined in profile-specific config files.
- Parameters specified in the [`modules.config`](../conf/modules.config).
- Parameters specified in the [`nextflow.config`](../conf/nextflow.config).

## 📁 Config Files

ClinBioNGS uses several layered configuration files:

- [`nextflow.config`](../nextflow.config): Global defaults and profile declarations.
- [`base.config`](../conf/base.config): Computational resources using specific labels.
- [`modules.config`](../conf/modules.config): Process-specific parameters.
- `profile-specific files` in [conf/](../conf/) directory.

## 🧩 Profile-Based Customization

Specific Nextflow parameters can be grouped into configuration files (found in [`conf/`](../conf/)) and can be called as profiles (`-profile <profile>` in the command line).

Profiles allow easy switching between environments and panel-specific settings using dedicated *.config* files :

- **Execution environments**: [`sge`](../conf/sge.config) and [`slurm`](../conf/slurm.config) execute the pipeline on HPC clusters with [Sun Grid Engine (SGE)](http://en.wikipedia.org/wiki/Oracle_Grid_Engine) or [SLURM](https://slurm.schedmd.com/documentation.html) schedulers, respectively.
- **NGS panels**: [`tso500`](../conf/tso500.config), [`opa`](../conf/opa.config), [`oca`](../conf/oca.config), and SEQC2 oncopanels ([`agl`](../conf/agl.config), [`brp`](../conf/brp.config), [`idt`](../conf/idt.config)[`igt`](../conf/igt.config), [`ilm`](../conf/ilm.config), [`tfs`](../conf/tfs.config)) execute the pipeline using the corresponding panel parameters.
- **Custom panels**: [`custom`](../conf/custom.config) profile for using parameters in the corresponding config file.

Multiple types of profiles can be defined by separating them with a comma (e.g., `-profile slurm,tso500`).

## 🚀 Initial Parameters

To launch the pipeline, the following parameters must be set (defaults are in [`nextflow.config`](../nextflow.config)):

### `--projectDir`

Output directory where results will be saved. Must be and ABSOLUTE path. ClinBioNGS pipeline directory by default.

### `--dataDir`

Storage location for processed data (FASTQ, BAM, or VCF). Must be and ABSOLUTE path. ClinBioNGS pipeline directory by default.

### `--runName`

Unique identifier for the analysis (used to name the output `--projectDir` and `--dataDir` directories, and the `--sampleSheet`).

### `--startingDataDir`

Path to the starting data files. ClinBioNGS pipeline directory by default. Supported formats:

- **FASTQ or BAM files**:
  - Files named `<sample>_<DNA/RNA>*.fastq*` or `<sample>_<DNA/RNA>*.bam`.
  - `<sample>` must correspond to the name defined in `--sampleSheet`.
  - Symbolic links are allowed.
  - Use `--startingDataType FASTQ` (default) or `--startingDataType BAM`.

- **Illumina BCL folder**:
  - Directory with raw BCL files and `SampleSheet.csv`.
  - Use `--startingDataType BCL`.

- **Ion Torrent results folder**:
  - Output folder of Ion Torrent platform with `Final_Results_Files/<...>.tar.xz` directory of results for each sample.
  - It contains the starting uBAM (`<...>_rawlib.basecaller.bam`) and other files for preparing the initial resources.
  - Use `--startingDataType BAM` and set `--prepareIontorrentBam`.

### `--sampleSheet`

Describes sample identifiers. By default, it is expected to be located in [`resources/sampleSheets/`](../resources/sampleSheets/) named by `SampleSheet_<runName>.csv`. Otherwise, it is searched in `--startingDataDir` and copied to the sample sheets directory with the default naming for being reused in future analyses.

File format described in [Input files](docs/input_files.md) section.

Other alternative locations:

- **Manual input** defined in `--sampleSheet`:

- **Illumina BCL folder**:
  - Use `SampleSheet.csv` from the run folder or copy it to default path following naming conventions.

- **Ion Torrent results folder**:
  - Can be auto-generated from `Info.csv` in `--startingDataDir`.
  - Use `--prepareIontorrentSamplesheet` (preconfigured in Ion Torrent panel profiles).

### `--manifestDir`, `--dnaManifest`, `--rnaManifest`

Manifest directory and the corresponding DNA and RNA input manifest files.

- The [profiled](#-profile-based-customization) NGS panels have their own directory in [`resources/manifests/`](../resources/manifests/).

- For non-profiled, they must be prepared by the user (see [Input files](docs/input_files.md) section).

### Other parameters

Additional parameters can be consulted in [`nextflow.config`](../nextflow.config) and [`modules.config`](../conf/modules.config) files.

## 🧠 Resource Allocation System

ClinBioNGS uses a tag-based system for allocating resources to each process.

- Tags: `min`, `low`, `med`, `high`, `extra`.
- Control CPU, memory, and time limits per process
- If a task exceeds its limits, it is automatically retried with a higher resource level.
- Each process' tag is defined in the corresponding script in [modules/process](../modules/process/) folder.

Default labels:

|   Label  |  Value  |   |   |   |   |   Label  | Value |
|:--------:|:-------:|:-:|---|---|:-:|:--------:|:-----:|
|  minCpu  |  1 CPU  |   |   |   |   |  minMem  |  6 GB |
|  lowCpu  |  4 CPUs |   |   |   |   |  lowMem  | 12 GB |
|  medCpu  |  8 CPUs |   |   |   |   |  medMem  | 24 GB |
|  highCpu | 16 CPUs |   |   |   |   |  highMem | 48 GB |
| extraCpu | 32 CPUs |   |   |   |   | extraMem | 96 GB |

Default maximum values:

```text
params {
  maxCpus   = 32
  maxMemory = 200.GB
  maxTime   = 24.h
}
```

Both labels and maximum values can be modified by the user in the [`base.config`](../conf/base.config) file.

---

For additional details, see configuration examples provided in the [`conf/`](../conf/) directory.
