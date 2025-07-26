#!/usr/bin/env nextflow 


/*
===========================================================
					ClinBioNGS pipeline
===========================================================
An integrated clinical bioinformatics pipeline 
for the analysis of somatic NGS cancer panels.
-----------------------------------------------------------
*/

nextflow.enable.dsl=2


/*
========================================
 VALIDATE INPUTS
========================================
*/

// Check required parameters
def checkPathParamList = [
	params.codeDir,
	params.dataDir,
	params.projectDir,
]
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)


/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/
params.runDir                   = "${params.projectDir}/${params.runName}"
params.dataRunDir               = "${params.dataDir}/${params.runName}/Data"
params.logsDir                  = "${params.runDir}/Logs"
params.analysisDir              = "${params.runDir}/Analysis"
params.resultsDir               = "${params.runDir}/Results"
params.dnaTargetBed             = "${params.manifestDir}/${params.seqPanel}_DNA_${params.genomeVersionHg}.bed"
params.dnaOffTargetBed          = "${params.manifestDir}/${params.seqPanel}_DNA_${params.genomeVersionHg}_complement.bed"
params.dnaTargetIntervalList    = "${params.manifestDir}/${params.seqPanel}_DNA_${params.genomeVersionHg}.interval_list"
params.extBases                 = params.variantIntervalPadding && params.variantIntervalPadding != '0' ? 
									"_ext${params.variantIntervalPadding}" : ''
params.dnaTargetBedExt          = "${params.manifestDir}/${params.seqPanel}_DNA_${params.genomeVersionHg}${params.extBases}.bed"
params.dnaTargetIntervalListExt = "${params.manifestDir}/${params.seqPanel}_DNA_${params.genomeVersionHg}${params.extBases}.interval_list"
params.dnaTargetGenes           = "${params.manifestDir}/${params.seqPanel}_DNA_genes.txt"
params.dnaTargetChr             = "${params.manifestDir}/${params.seqPanel}_DNA_chromosomes.txt"
params.rnaTargetBed             = "${params.manifestDir}/${params.seqPanel}_RNA_${params.genomeVersionHg}.bed"
params.rnaTargetGenes           = "${params.manifestDir}/${params.seqPanel}_RNA_genes.txt"


/*
========================================
 IMPORT SUBWORKFLOWS & PROCESSES
========================================
*/

// Workflows
include { PREPARE_TOOL_IMAGES                  } from './modules/subworkflow/prepareImages'
include { RESOURCES                            } from './modules/subworkflow/prepareResources'
include { 
	FASTQ_GENERATION as FASTQ_GENERATION_DNA; 
	FASTQ_GENERATION as FASTQ_GENERATION_RNA   } from './modules/subworkflow/fastqGeneration'
include { 
	FASTQ_PROCESSING as FASTQ_PROCESSING_DNA; 
	FASTQ_PROCESSING as FASTQ_PROCESSING_RNA   } from './modules/subworkflow/fastqProcessing'
include { 
	FASTQ_QC as FASTQ_QC_DNA; 
	FASTQ_QC as FASTQ_QC_RNA                   } from './modules/subworkflow/fastqQC'
include { ALIGNMENT_DNA                        } from './modules/subworkflow/alignmentDna'
include { ALIGNMENT_RNA                        } from './modules/subworkflow/alignmentRna'
include { 
	DEDUPLICATION as DEDUPLICATION_DNA; 
	DEDUPLICATION as DEDUPLICATION_RNA         } from './modules/subworkflow/deduplication'
include { REALIGNMENT_DNA                      } from './modules/subworkflow/realignmentDna'
include { REALIGNMENT_RNA                      } from './modules/subworkflow/realignmentRna'
include { VARIANT_CALLING                      } from './modules/subworkflow/variantCalling'
include { VARIANT_ANNOTATION                   } from './modules/subworkflow/variantAnnotation'
include { CNA                                  } from './modules/subworkflow/cna'
include { MSI                                  } from './modules/subworkflow/msi'
include { FUSION                               } from './modules/subworkflow/fusion'
include { SPLICING                             } from './modules/subworkflow/splicing'

// Processes
include { 
	PREPARE_SAMPLE_FILES as PREPARE_SAMPLE_FILES_DNA; 
	PREPARE_SAMPLE_FILES as PREPARE_SAMPLE_FILES_RNA } from './modules/process/prepareSampleFiles'
include { 
	MULTIQC_BCL_CONVERT; 
	MULTIQC_RUN as MULTIQC_RUN_DNA;
	MULTIQC_RUN as MULTIQC_RUN_RNA;                  } from './modules/process/multiQC'
include { TMB                                        } from './modules/process/tmb'
include { RUN_LIBRARY                                } from './modules/process/runLibrary'
include { UNIQUE_LIBRARY                             } from './modules/process/uniqueLibrary'
include { SAMPLE_REPORT                              } from './modules/process/sampleReport'


/*
========================================
 BUILD CHANNELS
========================================
*/

// Directories
chStartingDataDir = Channel.fromPath(params.startingDataDir, checkIfExists: true).collect()

// Tumor names
chTumorNamesFile = Channel.fromPath(params.tumorNamesFile, checkIfExists: true).collect()

// Sample info
chClinicalInfoFile = Channel.fromPath(params.clinicalInfoFile, checkIfExists: true).collect()

// Whitelist genes
chWhitelistGenes = Channel.fromPath(params.whitelistGenesCsv, checkIfExists: true).collect()

// multiQC config
chMultiqcConfig = params.multiqcConfig && file(params.multiqcConfig).exists() ? Channel.fromPath(params.multiqcConfig).collect() : []

// Variant libraries
chUniqueLibrary = file("${params.projectDir}/${params.libraryName}*.db").any{file(it).exists()} ? 
	Channel.fromPath("${params.projectDir}/${params.libraryName}*.db").collect(sort: true).reverse().flatten().take(1).collect() : []
chRunLibrary    = file("${params.resultsDir}/*${params.runName}*.db").any{file(it).exists()} ? 
	Channel.fromPath("${params.resultsDir}/*${params.runName}*.db").collect(sort: true).reverse().flatten().take(1).collect() : []

// Annotation resources
chMmrGenes          = params.mmrGenesFile && file(params.mmrGenesFile).exists() ? 
	Channel.fromPath(params.mmrGenesFile).collect() : []
chSopOncogenic      = params.sopOncogenicFile && file(params.sopOncogenicFile).exists() ? 
	Channel.fromPath(params.sopOncogenicFile).collect() : []
chGenieMutCounts    = params.genieMutCountsFile && file(params.genieMutCountsFile).exists() ? 
	Channel.fromPath(params.genieMutCountsFile).collect() : []
chGenieMutOncogenic = params.genieMutOncogenicFile && file(params.genieMutOncogenicFile).exists() ? 
	Channel.fromPath(params.genieMutOncogenicFile).collect() : []
chGenieCna          = params.genieCnaFile && file(params.genieCnaFile).exists() ? 
	Channel.fromPath(params.genieCnaFile).collect() : []
chGenieFusion       = params.genieFusionFile && file(params.genieFusionFile).exists() ? 
	Channel.fromPath(params.genieFusionFile).collect() : []

// SmallVariant files
chTorrentVariantCallerJson = file(params.torrentVariantCallerJson).exists() ? 
	Channel.fromPath(params.torrentVariantCallerJson).collect() : []
chVariantConsensusVcfHeader = file(params.variantConsensusVcfHeader).exists() ? 
	Channel.fromPath(params.variantConsensusVcfHeader).collect() : Channel.empty()
chVariantAnnotVcfHeader     = file(params.variantAnnotationVcfHeader).exists() ? 
	Channel.fromPath(params.variantAnnotationVcfHeader).collect() : Channel.empty()
chPanelRecurrentMutations   = params.variantRecurrentFile && file(params.variantRecurrentFile).exists() ? 
	Channel.fromPath(params.variantRecurrentFile).collect() : []

// CNA files
chCnvkitBaseline      = params.cnvkitBaseline && file(params.cnvkitBaseline).exists() ? 
	Channel.fromPath(params.cnvkitBaseline).collect() : Channel.empty()
chCnaVcfHeader        = params.cnaVcfHeader && file(params.cnaVcfHeader).exists() ? 
	Channel.fromPath(params.cnaVcfHeader).collect() : Channel.empty()

// MSI files
chMsiBaseline = params.msiBaseline && file(params.msiBaseline).exists() ? Channel.fromPath(params.msiBaseline).collect() : Channel.empty()

// Fusion files
chFusionVcfHeader = params.fusionVcfHeader && file(params.fusionVcfHeader).exists() ? Channel.fromPath(params.fusionVcfHeader).collect() : Channel.empty()
chFusionWhitelist = params.fusionWhitelistFile && file(params.fusionWhitelistFile).exists() ? 
	Channel.fromPath(params.fusionWhitelistFile).collect() : []

// Splicing files
chSplicingVcfHeader = params.splicingVcfHeader && file(params.splicingVcfHeader).exists() ? Channel.fromPath(params.splicingVcfHeader).collect() : Channel.empty()
chSplicingWhitelist = params.splicingWhitelistFile && file(params.splicingWhitelistFile).exists() ? 
	Channel.fromPath(params.splicingWhitelistFile).collect() : []

// FASTQ files
chDnaFastq = Channel.fromPath("${params.dataRunDir}/**_DNA*.fastq*")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chRnaFastq = Channel.fromPath("${params.dataRunDir}/**_RNA*.fastq*")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()

// BAM files
chDnaBamAligned = Channel.fromPath("${params.dataRunDir}/**_DNA${params.bamAligned_suffix}.bam")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaBaiAligned = Channel.fromPath("${params.dataRunDir}/**_DNA${params.bamAligned_suffix}*.bai")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaBamRealigned = Channel.fromPath("${params.dataRunDir}/**_DNA${params.bamRealigned_suffix}.bam")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaBaiRealigned = Channel.fromPath("${params.dataRunDir}/**_DNA${params.bamRealigned_suffix}*.bai")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chRnaBamAligned = Channel.fromPath("${params.dataRunDir}/**_RNA${params.bamAligned_suffix}.bam")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaBaiAligned = Channel.fromPath("${params.dataRunDir}/**_RNA${params.bamAligned_suffix}*.bai")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaBamRealigned = Channel.fromPath("${params.dataRunDir}/**_RNA${params.bamRealigned_suffix}.bam")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaBaiRealigned = Channel.fromPath("${params.dataRunDir}/**_RNA${params.bamRealigned_suffix}*.bai")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()					

// FastqGeneration files
chDnaBclConvertReport = Channel.fromPath("${params.analysisDir}/**BCL_CONVERT_DNA/Reports/*").collect()
chRnaBclConvertReport = Channel.fromPath("${params.analysisDir}/**BCL_CONVERT_RNA/Reports/*").collect()
chDnaBamHeaderInfo    = []

// FastqProcessing files
chDnaFastqProcessingFiles = Channel.fromPath("${params.analysisDir}/**FASTQ_PROCESSING/*_DNA*.{html,json}").collect()
chRnaFastqProcessingFiles = Channel.fromPath("${params.analysisDir}/**FASTQ_PROCESSING/*_RNA*.{html,json}").collect()

// FastQC files
chDnaFastqcFiles = Channel.fromPath("${params.analysisDir}/**FastQC_Reports/*_DNA*.{html,zip}").collect()
chRnaFastqcFiles = Channel.fromPath("${params.analysisDir}/**FastQC_Reports/*_RNA*.{html,zip}").collect()

// RNA alignment results
chRnaRealignmentSJ     = Channel.fromPath("${params.analysisDir}/**REALIGNMENT/*_RNA*SJ.out.tab")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaRealignmentJunction = Channel.fromPath("${params.analysisDir}/**REALIGNMENT/*_RNA*.out.junction")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()

// Deduplication files
chDnaDedupFiles = Channel.fromPath("${params.analysisDir}/**DEDUPLICATION/*_DNA*.{html,json,logs,tsv,txt}").collect()
chRnaDedupFiles = Channel.fromPath("${params.analysisDir}/**DEDUPLICATION/*_RNA*.{html,json,logs,tsv,txt}").collect()

// BAM QC results
chDnaBamQcRawResults   = Channel.fromPath("${params.analysisDir}/**_DNA_RawBamQC.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaPerBaseCovBed     = Channel.fromPath("${params.analysisDir}/**_DNA*per-base.bed.gz")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaBamQcFinalResults = Channel.fromPath("${params.analysisDir}/**_DNA_BamQC.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chDnaBamQcFiles        = Channel.fromPath("${params.analysisDir}/**_DNA_{BamQC,CoverageBy,GlobalCoverage}*.{png,txt}")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chRnaBamQcRawResults   = Channel.fromPath("${params.analysisDir}/**_RNA_RawBamQC.txt")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaPerBaseCovBed     = Channel.fromPath("${params.analysisDir}/**_RNA*per-base.bed.gz")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaBamQcFinalResults = Channel.fromPath("${params.analysisDir}/**_RNA_BamQC.txt")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chRnaBamQcFiles        = Channel.fromPath("${params.analysisDir}/**_RNA_{BamQC,CoverageBy,GlobalCoverage}*.{png,txt}")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()					
chDnaBamQcRunFiles     = Channel.fromPath("${params.analysisDir}/**REALIGNMENT/BamQC_Reports/*_DNA*.{gz,txt}").collect()
chRnaBamQcRunFiles     = Channel.fromPath("${params.analysisDir}/**REALIGNMENT/BamQC_Reports/*_RNA*.{gz,txt}").collect()

// Small variant results
chVariantConsensusVcf     = Channel.fromPath("${params.dataRunDir}/**_DNA_SmallVariant_consensus.vcf")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chVariantConsensusResults = Channel.fromPath("${params.analysisDir}/**_DNA_SmallVariant_consensus.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chVariantAnnotResults = Channel.fromPath("${params.analysisDir}/**_DNA_SmallVariant_annot.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chVariantAnnotFiles = Channel.fromPath("${params.analysisDir}/**_DNA_SmallVariant*.{png,txt}")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()

// TMB results
chTmbResults = Channel.fromPath("${params.analysisDir}/**_DNA_TMB_results.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chTmbTrace = Channel.fromPath("${params.analysisDir}/**_DNA_TMB_trace.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()

// CNA results
chCnaGeneResults = Channel.fromPath("${params.analysisDir}/**_DNA_CNA_gene.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chCnaArmResults = Channel.fromPath("${params.analysisDir}/**_DNA_CNA_arm.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()
chCnaFiles = Channel.fromPath("${params.analysisDir}/**_DNA_CNA*.{png,txt}")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()

// MSI results
chMsiResults = Channel.fromPath("${params.analysisDir}/**_DNA_MSI_results.txt")
	.map{ [it.simpleName.toString().split("_DNA")[0], it] }.groupTuple()

// Fusion results
chFusionResults = Channel.fromPath("${params.analysisDir}/**_RNA_Fusion_results.txt")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chFusionFiles   = Channel.fromPath("${params.analysisDir}/**_RNA_Fusion_*.{png,txt}")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()

// Splicing results
chSplicingResults = Channel.fromPath("${params.analysisDir}/**_RNA_Splicing_results.txt")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()
chSplicingFiles   = Channel.fromPath("${params.analysisDir}/**_RNA_Splicing_*.{png,txt}")
	.map{ [it.simpleName.toString().split("_RNA")[0], it] }.groupTuple()

// Sample files
chDnaSampleFile  = []
chDnaSampleSheet = Channel.empty()
chRnaSampleFile  = []
chRnaSampleSheet = Channel.empty()


/*
========================================
 RUN WORKFLOW
========================================
*/

workflow {

	if (params.prepareImages) {

		// Prepare singularity images (only this subworkflow is executed)
		PREPARE_TOOL_IMAGES()

	} else {

		// Prepare resources files
		RESOURCES(chStartingDataDir, chTumorNamesFile)
		
		// MANE resources
		chManeGtf                  = RESOURCES.out.maneGtf
		chManeGene                 = RESOURCES.out.maneGene
		chManeExon                 = RESOURCES.out.maneExon
		chManeIntron               = RESOURCES.out.maneIntron
		chManeCoding               = RESOURCES.out.maneCoding

		// DNA resources
		chDnaTargetBed             = RESOURCES.out.dnaBed
		chDnaTargetIntervalList    = RESOURCES.out.dnaIntervalList
		chDnaTargetBedExt          = RESOURCES.out.dnaBedExt
		chDnaTargetIntervalListExt = RESOURCES.out.dnaIntervalListExt
		chDnaTargetGenes           = RESOURCES.out.dnaGenes
		chDnaTargetChr             = RESOURCES.out.dnaChr
		chDnaGenomeFastaDir        = RESOURCES.out.dnaGenomeFastaDir
		chDnaGenomeFastaFiles      = RESOURCES.out.dnaGenomeFastaFiles
		chDnaGenomeIndexFiles      = RESOURCES.out.dnaGenomeIndexFiles
		chDnaPanelHotspotsBed      = RESOURCES.out.dnaPanelHotspotsBed
		chDnaPanelBlacklistBed     = RESOURCES.out.dnaPanelBlacklistBed
		chDnaOffTargetBed          = RESOURCES.out.dnaOffTargetBed
		chVepCacheDir              = RESOURCES.out.vepCacheDir
		chVepCaddSnv	           = RESOURCES.out.vepCaddSnv
		chVepCaddIndel	           = RESOURCES.out.vepCaddIndel
		chVepRevel	               = RESOURCES.out.vepRevel
		chVepAlphaMissense	       = RESOURCES.out.vepAlphaMissense
		chClinvarVcf	           = RESOURCES.out.clinvarVcf
		chCivicVcf	               = RESOURCES.out.civicVcf
		chHotspotsWhitelistBed	   = RESOURCES.out.hotspotsWhitelistBed
		chCancerhotspotsResults    = RESOURCES.out.cancerhotspotsResults
		chCgiOncogenic	           = RESOURCES.out.cgiOncogenic
		chProblematicRegionsBed    = RESOURCES.out.problematicRegionsBed
		chCtrRegionsBed            = RESOURCES.out.ctrRegionsBed

		chArmCoordinates           = RESOURCES.out.armCoordinates
		chAscetsResources          = RESOURCES.out.ascetsResources

		// RNA resources
		chRnaTargetBed             = RESOURCES.out.rnaBed
		chRnaTargetGenes           = RESOURCES.out.rnaGenes
		chRnaGenomeFastaFiles      = RESOURCES.out.rnaGenomeFastaFiles
		chRnaGenomeIndexDir        = RESOURCES.out.rnaGenomeIndexDir
		chCtatLibDir               = RESOURCES.out.ctatLibDir
		chMitelmanFusion           = RESOURCES.out.mitelmanFusion

		// Other resources
		chHg38ToHg19ChainFile      = RESOURCES.out.hg38ToHg19ChainFile
		chCancerdrivers            = RESOURCES.out.cancerdrivers
		chCivicOncogenic	       = RESOURCES.out.civicOncogenic
		chCivicClinical	           = RESOURCES.out.civicClinical
		chMouseAndHumanIndex	   = RESOURCES.out.mouseAndHumanIndex

		// Sample sheet
		chSampleSheet              = RESOURCES.out.sampleSheet
		
		// Generate only the resources or perform the analysis
		if (!params.resourcesOnly) {	

			// DNA analysis
			if (params.DNA) {

				// Prepare sample files
				PREPARE_SAMPLE_FILES_DNA(chSampleSheet, chClinicalInfoFile, chTumorNamesFile, 'DNA')
				chDnaSampleFile  = PREPARE_SAMPLE_FILES_DNA.out.txt
				chDnaSampleSheet = PREPARE_SAMPLE_FILES_DNA.out.csv

				// DNA samples 
				chDnaSamples              = chDnaSampleFile.splitCsv(sep: '\t', header: true).map(row -> "${row.Sample}")
				chDnaReads                = chDnaSamples.map{ [it] }.combine(chDnaFastq, by: 0)
				chDnaBamAligned           = chDnaSamples.map{ [it] }.combine(chDnaBamAligned, by: 0)
				chDnaBaiAligned           = chDnaSamples.map{ [it] }.combine(chDnaBaiAligned, by: 0)
				chDnaBamRealigned         = chDnaSamples.map{ [it] }.combine(chDnaBamRealigned, by: 0)
				chDnaBaiRealigned         = chDnaSamples.map{ [it] }.combine(chDnaBaiRealigned, by: 0)
				chDnaBamQcRawResults      = chDnaSamples.map{ [it] }.combine(chDnaBamQcRawResults, by: 0)
				chDnaPerBaseCovBed        = chDnaSamples.map{ [it] }.combine(chDnaPerBaseCovBed, by: 0)
				chDnaBamQcFinalResults    = chDnaSamples.map{ [it] }.combine(chDnaBamQcFinalResults, by: 0)
				chDnaBamQcFiles           = chDnaSamples.map{ [it] }.combine(chDnaBamQcFiles, by: 0)
				chVariantConsensusVcf     = chDnaSamples.map{ [it] }.combine(chVariantConsensusVcf, by: 0)
				chVariantConsensusResults = chDnaSamples.map{ [it] }.combine(chVariantConsensusResults, by: 0)
				chVariantAnnotResults     = chDnaSamples.map{ [it] }.combine(chVariantAnnotResults, by: 0)
				chVariantAnnotFiles       = chDnaSamples.map{ [it] }.combine(chVariantAnnotFiles, by: 0)
				chTmbResults              = chDnaSamples.map{ [it] }.combine(chTmbResults, by: 0)
				chTmbTrace                = chDnaSamples.map{ [it] }.combine(chTmbTrace, by: 0)
				chCnaGeneResults          = chDnaSamples.map{ [it] }.combine(chCnaGeneResults, by: 0)
				chCnaArmResults           = chDnaSamples.map{ [it] }.combine(chCnaArmResults, by: 0)
				chCnaFiles                = chDnaSamples.map{ [it] }.combine(chCnaFiles, by: 0)
				chMsiResults              = chDnaSamples.map{ [it] }.combine(chMsiResults, by: 0)

				// FASTQ generation
				if (params.fastqGeneration) {
					FASTQ_GENERATION_DNA(chStartingDataDir, chDnaSampleFile, chDnaSampleSheet, 'DNA')
					chDnaReads            = chDnaSamples.map{ [it] }.combine(FASTQ_GENERATION_DNA.out.fastq, by: 0)
					chDnaBclConvertReport = FASTQ_GENERATION_DNA.out.bclConvertReport
					chDnaBamHeaderInfo    = FASTQ_GENERATION_DNA.out.bamHeaderInfo
				}

				// FASTQ processing
				if (params.fastqProcessing) {
					FASTQ_PROCESSING_DNA(chDnaReads, chDnaSampleFile, chMouseAndHumanIndex, 'DNA')
					chDnaReads = chDnaSamples.map{ [it] }.combine(FASTQ_PROCESSING_DNA.out.fastq, by: 0)
					chDnaFastqProcessingFiles = FASTQ_PROCESSING_DNA.out.files
				}

				// FASTQ QC
				if (params.fastqQC) { 
					FASTQ_QC_DNA(chDnaReads, chMultiqcConfig, 'DNA')
					chDnaFastqcFiles = FASTQ_QC_DNA.out.files
				}

				// ALIGNMENT + DEDUPLICATION
				if (params.deduplicationDna) {

					if (params.alignment) {
						ALIGNMENT_DNA(chDnaReads, chDnaGenomeFastaFiles, chDnaGenomeIndexFiles, chDnaTargetBed, chDnaBamHeaderInfo, chMultiqcConfig)
						chDnaBamAligned      = chDnaSamples.map{ [it] }.combine(ALIGNMENT_DNA.out.bam, by: 0)
						chDnaBaiAligned      = chDnaSamples.map{ [it] }.combine(ALIGNMENT_DNA.out.bai, by: 0)
						chDnaBamQcRawResults = chDnaSamples.map{ [it] }.combine(ALIGNMENT_DNA.out.results, by: 0)
					}

					DEDUPLICATION_DNA(chDnaBamAligned, chDnaBaiAligned, chDnaGenomeFastaFiles, chDnaTargetBed, chMultiqcConfig, 'DNA')
					chDnaReads      = chDnaSamples.map{ [it] }.combine(DEDUPLICATION_DNA.out.bam, by: 0)
					chDnaDedupFiles = DEDUPLICATION_DNA.out.files
				}

				// REALIGNMENT
				if (params.realignment) {
					chDnaRealignBam = params.deduplicationDna ? true: false
					REALIGNMENT_DNA(
						chDnaReads, chDnaRealignBam, chDnaGenomeFastaFiles, chDnaGenomeIndexFiles, chDnaTargetBed, 
						chDnaBamHeaderInfo, chMultiqcConfig, chDnaSampleFile, chDnaTargetGenes, chWhitelistGenes, 
						chManeGtf, chManeGene, chManeExon, chManeCoding, chDnaBamQcRawResults
					)
					chDnaBamRealigned      = chDnaSamples.map{ [it] }.combine(REALIGNMENT_DNA.out.bam, by: 0)
					chDnaBaiRealigned      = chDnaSamples.map{ [it] }.combine(REALIGNMENT_DNA.out.bai, by: 0)
					chDnaPerBaseCovBed     = chDnaSamples.map{ [it] }.combine(REALIGNMENT_DNA.out.perBaseCovBed, by: 0)
					chDnaBamQcFinalResults = chDnaSamples.map{ [it] }.combine(REALIGNMENT_DNA.out.results, by: 0)
					chDnaBamQcFiles        = chDnaSamples.map{ [it] }.combine(REALIGNMENT_DNA.out.files, by: 0)
					chDnaBamQcRunFiles     = REALIGNMENT_DNA.out.runFiles
				}

				// SMALL VARIANT: CALLING
				if (params.smallVariantCalling) {
					VARIANT_CALLING(
						chDnaBamRealigned.join(chDnaBaiRealigned), chDnaGenomeFastaFiles, chDnaGenomeFastaDir, 
						chDnaTargetBed, chDnaTargetBedExt, chDnaTargetIntervalListExt, chDnaOffTargetBed, 
						chDnaTargetChr, chTorrentVariantCallerJson, chHg38ToHg19ChainFile, chDnaPanelHotspotsBed, 
						chDnaPanelBlacklistBed, chVariantConsensusVcfHeader
					)
					chVariantConsensusVcf     = chDnaSamples.map{ [it] }.combine(VARIANT_CALLING.out.vcf, by: 0)
					chVariantConsensusResults = chDnaSamples.map{ [it] }.combine(VARIANT_CALLING.out.results, by: 0)
				}

				// SMALL VARIANT: ANNOTATION
				if (params.smallVariantAnnotation) {
					VARIANT_ANNOTATION(
						chVariantConsensusVcf.flatMap{ [it[1]] }.collect(), chVariantConsensusResults,
						chVepCacheDir, chVepCaddSnv, chVepCaddIndel, chVepRevel, chVepAlphaMissense, 
						chClinvarVcf, chCivicVcf, chHotspotsWhitelistBed, chProblematicRegionsBed, 
						chCtrRegionsBed, chMmrGenes, chCancerdrivers, chCancerhotspotsResults, chSopOncogenic, 
						chGenieMutCounts, chGenieMutOncogenic, chCgiOncogenic, chCivicOncogenic, 
						chCivicClinical, chPanelRecurrentMutations, chDnaPanelHotspotsBed, chDnaTargetGenes, 
						chWhitelistGenes, chManeGtf, chManeGene, chManeExon, chVariantAnnotVcfHeader, chDnaSampleFile
					)
					chVariantAnnotResults = chDnaSamples.map{ [it] }.combine(VARIANT_ANNOTATION.out.results, by: 0)
					chVariantAnnotFiles   = chDnaSamples.map{ [it] }.combine(VARIANT_ANNOTATION.out.files, by: 0)
				}

				// TMB
				if (params.tmb) {
					TMB(
						chVariantAnnotResults.join(chDnaPerBaseCovBed, by: 0), chDnaTargetBed, 
						chManeCoding, chProblematicRegionsBed, '08_TMB')
					chTmbResults = chDnaSamples.map{ [it] }.combine(TMB.out.results, by: 0)
					chTmbTrace   = chDnaSamples.map{ [it] }.combine(TMB.out.trace, by: 0)
				}

				// CNA
				if (params.cna) {
					CNA(
						chDnaBamRealigned.join(chDnaBaiRealigned), chDnaTargetBed, chCnvkitBaseline, 
						chAscetsResources, chArmCoordinates, chDnaGenomeFastaFiles,
						chDnaTargetGenes, chWhitelistGenes, chManeGtf, chManeGene, chManeExon,
						chCancerdrivers, chGenieCna, chCivicClinical,
						chCnaVcfHeader, chDnaSampleFile, chDnaPerBaseCovBed
					)
					chCnaGeneResults = chDnaSamples.map{ [it] }.combine(CNA.out.geneResults, by: 0)
					chCnaArmResults  = chDnaSamples.map{ [it] }.combine(CNA.out.armResults, by: 0)
					chCnaFiles       = chDnaSamples.map{ [it] }.combine(CNA.out.files, by: 0)
				}

				// MSI
				if (params.msi) {
					MSI(chDnaBamRealigned.join(chDnaBaiRealigned), chDnaTargetBed, chMsiBaseline)
					chMsiResults = chDnaSamples.map{ [it] }.combine(MSI.out.results, by: 0)
				}

				// MULTIQC of results
				if (params.multiqcRun) {
					MULTIQC_RUN_DNA(
						chDnaBclConvertReport.concat(chDnaFastqProcessingFiles, chDnaFastqcFiles, chDnaDedupFiles, chDnaBamQcRunFiles).collect(), 
						chMultiqcConfig, '09_MULTIQC_RUN', 'DNA'
					)
				}				

			} else {
				chDnaBclConvertReport = []
			}

			// RNA analysis
			if (params.RNA) {

				// Sample files
				PREPARE_SAMPLE_FILES_RNA(chSampleSheet, chClinicalInfoFile, chTumorNamesFile, "RNA")
				chRnaSampleFile = PREPARE_SAMPLE_FILES_RNA.out.txt
				chRnaSampleSheet = PREPARE_SAMPLE_FILES_RNA.out.csv

				// RNA samples 
				chRnaSamples             = chRnaSampleFile.splitCsv(sep: '\t', header: true).map(row -> "${row.Sample}")
				chRnaReads               = chRnaSamples.map{ [it] }.combine(chRnaFastq, by: 0)
				chRnaBamAligned          = chRnaSamples.map{ [it] }.combine(chRnaBamAligned, by: 0)
				chRnaBaiAligned          = chRnaSamples.map{ [it] }.combine(chRnaBaiAligned, by: 0)
				chRnaBamRealigned        = chRnaSamples.map{ [it] }.combine(chRnaBamRealigned, by: 0)
				chRnaBaiRealigned        = chRnaSamples.map{ [it] }.combine(chRnaBaiRealigned, by: 0)
				chRnaRealignmentSJ       = chRnaSamples.map{ [it] }.combine(chRnaRealignmentSJ, by: 0)
				chRnaRealignmentJunction = chRnaSamples.map{ [it] }.combine(chRnaRealignmentJunction, by: 0)
				chRnaBamQcRawResults     = chRnaSamples.map{ [it] }.combine(chRnaBamQcRawResults, by: 0)
				chRnaPerBaseCovBed       = chRnaSamples.map{ [it] }.combine(chRnaPerBaseCovBed, by: 0)
				chRnaBamQcFinalResults   = chRnaSamples.map{ [it] }.combine(chRnaBamQcFinalResults, by: 0)
				chRnaBamQcFiles          = chRnaSamples.map{ [it] }.combine(chRnaBamQcFiles, by: 0)
				chFusionResults          = chRnaSamples.map{ [it] }.combine(chFusionResults, by: 0)
				chFusionFiles            = chRnaSamples.map{ [it] }.combine(chFusionFiles, by: 0)
				chSplicingResults        = chRnaSamples.map{ [it] }.combine(chSplicingResults, by: 0)
				chSplicingFiles          = chRnaSamples.map{ [it] }.combine(chSplicingFiles, by: 0)

				// FASTQ generation
				if (params.fastqGeneration) {
					FASTQ_GENERATION_RNA(chStartingDataDir, chRnaSampleFile, chRnaSampleSheet, 'RNA')
					chRnaReads            = chRnaSamples.map{ [it] }.combine(FASTQ_GENERATION_RNA.out.fastq, by: 0)
					chRnaBclConvertReport = FASTQ_GENERATION_RNA.out.bclConvertReport
				}

				// FASTQ processing
				if (params.fastqProcessing) {
					FASTQ_PROCESSING_RNA(chRnaReads, chRnaSampleFile, chMouseAndHumanIndex, 'RNA')
					chRnaReads = chRnaSamples.map{ [it] }.combine(FASTQ_PROCESSING_RNA.out.fastq, by: 0)
					chRnaFastqProcessingFiles = FASTQ_PROCESSING_RNA.out.files
				}

				// FASTQ QC
				if (params.fastqQC) { 
					FASTQ_QC_RNA(chRnaReads, chMultiqcConfig, 'RNA')
					chRnaFastqcFiles = FASTQ_QC_RNA.out.files
				}

				// ALIGNMENT + DEDUPLICATION
				if (params.deduplicationRna) {

					if (params.alignment) {
						ALIGNMENT_RNA(chRnaReads, chRnaGenomeFastaFiles, chRnaGenomeIndexDir, chRnaTargetBed, chMultiqcConfig)
						chRnaBamAligned      = chRnaSamples.map{ [it] }.combine(ALIGNMENT_RNA.out.bam, by: 0)
						chRnaBaiAligned      = chRnaSamples.map{ [it] }.combine(ALIGNMENT_RNA.out.bai, by: 0)
						chRnaBamQcRawResults = chRnaSamples.map{ [it] }.combine(ALIGNMENT_RNA.out.results, by: 0)
					}

					DEDUPLICATION_RNA(chRnaBamAligned, chRnaBaiAligned, chRnaGenomeFastaFiles, chRnaTargetBed, chMultiqcConfig, 'RNA')
					chRnaReads      = chRnaSamples.map{ [it] }.combine(DEDUPLICATION_RNA.out.bam, by: 0)
					chRnaDedupFiles = DEDUPLICATION_RNA.out.files
				}

				// REALIGNMENT
				if (params.realignment) {
					chRnaRealignBam = params.deduplicationRna ? true: false
					REALIGNMENT_RNA(
						chRnaReads, chRnaRealignBam, chRnaGenomeFastaFiles, chRnaGenomeIndexDir, chRnaTargetBed, 
						chMultiqcConfig, chRnaSampleFile, chRnaTargetGenes, chWhitelistGenes, 
						chManeGtf, chManeGene, chManeExon, chManeCoding, chRnaBamQcRawResults
					)
					chRnaBamRealigned        = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.bam, by: 0)
					chRnaBaiRealigned        = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.bai, by: 0)
					chRnaRealignmentSJ       = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.SJ, by: 0)
					chRnaRealignmentJunction = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.junction, by: 0)
					chRnaPerBaseCovBed       = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.perBaseCovBed, by: 0)
					chRnaBamQcFinalResults   = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.results, by: 0)
					chRnaBamQcFiles          = chRnaSamples.map{ [it] }.combine(REALIGNMENT_RNA.out.files, by: 0)
					chRnaBamQcRunFiles       = REALIGNMENT_RNA.out.runFiles
				}

				// FUSION
				if (params.fusion) {
					FUSION(
						chRnaBamRealigned.join(chRnaBaiRealigned), chRnaRealignmentJunction, chCtatLibDir, 
						chRnaTargetBed, chRnaTargetGenes, chWhitelistGenes, chFusionWhitelist, 
						chManeGtf, chManeGene, chManeExon, chManeIntron, chCancerdrivers, 
						chMitelmanFusion, chGenieFusion, chCivicClinical, chHg38ToHg19ChainFile, 
						chRnaGenomeFastaFiles, chFusionVcfHeader, chRnaSampleFile, chRnaPerBaseCovBed
					)
					chFusionResults = chRnaSamples.map{ [it] }.combine(FUSION.out.results, by: 0)
					chFusionFiles   = chRnaSamples.map{ [it] }.combine(FUSION.out.files, by: 0)
				}

				// SPLICING
				if (params.splicing) {
					SPLICING(
						chRnaRealignmentSJ.join(chRnaRealignmentJunction), chCtatLibDir, 
						chRnaTargetBed, chRnaTargetGenes, chWhitelistGenes, chSplicingWhitelist, 
						chManeGtf, chManeGene, chManeExon, chManeIntron, 
						chCancerdrivers, chCivicClinical, chHg38ToHg19ChainFile, 
						chRnaGenomeFastaFiles, chSplicingVcfHeader, chRnaSampleFile, chRnaPerBaseCovBed, 
						chRnaSamples.map{ [it] }.combine(chVariantAnnotResults, by: 0)
					)
					chSplicingResults = chRnaSamples.map{ [it] }.combine(SPLICING.out.results, by: 0)
					chSplicingFiles   = chRnaSamples.map{ [it] }.combine(SPLICING.out.files, by: 0)
				}

				// MULTIQC of results
				if (params.multiqcRun) {
					MULTIQC_RUN_RNA(
						chRnaBclConvertReport.concat(chRnaFastqProcessingFiles, chRnaFastqcFiles, chRnaDedupFiles, chRnaBamQcRunFiles).collect(),
						chMultiqcConfig, '09_MULTIQC_RUN', 'RNA'
					)
				}				
				
			} else {
				chRnaBclConvertReport = []
			}

			// MULTIQC of BCL convert
			if (params.startingDataType.toUpperCase() == 'BCL' && params.fastqGeneration && (params.DNA || params.RNA)) {
				MULTIQC_BCL_CONVERT(chDnaBclConvertReport, chRnaBclConvertReport, chMultiqcConfig, '00_FASTQ_GENERATION', '02_MULTIQC_BCL_CONVERT')
			}

			// RUN LIBRARY
			if (params.runLibrary) {
				chResults = Channel.empty()
					.concat(
						chDnaBamQcFinalResults.transpose().map{ [it[1]] }, 
						chVariantConsensusResults.transpose().map{ [it[1]] },
						chVariantAnnotResults.transpose().map{ [it[1]] }, 
						chTmbResults.transpose().map{ [it[1]] },
						chCnaGeneResults.transpose().map{ [it[1]] }, 
						chCnaArmResults.transpose().map{ [it[1]] }, 
						chMsiResults.transpose().map{ [it[1]] }, 
						chRnaBamQcFinalResults.transpose().map{ [it[1]] }, 
						chFusionResults.transpose().map{ [it[1]] }, 
						chSplicingResults.transpose().map{ [it[1]] })
					.collect()
				
				RUN_LIBRARY(chResults, chDnaSampleFile, chRnaSampleFile, '09_VARIANT_LIBRARY', '01_RUN_LIBRARY')
				chRunLibrary = RUN_LIBRARY.out.db
			}

			// UNIQUE LIBRARY
			if (params.uniqueLibrary) {
				UNIQUE_LIBRARY(chUniqueLibrary, chRunLibrary, Channel.fromPath(params.projectDir), '09_VARIANT_LIBRARY', '02_UNIQUE_LIBRARY')
				chUniqueLibrary = UNIQUE_LIBRARY.out.db.collect()
			}

			// SAMPLE REPORT
			if (params.sampleReport) {
				chResults = Channel.empty()
					.concat(
						chDnaBamQcFiles.transpose(), chVariantAnnotFiles.transpose(), 
						chCnaFiles.transpose(), chTmbResults.transpose(), chTmbTrace.transpose(), chMsiResults.transpose(), 
						chRnaBamQcFiles.transpose(), chFusionFiles.transpose(), chSplicingFiles.transpose())
					.groupTuple()
				chReportRmdFiles = Channel.fromPath("${params.codeDir}/*.Rmd").collect()
				SAMPLE_REPORT(
					chResults, chDnaSampleFile, chRnaSampleFile, chDnaPanelHotspotsBed, 
					chCivicClinical, chUniqueLibrary, chReportRmdFiles, '10_SAMPLE_REPORT')
			}
		}
	}
}
