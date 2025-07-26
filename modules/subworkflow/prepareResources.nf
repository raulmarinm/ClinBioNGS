/*
 * Prepare resources workflow
 */

include { 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_MANE_GTF;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_MANE_SUMMARY; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CYTOBAND_BED;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CADD_SNV;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CADD_SNV_TBI;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CADD_INDEL;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CADD_INDEL_TBI;
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_REVEL_FILE; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_ALPHAMISSENSE_FILE; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CLINVAR_VCF; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CLINVAR_VCF_TBI; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CIVIC_VCF; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_HOTSPOTS_BED; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CANCERHOTSPOTS; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CGI_MUTATIONS; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_UCSC_UNUSUAL; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_ENCODE_BLACKLIST; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GRC_EXCLUSIONS; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GIAB_REPEATS; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GIAB_SEGDUPS; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GIAB_LOWMAP; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GIAB_MHC; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_GIAB_VDJ; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CTR_BED; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CIVIC_VARIANT; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CIVIC_MOLECULAR; 
	DOWNLOAD_ANNOT_FILE as DOWNLOAD_CIVIC_EVIDENCE;     } from '../process/downloadFiles'

include { PREPARE_MANE_FILES                            } from '../process/prepareManeFiles'
include { 
	DOWNLOAD_GENOME as DOWNLOAD_GENOME_DNA; 
	INDEX_GENOME_BWAMEM2; INDEX_GENOME_TMAP;
	DOWNLOAD_GENOME as DOWNLOAD_GENOME_MOUSE;           } from '../process/prepareGenomeFiles'
include { 
	DOWNLOAD_CHAIN_FILE as DOWNLOAD_CHAIN_HG38TOHG19;
	DOWNLOAD_CHAIN_FILE as DOWNLOAD_CHAIN_HG19TOHG38;   } from '../process/downloadChainFile'
include { 
	MANIFEST_TO_BED as MANIFEST_TO_BED_DNA; 
	MANIFEST_TO_BED as MANIFEST_TO_BED_RNA;             } from '../process/manifestToBed'
include { 
	PREPARE_MANIFEST_IONTORRENT_DNA; 
	PREPARE_HOTSPOTS_IONTORRENT; 
	PREPARE_MANIFEST_IONTORRENT_RNA;                    } from '../process/prepareIonTorrentFiles.nf'
include { 
	BED_LIFTOVER_MANIFEST as BED_LIFTOVER_MANIFEST_DNA; 
	BED_LIFTOVER_MANIFEST as BED_LIFTOVER_MANIFEST_RNA;
	BED_LIFTOVER_VARIANT as BED_LIFTOVER_HOTSPOT_DNA; 
	BED_LIFTOVER_VARIANT as BED_LIFTOVER_BLACKLIST_DNA; } from '../process/bedLiftover'
include { 
	BED_ANNOTATE_GENES as BED_ANNOTATE_GENES_DNA; 
	BED_ANNOTATE_GENES as BED_ANNOTATE_GENES_RNA;       } from '../process/bedAnnotateGenes'
include { 
	PREPARE_TARGET_BED as PREPARE_TARGET_BED_DNA; 
	PREPARE_TARGET_BED as PREPARE_TARGET_BED_RNA; 
	PREPARE_TARGET_GENES as PREPARE_TARGET_GENES_DNA; 
	PREPARE_TARGET_GENES as PREPARE_TARGET_GENES_RNA;
	PREPARE_TARGET_CHROM;                               } from '../process/prepareTargetFiles'
include { BED_ADD_PADDING                               } from '../process/bedAddPadding'
include { 
	BED_TO_INTERVAL_LIST as BED_TO_INTERVAL_LIST;
	BED_TO_INTERVAL_LIST as BED_TO_INTERVAL_LIST_EXT    } from '../process/bedToIntervalList'
include { BED_COMPLEMENT                                } from '../process/bedComplement'
include { DOWNLOAD_VEP_CACHE                            } from '../process/vep'
include { PREPARE_REVEL_FILE; INDEX_ALPHAMISSENSE_FILE  } from '../process/prepareVepFiles'
include { PREPARE_HOTSPOTS_BED; PREPARE_CANCERHOTSPOTS  } from '../process/prepareCancerhotspots'
include { 
	BIGBED_TO_BED as BIGBED_TO_BED_UCSC_UNUSUAL;
	BIGBED_TO_BED as BIGBED_TO_BED_ENCODE_BLACKLIST;
	BIGBED_TO_BED as BIGBED_TO_BED_GRC_EXCLUSIONS;      } from '../process/bigBedToBed'
include { BED_MERGE as BED_MERGE_PROBLEMATIC_REGIONS    } from '../process/bedMerge'
include { PREPARE_ARM_COORDINATES                       } from '../process/prepareArmCoordinates'
include { DOWNLOAD_ASCETS                               } from '../process/downloadAscets'
include { PREPARE_CTAT_LIB                              } from '../process/prepareCtatResource'
include { PREPARE_MITELMAN_FUSION                       } from '../process/prepareMitelmanFile'
include { PREPARE_CANCERDRIVERS                         } from '../process/prepareCancerdrivers'
include { PREPARE_CGI_MUTATIONS                         } from '../process/prepareCgiMutations'
include { PREPARE_CIVIC_VCF; PREPARE_CIVIC_EVIDENCE     } from '../process/prepareCivicFiles'
include { XENGSORT_INDEX                                } from '../process/xengsort'
include { PREPARE_SAMPLESHEET                           } from '../process/prepareSampleSheet'

// Prepare resources files if required

workflow RESOURCES {

	take:
		startingDir
		tumorNames

	main:

		// MANE files
		def maneGtfName     = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.ensembl_genomic.gtf.gz"
		def maneSummaryName = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.summary.txt.gz"
		def maneGeneName    = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.genes.txt"
		def maneExonName    = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.exons.txt"
		def maneIntronName  = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.introns.txt"
		def maneCodingName  = "MANE.${params.genomeVersionGrc}.v${params.maneVersion}.coding.txt"
		chManeGtf = file("${params.annotationDir}/${maneGtfName}").exists() ? 
			Channel.fromPath("${params.annotationDir}/${maneGtfName}").collect() : 
			DOWNLOAD_MANE_GTF(params.maneGtfLink, maneGtfName, 'DOWNLOAD_MANE_GTF').file
		chManeSummary = file("${params.annotationDir}/${maneSummaryName}").exists() ? 
			Channel.fromPath("${params.annotationDir}/${maneSummaryName}").collect() : 
			DOWNLOAD_MANE_SUMMARY(params.maneSummaryLink, maneSummaryName, 'DOWNLOAD_MANE_SUMMARY').file
		if ([
			"${params.annotationDir}/${maneGeneName}.gz", "${params.annotationDir}/${maneExonName}.gz", 
			"${params.annotationDir}/${maneIntronName}.gz", "${params.annotationDir}/${maneCodingName}.gz"].any{!(file(it).exists())}) {
				PREPARE_MANE_FILES(chManeGtf, chManeSummary, maneGeneName, maneExonName, maneIntronName, maneCodingName)
				chManeGene   = PREPARE_MANE_FILES.out.genes
				chManeExon   = PREPARE_MANE_FILES.out.exons
				chManeIntron = PREPARE_MANE_FILES.out.introns
				chManeCoding = PREPARE_MANE_FILES.out.coding
		} else {
			chManeGene   = Channel.fromPath("${params.annotationDir}/${maneGeneName}.gz").collect()
			chManeExon   = Channel.fromPath("${params.annotationDir}/${maneExonName}.gz").collect()
			chManeIntron = Channel.fromPath("${params.annotationDir}/${maneIntronName}.gz").collect()
			chManeCoding = Channel.fromPath("${params.annotationDir}/${maneCodingName}.gz").collect()
		}

		// UCSC cytoband file
		def cytobandName = "UCSC_${params.genomeVersionHg}_cytoBand.bed.gz"
		chCytobandBed = file("${params.annotationDir}/${cytobandName}").exists() ? 
			Channel.fromPath("${params.annotationDir}/${cytobandName}").collect() : 
			DOWNLOAD_CYTOBAND_BED(params.cytobandLink, cytobandName, 'DOWNLOAD_CYTOBAND_BED').file

		// hg38ToHg19 and hg19ToHg38 chain files
		def liftOverChainDir = "${params.genomeDir}/Homo_sapiens/UCSC/LiftOver_Chain_Files"
		chHg38ToHg19ChainFile = file("${liftOverChainDir}/hg38ToHg19.over.chain").exists() ? 
			Channel.fromPath("${liftOverChainDir}/hg38ToHg19.over.chain").collect() : 
			DOWNLOAD_CHAIN_HG38TOHG19('hg38', 'hg38ToHg19', liftOverChainDir).chain
		chHg19ToHg38ChainFile = file("${liftOverChainDir}/hg19ToHg38.over.chain").exists() ? 
			Channel.fromPath("${liftOverChainDir}/hg19ToHg38.over.chain").collect() : 
			DOWNLOAD_CHAIN_HG19TOHG38('hg19', 'hg19ToHg38', liftOverChainDir).chain

		// DNA resources
		chDnaTargetBed             = Channel.empty()
		chDnaTargetIntervalList    = Channel.empty()
		chDnaTargetBedExt          = Channel.empty()
		chDnaTargetIntervalListExt = Channel.empty()
		chDnaTargetGenes           = Channel.empty()
		chDnaTargetChr             = Channel.empty()
		chDnaGenomeFastaDir        = Channel.empty()
		chDnaGenomeFastaFiles      = Channel.empty()
		chDnaGenomeIndexFiles      = Channel.empty()
		chDnaPanelHotspotsBed      = []
		chDnaPanelBlacklistBed     = []
		chDnaOffTargetBed          = Channel.empty()
		chVepCacheDir              = Channel.empty()
		chVepCaddSnv               = Channel.empty()
		chVepCaddIndel             = Channel.empty()
		chVepRevel                 = Channel.empty()
		chVepAlphaMissense         = Channel.empty()
		chClinvarVcf               = Channel.empty()
		chCivicVcf                 = Channel.empty()
		chHotspotsWhitelistBed     = Channel.empty()
		chCancerhotspotsResults    = Channel.empty()
		chCgiOncogenic             = Channel.empty()
		chProblematicRegionsBed    = Channel.empty()
		chCtrRegionsBed            = Channel.empty()
		chArmCoordinates           = Channel.empty()
		chAscetsResources          = Channel.empty()
		if (params.DNA) {

			// Genome FASTA files
			def dnaGenomeFastaDir = "Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/WholeGenomeFasta"
			def dnaGenomeStepsList = [
				params.alignment, params.deduplicationDna, params.realignment, 
				params.smallVariantCalling, params.smallVariantAnnotation, params.cna, params.pdx
			]
			if (dnaGenomeStepsList.any{it}) {
				if ([
					"${params.genomeDir}/${dnaGenomeFastaDir}/genome.fa", 
					"${params.genomeDir}/${dnaGenomeFastaDir}/genome.fa.fai", 
					"${params.genomeDir}/${dnaGenomeFastaDir}/genome.dict"].any{!(file(it).exists())}) {
					DOWNLOAD_GENOME_DNA(params.dnaGenomeLink, dnaGenomeFastaDir, "DOWNLOAD_GENOME_DNA")
					chDnaGenomeFastaDir   = DOWNLOAD_GENOME_DNA.out.fastaDir
					chDnaGenomeFastaFiles = DOWNLOAD_GENOME_DNA.out.fastaFiles
				} else {
					chDnaGenomeFastaDir   = Channel.fromPath("${params.genomeDir}/${dnaGenomeFastaDir}").collect()
					chDnaGenomeFastaFiles = Channel.fromPath("${params.genomeDir}/${dnaGenomeFastaDir}/*.{dict,fa,fai,xml}").collect()
				}
			}

			// Genome index files
			if (params.seqPlatform.toLowerCase() == 'iontorrent'){
				def dnaGenomeIndexDir = "${params.genomeDir}/Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/TMAPIndex/v${params.tmapTvcVersion}"
				if ([params.alignment, params.realignment].any{it}) {
					chDnaGenomeIndexFiles = file("${dnaGenomeIndexDir}/genome.fa").exists() ? 
						Channel.fromPath("${dnaGenomeIndexDir}/genome.fa*").collect() : 
						INDEX_GENOME_TMAP(chDnaGenomeFastaFiles, dnaGenomeIndexDir).indexFiles
				}
			} else {
				def dnaGenomeIndexDir = "${params.genomeDir}/Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/BwaMem2Index/v${params.bwamem2Version}"
				if ([params.alignment, params.realignment].any{it}) {
					chDnaGenomeIndexFiles = file("${dnaGenomeIndexDir}/genome.fa").exists() ? 
						Channel.fromPath("${dnaGenomeIndexDir}/genome.fa*").collect() : 
						INDEX_GENOME_BWAMEM2(chDnaGenomeFastaFiles, dnaGenomeIndexDir).indexFiles
				}
			}

			// Target files from panel manifests
			if (params.seqPanel) {
				if (file(params.dnaTargetBed).exists()) {
					chDnaTargetBed = Channel.fromPath(params.dnaTargetBed).collect()
				} else {
					
					// BED file from manifest
					if (params.manifestToBed) {
						MANIFEST_TO_BED_DNA(file(params.dnaManifest, checkIfExists: true), 'MANIFEST_TO_BED_DNA')
						chDnaTargetBed = MANIFEST_TO_BED_DNA.out.bed
					} else if (params.prepareIontorrentManifest) {
						PREPARE_MANIFEST_IONTORRENT_DNA(startingDir)
						chDnaTargetBed = PREPARE_MANIFEST_IONTORRENT_DNA.out.bed
					} else {
						chDnaTargetBed = Channel.fromPath(params.dnaManifest, checkIfExists: true).collect()
					}

					// Liftover BED file
					if (params.liftoverManifest) {
						BED_LIFTOVER_MANIFEST_DNA(chDnaTargetBed, chHg19ToHg38ChainFile, 'BED_LIFTOVER_MANIFEST_DNA')
						chDnaTargetBed = BED_LIFTOVER_MANIFEST_DNA.out.bed
					}

					// Annotate genes to the BED file
					if (params.annotateGenesBed) {
						BED_ANNOTATE_GENES_DNA(chDnaTargetBed, chManeGene, 'BED_ANNOTATE_GENES_DNA')
						chDnaTargetBed = BED_ANNOTATE_GENES_DNA.out.bed
					}

					// Prepare BED file
					PREPARE_TARGET_BED_DNA(chDnaTargetBed, 'DNA', 'PREPARE_TARGET_BED_DNA')
					chDnaTargetBed = PREPARE_TARGET_BED_DNA.out.bed

				}			
				
				chDnaTargetIntervalList = file(params.dnaTargetIntervalList).exists() ? Channel.fromPath(params.dnaTargetIntervalList).collect() : 
					BED_TO_INTERVAL_LIST(chDnaTargetBed, chDnaGenomeFastaFiles, 'BED_TO_INTERVAL_LIST_DNA').interval_list
				
				if (params.variantIntervalPadding && params.variantIntervalPadding != '0') {
					chDnaTargetBedExt = file(params.dnaTargetBedExt).exists() ? Channel.fromPath(params.dnaTargetBedExt).collect() : 
						BED_ADD_PADDING(chDnaTargetBed, params.variantIntervalPadding, 'BED_ADD_PADDING_DNA').bed
					chDnaTargetIntervalListExt = file(params.dnaTargetIntervalListExt).exists() ? Channel.fromPath(params.dnaTargetIntervalListExt).collect() : 
						BED_TO_INTERVAL_LIST_EXT(chDnaTargetBedExt, chDnaGenomeFastaFiles, 'BED_TO_INTERVAL_LIST_EXT_DNA').interval_list
				} else {
					chDnaTargetBedExt = chDnaTargetBed
					chDnaTargetIntervalListExt = chDnaTargetIntervalList
				}			
				
				// Target genes from panel manifests
				chDnaTargetGenes = file(params.dnaTargetGenes).exists() ? Channel.fromPath(params.dnaTargetGenes).collect() : 
					PREPARE_TARGET_GENES_DNA(chManeGene, chCytobandBed, chDnaTargetBed, 'DNA', 'PREPARE_TARGET_GENES_DNA').genes
				
				// Target chromosomes from panel manifests
				chDnaTargetChr = file(params.dnaTargetChr).exists() ? Channel.fromPath(params.dnaTargetChr).collect() : 
					PREPARE_TARGET_CHROM(chDnaTargetBed, 'DNA', 'PREPARE_TARGET_CHROM_DNA').chrom
			}

			// Variant hotspots
			if (params.variantHotspotBed && file(params.variantHotspotBed).exists()) {
				chDnaPanelHotspotsBed = Channel.fromPath(params.variantHotspotBed).collect()
				if (params.liftoverVariantHotspotBed) {
					BED_LIFTOVER_HOTSPOT_DNA(chDnaPanelHotspotsBed, chHg19ToHg38ChainFile, 'BED_LIFTOVER_HOTSPOT_DNA')
					chDnaPanelHotspotsBed = BED_LIFTOVER_HOTSPOT_DNA.out.bed
				}
			} else if (params.prepareIontorrentVariantHotspots) {
				PREPARE_HOTSPOTS_IONTORRENT(startingDir, chHg19ToHg38ChainFile)
				chDnaPanelHotspotsBed = PREPARE_HOTSPOTS_IONTORRENT.out.bed	
			}

			// Variant blacklist
			if (params.variantBlacklistBed && file(params.variantBlacklistBed).exists()) {
				chDnaPanelBlacklistBed = Channel.fromPath(params.variantBlacklistBed).collect()
				if (params.liftoverVariantBlacklistBed) {
					BED_LIFTOVER_BLACKLIST_DNA(chDnaPanelBlacklistBed, chHg19ToHg38ChainFile, 'BED_LIFTOVER_BLACKLIST_DNA')
					chDnaPanelBlacklistBed = BED_LIFTOVER_BLACKLIST_DNA.out.bed
				}
			}

			// BAM clipping: off-target BED
			if (params.smallVariantCalling && (params.bamclipping)) {

				// Complementary target regions (not overlapping)
				chDnaOffTargetBed = file(params.dnaOffTargetBed).exists() ? Channel.fromPath(params.dnaOffTargetBed).collect() : 
					BED_COMPLEMENT(chDnaTargetBed, chDnaGenomeFastaFiles, 'BED_COMPLEMENT_DNA').bed
			}

			// Small variant annotation files
			if (params.smallVariantAnnotation) {

				// VEP cache
				if (params.downloadVepCache | !file("${params.annotationDir}/${params.vepSpecie}/${params.vepVersion}_${params.genomeVersionGrc}").exists()) {
					chVepCacheDir = DOWNLOAD_VEP_CACHE(params.vepSpecie, params.genomeVersionGrc, params.vepVersion).dir
				} else { chVepCacheDir = Channel.fromPath("${params.annotationDir}/${params.vepSpecie}").collect() }

				// CADD
				def caddSnvFile   = "CADD_v${params.caddVersion}_${params.genomeVersionGrc}_SNVs.tsv.gz"
				def caddIndelFile = "CADD_v${params.caddVersion}_${params.genomeVersionGrc}_InDels.tsv.gz"
				chVepCaddSnv = file("${params.annotationDir}/${caddSnvFile}").exists() ?
					Channel.fromPath(["${params.annotationDir}/${caddSnvFile}"]).collect() :
					DOWNLOAD_CADD_SNV(params.caddSnvLink, caddSnvFile, 'DOWNLOAD_CADD_SNV').file
				chVepCaddSnvTbi = file("${params.annotationDir}/${caddSnvFile}.tbi").exists() ?
					Channel.fromPath(["${params.annotationDir}/${caddSnvFile}.tbi"]).collect() :
					DOWNLOAD_CADD_SNV_TBI("${params.caddSnvLink}.tbi", "${caddSnvFile}.tbi", 'DOWNLOAD_CADD_SNV_TBI').file
				chVepCaddIndel = file("${params.annotationDir}/${caddIndelFile}").exists() ?
					Channel.fromPath(["${params.annotationDir}/${caddIndelFile}"]).collect() :
					DOWNLOAD_CADD_INDEL(params.caddIndelLink, caddIndelFile, 'DOWNLOAD_CADD_INDEL').file
				chVepCaddIndelTbi = file("${params.annotationDir}/${caddIndelFile}.tbi").exists() ?
					Channel.fromPath(["${params.annotationDir}/${caddIndelFile}.tbi"]).collect() :
					DOWNLOAD_CADD_INDEL_TBI("${params.caddIndelLink}.tbi", "${caddIndelFile}.tbi", 'DOWNLOAD_CADD_INDEL_TBI').file
				chVepCaddSnv   = chVepCaddSnv.concat(chVepCaddSnvTbi).collect()
				chVepCaddIndel = chVepCaddIndel.concat(chVepCaddIndelTbi).collect()

				// REVEL
				def revelFile = "revel-v${params.revelVersion}_${params.genomeVersionGrc}.tsv.gz"
				if (["${params.annotationDir}/${revelFile}", "${params.annotationDir}/${revelFile}.tbi"].any{!(file(it).exists())}) {
					def revelZip  = "revel-v${params.revelVersion}_all_chromosomes.zip"
					chRevelZip = file("${params.annotationDir}/${revelZip}").exists() ? 
						Channel.fromPath("${params.annotationDir}/${revelZip}").collect() : 
						DOWNLOAD_REVEL_FILE(params.revelLink, revelZip, 'DOWNLOAD_REVEL_FILE').file
					chVepRevel = PREPARE_REVEL_FILE(chRevelZip, revelFile).files
				} else {
					chVepRevel = Channel.fromPath(["${params.annotationDir}/${revelFile}", "${params.annotationDir}/${revelFile}.tbi"]).collect()
				}
				
				// AlphaMissense
				def alphaMissenseFile = "AlphaMissense_${params.genomeVersionHg}.tsv.gz"
				chVepAlphaMissense = file("${params.annotationDir}/${alphaMissenseFile}").exists() ? 
					Channel.fromPath("${params.annotationDir}/${alphaMissenseFile}").collect() : 
					DOWNLOAD_ALPHAMISSENSE_FILE(params.alphaMissenseLink, alphaMissenseFile, 'DOWNLOAD_ALPHAMISSENSE_FILE').file
				chVepAlphaMissenseTbi = file("${params.annotationDir}/${alphaMissenseFile}.tbi").exists() ?
					Channel.fromPath("${params.annotationDir}/${alphaMissenseFile}.tbi").collect() :
					INDEX_ALPHAMISSENSE_FILE(chVepAlphaMissense).file
				chVepAlphaMissense = chVepAlphaMissense.concat(chVepAlphaMissenseTbi).collect()
				
				// ClinVar
				def clinvarVcfFile = "clinvar_${params.clinvarDate}_${params.genomeVersionGrc}.vcf.gz"
				chClinvarVcf = file("${params.annotationDir}/${clinvarVcfFile}").exists() ?
					Channel.fromPath("${params.annotationDir}/${clinvarVcfFile}").collect() : 
					DOWNLOAD_CLINVAR_VCF(params.clinvarLink, clinvarVcfFile, 'DOWNLOAD_CLINVAR_VCF').file
				chClinvarVcfTbi = file("${params.annotationDir}/${clinvarVcfFile}.tbi").exists() ?
					Channel.fromPath("${params.annotationDir}/${clinvarVcfFile}.tbi").collect() : 
					DOWNLOAD_CLINVAR_VCF_TBI("${params.clinvarLink}.tbi", "${clinvarVcfFile}.tbi", 'DOWNLOAD_CLINVAR_VCF_TBI').file
				chClinvarVcf = chClinvarVcf.concat(chClinvarVcfTbi).collect()
				
				// CIViC
				def civicVcfFile = "CIViC_${params.civicDate}_${params.genomeVersionGrc}_accepted.vcf.gz"
				chCivicVcf = ["${params.annotationDir}/${civicVcfFile}", "${params.annotationDir}/${civicVcfFile}.tbi"].any{!(file(it).exists())} ?
					PREPARE_CIVIC_VCF(
						DOWNLOAD_CIVIC_VCF(params.civicVcfLink, "CIViC_${params.civicDate}_accepted.vcf", 'DOWNLOAD_CIVIC_VCF').file, 
						chHg19ToHg38ChainFile, chDnaGenomeFastaFiles, civicVcfFile).vcf :
					Channel.fromPath(["${params.annotationDir}/${civicVcfFile}", "${params.annotationDir}/${civicVcfFile}.tbi"]).collect()

				// Cancer hotspot regions (BED)
				def hotspotsFile = "known_somatic_sites_${params.genomeVersionHg}.bed"
				chHotspotsWhitelistBed = file("${params.annotationDir}/${hotspotsFile}").exists() ?
					Channel.fromPath("${params.annotationDir}/${hotspotsFile}").collect() :
					PREPARE_HOTSPOTS_BED(
						file("${params.annotationDir}/known_somatic_sites_GRCh37.bed").exists() ? 
							Channel.fromPath("${params.annotationDir}/known_somatic_sites_GRCh37.bed").collect() :
							DOWNLOAD_HOTSPOTS_BED(params.hotspotsWhitelistBedLink, "known_somatic_sites_GRCh37.bed", 'DOWNLOAD_HOTSPOTS_BED').file,
						chHg19ToHg38ChainFile, hotspotsFile).bed
				
				// Cancerhotspots resource (excel results)
				def cancerhotspotsFile  = "cancerhotspots_v${params.cancerhotspotsVersion}"
				chCancerhotspotsResults = file("${params.annotationDir}/${cancerhotspotsFile}.processed.txt").exists() ?
					Channel.fromPath("${params.annotationDir}/${cancerhotspotsFile}.processed.txt").collect() :
					PREPARE_CANCERHOTSPOTS(
						file("${params.annotationDir}/${cancerhotspotsFile}.xls").exists() ? 
							Channel.fromPath("${params.annotationDir}/${cancerhotspotsFile}.xls").collect() :
							DOWNLOAD_CANCERHOTSPOTS(params.cancerhotspotsLink, "${cancerhotspotsFile}.xls", 'DOWNLOAD_CANCERHOTSPOTS').file, 
						"${cancerhotspotsFile}.processed.txt").file
				
				// CGI oncogenic mutations
				def cgiMutationsFile = "CGI_mutations_${params.cgiMutationsDate}"
				chCgiOncogenic = file("${params.annotationDir}/${cgiMutationsFile}.processed.txt").exists() ? 
					Channel.fromPath("${params.annotationDir}/${cgiMutationsFile}.processed.txt").collect() : 
					PREPARE_CGI_MUTATIONS(
						chHg19ToHg38ChainFile,
						DOWNLOAD_CGI_MUTATIONS(params.cgiMutationsLink, "${cgiMutationsFile}.tsv", 'DOWNLOAD_CGI_MUTATIONS').file, 
						"${cgiMutationsFile}.processed.txt").file


				// UCSC Problematic regions
				def ucscUnusualRegions  = "UCSC_${params.genomeVersionHg}_unusualRegions_${params.ucscProblematicAccessDate}"
				def ucscEncodeBlacklist = "UCSC_${params.genomeVersionHg}_encodeBlacklist_${params.ucscProblematicAccessDate}"
				def ucscGrcExclusions   = "UCSC_${params.genomeVersionHg}_grcExclusions_${params.ucscProblematicAccessDate}"		
				chUcscUnusualRegions = file("${params.annotationDir}/${ucscUnusualRegions}.bed.gz").exists() ?
					Channel.fromPath("${params.annotationDir}/${ucscUnusualRegions}.bed.gz").collect() :
					BIGBED_TO_BED_UCSC_UNUSUAL(
						file("${params.annotationDir}/${ucscUnusualRegions}.bb").exists() ? 
							Channel.fromPath("${params.annotationDir}/${ucscUnusualRegions}.bb").collect() :
							DOWNLOAD_UCSC_UNUSUAL(params.ucscUnusualRegionsLink, "${ucscUnusualRegions}.bb", 'DOWNLOAD_UCSC_UNUSUAL').file, 
						"${ucscUnusualRegions}.bed", 'BIGBED_TO_BED_UCSC_UNUSUAL').bed
				chUcscEncodeBlacklist = file("${params.annotationDir}/${ucscEncodeBlacklist}.bed.gz").exists() ?
					Channel.fromPath("${params.annotationDir}/${ucscEncodeBlacklist}.bed.gz").collect() :
					BIGBED_TO_BED_ENCODE_BLACKLIST(
						file("${params.annotationDir}/${ucscEncodeBlacklist}.bb").exists() ? 
							Channel.fromPath("${params.annotationDir}/${ucscEncodeBlacklist}.bb").collect() :
							DOWNLOAD_ENCODE_BLACKLIST(params.ucscEncodeBlacklistLink, "${ucscEncodeBlacklist}.bb", 'DOWNLOAD_ENCODE_BLACKLIST').file, 
						"${ucscEncodeBlacklist}.bed", 'BIGBED_TO_BED_ENCODE_BLACKLIST').bed				
				chUcscGrcExclusions = file("${params.annotationDir}/${ucscGrcExclusions}.bed.gz").exists() ?
					Channel.fromPath("${params.annotationDir}/${ucscGrcExclusions}.bed.gz").collect() :
					BIGBED_TO_BED_GRC_EXCLUSIONS(
						file("${params.annotationDir}/${ucscGrcExclusions}.bb").exists() ? 
							Channel.fromPath("${params.annotationDir}/${ucscGrcExclusions}.bb").collect() :
							DOWNLOAD_GRC_EXCLUSIONS(params.ucscGrcExclusionsLink, "${ucscGrcExclusions}.bb", 'DOWNLOAD_GRC_EXCLUSIONS').file, 
						"${ucscGrcExclusions}.bed", 'BIGBED_TO_BED_GRC_EXCLUSIONS').bed
				chUcscProblematicRegions = Channel.empty().concat(chUcscUnusualRegions, chUcscEncodeBlacklist, chUcscGrcExclusions)
				
				// GIAB stratifications
				def giabRepeats = "GIAB_v${params.giabVersion}_${params.genomeVersionGrc}_AllTandemRepeatsandHomopolymers_slop5.bed.gz"
				def giabSegDup  = "GIAB_v${params.giabVersion}_${params.genomeVersionGrc}_segdups.bed.gz"
				def giabLowMap  = "GIAB_v${params.giabVersion}_${params.genomeVersionGrc}_lowmappabilityall.bed.gz"
				def giabMhc     = "GIAB_v${params.giabVersion}_${params.genomeVersionGrc}_MHC.bed.gz"
				def giabVdj     = "GIAB_v${params.giabVersion}_${params.genomeVersionGrc}_VDJ.bed.gz"
				chGiabRepeats = file("${params.annotationDir}/${giabRepeats}").exists() ?
					Channel.fromPath("${params.annotationDir}/${giabRepeats}").collect() : 
					DOWNLOAD_GIAB_REPEATS(params.giabRepeatsLink, giabRepeats, 'DOWNLOAD_GIAB_REPEATS').file
				chGiabSegDup = file("${params.annotationDir}/${giabSegDup}").exists() ?
					Channel.fromPath("${params.annotationDir}/${giabSegDup}").collect() : 
					DOWNLOAD_GIAB_SEGDUPS(params.giabSegDupLink, giabSegDup, 'DOWNLOAD_GIAB_SEGDUPS').file
				chGiabLowMap = file("${params.annotationDir}/${giabLowMap}").exists() ?
					Channel.fromPath("${params.annotationDir}/${giabLowMap}").collect() : 
					DOWNLOAD_GIAB_LOWMAP(params.giabLowMapLink, giabLowMap, 'DOWNLOAD_GIAB_LOWMAP').file
				chGiabMhc = file("${params.annotationDir}/${giabMhc}").exists() ?
					Channel.fromPath("${params.annotationDir}/${giabMhc}").collect() : 
					DOWNLOAD_GIAB_MHC(params.giabMhcLink, giabMhc, 'DOWNLOAD_GIAB_MHC').file
				chGiabVdj = file("${params.annotationDir}/${giabVdj}").exists() ?
					Channel.fromPath("${params.annotationDir}/${giabVdj}").collect() : 
					DOWNLOAD_GIAB_VDJ(params.giabVdjLink, giabVdj, 'DOWNLOAD_GIAB_VDJ').file
				chGiabStratifications = Channel.empty().concat(chGiabRepeats, chGiabSegDup, chGiabLowMap, chGiabMhc, chGiabVdj)
				
				// Merge BED files of problematic regions
				def problematicRegions = "problematicRegions_${params.genomeVersionGrc}_GIAB-v${params.giabVersion}_UCSC-${params.ucscProblematicAccessDate}.bed"
				chProblematicRegionsBed = file("${params.annotationDir}/${problematicRegions}.gz").exists() ?
					Channel.fromPath("${params.annotationDir}/${problematicRegions}.gz").collect() :
					BED_MERGE_PROBLEMATIC_REGIONS(
						chUcscProblematicRegions.concat(chGiabStratifications).collect(), 
						problematicRegions, 'BED_MERGE_PROBLEMATIC_REGIONS').bed
				
				// Consensus Targeted Regions (CTR) BED file
				def ctrRegionsFile = "CTR_${params.genomeVersionHg}.bed.gz"
				chCtrRegionsBed = file("${params.annotationDir}/${ctrRegionsFile}").exists() ?
					Channel.fromPath("${params.annotationDir}/${ctrRegionsFile}").collect() : 
					DOWNLOAD_CTR_BED(params.ctrRegionsLink, ctrRegionsFile, 'DOWNLOAD_CTR_BED').file

			}
			
			// CNA files
			if (params.cna) {
				def armCoordinatesName = "UCSC_${params.genomeVersionHg}_arm_coordinates.txt"
				chArmCoordinates = file("${params.annotationDir}/${armCoordinatesName}").exists() ? 
					Channel.fromPath("${params.annotationDir}/${armCoordinatesName}").collect() : 
					PREPARE_ARM_COORDINATES(chCytobandBed, armCoordinatesName).txt
				
				if (!params.amplicon) {
					chAscetsResources = file(params.ascetsResources).exists() ? 
						Channel.fromPath(params.ascetsResources).collect() : DOWNLOAD_ASCETS().resources
				}
				
			}

		}
		
		// RNA resources
		chRnaTargetBed        = Channel.empty()
		chRnaTargetGenes      = Channel.empty()
		chRnaGenomeFastaFiles = Channel.empty()
		chRnaGenomeIndexDir   = Channel.empty()
		chCtatLibDir          = Channel.empty()
		chMitelmanFusion      = Channel.empty()
		if (params.RNA) {

			// Target files from panel manifests
			if (params.seqPanel) {
				if (file(params.rnaTargetBed).exists()) {
					chRnaTargetBed = Channel.fromPath(params.rnaTargetBed).collect()
				} else {
					
					// BED file from manifest
					if (params.manifestToBed) {
						MANIFEST_TO_BED_RNA(file(params.rnaManifest, checkIfExists: true), 'MANIFEST_TO_BED_RNA')
						chRnaTargetBed = MANIFEST_TO_BED_RNA.out.bed
					} else if (params.prepareIontorrentManifest) {
						PREPARE_MANIFEST_IONTORRENT_RNA(startingDir)
						chRnaTargetBed = PREPARE_MANIFEST_IONTORRENT_RNA.out.bed
					} else {
						chRnaTargetBed = Channel.fromPath(params.rnaManifest, checkIfExists: true).collect()
					}

					// Liftover BED file
					if (params.liftoverManifest) {
						BED_LIFTOVER_MANIFEST_RNA(chRnaTargetBed, chHg19ToHg38ChainFile, 'BED_LIFTOVER_MANIFEST_RNA')
						chRnaTargetBed = BED_LIFTOVER_MANIFEST_RNA.out.bed
					}

					// Annotate genes to the BED file
					if (params.annotateGenesBed) {
						BED_ANNOTATE_GENES_RNA(chRnaTargetBed, chManeGene, 'BED_ANNOTATE_GENES_RNA')
						chRnaTargetBed = BED_ANNOTATE_GENES_RNA.out.bed
					}

					// Prepare BED file
					PREPARE_TARGET_BED_RNA(chRnaTargetBed, 'RNA', 'PREPARE_TARGET_BED_RNA')
					chRnaTargetBed = PREPARE_TARGET_BED_RNA.out.bed

				}
				
				chRnaTargetGenes = file(params.rnaTargetGenes).exists() ? Channel.fromPath(params.rnaTargetGenes).collect() : 
					PREPARE_TARGET_GENES_RNA(chManeGene, chCytobandBed, chRnaTargetBed, 'RNA', 'PREPARE_TARGET_GENES_RNA').genes
			}
			

			// Genome files & CTAT library
			def ctatLib = "${params.ctatLibDir}/${params.ctatLibName}"
			def rnaGenomeStepsList = [params.alignment, params.deduplicationRna, params.realignment, params.fusion, params.splicing]
			if (rnaGenomeStepsList.any{it}) {
				def rnaCtatLibFilesList = [
					"${ctatLib}/ref_genome.fa", "${ctatLib}/ref_genome.fa.fai", "${ctatLib}/ref_genome.fa.star.idx", 
					"${ctatLib}/fusion_annot_lib.idx", "${ctatLib}/cancer_splicing_lib/cancer_splicing.idx"
				]
				if (rnaCtatLibFilesList.any{!(file(it).exists())}) {
					PREPARE_CTAT_LIB(params.ctatLibName, params.ctatLibLink, params.ctatAnnotFilterLink, params.ctatLibSplicingLink)
					chRnaGenomeFastaFiles = PREPARE_CTAT_LIB.out.fastaFiles
					chRnaGenomeIndexDir   = PREPARE_CTAT_LIB.out.indexDir
					chCtatLibDir          = PREPARE_CTAT_LIB.out.libDir
				} else {
					chRnaGenomeFastaFiles = Channel.fromPath("${ctatLib}/ref_genome*.{fa,fai}").collect()
					chRnaGenomeIndexDir   = Channel.fromPath("${ctatLib}/ref_genome.fa.star.idx").collect()
					chCtatLibDir          = Channel.fromPath("${ctatLib}").collect()
				}
			}

			// Fusion files
			if (params.fusion) {				
				def mitelmanFile = "MitelmanDB_GeneFusions_${params.mitelmanSearchDate}.txt"
				chMitelmanFusion = file("${params.annotationDir}/${mitelmanFile}").exists() ? 
					Channel.fromPath("${params.annotationDir}/${mitelmanFile}").collect() : 
					PREPARE_MITELMAN_FUSION(Channel.fromPath(params.mitelmanSearchResult, checkIfExists: true).collect(), mitelmanFile).file
			}

		}

		// Network of Cancer Genes (NCG)
		def ncgFile = "NCG_v${params.ncgVersion}_cancerdrivers"
		chCancerdrivers = file("${params.annotationDir}/${ncgFile}.processed.txt").exists() ? 
			Channel.fromPath("${params.annotationDir}/${ncgFile}.processed.txt").collect() : 
			PREPARE_CANCERDRIVERS(
				Channel.fromPath("${params.annotationDir}/${ncgFile}.tsv", checkIfExists: true).collect(), 
				"${ncgFile}.processed.txt").file
		
		// CIViC
		def civicFile      = "CIViC_${params.civicDate}"
		def civicVariant   = "CIViC_${params.civicDate}_VariantSummaries.tsv"
		def civicMolecular = "CIViC_${params.civicDate}_MolecularProfileSummaries.tsv"
		def civicEvidence  = "CIViC_${params.civicDate}_ClinicalEvidenceSummaries.tsv"
		if ([
			"${params.annotationDir}/${civicFile}_Oncogenic.txt", 
			"${params.annotationDir}/${civicFile}_Clinical.txt"].any{!(file(it).exists())}) {
			PREPARE_CIVIC_EVIDENCE(
				tumorNames, 
				DOWNLOAD_CIVIC_VARIANT(params.civicVariantLink, civicVariant, 'DOWNLOAD_CIVIC_VARIANT').file, 
				DOWNLOAD_CIVIC_MOLECULAR(params.civicMolecularLink, civicMolecular, 'DOWNLOAD_CIVIC_MOLECULAR').file, 
				DOWNLOAD_CIVIC_EVIDENCE(params.civicEvidenceLink, civicEvidence, 'DOWNLOAD_CIVIC_EVIDENCE').file)
			chCivicOncogenic  = PREPARE_CIVIC_EVIDENCE.out.oncogenic
			chCivicClinical = PREPARE_CIVIC_EVIDENCE.out.clinical
		} else {
			chCivicOncogenic  = Channel.fromPath("${params.annotationDir}/${civicFile}_Oncogenic.txt").collect()
			chCivicClinical = Channel.fromPath("${params.annotationDir}/${civicFile}_Clinical.txt").collect()
		}
			

		// PDX processing: filter mouse reads from FASTQs
		chMouseAndHumanIndex  = Channel.empty()
		if (params.pdx) {

			// Mouse FASTA files
			def mouseGenomeFastaDir = "Mus_musculus/NCBI/GRCm38/Sequence/WholeGenomeFasta"
			if ([
					"${params.genomeDir}/${mouseGenomeFastaDir}/genome.fa", 
					"${params.genomeDir}/${mouseGenomeFastaDir}/genome.fa.fai", 
					"${params.genomeDir}/${mouseGenomeFastaDir}/genome.dict"].any{!(file(it).exists())}) {
					chMouseGenomeFastaFiles = DOWNLOAD_GENOME_MOUSE(params.mouseGenomeLink, mouseGenomeFastaDir, "DOWNLOAD_GENOME_MOUSE").fastaFiles
				} else {
					chMouseGenomeFastaFiles = Channel.fromPath("${params.genomeDir}/${mouseGenomeFastaDir}/*.{dict,fa,fai,xml}").collect()
				}
			
			// Mouse index file
			def mouseGenomeIndexDir = "${params.genomeDir}/Mus_musculus/NCBI/GRCm38/Sequence/xengsortIndex/v${params.xengsortVersion}/MouseAndHuman"
			def mouseAndHumanFile = 'GRC38_MouseHuman.h5'
			chMouseAndHumanIndex = file("${mouseGenomeIndexDir}/${mouseAndHumanFile}").exists() ? 
				Channel.fromPath("${mouseGenomeIndexDir}/${mouseAndHumanFile}").collect() : 
				XENGSORT_INDEX(chMouseGenomeFastaFiles, chDnaGenomeFastaFiles, mouseGenomeIndexDir, mouseAndHumanFile).index
		
		}

		// Sample sheet (if the analysis is performed)
		chSampleSheet = Channel.empty()
		if (!params.resourcesOnly) {
			if (file(params.sampleSheet).exists()) {
				chSampleSheet = Channel.fromPath(params.sampleSheet).collect()
			} else if (file("${params.sampleSheetDir}/SampleSheet_${params.runName}.csv").exists()) {
				chSampleSheet = Channel.fromPath("${params.sampleSheetDir}/SampleSheet_${params.runName}.csv").collect()
			} else {
				chSampleSheet = PREPARE_SAMPLESHEET(startingDir).csv
			}
		}


	

	emit:
		maneGtf             = chManeGtf
		maneGene            = chManeGene
		maneExon            = chManeExon
		maneIntron          = chManeIntron
		maneCoding          = chManeCoding
		hg38ToHg19ChainFile = chHg38ToHg19ChainFile

		dnaBed                = chDnaTargetBed
		dnaIntervalList       = chDnaTargetIntervalList
		dnaBedExt             = chDnaTargetBedExt
		dnaIntervalListExt    = chDnaTargetIntervalListExt
		dnaGenes              = chDnaTargetGenes
		dnaChr                = chDnaTargetChr
		dnaGenomeFastaDir     = chDnaGenomeFastaDir
		dnaGenomeFastaFiles   = chDnaGenomeFastaFiles
		dnaGenomeIndexFiles   = chDnaGenomeIndexFiles
		dnaPanelHotspotsBed   = chDnaPanelHotspotsBed
		dnaPanelBlacklistBed  = chDnaPanelBlacklistBed
		dnaOffTargetBed       = chDnaOffTargetBed
		vepCacheDir           = chVepCacheDir
		vepCaddSnv            = chVepCaddSnv
		vepCaddIndel          = chVepCaddIndel
		vepRevel              = chVepRevel
		vepAlphaMissense      = chVepAlphaMissense
		clinvarVcf            = chClinvarVcf
		civicVcf              = chCivicVcf
		hotspotsWhitelistBed  = chHotspotsWhitelistBed
		cancerhotspotsResults = chCancerhotspotsResults
		cgiOncogenic          = chCgiOncogenic
		problematicRegionsBed = chProblematicRegionsBed
		ctrRegionsBed         = chCtrRegionsBed
		
		armCoordinates        = chArmCoordinates
		ascetsResources       = chAscetsResources

		rnaBed                = chRnaTargetBed
		rnaGenes              = chRnaTargetGenes
		rnaGenomeFastaFiles   = chRnaGenomeFastaFiles
		rnaGenomeIndexDir     = chRnaGenomeIndexDir
		ctatLibDir            = chCtatLibDir
		mitelmanFusion        = chMitelmanFusion

		cancerdrivers         = chCancerdrivers
		civicOncogenic        = chCivicOncogenic
		civicClinical         = chCivicClinical
		mouseAndHumanIndex    = chMouseAndHumanIndex

		sampleSheet           = chSampleSheet
	
}

