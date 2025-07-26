/*
 * Small variant analysis - Annotation
 */

include { VEP_RUN; VEP_PROCESSING } from '../process/vep'
include { RESULTS_SMALL_VARIANT   } from '../process/resultsSmallVariant'


workflow VARIANT_ANNOTATION {

	take:
		vcf
		consensusResults
		vepCacheDir
		caddSnv
		caddIndel
		vepRevel
		vepAlphaMissense
		clinvarVcf
		civicVcf
		hotspotsWhitelistBed
		problematicRegionsBed
		ctrRegionsBed
		mmrGenes
		cancerdrivers
		cancerHotspotsResults
		sopOncogenic
		genieCounts
		genieOncogenic
		cgiOncogenic
		civicOncogenic
		civicClinical
		panelRecurrentMutations
		panelHotspotsBed
		targetGenes
		whitelistGenes
		maneGtf
		maneGene
		maneExon
		vcfHeader
		samplesFile
	
	main:

		// 1. VEP
		def outVcf = "VEP${params.vepVersion}_${params.runName}_DNA_SmallVariant.vcf"
		VEP_RUN(
			vcf, outVcf, vepCacheDir, caddSnv, caddIndel, vepRevel, vepAlphaMissense, clinvarVcf, civicVcf, 
			'07_VARIANT_ANNOTATION', '01_VEP_RUN'
		)		

		// 2. VEP processing
		VEP_PROCESSING(
			VEP_RUN.out.vcf, hotspotsWhitelistBed, problematicRegionsBed, ctrRegionsBed, mmrGenes, 
			cancerdrivers, cancerHotspotsResults, sopOncogenic, genieCounts, genieOncogenic, 
			cgiOncogenic, civicOncogenic, panelRecurrentMutations,
			'07_VARIANT_ANNOTATION', '02_VEP_PROCESSING'
		)

		// 3. Small variant results
		RESULTS_SMALL_VARIANT(
			consensusResults, VEP_PROCESSING.out.results, targetGenes, whitelistGenes, 
			maneGtf, maneGene, maneExon, genieOncogenic, civicClinical, 
			panelHotspotsBed, vcfHeader, samplesFile,
			'07_VARIANT_ANNOTATION', '03_RESULTS_SMALL_VARIANT'
		)
	
	emit:
		results = RESULTS_SMALL_VARIANT.out.results
		files   = RESULTS_SMALL_VARIANT.out.files

}