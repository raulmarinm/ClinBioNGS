/*
 * Fusion analysis
 */

include { STAR_FUSION    } from '../process/starFusion'
include { RESULTS_FUSION } from '../process/resultsFusion'


workflow FUSION {

	take:
		bam
		junctions
		ctatLib
		targetBed
		targetGenes
		whitelistGenes
		whitelistFusions
		maneGtf
		maneGene
		maneExon
		maneIntron
		cancerdrivers
		mitelmanFusion
		genieFusion
		civicClinical
		chainFile
		genomeFastaFiles
		vcfHeader
		samplesFile
		perBaseCovBed
	
	main:

		// 1. STAR-fusion
		STAR_FUSION(bam.join(junctions), ctatLib, '06_FUSION', '01_STAR_FUSION')
		// 2. Fusion results
		RESULTS_FUSION(
			STAR_FUSION.out.results.transpose().concat(perBaseCovBed.transpose()).groupTuple(), 
			targetBed, targetGenes, whitelistGenes, whitelistFusions, maneGtf, maneGene, maneExon, maneIntron, 
			cancerdrivers, mitelmanFusion, genieFusion, civicClinical, chainFile, genomeFastaFiles, vcfHeader, 
			samplesFile, '06_FUSION', '02_RESULTS_FUSION')

	
	emit:
		results = RESULTS_FUSION.out.results
		files   = RESULTS_FUSION.out.files
		
}
