/*
 * Splicing analysis
 */

include { CTAT_SPLICING    } from '../process/ctatSplicing'
include { RESULTS_SPLICING } from '../process/resultsSplicing'


workflow SPLICING {

	take:
		juntions
		ctatLib
		targetBed
		targetGenes
		whitelistGenes
		whitelistSplicing
		maneGtf
		maneGene
		maneExon
		maneIntron
		cancerdrivers
		civicClinical
		chainFile
		genomeFastaFiles
		vcfHeader
		samplesFile
		perBaseCovBed
		variants
	
	main:

		// 1. STAR-fusion
		CTAT_SPLICING(juntions, ctatLib, '06_SPLICING', '01_CTAT_SPLICING')

		// 2. Splicing results
		RESULTS_SPLICING(
			CTAT_SPLICING.out.results.transpose().concat(variants, perBaseCovBed.transpose()).groupTuple(), 
			targetBed, targetGenes, whitelistGenes, whitelistSplicing, maneGtf, maneGene, maneExon, maneIntron, 
			cancerdrivers, civicClinical, chainFile, genomeFastaFiles, vcfHeader, 
			samplesFile, '06_SPLICING', '02_RESULTS_SPLICING')

	emit:
		results = RESULTS_SPLICING.out.results
		files   = RESULTS_SPLICING.out.files
		
}
