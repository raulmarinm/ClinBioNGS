/*
 * CNA analysis
 */

include { CNVKIT      } from '../process/cnvkit'
include { ASCETS      } from '../process/ascets'
include { RESULTS_CNA } from '../process/resultsCna'


workflow CNA {

	take:
		bam
		targetBed
		baseline
		ascetsResources
		armCoordinates
		genomeFastaFiles
		targetGenes
		whitelistGenes
		maneGtf
		maneGene
		maneExon
		cancerdrivers
		genie
		civicClinical
		vcfHeader
		samplesFile
		perBaseCovBed
	
	main:

		// 1. CNVKIT
		CNVKIT(bam, baseline, '06_CNA', '01_CNVKIT')

		// 2. ASCETS
		if (params.amplicon) { 
			cnaResults = Channel.empty()
		} else {
			ASCETS(CNVKIT.out.segments, ascetsResources, armCoordinates, samplesFile, '06_CNA', '02_ASCETS')
			cnaResults = ASCETS.out.results
		}
		
		// 3. CNA results
		RESULTS_CNA(
			cnaResults.concat(CNVKIT.out.bins, CNVKIT.out.sex, perBaseCovBed.transpose()).groupTuple(), 
			targetBed, armCoordinates, genomeFastaFiles,
			targetGenes, whitelistGenes, maneGtf, maneGene, maneExon, 
			cancerdrivers, genie, civicClinical, vcfHeader, samplesFile, 
			'06_CNA', params.amplicon ? '02_RESULTS_CNA' : '03_RESULTS_CNA'
		)

	
	emit:
		geneResults = RESULTS_CNA.out.results
		armResults  = cnaResults
		files       = RESULTS_CNA.out.files.concat(CNVKIT.out.sex)
		
}
