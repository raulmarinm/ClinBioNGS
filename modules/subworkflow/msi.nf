/*
 * MSI analysis
 */

include { MSISENSORPRO } from '../process/msisensorpro'


workflow MSI {

	take:
		bam
		targetBed
		baseline
	
	main:

		// 1. CNVKIT
		MSISENSORPRO(bam, targetBed, baseline, '06_MSI', '01_MSISENSORPRO')

	emit:
		results = MSISENSORPRO.out.results
		
}
