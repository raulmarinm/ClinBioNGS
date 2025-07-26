/*
 * Read deduplication
 */

include { GENCORE             } from '../process/gencore'
include { GATK_MARKDUPLICATES } from '../process/markDuplicates'
include { UMITOOLS            } from '../process/umitools'
include { 
	MULTIQC as MULTIQC_UMITOOLS;
	MULTIQC as MULTIQC_MARKDUP } from '../process/multiQC'


workflow DEDUPLICATION {

	take:
		bam
		bai
		genomeFastaFiles
		targetBed
		multiqcConfig
		type
	
	main:

		chDedupFiles = []
		isUmi = type.toUpperCase() == 'RNA' ? params.umiRna: params.umiDna

		// UMI deduplication
		//	- Single-end: UMI-tools
		//	- Paired-end: Gencore
		if (isUmi) {
			if (params.singleEnd | params.umitools) {
				UMITOOLS(bam.join(bai), '04_DEDUPLICATION', '01_UMITOOLS', type)
				MULTIQC_UMITOOLS(UMITOOLS.out.metrics, multiqcConfig, '04_DEDUPLICATION', '02_MULTIQC_UMITOOLS', type)
				
				chBam        = UMITOOLS.out.bam
				chDedupFiles = UMITOOLS.out.files
			} else {
				GENCORE(bam, genomeFastaFiles, targetBed, '04_DEDUPLICATION', '01_GENCORE', type)
				
				chBam = GENCORE.out.bam
			}
		} else {
			// Non-UMI deduplication: GATK Markduplicates
			GATK_MARKDUPLICATES(bam, '04_DEDUPLICATION', '01_GATK_MARKDUPLICATES', type)
			MULTIQC_MARKDUP(GATK_MARKDUPLICATES.out.metrics, multiqcConfig, '04_DEDUPLICATION', '02_MULTIQC_MARKDUP', type)
			
			chBam = GATK_MARKDUPLICATES.out.bam
			chDedupFiles = GATK_MARKDUPLICATES.out.files
		}

	emit:
		bam   = chBam
		files = chDedupFiles
}
