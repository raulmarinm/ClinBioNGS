/*
 * Make FASTQ quality control
 */

include { FASTQC  } from '../process/fastqc'
include { MULTIQC } from '../process/multiQC'


workflow FASTQ_QC {

	take:
		reads
		multiqcConfig
		type
	
	main:
		FASTQC(reads, "02_FASTQ_QC", "01_FASTQC", type)
		MULTIQC(FASTQC.out.report, multiqcConfig, "02_FASTQ_QC", '02_MULTIQC_FASTQC', type)

	emit:
		files = FASTQC.out.files.collect()

}
