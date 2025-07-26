/*
 * Process FASTQ files
 */

include { FASTQ_MERGE_LANES    } from '../process/fastqMergeLanes'
include { FASTQ_UMI_TRANSFER   } from '../process/fastqUmiTransfer'
include { XENGSORT_CLASSIFY    } from '../process/xengsort'
include { FASTQ_CHANGE_UMI_SEP } from '../process/fastqChangeUmiSep'
include { FASTP                } from '../process/fastp'


workflow FASTQ_PROCESSING {

	take:
		reads
		samplesFile
		mouseAndHumanIndex
		type
	
	main:

		chFastq   = reads
		chSamples = samplesFile.splitCsv(sep: '\t', header: true).map(row -> "${row.Sample}")
		isUmi     = type.toUpperCase() == "RNA" ? params.umiRna: params.umiDna

		// FASTQ_MERGE_LANES: merge FASTQ files
		if (params.fastqMergeLanes) {
			FASTQ_MERGE_LANES(chFastq, '01_FASTQ_PROCESSING', 'FASTQ_MERGE_LANES', type)
			chFastq = chSamples.map{ [it] }.combine(FASTQ_MERGE_LANES.out.fastq, by: 0)
		}

		// FASTQ_UMI_TRANSFER: transfer UMI from a separate FASTQ to the new FASTQ header
		if (isUmi && params.fastqUmiTransfer) {
			FASTQ_UMI_TRANSFER(chFastq, '01_FASTQ_PROCESSING', 'FASTQ_UMI_TRANSFER', type)
			chFastq = chSamples.map{ [it] }.combine(FASTQ_UMI_TRANSFER.out.fastq, by: 0)
		}

		// PDX processing: filter mouse reads from FASTQs with XENGSORT_CLASSIFY
		if (params.pdx) {
			XENGSORT_CLASSIFY(chFastq, mouseAndHumanIndex, '01_FASTQ_PROCESSING', 'FILTER_MOUSE_READS', type)
			chFastq = chSamples.map{ [it] }.combine(XENGSORT_CLASSIFY.out.fastq, by: 0)
		}

		// FASTQ_CHANGE_UMI_SEP: Change UMI separator from the FASTQ header
		if (isUmi && params.fastqChangeUmiSep) {
			FASTQ_CHANGE_UMI_SEP(chFastq.transpose(), '01_FASTQ_PROCESSING', 'FASTQ_CHANGE_UMI_SEP', type)
			chFastq = FASTQ_CHANGE_UMI_SEP.out.fastq
				.map{ a, b -> [a, b.simpleName] }.groupTuple(sort: true).transpose()
				.join(FASTQ_CHANGE_UMI_SEP.out.fastq.map{ a, b -> [a, b.simpleName, b] }, by: [0,1])
				.map{ a, b, c -> [a, c] }.groupTuple()
		}

		// FASTP: preprocessing for FASTQ files
		if (params.fastp) {
			FASTP(chFastq, '01_FASTQ_PROCESSING', 'FASTP', type)
			chFastq = chSamples.map{ [it] }.combine(FASTP.out.fastq, by: 0)
		}

		
	emit:
		fastq = chFastq
		files = FASTP.out.files
	
}
