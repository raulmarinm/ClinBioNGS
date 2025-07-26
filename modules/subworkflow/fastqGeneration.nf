/*
 * Generate FASTQ files
 */

include { BCL_CONVERT       } from '../process/bclConvert'
include { 
	EXTRACT_UBAM_IONTORRENT; 
	PREPARE_UBAM_IONTORRENT } from '../process/prepareIonTorrentFiles.nf'
include { BAM_ADD_UMI_QNAME } from '../process/bamAddUmiQname'
include { UBAM_TO_FASTQ     } from '../process/uBamToFastq'


workflow FASTQ_GENERATION {

	take:
		startingDir
		samplesFile
		sampleSheet
		type
	
	main:

		chSamples = samplesFile.splitCsv(sep: '\t', header: true).map(row -> "${row.Sample}")
		isUmi     = type.toUpperCase() == "RNA" ? params.umiRna: params.umiDna

		chBclConvertReport = []
		chBamHeaderInfo    = Channel.empty()

		// Start from FASTQ files
		if (params.startingDataType.toUpperCase() == 'FASTQ') {
			chStartingFastq = Channel.fromPath("${params.startingDataDir}/**_${type}*.fastq*")
				.map{ [it.simpleName.toString().split("_${type}")[0], it] }
				.groupTuple()
			chFastq = chSamples.map{ [it] }.combine(chStartingFastq, by: 0)
			
		// Convert BCL to FASTQ files
		} else if (params.startingDataType.toUpperCase() == 'BCL') {
			// 1. BCL_CONVERT
			BCL_CONVERT(startingDir, sampleSheet, '00_FASTQ_GENERATION', '01_BCL_CONVERT', type)
			chBclConvertReport = BCL_CONVERT.out.report.collect()
			chBclConvertFastq  = BCL_CONVERT.out.fastq
				.map { [it.simpleName, it] }
				.transpose()
				.map { [it[0].toString().split("_${type}")[0], it[1]] }
				.groupTuple()
			chFastq = chSamples.map{ [it] }.combine(chBclConvertFastq, by: 0)		
		
		// Convert BAM to FASTQ files
		} else if (params.startingDataType.toUpperCase() == 'BAM') {
			// 0. PREPARE UBAM
			if(params.prepareIontorrentBam) {
				chSamplesBarcode  = samplesFile.splitCsv(sep: '\t', header: true).map(row -> ["${row.Sample}", "${row.Barcode}"])
				EXTRACT_UBAM_IONTORRENT(startingDir, '00_FASTQ_GENERATION', '00_EXTRACT_UBAM_IONTORRENT', type)
				chStartingBam = EXTRACT_UBAM_IONTORRENT.out.bam.collect()
				PREPARE_UBAM_IONTORRENT(chSamplesBarcode, chStartingBam, '00_FASTQ_GENERATION', '01_PREPARE_UBAM_IONTORRENT', type)
				chBam = PREPARE_UBAM_IONTORRENT.out.bam
			} else { 
				chStartingBam = Channel.fromPath("${params.startingDataDir}/**_${type}*.bam")
					.map{ [it.simpleName.toString().split("_${type}")[0], it] }
					.groupTuple()
				chBam = chSamples.map{ [it] }.combine(chStartingBam, by: 0)
			}
			// 1. BAM_ADD_UMI_QNAME (if Umi) + UBAM_TO_FASTQ
			if (isUmi && params.bamAddUmiQname) {
				BAM_ADD_UMI_QNAME(chBam, '00_FASTQ_GENERATION', '01_BAM_ADD_UMI_QNAME', type)
				UBAM_TO_FASTQ(BAM_ADD_UMI_QNAME.out.bam, '00_FASTQ_GENERATION', '02_UBAM_TO_FASTQ', type)
			} else { 
				UBAM_TO_FASTQ(chBam, '00_FASTQ_GENERATION', '01_UBAM_TO_FASTQ', type)
			}
			chFastq         = UBAM_TO_FASTQ.out.fastq
			chBamHeaderInfo = UBAM_TO_FASTQ.out.headerInfo
		
		} else { chFastq = Channel.empty() }

	
	emit:
		fastq            = chFastq
		bclConvertReport = chBclConvertReport
		bamHeaderInfo    = chBamHeaderInfo
}

