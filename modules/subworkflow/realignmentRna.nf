/*
 * RNA FASTQ realignment
 */

include { BAM_TO_FASTQ               } from '../process/bamToFastq'
include { STAR_REALIGNMENT as STAR   } from '../process/star'
include { COLLECT_ALIGNMENT_METRICS  } from '../process/collectAlignmentMetrics'
include { COLLECT_INSERTSIZE_METRICS } from '../process/collectInsertSizeMetrics'
include { COUNT_ONTARGET_READS       } from '../process/countOnTargetReads'
include { MOSDEPTH                   } from '../process/mosdepth'
include { MULTIQC                    } from '../process/multiQC'
include { RESULTS_BAM_QC             } from '../process/resultsBamQC'


workflow REALIGNMENT_RNA {

	take:
		reads
		isInputBam
		genomeFastaFiles
		genomeIndexDir
		targetBed
		multiqcConfig
		samplesFile
		targetGenes
		whitelistGenes
		maneGtf
		maneGene
		maneExon
		maneCoding
		bamQcRawResults
	
	main:

		chRunFiles = []

		if (isInputBam) {
			// 0. BAM to FASTQ
			BAM_TO_FASTQ(reads, '05_REALIGNMENT', '00_BAM_TO_FASTQ', 'RNA')
			chFastq = BAM_TO_FASTQ.out.fastq
		} else {
			chFastq = reads
		}

		// 1. STAR
		STAR(chFastq, genomeIndexDir, '05_REALIGNMENT', '01_STAR')

		// 2. BAM QC
		COLLECT_ALIGNMENT_METRICS(STAR.out.bam, genomeFastaFiles, '05_REALIGNMENT', '02_COLLECT_ALIGNMENT_METRICS', 'RNA')
		MOSDEPTH(STAR.out.bam.join(STAR.out.bai), targetBed, '05_REALIGNMENT', '02_MOSDEPTH', 'RNA')
		if (!params.singleEnd) {
			COLLECT_INSERTSIZE_METRICS(STAR.out.bam, genomeFastaFiles, '05_REALIGNMENT', '02_COLLECT_INSERTSIZE_METRICS', 'RNA')
			chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics.concat(COLLECT_INSERTSIZE_METRICS.out.metrics)
			chRunFiles = COLLECT_ALIGNMENT_METRICS.out.files.concat(COLLECT_INSERTSIZE_METRICS.out.files, MOSDEPTH.out.files)
		} else {
			chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics
			chRunFiles = COLLECT_ALIGNMENT_METRICS.out.files.concat(MOSDEPTH.out.files)
		}
		COUNT_ONTARGET_READS(STAR.out.bam, targetBed, '05_REALIGNMENT', '02_COUNT_ONTARGET_READS', 'RNA')
		MULTIQC(
			chAlignmentMetrics.concat(MOSDEPTH.out.metrics.transpose()).groupTuple(), 
			multiqcConfig, '05_REALIGNMENT', '03_MULTIQC_PICARD_MOSDEPTH', 'RNA')
		
		chAlignmentMetrics = COUNT_ONTARGET_READS.out.metrics
			.concat(chAlignmentMetrics, MOSDEPTH.out.perBaseBed.transpose(), bamQcRawResults)
			.groupTuple()

		RESULTS_BAM_QC(
			chAlignmentMetrics, samplesFile, targetBed, targetGenes, whitelistGenes, 
			maneGtf, maneGene, maneExon, maneCoding, '05_REALIGNMENT', '04_RESULTS_BAM_QC', 'RNA'
		)

	
	emit:
		bam           = STAR.out.bam
		bai           = STAR.out.bai
		SJ            = STAR.out.SJ
		junction      = STAR.out.junction
		perBaseCovBed = MOSDEPTH.out.perBaseBed
		results       = RESULTS_BAM_QC.out.results
		files         = RESULTS_BAM_QC.out.files
		runFiles      = chRunFiles
		
}
