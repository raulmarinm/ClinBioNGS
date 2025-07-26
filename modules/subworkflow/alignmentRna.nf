/*
 * RNA FASTQ alignment
 */

include { STAR_ALIGNMENT as STAR     } from '../process/star'
include { COLLECT_ALIGNMENT_METRICS  } from '../process/collectAlignmentMetrics'
include { COLLECT_INSERTSIZE_METRICS } from '../process/collectInsertSizeMetrics'
include { COUNT_ONTARGET_READS       } from '../process/countOnTargetReads'
include { MULTIQC                    } from '../process/multiQC'
include { RESULTS_BAM_QC_RAW         } from '../process/resultsBamQC'


workflow ALIGNMENT_RNA {

	take:
		reads
		genomeFastaFiles
		genomeIndexDir
		targetBed
		multiqcConfig
	
	main:

		// 1. STAR
		STAR(reads, genomeIndexDir, '03_ALIGNMENT', '01_STAR')
		
		// 2. BAM QC
		COLLECT_ALIGNMENT_METRICS(STAR.out.bam, genomeFastaFiles, '03_ALIGNMENT', '02_COLLECT_ALIGNMENT_METRICS', 'RNA')
		chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics
		if(!params.singleEnd){
			COLLECT_INSERTSIZE_METRICS(STAR.out.bam, genomeFastaFiles, '03_ALIGNMENT', '02_COLLECT_INSERTSIZE_METRICS', 'RNA')
			chAlignmentMetrics = chAlignmentMetrics.concat(COLLECT_INSERTSIZE_METRICS.out.metrics).groupTuple()
		}
		COUNT_ONTARGET_READS(STAR.out.bam, targetBed, '03_ALIGNMENT', '02_COUNT_ONTARGET_READS', 'RNA')
		MULTIQC(chAlignmentMetrics, multiqcConfig, '03_ALIGNMENT', '03_MULTIQC_PICARD', 'RNA')
		RESULTS_BAM_QC_RAW(
			COUNT_ONTARGET_READS.out.metrics.concat(chAlignmentMetrics.transpose()).groupTuple(), 
			'03_ALIGNMENT', '04_RESULTS_BAM_QC_RAW', 'RNA')


	emit:
		bam      = STAR.out.bam
		bai      = STAR.out.bai
		results  = RESULTS_BAM_QC_RAW.out.results
		
}
