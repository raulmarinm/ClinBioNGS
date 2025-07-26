/*
 * DNA FASTQ alignment
 */

include { BWAMEM2_ALIGNMENT as BWAMEM2 } from '../process/bwamem2'
include { FASTQ_TO_BAM                 } from '../process/fastqToBam'
include { TMAP_ALIGNMENT as TMAP       } from '../process/tmap'
include { COLLECT_ALIGNMENT_METRICS    } from '../process/collectAlignmentMetrics'
include { COLLECT_INSERTSIZE_METRICS   } from '../process/collectInsertSizeMetrics'
include { COUNT_ONTARGET_READS         } from '../process/countOnTargetReads'
include { MULTIQC                      } from '../process/multiQC'
include { RESULTS_BAM_QC_RAW           } from '../process/resultsBamQC'


workflow ALIGNMENT_DNA {

	take:
		reads
		genomeFastaFiles
		genomeIndexFiles
		targetBed
		bamHeaderInfo
		multiqcConfig
	
	main:

		// 1. BWAMEM2 or TMAP
		if (params.seqPlatform.toLowerCase() == 'iontorrent') {
			FASTQ_TO_BAM(reads.combine(bamHeaderInfo, by: 0), '03_ALIGNMENT', '00_FASTQ_TO_BAM', 'DNA')
			TMAP(FASTQ_TO_BAM.out.bam, genomeIndexFiles, targetBed, '03_ALIGNMENT', '01_TMAP')
			chAlignedBam = TMAP.out.bam
			chAlignedBai = TMAP.out.bai
		} else {
			BWAMEM2(reads, genomeIndexFiles, '03_ALIGNMENT', '01_BWAMEM2')
			chAlignedBam = BWAMEM2.out.bam
			chAlignedBai = BWAMEM2.out.bai
		}
				

		// 2. BAM QC
		COLLECT_ALIGNMENT_METRICS(chAlignedBam, genomeFastaFiles, '03_ALIGNMENT', '02_COLLECT_ALIGNMENT_METRICS', 'DNA')
		chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics
		if (!params.singleEnd) {
			COLLECT_INSERTSIZE_METRICS(chAlignedBam, genomeFastaFiles, '03_ALIGNMENT', '02_COLLECT_INSERTSIZE_METRICS', 'DNA')
			chAlignmentMetrics = chAlignmentMetrics.concat(COLLECT_INSERTSIZE_METRICS.out.metrics).groupTuple()
		}
		COUNT_ONTARGET_READS(chAlignedBam, targetBed, '03_ALIGNMENT', '02_COUNT_ONTARGET_READS', 'DNA')
		MULTIQC(chAlignmentMetrics, multiqcConfig, '03_ALIGNMENT', '03_MULTIQC_PICARD', 'DNA')
		RESULTS_BAM_QC_RAW(
			COUNT_ONTARGET_READS.out.metrics.concat(chAlignmentMetrics.transpose()).groupTuple(), 
			'03_ALIGNMENT', '04_RESULTS_BAM_QC_RAW', 'DNA')

	
	emit:
		bam     = chAlignedBam
		bai     = chAlignedBai
		results = RESULTS_BAM_QC_RAW.out.results
		
}
