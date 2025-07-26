/*
 * DNA FASTQ realignment
 */

include { BAM_TO_FASTQ                       } from '../process/bamToFastq'
include { BAM_DEDUP_INDEX; BAM_REALIGN_INDEX } from '../process/bamIndex'
include { FASTQ_TO_BAM                       } from '../process/fastqToBam'
include { BWAMEM2_REALIGNMENT as BWAMEM2     } from '../process/bwamem2'
include { ABRA2                              } from '../process/abra2'
include { TMAP_REALIGNMENT as TMAP           } from '../process/tmap'
include { COLLECT_ALIGNMENT_METRICS          } from '../process/collectAlignmentMetrics'
include { COLLECT_INSERTSIZE_METRICS         } from '../process/collectInsertSizeMetrics'
include { COUNT_ONTARGET_READS               } from '../process/countOnTargetReads'
include { MOSDEPTH                           } from '../process/mosdepth'
include { MULTIQC                            } from '../process/multiQC'
include { RESULTS_BAM_QC                     } from '../process/resultsBamQC'


workflow REALIGNMENT_DNA {

	take:
		reads
		isInputBam
		genomeFastaFiles
		genomeIndexFiles
		targetBed
		bamHeaderInfo
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

		// 1. BWAMEM2 + ABRA2 or TMAP
		if (params.seqPlatform.toLowerCase() == 'iontorrent') {
			if (isInputBam) {
				BAM_DEDUP_INDEX(reads, '04_DEDUPLICATION', '02_BAM_INDEX', 'DNA')
				chInBam = reads.join(BAM_DEDUP_INDEX.out.bai)
			} else { 
				chInBam = FASTQ_TO_BAM(reads.combine(bamHeaderInfo, by: 0), '05_REALIGNMENT', '00_FASTQ_TO_BAM', 'DNA').bam
			}
			
			TMAP(chInBam, genomeIndexFiles, targetBed, '05_REALIGNMENT', '01_TMAP')
			chAlignedBam = TMAP.out.bam
			chAlignedBai = TMAP.out.bai

		} else {
			if (isInputBam) {
				// 0. BAM to FASTQ
				BAM_TO_FASTQ(reads, '05_REALIGNMENT', '00_BAM_TO_FASTQ', 'DNA')
				chFastq = BAM_TO_FASTQ.out.fastq
			} else { chFastq = reads }

			BWAMEM2(chFastq, genomeIndexFiles, '05_REALIGNMENT', '01_1_BWAMEM2')
			ABRA2(BWAMEM2.out.bam.join(BWAMEM2.out.bai), genomeFastaFiles, targetBed, '05_REALIGNMENT', '01_2_ABRA2', 'DNA')
			BAM_REALIGN_INDEX(ABRA2.out.bam, '05_REALIGNMENT', '01_3_BAM_INDEX', 'DNA')

			chAlignedBam = ABRA2.out.bam
			chAlignedBai = BAM_REALIGN_INDEX.out.bai
		}				

		// 2. BAM QC
		COLLECT_ALIGNMENT_METRICS(chAlignedBam, genomeFastaFiles, '05_REALIGNMENT', '02_COLLECT_ALIGNMENT_METRICS', 'DNA')
		MOSDEPTH(chAlignedBam.join(chAlignedBai), targetBed, '05_REALIGNMENT', '02_MOSDEPTH', 'DNA')
		if (!params.singleEnd) {
			COLLECT_INSERTSIZE_METRICS(chAlignedBam, genomeFastaFiles, '05_REALIGNMENT', '02_COLLECT_INSERTSIZE_METRICS', 'DNA')
			chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics.concat(COLLECT_INSERTSIZE_METRICS.out.metrics)
			chRunFiles = COLLECT_ALIGNMENT_METRICS.out.files.concat(COLLECT_INSERTSIZE_METRICS.out.files, MOSDEPTH.out.files)
		} else {
			chAlignmentMetrics = COLLECT_ALIGNMENT_METRICS.out.metrics
			chRunFiles = COLLECT_ALIGNMENT_METRICS.out.files.concat(MOSDEPTH.out.files)
		}
		COUNT_ONTARGET_READS(chAlignedBam, targetBed, '05_REALIGNMENT', '02_COUNT_ONTARGET_READS', 'DNA')
		MULTIQC(
			chAlignmentMetrics.concat(MOSDEPTH.out.metrics.transpose()).groupTuple(), 
			multiqcConfig, '05_REALIGNMENT', '03_MULTIQC_PICARD_MOSDEPTH', 'DNA')

		chAlignmentMetrics = COUNT_ONTARGET_READS.out.metrics
			.concat(chAlignmentMetrics, MOSDEPTH.out.perBaseBed.transpose(), bamQcRawResults)
			.groupTuple()
		
		RESULTS_BAM_QC(
			chAlignmentMetrics, samplesFile, targetBed, targetGenes, whitelistGenes, 
			maneGtf, maneGene, maneExon, maneCoding, '05_REALIGNMENT', '04_RESULTS_BAM_QC', 'DNA'
		)


	emit:
		bam           = chAlignedBam
		bai           = chAlignedBai
		perBaseCovBed = MOSDEPTH.out.perBaseBed
		results       = RESULTS_BAM_QC.out.results
		files         = RESULTS_BAM_QC.out.files
		runFiles      = chRunFiles
		
}
