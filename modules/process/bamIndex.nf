/*
 * Index BAM file
 */

process BAM_DEDUP_INDEX {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.bamDedupIndex_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bamDedupIndex_publishMode,
		pattern: "${sample}_${type}*.bai"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam)
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools index ${bam} ${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process BAM_REALIGN_INDEX {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.bamRealignIndex_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bamRealignIndex_publishMode,
		pattern: "${sample}_${type}*.bai"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam)
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools index ${bam} ${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
