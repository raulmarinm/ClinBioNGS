/*
 * Read deduplication - UMI-tools dedup
 */

process UMITOOLS {

	tag "${sample}"

	label 'umitools'
	label 'minCpu'
	label 'medMem'

	publishDir(
		enabled: params.umitools_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.umitools_publishMode,
		pattern: "${sample}_${type}*.bam"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_umitools*",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('input/*'), path('input/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		tuple val(sample), path("${sample}_${type}_umitools*"), emit: metrics
		path "${sample}_${type}_umitools*", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def paired    = params.singleEnd ? "" : '--paired'
		def suffix    = params.bamDedup_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls input/*.bam)

		umi_tools dedup \\
			-I \${bam} \\
			-S ${sample}_${type}${suffix}.bam \\
			--output-stats ${sample}_${type}_umitools \\
			--log ${sample}_${type}_umitools.logs \\
			--umi-separator="${params.umitools_umiDelim}" \\
			--edit-distance-threshold ${params.umitools_distanceThreshold} \\
			${paired} \\
			${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
