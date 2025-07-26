/*
 * Read deduplication - Gencore
 */

process GENCORE {

	tag "${sample}"

	label 'gencore'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		enabled: params.gencore_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.gencore_publishMode,
		pattern: "${sample}_${type}*.bam"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_gencore.{html,json}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('input/*')
		path('genome/*')
		path bed
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		path "${sample}_${type}_gencore.{html,json}", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix    = params.bamDedup_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls input/*.bam)
		fasta=\$(ls genome/*genome.fa)

		gencore -i \${bam} \\
			-o ${sample}_${type}${suffix}.bam \\
			-r \${fasta} \\
			-b ${bed} \\
			--supporting_reads ${params.gencore_supportingReads} \\
			--ratio_threshold ${params.gencore_ratioThreshold} \\
			--score_threshold ${params.gencore_scoreThreshold} \\
			--umi_prefix "${params.gencore_umiDelim}" \\
			--umi_diff_threshold ${params.gencore_umiDiffThreshold} \\
			--json ${sample}_${type}_gencore.json \\
			--html ${sample}_${type}_gencore.html \\
			${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
