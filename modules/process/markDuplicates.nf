/*
 * Read deduplication - GATK MarkDuplicates
 */

process GATK_MARKDUPLICATES {

	tag "${sample}"

	label 'gatk4'
	label 'minCpu'
	label 'highMem'

	publishDir(
		enabled: params.gatkMarkduplicates_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.gatkMarkduplicates_publishMode,
		pattern: "${sample}_${type}*.bam"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_MarkDuplicates.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('input/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		tuple val(sample), path("${sample}_${type}_MarkDuplicates.txt"), emit: metrics
		path "${sample}_${type}_MarkDuplicates.txt", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix    = params.bamDedup_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls input/*.bam)
		gatk MarkDuplicates \\
			-I \${bam} \\
			-O ${sample}_${type}${suffix}.bam \\
			--REMOVE_DUPLICATES true \\
			--METRICS_FILE ${sample}_${type}_MarkDuplicates.txt \\
			--java-options "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=1" \\
			${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
