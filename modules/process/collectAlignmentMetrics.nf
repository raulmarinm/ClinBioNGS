/*
 * BAM quality control - CollectAlignmentSummaryMetrics
 */

process COLLECT_ALIGNMENT_METRICS {

	tag "${sample}"

	label 'gatk4'
	label 'minCpu'
	label 'medMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/BamQC_Reports",
		mode: 'copy',
		pattern: "${sample}_${type}_AlignmentSummaryMetrics.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam)
		path('genome/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}_AlignmentSummaryMetrics.txt"), emit: metrics
		path "${sample}_${type}_AlignmentSummaryMetrics.txt", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		fasta=\$(ls genome/*genome.fa)
		gatk CollectAlignmentSummaryMetrics \\
			-I ${bam} \\
			-R \${fasta} \\
			-O ${sample}_${type}_AlignmentSummaryMetrics.txt \\
			--java-options "-XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=1" \\
			${args}
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
