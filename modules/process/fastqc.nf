/*
 * FASTQ quality control - FastQC
 */

process FASTQC {

	tag "${sample}"

	label 'fastqc'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/FastQC_Reports",
		mode: 'copy',
		pattern: "${sample}_${type}*.{html,zip}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('data/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.{html,zip}"), emit: report
		path "${sample}_${type}*.{html,zip}", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		for fq in \$(ls data/*.fastq*)
		do
		fastqc \${fq} --outdir . ${args}
		done
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
