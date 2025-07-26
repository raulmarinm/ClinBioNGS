/*
 * MultiQC scripts
 */

// Collect BCL Convert metrics and make a report
process MULTIQC_BCL_CONVERT {

	tag "${params.runName}"

	label 'multiqc'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${subworkflow}",
		mode: 'copy',
		pattern: '*multiqc*',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path 'DNA/*'
		path 'RNA/*'
		path multiqc_config
		val subworkflow
		val step

	output:
		path '*multiqc*'
		path '*.{err,log,out,run,sh}'
	
	script:
		def args       = task.ext.args ?: ''
		def config     = multiqc_config ? "--config ${multiqc_config}" : ''
		def logPrefix  = "MULTIQC_BCL_CONVERT"

		"""
		multiqc --module bclconvert --force --title "${params.runName}" ${config} ${args} .

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

// Collect run metrics and make a report
process MULTIQC_RUN {

	tag "${params.runName}"

	label 'multiqc'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: params.resultsDir,
		mode: 'copy',
		pattern: '*multiqc*',
	)
	publishDir(
		path: { "${params.logsDir}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path 'input/*'
		path multiqc_config
		val step
		val type

	output:
		path '*multiqc*'
		path '*.{err,log,out,run,sh}'
	
	script:
		def args       = task.ext.args ?: ''
		def config     = multiqc_config ? "--config ${multiqc_config}" : ''
		def run        = params.runName
		def logPrefix  = "${step}_${type}"

		"""
		multiqc --force --title "${run}_${type}" ${config} ${args} .

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

// Collect sample metrics and make a report
process MULTIQC {

	tag "${sample}"

	label 'multiqc'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: '*multiqc*',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('reports/*')
		path multiqc_config
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}_multiqc_report_data"), emit: data
		path '*multiqc*'
		path '*.{err,log,out,run,sh}'
	
	script:
		def args       = task.ext.args ?: ''
		def config     = multiqc_config ? "--config ${multiqc_config}" : ''
		def logPrefix  = "${step}_${type}_${sample}"

		"""
		multiqc --force --title "${sample}_${type}" ${config} ${args} .

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
