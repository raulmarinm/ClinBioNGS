/*
 * MSI analysis - MSIsensor-PRO
 */

process MSISENSORPRO {

	tag "${sample}"

	label 'msisensorpro'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_DNA_MSI*",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		path baseline
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA_MSI_results.txt"), emit: results
		path "${sample}_DNA_MSI*"
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		msisensor-pro pro \\
			-b ${task.cpus} \\
			-d ${baseline} \\
			-t ${bam} \\
			-e ${bed} \\
			-o ${sample}_DNA_MSI \\
			-c ${params.msisensorpro_minCov} \\
			${args}

		mv ${sample}_DNA_MSI ${sample}_DNA_MSI_results.txt
		mv ${sample}_DNA_MSI_all ${sample}_DNA_MSI_allSites.txt
		mv ${sample}_DNA_MSI_unstable ${sample}_DNA_MSI_unstableSites.txt
		mv ${sample}_DNA_MSI_dis ${sample}_DNA_MSI.dist

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
