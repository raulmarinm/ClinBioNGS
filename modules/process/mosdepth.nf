/*
 * BAM quality control - Mosdepth
 */

process MOSDEPTH {

	tag "${sample}"

	label 'mosdepth'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/BamQC_Reports",
		mode: 'copy',
		pattern: "${sample}_${type}*.{bed,txt}*",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*per-base.bed.gz"), emit: perBaseBed
		tuple val(sample), path("${sample}_${type}*.{bed,txt}*"), emit: metrics
		path "${sample}_${type}*.{bed,txt}*", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def useMedian = params.mosdepth_useMedian ? '--use-median' : ''
		def fileName  = bam.getBaseName()
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		mosdepth \\
			--by ${bed} \\
			--threads ${task.cpus} \\
			${fileName} \\
			${bam} \\
			--mapq ${params.mosdepth_mapq} \\
			${useMedian} \\
			${args}
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
