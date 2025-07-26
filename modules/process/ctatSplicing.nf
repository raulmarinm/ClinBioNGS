/*
 * Splicing analysis - CTAT-splicing
 */

process CTAT_SPLICING {

	tag "${sample}"

	label 'ctatSplicing'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/CTAT-Splicing",
		mode: 'copy',
		pattern: "${sample}_RNA*.{introns,prelim}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.command.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(SJ), path(chimJ)
		path(ctatLib)
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_RNA*.{introns,prelim}"), emit: results
		path '*.command.{err,log,out,run,sh}'
	
	script:
		def args             = task.ext.args ?: ''
		def logPrefix        = "${step}_RNA_${sample}"

		"""
		/usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py --ctat_genome_lib ${ctatLib} --sample_name ${sample}_RNA \\
			--SJ_tab_file ${SJ} --chimJ_file ${chimJ} --output_prefix ${sample}_RNA ${args}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
