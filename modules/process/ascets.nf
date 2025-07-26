/*
 * CNA analysis - Arm-level Somatic Copy-number Events in Targeted Sequencing (ASCETS)
 */

process ASCETS {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "{Ascets_Reports,${sample}_DNA_CNA_arm.txt}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(segments), path(sex)
		path resources
		path armCoordinates
		path samplesFile
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA_CNA_arm.txt"), emit: results
		path 'Ascets_Reports'
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix = "${step}_DNA_${sample}"

		"""
		ascets.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=DNA \\
			--samplesFile=${samplesFile} \\
			--resources=${resources} \\
			--segments=${segments} \\
			--sex=${sex} \\
			--armCoordinates=${armCoordinates} \\
			--minBoc=${params.ascets_minBoc} \\
			--callingThreshold=${params.ascets_callingThreshold} \\
			--fractionThreshold=${params.ascets_fractionThreshold}
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

