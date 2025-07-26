/*
 * Prepare resources - UCSC arm genomic coordinates
 */

process PREPARE_ARM_COORDINATES {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path cytoband
		val fileName

	output:
		path '*.txt', emit: txt
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_ARM_COORDINATES"

		"""
		prepareArmCoordinates.R --cytoband=${cytoband} --fileName=${fileName}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

