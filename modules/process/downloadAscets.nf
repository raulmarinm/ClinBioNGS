/*
 * Prepare resources - Download ASCETS resources
 */

process DOWNLOAD_ASCETS {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outCnaDir'

	output:
		path 'ascets_resources.R', emit: resources
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "DOWNLOAD_ASCETS"

		"""
		wget --no-check-certificate https://raw.githubusercontent.com/beroukhim-lab/ascets/master/ascets_resources.R

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
