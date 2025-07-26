/*
 * Prepare resources - Mitelman DB files
 */

process PREPARE_MITELMAN_FUSION{

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path inFile
		val outName

	output:
		path "${outName}", emit: file
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_MITELMAN_FUSION"

		"""
		prepareMitelmanFusion.R --in=${inFile} --out=${outName}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

