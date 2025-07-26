/*
 * Prepare resources - CGI oncogenic mutations
 */

process PREPARE_CGI_MUTATIONS{

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path chain
		path inFile
		val outName

	output:
		path "${outName}", emit: file
		path '*.{err,log,out,run,sh,tsv}'

	script:
		def logPrefix = "PREPARE_CGI_MUTATIONS"

		"""
		prepareCgiMutations.R --chain=${chain} --in=${inFile} --out=${outName}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}
