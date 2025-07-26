/*
 * Prepare resources - Download files
 */

// Annotation files
process DOWNLOAD_ANNOT_FILE {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		val link
		val fileName
		val step

	output:
		path "${fileName}", emit: file
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "${step}"

		"""
		wget --no-check-certificate -O ${fileName} ${link}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
}
