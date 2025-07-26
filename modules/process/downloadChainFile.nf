/*
 * Prepare resources - Download UCSC LiftOver chain file
 */

process DOWNLOAD_CHAIN_FILE {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${dir}",
		mode: 'copy',
		pattern: "${name}.over.chain"
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		val genome
		val name
		val dir

	output:
		path "${name}.over.chain", emit: chain
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "DOWNLOAD_CHAIN_FILE_${name.toUpperCase()}"

		"""
		wget --no-check-certificate https://hgdownload.cse.ucsc.edu/goldenpath/${genome}/liftOver/${name}.over.chain.gz
		gunzip ${name}.over.chain.gz

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
