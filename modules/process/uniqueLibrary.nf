/*
 * Create or update a unique varaint library 
 */

process UNIQUE_LIBRARY {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.projectDir}",
		mode: 'copy',
		pattern: '*.db',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path 'input/*'
		path runLibrary
		path projectDir
		val subworkflow
		val step

	output:
		path '*.db', emit: db
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix  = "${step}"

		"""
		uniqueLibrary.R \\
			--run=${params.runName} \\
			--libraryName=${params.libraryName} \\
			--runLibrary=${runLibrary} \\
			--resultsDir=${params.resultsDir}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

