/*
 * Prepare sample files
 * 	- Split the sample sheet into DNA/RNA Sample files
 */

process PREPARE_SAMPLE_FILES {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: { params.runDir },
		mode: 'copy',
		pattern: '*.{csv,txt}'
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path sampleSheet
		path sampleInfo
		path tumorNames
		val type

	output:
		path '*.txt', emit: txt
		path '*.csv', emit: csv
		path '*.{err,log,out,run,sh}'

	script:
		def panel     = params.seqPanel.toLowerCase()
		def isBcl     = params.startingDataType.toUpperCase() == 'BCL'
		def isUmi     = type.toUpperCase() == 'RNA' ? params.umiRna: params.umiDna
		def logPrefix = "PREPARE_SAMPLE_FILES_${type}"

		"""
		prepareSampleFiles.R \\
			--sampleSheet=${sampleSheet} \\
			--sampleInfo=${sampleInfo} \\
			--tumorNames=${tumorNames} \\
			--panel=${panel} \\
			--run=${params.runName} \\
			--type=${type} \\
			--isIontorrentSamplesheet=${params.prepareIontorrentSamplesheet} \\
			--isBcl=${isBcl} \\
			--isUmi=${isUmi}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

