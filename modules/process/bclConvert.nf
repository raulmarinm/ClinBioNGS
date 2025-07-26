/*
 * FASTQ generation - BCL Convert
 * 	- Convert BCL to FASTQ files
 * 	- Change FASTQ names and report them
 */

process BCL_CONVERT {

	tag "${params.runName}"

	label 'bclconvert'
	label 'highCpu'
	label 'lowMem'

	publishDir(
		enabled: params.bclConvert_publish,
		path: "${params.dataRunDir}/BCL_CONVERT/",
		mode: params.bclConvert_publishMode,
		pattern: '*.fastq*'
	)
	publishDir(
		path: "${params.analysisDir}/${subworkflow}/${step}_${type}",
		mode: 'copy',
		pattern: '{Logs,Reports}/*',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path bclDir
		path sampleSheet
		val subworkflow
		val step
		val type

	output:
		path '*.fastq*', emit: fastq
		path 'Reports/*', emit: report
		path 'Logs/*'
		path '*.{err,log,out,run,sh}'
	
	script:
		def threads    = task.cpus < 3 ? '1' : task.cpus.intdiv(3)
		def args       = task.ext.args ?: ''
		def logPrefix  = "${step}_${type}"

		"""
		bcl-convert --bcl-input-directory ${bclDir} --output-directory . --sample-sheet ${sampleSheet} --no-lane-splitting true \\
		--bcl-num-parallel-tiles 1 --bcl-num-conversion-threads ${threads} --bcl-num-compression-threads ${threads} \\
		--bcl-num-decompression-threads ${threads} --force --strict-mode true ${args}

		rm -f ./core*
		rm -f ./Undetermined*.fastq*

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
