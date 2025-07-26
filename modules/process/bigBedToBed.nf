/*
 * Prepare resources - Convert big BED (bb) to BED files
 */

process BIGBED_TO_BED {

	tag "${params.runName}"

	label 'ucscBigbedtobed'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path bigBed
		val outBed
		val step

	output:
		path "${outBed}.gz", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "${step}"

		"""
		bigBedToBed ${bigBed} ${outBed}
		gzip ${outBed}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
