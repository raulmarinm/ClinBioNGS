/*
 * Merge BED files
 */

process BED_MERGE {

	tag "${params.runName}"

	label 'bedtools'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path inBed
		val outBed
		val step

	output:
		path "${outBed}.gz", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "${step}"

		"""
		zcat ${inBed} | cut -f1-3 | sort -V -k1,1 -k2,2 | uniq | bedtools merge -i stdin > ${outBed}
		gzip ${outBed}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
