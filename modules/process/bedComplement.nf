/*
 * Prepare resources - Return complementary regions (not overlapping) from BED file
 */

process BED_COMPLEMENT {

	tag "${params.runName}"

	label 'bedtools'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		path('genome/*')
		val step

	output:
		path "*.bed", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def fileName  = bed.getSimpleName()
		def logPrefix = "${step}"

		"""
		fai=\$(ls genome/*.fai)
		bedtools complement -i ${bed} -g \${fai} > ${fileName}_complement.bed

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
