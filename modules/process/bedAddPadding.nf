/*
 * Prepare resources - Add BED padding
 */

process BED_ADD_PADDING {

	tag "${params.runName}"

	label 'bioawk'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		val padding
		val step

	output:
		path "*.bed", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def fileName  = bed.getSimpleName()
		def logPrefix = "${step}"

		"""
		bioawk -c bed '{print \$chrom,\$start-${padding},\$end+${padding},\$name}' ${bed} > ${fileName}_ext${padding}.bed


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
