/*
 * Prepare resources - BED to interval list
 */

process BED_TO_INTERVAL_LIST {

	tag "${params.runName}"

	label 'gatk4'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		path('genome/*')
		val step

	output:
		path "*.interval_list", emit: interval_list
		path '*.{err,log,out,run,sh}'

	script:
		def panel     = params.seqPanel.toUpperCase()
		def genome    = params.genomeVersionHg
		def fileName  = bed.getSimpleName()
		def args = task.ext.args ?: ''
		def logPrefix = "${step}"

		"""
		dict=\$(ls genome/*.dict)
		gatk BedToIntervalList -I ${bed} -O ${fileName}.interval_list -SD \${dict} ${args}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
