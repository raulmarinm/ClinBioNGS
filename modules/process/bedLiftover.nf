/*
 * Prepare resources - BED liftover
 */

process BED_LIFTOVER_MANIFEST {

	tag "${params.runName}"

	label 'ucscLiftover'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		path chain
		val step

	output:
		path "*.lifted.bed", emit: bed
		path '*.{bed,err,log,out,run,sh}'

	script:
		def fileName  = bed.getSimpleName()
		def logPrefix = "${step}"

		"""
		liftOver ${bed} ${chain} ${fileName}.lifted.bed ${fileName}.unlifted.bed
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process BED_LIFTOVER_VARIANT {

	tag "${params.runName}"

	label 'ucscLiftover'
	label 'minCpu'
	label 'minMem'
	label 'outSmallVariantDir'

	input:
		path bed
		path chain
		val step

	output:
		path "*.lifted.bed", emit: bed
		path '*.{bed,err,log,out,run,sh}'

	script:
		def fileName  = bed.getSimpleName()
		def logPrefix = "${step}"

		"""
		liftOver ${bed} ${chain} ${fileName}.lifted.bed ${fileName}.unlifted.bed
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
