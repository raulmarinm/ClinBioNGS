/*
 * Prepare singularity images - Download or Build them
 */

process SIF_DOWNLOAD {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outSingularityDir'

	input:
		val link
		val tool
		val version

	output:
		path '*.sif'
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "SIF_DOWNLOAD_${tool}_${version}"

		"""
		wget --no-check-certificate -O ${tool}_${version}.sif ${link}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process SIF_BUILD {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outSingularityDir'

	input:
		val container
		val tool
		val version

	output:
		path '*.sif'
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "SIF_BUILD_${tool}_${version}"

		if( params.apptainer )
		"""
		apptainer build ${tool}_${version}.sif ${container}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

		else
		"""
		singularity build ${tool}_${version}.sif ${container}
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
