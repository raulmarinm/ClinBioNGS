/*
 * Prepare resources - CTAT library
 */

process PREPARE_CTAT_LIB {

	tag "${params.runName}"

	label 'ctatSplicing'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: params.ctatLibDir,
		mode: params.prepareCtatLib_publishMode,
		pattern: "${libName}"
	)

	input:
		val libName
		val libLink
		val annotFilterLink
		val libSplicingLink

	output:
		path "${libName}", emit: libDir
		path "${libName}/*genome*.{fa,fai}", emit: fastaFiles
		path "${libName}/*star.idx", emit: indexDir
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_CTAT_LIB"

		"""
		wget --no-check-certificate ${libLink}

		libFile=\$(basename ${libLink})
		tar -zxvf \${libFile}

		libDir=\$(echo "\${libFile%.tar.gz}")
		mv \${libDir}/ctat_genome_lib_build_dir ${libName}
		rmdir \${libDir}

		wget --no-check-certificate ${annotFilterLink}
		annotFilterFile=\$(basename ${annotFilterLink})
		mv \${annotFilterFile} ${libName}/

		wget --no-check-certificate ${libSplicingLink}
		intronsFile=\$(basename ${libSplicingLink})

		/usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py --cancer_introns_tsv \${intronsFile} --genome_lib_dir ${libName}

		rm \${libFile}*
		rm \${intronsFile}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
