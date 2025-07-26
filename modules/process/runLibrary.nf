/*
 * Collect run results and create a variant library
 */

process RUN_LIBRARY {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.resultsDir}",
		mode: 'copy',
		pattern: '*.db',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path './*'
		path dnaFile
		path rnaFile
		val subworkflow
		val step

	output:
		path '*.db', emit: db
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix  = "${step}"

		"""
		runLibrary.R \\
			--run=${params.runName} \\
			--dnaSampleFile=${dnaFile} \\
			--rnaSampleFile=${rnaFile} \\
			--variantMinCallers=${params.variantMinCallers} \\
			--variantMinAD=${params.variantMinAD} \\
			--variantMinDP=${params.variantMinDP} \\
			--variantMinVAF=${params.variantMinVAF} \\
			--variantSomaticMaxpVAF=${params.variantSomaticMaxpVAF} \\
			--variantSomaticMaxVAF=${params.variantSomaticMaxVAF} \\
			--cnaMinTargetBins=${params.cnaMinTargetBins} \\
			--cnaHighAmpCN=${params.cnaHighAmpCN} \\
			--cnaHighDelCN=${params.cnaHighDelCN} \\
			--cnaArmCallingThreshold=${params.ascets_callingThreshold} \\
			--fusionCallingMinAD=${params.fusionCallingMinAD} \\
			--fusionCallingMinFFPM=${params.fusionCallingMinFFPM} \\
			--fusionFlagMinAD=${params.fusionFlagMinAD} \\
			--fusionFlagMinNonFused=${params.fusionFlagMinNonFused} \\
			--fusionFlagMinDP=${params.fusionFlagMinDP} \\
			--fusionFlagMinVAF=${params.fusionFlagMinVAF} \\
			--splicingCallingMinAD=${params.splicingCallingMinAD} \\
			--splicingFlagMinAD=${params.splicingFlagMinAD} \\
			--splicingFlagMinDP=${params.splicingFlagMinDP} \\
			--splicingFlagMinVAF=${params.splicingFlagMinVAF}
		 

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

