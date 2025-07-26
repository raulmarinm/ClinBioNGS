/*
 * TMB analysis
 */

process TMB {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${sample}_DNA_TMB*.xlsx",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${step}",
		mode: 'copy',
		pattern: "${sample}_DNA_TMB*.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(variants), path(perBaseCovBed)
		path targetBed
		path maneCoding
		path problematicRegionsBed
		val step

	output:
		tuple val(sample), path("${sample}_DNA_TMB_results.txt"), emit: results
		tuple val(sample), path("${sample}_DNA_TMB_trace.txt"), emit: trace
		path '*.{err,log,out,run,sh,xlsx}'
	
	script:
		def logPrefix  = "${step}_DNA_${sample}"

		"""
		tmb.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=DNA \\
			--variants=${variants} \\
			--perBaseCovBed=${perBaseCovBed} \\
			--targetBed=${targetBed} \\
			--maneCoding=${maneCoding} \\
			--problematicRegionsBed=${problematicRegionsBed} \\
			--minAD=${params.tmbMinAD} \\
			--minDP=${params.tmbMinDP} --minVAF=${params.tmbMinVAF} \\
			--somaticMaxpVAF=${params.tmbSomaticMaxpVAF} \\
			--somaticMaxVAF=${params.tmbSomaticMaxVAF}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

