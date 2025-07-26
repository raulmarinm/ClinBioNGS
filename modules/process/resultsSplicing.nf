/*
 * Splicing analysis - Collect metrics and report the results
 */

process RESULTS_SPLICING {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_RNA_Splicing*.xlsx"
	)
	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_RNA_Splicing.vcf"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_RNA_Splicing*.{png,txt}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		path targetBed
		path targetGenes
		path whitelistGenes
		path whitelistSplicing
		path maneGtf
		path maneGene
		path maneExon
		path maneIntron
		path cancerdrivers
		path civicClinical
		path chain
		path('genome/*')
		path vcfHeader
		path samplesFile
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_RNA_Splicing_results.txt"), emit: results
		tuple val(sample), path("${sample}_RNA_Splicing*.{png,txt}"), emit: files
		path "${sample}_RNA_Splicing.vcf"
		path "${params.runName}_${sample}_RNA_Splicing*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastaName = "${params.ctatLibDir}/${params.ctatLibName}/ref_genome.fa"
		def logPrefix  = "${step}_RNA_${sample}"

		"""
		fasta=\$(ls genome/*genome.fa)

		resultsSplicing.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=RNA \\
			--samplesFile=${samplesFile} \\
			--targetBed=${targetBed} \\
			--targetGenes=${targetGenes} \\
			--whitelistGenes=${whitelistGenes} \\
			--whitelistSplicing=${whitelistSplicing} \\
			--maneGtf=${maneGtf} \\
			--maneGene=${maneGene} \\
			--maneExon=${maneExon} \\
			--maneIntron=${maneIntron} \\
			--cancerdrivers=${cancerdrivers} \\
			--civicClinical=${civicClinical} \\
			--chain=${chain} \\
			--fasta=\${fasta} \\
			--fastaName=${fastaName} \\
			--vcfHeader=${vcfHeader} \\
			--callingMinAD=${params.splicingCallingMinAD} \\
			--flagMinAD=${params.splicingFlagMinAD} \\
			--flagMinDP=${params.splicingFlagMinDP} \\
			--flagMinVAF=${params.splicingFlagMinVAF}


		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
