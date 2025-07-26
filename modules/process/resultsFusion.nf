/*
 * Fusion analysis - Collect metrics and report the results
 */

process RESULTS_FUSION {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_RNA_Fusion*.xlsx"
	)
	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_RNA_Fusion.vcf"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_RNA_Fusion*.{png,txt}",
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
		path whitelistFusions
		path maneGtf
		path maneGene
		path maneExon
		path maneIntron
		path cancerdrivers
		path mitelmanFusion
		path genieFusion
		path civicClinical
		path chain
		path('genome/*')
		path vcfHeader
		path samplesFile
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_RNA_Fusion_results.txt"), emit: results
		tuple val(sample), path("${sample}_RNA_Fusion*.{png,txt}"), emit: files
		path "${sample}_RNA_Fusion.vcf"
		path "${params.runName}_${sample}_RNA_Fusion*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastaName = "${params.ctatLibDir}/${params.ctatLibName}/ref_genome.fa"
		def logPrefix  = "${step}_RNA_${sample}"

		"""
		fasta=\$(ls genome/*genome.fa)

		resultsFusion.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=RNA \\
			--samplesFile=${samplesFile} \\
			--targetBed=${targetBed} \\
			--targetGenes=${targetGenes} \\
			--whitelistGenes=${whitelistGenes} \\
			--whitelistFusions=${whitelistFusions} \\
			--maneGtf=${maneGtf} \\
			--maneGene=${maneGene} \\
			--maneExon=${maneExon} \\
			--maneIntron=${maneIntron} \\
			--cancerdrivers=${cancerdrivers} \\
			--mitelmanFusion=${mitelmanFusion} \\
			--genieFusion=${genieFusion} \\
			--civicClinical=${civicClinical} \\
			--chain=${chain} \\
			--fasta=\${fasta} \\
			--fastaName=${fastaName} \\
			--vcfHeader=${vcfHeader} \\
			--callingMinJunctionReads=${params.fusionCallingMinJunctionReads} \\
			--callingMinSpanningReads=${params.fusionCallingMinSpanningReads} \\
			--callingMinAD=${params.fusionCallingMinAD} \\
			--callingMinFFPM=${params.fusionCallingMinFFPM} \\
			--flagMinAD=${params.fusionFlagMinAD} \\
			--flagMinADNonFused=${params.fusionFlagMinNonFused} \\
			--flagMinDP=${params.fusionFlagMinDP} \\
			--flagMinVAF=${params.fusionFlagMinVAF}

		rm -f Rplots*

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
