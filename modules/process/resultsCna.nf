/*
 * CNA analysis - Collect metrics and report the results
 */

process RESULTS_CNA {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'medMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_DNA_CNA*.xlsx"
	)
	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_DNA_CNA.vcf"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/CNA_ByGene",
		mode: 'copy',
		pattern: "${sample}_DNA_CNA_ByGene*.png",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/CNA_ByChr",
		mode: 'copy',
		pattern: "${sample}_DNA_CNA_ByChr*.png",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_DNA*{CNA.png,.txt}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		path targetBed
		path armCoordinates
		path('genome/*')
		path targetGenes
		path whitelistGenes
		path maneGtf
		path maneGene
		path maneExon
		path cancerdrivers
		path genie
		path civicClinical
		path vcfHeader
		path samplesFile
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA_CNA_gene.txt"), emit: results
		tuple val(sample), path("${sample}_DNA_CNA*.{png,txt}"), emit: files
		path "${sample}_DNA_CNA.vcf"
		path "${params.runName}_${sample}_DNA_CNA*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastaName = "${params.genomeDir}/Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/WholeGenomeFasta/genome.fa"
		def logPrefix = "${step}_DNA_${sample}"

		"""
		fasta=\$(ls genome/*genome.fa)

		resultsCna.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=DNA \\
			--samplesFile=${samplesFile} \\
			--targetBed=${targetBed} \\
			--fasta=\${fasta} \\
			--armCoordinates=${armCoordinates} \\
			--fastaName=${fastaName} \\
			--targetGenes=${targetGenes} \\
			--whitelistGenes=${whitelistGenes} \\
			--maneGtf=${maneGtf} \\
			--maneGene=${maneGene} \\
			--maneExon=${maneExon} \\
			--cancerdrivers=${cancerdrivers} \\
			--genie=${genie} \\
			--civicClinical=${civicClinical} \\
			--vcfHeader=${vcfHeader} \\
			--dropLowCoverage=${params.cnvkit_dropLowCoverage} \\
			--callingThresholds=${params.cnaCallingThresholds} \\
			--minBins=${params.cnaMinTargetBins} \\
			--highAmpCN=${params.cnaHighAmpCN} \\
			--highDelCN=${params.cnaHighDelCN} \\
			--plotGlobalLabelAllGenes=${params.cna_plotGlobalLabelAllGenes} \\
			--plotByGeneAll=${params.cna_plotByGeneAll}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
