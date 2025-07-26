/*
 * Small variant annotation - Collect metrics and report the results
 */

process RESULTS_SMALL_VARIANT {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_DNA_SmallVariant*.xlsx"
	)
	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_DNA_SmallVariant_annot.vcf"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/SmallVariant_ByGene",
		mode: 'copy',
		pattern: "${sample}_DNA_SmallVariant*.png",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_DNA_SmallVariant_annot.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(consensus)
		path annotation
		path targetGenes
		path whitelistGenes
		path maneGtf
		path maneGene
		path maneExon
		path genieOncogenic
		path civicClinical
		path panelHotspots
		path vcfHeader
		path samplesFile
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA_SmallVariant_annot.txt"), emit: results
		tuple val(sample), path("${sample}_DNA_SmallVariant*.{png,txt}"), emit: files
		path "${sample}_DNA_SmallVariant_annot.vcf"
		path "${params.runName}_${sample}_DNA_SmallVariant*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastaName = "${params.genomeDir}/Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/WholeGenomeFasta/genome.fa"
		def logPrefix = "${step}_DNA_${sample}"

		"""
		resultsSmallVariant.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=DNA \\
			--samplesFile=${samplesFile} \\
			--consensus=${consensus} \\
			--annotation=${annotation} \\
			--targetGenes=${targetGenes} \\
			--whitelistGenes=${whitelistGenes} \\
			--maneGtf=${maneGtf} \\
			--maneGene=${maneGene} \\
			--maneExon=${maneExon} \\
			--genieOncogenic=${genieOncogenic} \\
			--civicClinical=${civicClinical} \\
			--panelHotspots=${panelHotspots} \\
			--vcfHeader=${vcfHeader} \\
			--fastaName=${fastaName} \\
			--minVAF=${params.variantMinVAF} \\
			--somaticMaxpVAF=${params.variantSomaticMaxpVAF} \\
			--somaticMaxVAF=${params.variantSomaticMaxVAF}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

