/*
 * Small variant calling - VCF consensus
 */

process VCF_CONSENSUS {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'medMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_DNA_SmallVariant*.xlsx"
	)
	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_DNA_SmallVariant_consensus.vcf"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_DNA*.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		val callers
		path chain
		path panelHotspotsBed
		path panelBlacklistBed
		path vcfHeader
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA_SmallVariant_consensus.vcf"), emit: vcf
		tuple val(sample), path("${sample}_DNA_SmallVariant_consensus.txt"), emit: results
		path "${sample}_DNA*.txt"
		path "${params.runName}_${sample}_DNA_SmallVariant*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def genome     = "${params.genomeDir}/Homo_sapiens/NCBI/${params.genomeVersionGrc}/Sequence/WholeGenomeFasta/genome.fa"
		def logPrefix  = "${step}_DNA_${sample}"

		"""
		vcfConsensus.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=DNA \\
			--chain=${chain} \\
			--panelHotspotsBed=${panelHotspotsBed} \\
			--panelBlacklistBed=${panelBlacklistBed} \\
			--genome=${genome} \\
			--vcfHeader=${vcfHeader} \\
			--callers=${callers} \\
			--minCallers=${params.variantMinCallers} \\
			--minAD=${params.variantMinAD} \\
			--minDP=${params.variantMinDP} \\
			--minVAF=${params.variantMinVAF} \\
			--minADHotspot=${params.variantMinAD_PanelHotspot} \\
			--minDPHotspot=${params.variantMinDP_PanelHotspot} \\
			--minVAFHotspot=${params.variantMinVAF_PanelHotspot}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
