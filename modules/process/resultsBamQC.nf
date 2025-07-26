/*
 * BAM quality control - Collect metrics and report the results
 */

// Results from the realignment step
process RESULTS_BAM_QC {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'medMem'

	publishDir(
		path: "${params.resultsDir}/${sample}",
		mode: 'copy',
		pattern: "${params.runName}_${sample}_${type}*.xlsx",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/BamQC_CoverageByGene",
		mode: 'copy',
		pattern: "${sample}_${type}_CoverageByGene*.png",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/BamQC_CoverageByChr",
		mode: 'copy',
		pattern: "${sample}_${type}_CoverageByChr*.png",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_{BamQC.txt,GlobalCoverage.png}",
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_CoverageBy*.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		path samplesFile
		path targetBed
		path targetGenes
		path whitelistGenes
		path maneGtf
		path maneGene
		path maneExon
		path maneCoding
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}_BamQC.txt"), emit: results
		tuple val(sample), path("${sample}_${type}*.{png,txt}"), emit: files
		path "${params.runName}_${sample}_${type}*.xlsx"
		path '*.{err,log,out,run,sh}'
	
	script:
		def pctXCoverage            = type.toUpperCase() == "RNA" ? params.pctXCoverageRna: params.pctXCoverageDna
		def lowCoverage             = type.toUpperCase() == "RNA" ? params.lowCoverageRna: params.lowCoverageDna
		def plotMarkChrCoverage     = type.toUpperCase() == "RNA" ? params.bamqc_plotMarkChrCoverageRna: params.bamqc_plotMarkChrCoverageDna
		def plotGlobalLabelAllGenes = type.toUpperCase() == "RNA" ? params.bamqc_plotGlobalLabelAllGenesRna: params.bamqc_plotGlobalLabelAllGenesDna
		def plotByGeneAll           = type.toUpperCase() == "RNA" ? params.bamqc_plotByGeneAllRna: params.bamqc_plotByGeneAllDna
		def logPrefix               = "${step}_${type}_${sample}"

		"""
		resultsBamQc.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--type=${type} \\
			--singleEnd=${params.singleEnd} \\
			--samplesFile=${samplesFile} \\
			--targetBed=${targetBed} \\
			--targetGenes=${targetGenes} \\
			--whitelistGenes=${whitelistGenes} \\
			--maneGtf=${maneGtf} \\
			--maneGene=${maneGene} \\
			--maneExon=${maneExon} \\
			--maneCoding=${maneCoding} \\
			--pctXCoverage=${pctXCoverage} \\
			--lowCoverage=${lowCoverage} \\
			--plotMarkChrCoverage=${plotMarkChrCoverage} \\
			--plotGlobalLabelAllGenes=${plotGlobalLabelAllGenes} \\
			--plotByGeneAll=${plotByGeneAll}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}


// Results from the alignment step
process RESULTS_BAM_QC_RAW {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_RawBamQC.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}_RawBamQC.txt"), emit: results
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix = "${step}_${type}_${sample}"

		"""
		resultsBamQcRaw.R --run=${params.runName} --sample=${sample} --type=${type} --singleEnd=${params.singleEnd}
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
