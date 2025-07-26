/*
 * Collect results and create a report for each sample
 */

process SAMPLE_REPORT {

	tag "${sample}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.resultsDir}",
		mode: 'copy',
		pattern: '*.html',
	)
	publishDir(
		path: { "${params.logsDir}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('./*')
		path dnaFile
		path rnaFile
		path panelHotspots
		path civicEvidence
		path libraryDb
		path rmdFiles
		val step

	output:
		path '*.html', emit: html
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix  = "${step}_${sample}"

		"""
		sampleReport.R \\
			--run=${params.runName} \\
			--sample=${sample} \\
			--panel=${params.seqPanel} \\
			--platform=${params.seqPlatform} \\
			--dnaSampleFile=${dnaFile} \\
			--rnaSampleFile=${rnaFile} \\
			--panelHotspots=${panelHotspots} \\
			--civicEvidence=${civicEvidence} \\
			--libraryDb=${libraryDb} \\
			--variantMinCallers=${params.variantMinCallers} \\
			--variantMinAD=${params.variantMinAD} \\
			--variantMinDP=${params.variantMinDP} \\
			--variantMinVAF=${params.variantMinVAF} \\
			--variantSomaticMaxpVAF=${params.variantSomaticMaxpVAF} \\
			--variantSomaticMaxVAF=${params.variantSomaticMaxVAF} \\
			--cnaMinTargetBins=${params.cnaMinTargetBins} \\
			--cnaHighAmpCN=${params.cnaHighAmpCN} \\
			--cnaHighDelCN=${params.cnaHighDelCN} \\
			--cnaPlotByGeneAll=${params.cna_plotByGeneAll} \\
			--fusionCallingMinAD=${params.fusionCallingMinAD} \\
			--fusionCallingMinFFPM=${params.fusionCallingMinFFPM} \\
			--fusionFlagMinAD=${params.fusionFlagMinAD} \\
			--fusionFlagMinNonFused=${params.fusionFlagMinNonFused} \\
			--fusionFlagMinDP=${params.fusionFlagMinDP} \\
			--fusionFlagMinVAF=${params.fusionFlagMinVAF} \\
			--splicingCallingMinAD=${params.splicingCallingMinAD} \\
			--splicingFlagMinAD=${params.splicingFlagMinAD} \\
			--splicingFlagMinDP=${params.splicingFlagMinDP} \\
			--splicingFlagMinVAF=${params.splicingFlagMinVAF} \\
			--covByExonAllGenesDna=${params.qc_covByExonAllGenesDna} \\
			--covByExonAllGenesRna=${params.qc_covByExonAllGenesRna} \\
			--bamqcPlotByGeneAllDna=${params.bamqc_plotByGeneAllDna} \\
			--bamqcPlotByGeneAllRna=${params.bamqc_plotByGeneAllRna} \\
			--totalReadsDna=${params.qcRange_totalReadsDna} \\
			--alignedReadsDna=${params.qcRange_alignedReadsDna} \\
			--alignedReadsPctDna=${params.qcRange_alignedReadsPctDna} \\
			--ontargetReadsDna=${params.qcRange_ontargetReadsDna} \\
			--ontargetReadsPctDna=${params.qcRange_ontargetReadsPctDna} \\
			--hqAlignedReadsDna=${params.qcRange_hqAlignedReadsDna} \\
			--hqAlignedReadsPctDna=${params.qcRange_hqAlignedReadsPctDna} \\
			--medianReadLengthDna=${params.qcRange_medianReadLengthDna} \\
			--medianInsertSizeDna=${params.qcRange_medianInsertSizeDna} \\
			--uniqueReadsDna=${params.qcRange_uniqueReadsDna} \\
			--duplicateReadsPctDna=${params.qcRange_duplicateReadsPctDna} \\
			--medianCoverageDna=${params.qcRange_medianCoverageDna} \\
			--meanCoverageDna=${params.qcRange_meanCoverageDna} \\
			--pct04xMeanCoverageDna=${params.qcRange_pct04xMeanCoverageDna} \\
			--pct100xCoverageDna=${params.qcRange_pct100xCoverageDna} \\
			--pct1000xCoverageDna=${params.qcRange_pct1000xCoverageDna} \\
			--totalReadsRna=${params.qcRange_totalReadsRna} \\
			--alignedReadsRna=${params.qcRange_alignedReadsRna} \\
			--alignedReadsPctRna=${params.qcRange_alignedReadsPctRna} \\
			--ontargetReadsRna=${params.qcRange_ontargetReadsRna} \\
			--ontargetReadsPctRna=${params.qcRange_ontargetReadsPctRna} \\
			--hqAlignedReadsRna=${params.qcRange_hqAlignedReadsRna} \\
			--hqAlignedReadsPctRna=${params.qcRange_hqAlignedReadsPctRna} \\
			--medianReadLengthRna=${params.qcRange_medianReadLengthRna} \\
			--medianInsertSizeRna=${params.qcRange_medianInsertSizeRna} \\
			--uniqueReadsRna=${params.qcRange_uniqueReadsRna} \\
			--duplicateReadsPctRna=${params.qcRange_duplicateReadsPctRna} \\
			--medianCoverageRna=${params.qcRange_medianCoverageRna} \\
			--meanCoverageRna=${params.qcRange_meanCoverageRna} \\
			--pct04xMeanCoverageRna=${params.qcRange_pct04xMeanCoverageRna} \\
			--pct100xCoverageRna=${params.qcRange_pct100xCoverageRna} \\
			--pct1000xCoverageRna=${params.qcRange_pct1000xCoverageRna}
		 

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
