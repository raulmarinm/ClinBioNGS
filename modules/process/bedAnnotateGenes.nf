/*
 * Prepare resources - Annotate genes to the BED file
 */

process BED_ANNOTATE_GENES {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		path genes
		val step

	output:
		path "*.bed", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def fileName  = bed.getBaseName()
		def logPrefix = "${step}"
		
		"""
		bedAnnotateGenes.R --bed=${bed} --genes=${genes} --outFile=${fileName}.annot.bed


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
