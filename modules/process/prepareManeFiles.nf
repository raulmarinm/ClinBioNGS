/*
 * Prepare resources - MANE files
 */

process PREPARE_MANE_FILES {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path gtf
		path summary
		val genes
		val exons
		val introns
		val coding

	output:
		path "${genes}*", emit: genes
		path "${exons}*", emit: exons
		path "${introns}*", emit: introns
		path "${coding}*", emit: coding
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_MANE_FILES"

		"""
		prepareManeFiles.R --gtf=${gtf} --summary=${summary} --genes=${genes} --exons=${exons} --introns=${introns} --coding=${coding}
		gzip --best ${genes}
		gzip --best ${exons}
		gzip --best ${introns}
		gzip --best ${coding}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
