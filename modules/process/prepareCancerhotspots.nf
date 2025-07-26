/*
 * Prepare resources - Cancer Hotspots
 */

process PREPARE_HOTSPOTS_BED{

	tag "${params.runName}"

	label 'ucscLiftover'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path bed
		path chain
		val outName

	output:
		path "${outName}", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def fileName  = bed.getSimpleName()
		def logPrefix = "PREPARE_HOTSPOTS_BED"

		"""
		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^#/ {next}; {\$1="chr"\$1; print \$1,\$2,\$3}' ${bed} > ${fileName}_ucsc.bed
		liftOver ${fileName}_ucsc.bed ${chain} ${outName} ${fileName}.unlifted.bed

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}


process PREPARE_CANCERHOTSPOTS {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'lowMem'
	label 'outAnnotDir'

	input:
		path excel
		val outFile

	output:
		path "${outFile}", emit: file
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_CANCERHOTSPOTS"

		"""
		prepareCancerHotspots.R --excel=${excel} --outFile=${outFile}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

