/*
 * Prepare resources - Manifest txt to BED file
 */

process MANIFEST_TO_BED {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path manifest
		val step

	output:
		path "*.bed", emit: bed
		path '*.{err,log,out,run,sh}'
		
	script:
		def fileName  = manifest.getSimpleName()
		def logPrefix = "${step}"

		"""
		sed -n '/\\[Regions\\]/,\$p' ${manifest} | \\
			tail -n+3 | cut -f2-5 | \\
			awk 'BEGIN {FS="\t";OFS="\t"}{print \$1,\$2-1,\$3,\$4}' | sort -V -k1,1 -k2,2  > ${fileName}.bed
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
