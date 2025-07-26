/*
 * Prepare resources - Target files
 * 	- Discard any non-primary chromosome
 * 	- Process de gene column and keep only with the gene name (if required)
 * 	- Merge overlapping regions into a single interval and collapse the gene name
 */

process PREPARE_TARGET_BED {

	tag "${params.runName}"

	label 'bedtools'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		val type
		val step

	output:
		path "*.bed", emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def sep       = params.geneColumnSep
		def panel     = params.seqPanel
		def genome    = params.genomeVersionHg
		def fileName  = "${panel}_${type}_${genome}.bed"
		def logPrefix = "${step}"

		if( panel.toUpperCase() == 'TSO500' )
		// Remove "MSI_" and "LoH_" from gene column
		"""
		sort -V -k1,1 -k2,2 ${bed} | \\
			awk '/^chr[0-9XY]*\\t/ {print \$0}' | \\
			awk 'BEGIN {FS="\\t";OFS="\\t"}{gsub(/MSI_/,"",\$4); gsub(/LoH_/,"",\$4); split(\$4, x, "${sep}"); print \$1,\$2,\$3,x[1]}' | \\
			bedtools merge -i stdin -c 4 -o distinct > ${fileName}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
		else
		// By default, the 4th column corresponds to the gene name
		"""
		sort -V -k1,1 -k2,2 ${bed} | \\
			awk '/^chr[0-9XY]*\\t/ {print \$0}' | \\
			awk 'BEGIN {FS="\\t";OFS="\\t"}{split(\$4, x, "${sep}"); print \$1,\$2,\$3,x[1]}' | \\
			bedtools merge -i stdin -c 4 -o distinct > ${fileName}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}


process PREPARE_TARGET_GENES {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path mane
		path cytoband
		path target
		val type
		val step

	output:
		path "*.txt", emit: genes
		path '*.{err,log,out,run,sh}'

	script:
		def panel     = params.seqPanel
		def logPrefix = "${step}"
		
		"""
		prepareTargetGenes.R --panel=${panel} --mane=${mane} --cytoband=${cytoband} --target=${target} --type=${type}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process PREPARE_TARGET_CHROM {

	tag "${params.runName}"

	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path bed
		val type
		val step

	output:
		path "*.txt", emit: chrom

	script:
		def logPrefix = "${step}"
		
		"""
		cut -f1 ${bed} | sort -V -k1,1 | uniq > ${params.seqPanel}_${type}_chromosomes.txt

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
