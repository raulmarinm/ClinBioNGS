/*
 * Prepare resources - genome files
 */

process DOWNLOAD_GENOME {

	tag "${params.runName}"

	label 'awscli'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.genomeDir}/${dir}",
		mode: params.genomeFastaDir_publishMode,
		pattern: '*.{dict,fa,fai,xml}'
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		val link
		val dir
		val step

	output:
		path "${dir}", emit: fastaDir
		path '*.{dict,fa,fai,xml}', emit: fastaFiles
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "${step}"

		"""
		aws s3 --no-sign-request sync ${link} ./${dir}

		for file in \$(ls ${dir}); do ln -s ${dir}/\${file} \${file}; done

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process INDEX_GENOME_BWAMEM2 {

	tag "${params.runName}"

	label 'bwamem2'
	label 'minCpu'
	label 'extraMem'

	publishDir(
		path: "${indexDir}",
		mode: params.dnaGenomeIndexDir_publishMode,
		pattern: 'genome.fa*'
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path('data/*')
		val indexDir

	output:
		path 'genome.fa*', emit: indexFiles
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "INDEX_GENOME_BWAMEM2"

		"""
		ln -s data/genome.fa genome.fa
		bwa-mem2 index genome.fa

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process INDEX_GENOME_TMAP {

	tag "${params.runName}"

	label 'tmap'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${indexDir}",
		mode: params.dnaGenomeIndexDir_publishMode,
		pattern: 'genome.fa*'
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path('data/*')
		val indexDir

	output:
		path 'genome.fa*', emit: indexFiles
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "INDEX_GENOME_TMAP"

		"""
		ln -s data/genome.fa genome.fa
		tmap index -f genome.fa

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

