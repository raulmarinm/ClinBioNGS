/*
 * PDX processing - Xengsort: filter mouse reads from FASTQs
 */

process XENGSORT_CLASSIFY {

	tag "${sample}"

	label 'xengsort'
	label 'medCpu'
	label 'lowMem'

	publishDir(
		enabled: params.xengsort_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.xengsort_publishMode,
		pattern: "${sample}_${type}*graft*.fastq.gz"
	)
	publishDir(
		path: { "${params.logsDir}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads)
		path index
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*graft*.fastq.gz"), emit: fastq
		path '*.{err,log,out,run,sh}'

	script:
		def fastqs    = params.singleEnd ? "--fastq <(zcat ${reads})" : "--fastq <(zcat ${reads[0]}) --pairs <(zcat ${reads[1]})"
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}"

		"""
		xengsort classify \\
			--index ${index} \\
			--threads ${task.cpus} \\
			--prefix ${sample}_${type} \\
			${args} \\
			${fastqs}
		
		for f in *.fq
		do
		out=\$(echo "\${f}" | sed -e 's/.1.fq/_R1.fastq/' -e 's/.2.fq/_R2.fastq/' -e 's/.fq/.fastq/')
		mv -- "\${f}" "\${out}"
		done
		gzip *.fastq

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}


process XENGSORT_INDEX {

	tag "${params.runName}"

	label 'xengsort'
	label 'highCpu'
	label 'highMem'

	publishDir(
		path: "${indexDir}",
		mode: params.xengsortIndex_publishMode,
		pattern: "${indexName}"
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path('mouse/*')
		path('human/*')
		val indexDir
		val indexName

	output:
		path "${indexName}", emit: index
		path '*.{err,log,out,run,sh}'

	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "XENGSORT_INDEX"

		"""
		ln -s mouse/genome.fa mouse.fa
		ln -s human/genome.fa human.fa

		xengsort index \\
			${indexName} \\
			-H mouse.fa \\
			-G human.fa \\
			--threads ${task.cpus} \\
			--nobjects ${params.xengsortIndex_nobjects} \\
			${args}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

