/*
 * Remove specific tags from the BAM file
 */

process BAM_REMOVE_TAGS {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('data/*'), path('data/*')
		val tags
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), path("${sample}_${type}*.bai"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls data/*.bam)
		filename=\$(basename \${bam})

		samtools view -h --remove-tag '${tags}' -o \${filename} \${bam} ${args}

		samtools index \${filename}


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
