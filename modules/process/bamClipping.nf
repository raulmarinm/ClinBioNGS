/*
 * Variant calling - Clip the ends of read alignments if they intersect with off-target regions
 */

process BAM_CLIPPING {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		enabled: params.bamclipping_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bamclipping_publishMode,
		pattern: "${sample}_${type}*.{bam,bai}"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: '*clipstats.txt',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), path("${sample}_${type}*.bai"), emit: bam
		path '*clipstats.txt'
		path '*.{err,log,out,run,sh}'
	
	script:
		def bothEnds  = params.bamclipping_bothEnds ? '--both-ends' : ''
		def suffix    = params.bamclipping_suffix 
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools ampliconclip \\
			${bothEnds} \\
			${args} \\
			-b ${bed} \\
			-f ${sample}_${type}_clipstats.txt \\
			-o tmp.bam \\
			${bam}
		samtools sort ${args2} \\
			-@ ${task.cpus} \\
			-o ${sample}_${type}${suffix}.bam \\
			tmp.bam
		
		samtools index -@ ${task.cpus} ${sample}_${type}${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
