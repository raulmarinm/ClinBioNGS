/*
 * Read alignment - TMAP
 */

process TMAP_ALIGNMENT {

	tag "${sample}"

	label 'tmap'
	label 'highCpu'
	label 'highMem'

	publishDir(
		enabled: params.tmapAlignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.tmapAlignment_publishMode,
		pattern: "${sample}_DNA*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/DNA/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('data/*'), path('data/*')
		path('index/*')
		path(bed)
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_DNA*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def paired    = params.singleEnd ? '' : '--pairing 2'
		def suffix    = params.bamAligned_suffix ?: ''
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		bam=\$(ls data/*.bam)
		tmap mapall \\
			-f index/genome.fa \\
			--reads-format bam \\
			-r \${bam} \\
			--bed-file ${bed} \\
			--context \\
			--num-threads ${task.cpus} \\
			-v \\
			-Y \\
			--rand-read-name \\
			--prefix-exclude 5 \\
			${paired} \\
			${args} \\
			-o 2 stage1 map4 \\
			| samtools sort - -@ ${task.cpus} -f ${sample}_DNA${suffix}.bam ${args2}
		
		samtools index ${sample}_DNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process TMAP_REALIGNMENT {

	tag "${sample}"

	label 'tmap'
	label 'highCpu'
	label 'highMem'

	publishDir(
		enabled: params.tmapRealignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.tmapRealignment_publishMode,
		pattern: "${sample}_DNA*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/DNA/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('data/*'), path('data/*')
		path('index/*')
		path(bed)
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_DNA*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def paired    = params.singleEnd ? '' : '--pairing 2'
		def suffix    = params.bamRealigned_suffix ?: ''
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		bam=\$(ls data/*.bam)
		tmap mapall \\
			-f index/genome.fa \\
			--reads-format bam \\
			-r \${bam} \\
			--bed-file ${bed} \\
			--context \\
			--num-threads ${task.cpus} \\
			-v \\
			-Y \\
			--rand-read-name \\
			--prefix-exclude 5 \\
			${paired} \\
			${args} \\
			-o 2 stage1 map4 \\
			| samtools sort - -@ ${task.cpus} -f ${sample}_DNA${suffix}.bam ${args2}
		
		samtools index ${sample}_DNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
