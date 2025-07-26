/*
 * FASTQ alignment - BWAMEM2
 */

process BWAMEM2_ALIGNMENT {

	tag "${sample}"

	label 'bwamem2'
	label 'highCpu'
	label 'highMem'

	publishDir(
		enabled: params.bwamem2Alignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bwamem2Alignment_publishMode,
		pattern: "${sample}_DNA*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/DNA/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads)
		path('index/*')
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_DNA*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def library   = params.runName
		def platform  = params.seqPlatform
		def suffix    = params.bamAligned_suffix ?: ''
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		bwa-mem2 mem \\
			index/genome.fa \\
			${reads} \\
			-t ${task.cpus} \\
			-M \\
			-R "@RG\\tID:${sample}_DNA\\tSM:${sample}_DNA\\tLB:${library}\\tPL:${platform}" \\
			-K ${params.bwamem2Alignment_inBases} \\
			${args} \\
			| samtools sort ${args2} -@ ${task.cpus} -o ${sample}_DNA${suffix}.bam -
		
		samtools index -@ ${task.cpus} ${sample}_DNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process BWAMEM2_REALIGNMENT {

	tag "${sample}"

	label 'bwamem2'
	label 'highCpu'
	label 'highMem'

	publishDir(
		enabled: params.bwamem2Realignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bwamem2Realignment_publishMode,
		pattern: "${sample}_DNA*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/DNA/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads)
		path('index/*')
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_DNA*.bai"), emit: bai
		path '*.{err,log,out,run,sh}'
	
	script:
		def library   = params.runName
		def platform  = params.seqPlatform
		def suffix    = params.bamRealigned_suffix ?: ''
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		bwa-mem2 mem \\
			index/genome.fa \\
			${reads} \\
			-t ${task.cpus} \\
			-M \\
			-R "@RG\\tID:${sample}_DNA\\tSM:${sample}_DNA\\tLB:${library}\\tPL:${platform}" \\
			-K ${params.bwamem2Realignment_inBases} \\
			${args} \\
			| samtools sort ${args2} -@ ${task.cpus} -o ${sample}_DNA${suffix}.bam -
		
		samtools index -@ ${task.cpus} ${sample}_DNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}



process BWAMEM2_IONTORRENT {

	tag "${sample}"

	label 'bwamem2'
	label 'medCpu'
	label 'highMem'

	publishDir(
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: "${publishMode}",
		pattern: "${sample}_DNA*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads), path(bamRG)
		path('index/*')
		val suffix
		val publishMode
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.bam"), path("${sample}_DNA*.bai"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def library   = params.runName
		def platform  = params.seqPlatform
		def args      = task.ext.args ?: ''
		def args2     = task.ext.args2 ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		cp ${bamRG} ${sample}_DNA.sam

		bwa-mem2 mem index/genome.fa ${reads} -t ${task.cpus} -M -C ${args} >> ${sample}_DNA.sam

		samtools sort ${args2} -@ ${task.cpus} -o ${sample}_DNA${suffix}.bam ${sample}_DNA.sam
		
		samtools index -@ ${task.cpus} ${sample}_DNA${suffix}.bam

		rm ${sample}_DNA.sam

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
