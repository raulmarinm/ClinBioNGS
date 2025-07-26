/*
 * FASTQ alignment - STAR
 */

process STAR_ALIGNMENT {

	tag "${sample}"

	label 'starFusion'
	label 'medCpu'
	label 'highMem'

	publishDir(
		enabled: params.starAlignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.starAlignment_publishMode,
		pattern: "${sample}_RNA*.{bam,bai}"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_RNA*.{junction,tab}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/RNA/${step}" },
		mode: 'copy',
		pattern: '*.command.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads)
		path(index)
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_RNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_RNA*.bai"), emit: bai
		tuple val(sample), path('*SJ.out.tab'), emit: SJ
		tuple val(sample), path('*.out.junction'), emit: junction
		tuple val(sample), path('*ReadsPerGene.out.tab')
		path '*.command.{err,log,out,run,sh}'
	
	script:
		def library   = params.runName
		def platform  = params.seqPlatform
		def suffix    = params.bamAligned_suffix ?: ''
		def mainArgs  = task.ext.mainArgs ?: ''
		def otherArgs = task.ext.otherArgs ?: ''
		def logPrefix = "${step}_RNA_${sample}"

		"""
		STAR \\
			--runThreadN ${task.cpus} \\
			--genomeDir ${index} \\
			--readFilesIn ${reads} \\
			--readFilesCommand zcat \\
			--outFileNamePrefix ${sample}_RNA${suffix}. \\
			--outSAMtype BAM SortedByCoordinate \\
			--outSAMattrRGline ID:${sample}_RNA SM:${sample}_RNA LB:${library} PL:${platform} \\
			${mainArgs} \\
			${otherArgs}

		mv ${sample}_RNA${suffix}.Aligned.sortedByCoord.out.bam ${sample}_RNA${suffix}.bam

		samtools index -@ ${task.cpus} ${sample}_RNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

process STAR_REALIGNMENT {

	tag "${sample}"

	label 'starFusion'
	label 'medCpu'
	label 'highMem'

	publishDir(
		enabled: params.starRealignment_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.starRealignment_publishMode,
		pattern: "${sample}_RNA*.{bam,bai}"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_RNA*.{junction,tab}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/RNA/${step}" },
		mode: 'copy',
		pattern: '*.command.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads)
		path(index)
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_RNA*.bam"), emit: bam
		tuple val(sample), path("${sample}_RNA*.bai"), emit: bai
		tuple val(sample), path('*SJ.out.tab'), emit: SJ
		tuple val(sample), path('*.out.junction'), emit: junction
		tuple val(sample), path('*ReadsPerGene.out.tab')
		path '*.command.{err,log,out,run,sh}'
	
	script:
		def library   = params.runName
		def platform  = params.seqPlatform
		def suffix    = params.bamRealigned_suffix ?: ''
		def mainArgs  = task.ext.mainArgs ?: ''
		def otherArgs = task.ext.otherArgs ?: ''
		def logPrefix = "${step}_RNA_${sample}"

		"""
		STAR \\
			--runThreadN ${task.cpus} \\
			--genomeDir ${index} \\
			--readFilesIn ${reads} \\
			--readFilesCommand zcat \\
			--outFileNamePrefix ${sample}_RNA${suffix}. \\
			--outSAMtype BAM SortedByCoordinate \\
			--outSAMattrRGline ID:${sample}_RNA SM:${sample}_RNA LB:${library} PL:${platform} \\
			${mainArgs} \\
			${otherArgs}

		mv ${sample}_RNA${suffix}.Aligned.sortedByCoord.out.bam ${sample}_RNA${suffix}.bam

		samtools index -@ ${task.cpus} ${sample}_RNA${suffix}.bam

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
