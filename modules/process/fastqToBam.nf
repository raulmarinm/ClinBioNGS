/*
 * Read alignment - FASTQ to BAM
 */

process FASTQ_TO_BAM {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.fastqToBam_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.fastqToBam_publishMode,
		pattern: "${sample}_${type}*.{bam,bai}"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(reads), path(headerInfo)
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), path("${sample}_${type}*.bai"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastq     = params.singleEnd ? "-0 ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
		def tags      = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.bamTags_IonTorrent : params.bamTags
		def suffix    = params.fastqToBam_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools import ${fastq} -o tmp.bam -T '${tags}' ${args}
		samtools view -H tmp.bam | grep "@HD" > tmp_header.txt
		cat ${headerInfo} >> tmp_header.txt

		cat tmp_header.txt <(samtools view tmp.bam) | samtools view -bh -o ${sample}_${type}${suffix}.bam
		samtools index -@ ${task.cpus} ${sample}_${type}${suffix}.bam

		rm tmp.bam
		rm tmp_header.txt
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
