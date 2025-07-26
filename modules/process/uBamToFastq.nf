/*
 * FASTQ generation - uBAM to FASTQ
 * 	- Shuffle the BAM + Convert to FASTQ (add BAM tags to FASTQ header)
 */

process UBAM_TO_FASTQ {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.ubamToFastq_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.ubamToFastq_publishMode,
		pattern: "${sample}_${type}*.fastq*"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam)
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.fastq*"), emit: fastq
		tuple val(sample), path("${sample}_${type}_header_info.txt"), emit: headerInfo
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix    = params.ubamToFastq_suffix ?: ''
		def fastqs    = params.singleEnd ? "-0 ${sample}_${type}${suffix}.fastq.gz" : 
			"-1${sample}_${type}${suffix}.R1.fastq.gz -2 ${sample}_${type}${suffix}.R2.fastq.gz -s tmp_singletons.fastq -0 tmp_other.fastq"
		def tags      = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.bamTags_IonTorrent : params.bamTags
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools collate ${bam} -O | samtools fastq ${fastqs} -T '${tags}' ${args}
		samtools view -H ${bam} | grep -e "^@[CO,PG,RG]" > ${sample}_${type}_header_info.txt
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
