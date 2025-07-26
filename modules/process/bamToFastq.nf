/*
 * BAM to FASTQ
 * 	- Shuffle the BAM + Convert to FASTQ (add BAM tags to FASTQ header)
 */

process BAM_TO_FASTQ {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'minMem'

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
		tuple val(sample), path("${sample}_${type}_RG.txt"), emit: RG
		path '*.{err,log,out,run,sh}'
	
	script:
		def fastqs    = params.singleEnd ? "-0 ${sample}_${type}.fastq.gz" : 
			"-1 ${sample}_${type}_R1.fastq.gz -2 ${sample}_${type}_R2.fastq.gz -s tmp_singletons.fastq -0 tmp_other.fastq"
		def tags      = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.bamTags_IonTorrent : params.bamTags
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		samtools collate ${bam} -O | samtools fastq ${fastqs} -T '${tags}' ${args}
		rm -f tmp_*

		samtools view ${bam} -H | grep "@RG" > ${sample}_${type}_RG.txt


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
