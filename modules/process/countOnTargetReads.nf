/*
 * BAM quality control - Count on-target reads
 */

process COUNT_ONTARGET_READS {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/BamQC_Reports",
		mode: 'copy',
		pattern: "${sample}_${type}_OnTargetReads.txt",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam)
		path bed
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}_OnTargetReads.txt"), emit: metrics
		path '*.{err,log,out,run,sh}'
	
	script:
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		target_reads=\$(samtools view -L ${bed} -c -F 0x900 ${args} ${bam})

		echo -e "SAMPLE\\tONTARGET_READS" > ${sample}_${type}_OnTargetReads.txt
		echo -e "${sample}_${type}\\t\${target_reads}" >> ${sample}_${type}_OnTargetReads.txt

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
