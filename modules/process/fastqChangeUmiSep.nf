/*
 * FASTQ processing - Change UMI separator from the FASTQ header
 */

process FASTQ_CHANGE_UMI_SEP {

	tag "${sample}"

	label 'bioawk'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.fastqChangeUmiSep_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.fastqChangeUmiSep_publishMode,
		pattern: "${sample}_${type}*.fastq*"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('data/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.fastq*"), emit: fastq
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix    = params.fastqChangeUmiSep_suffix ?: ''
		def from      = params.fastqChangeUmiSep_from
		def to        = params.fastqChangeUmiSep_to
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		fastq=\$(ls data/*.fastq*)
		filename=\$(basename \${fastq} | sed -e 's/.fastq/${suffix}.fastq/')
		
		bioawk -c fastx '{gsub(/\\${from}/,"${to}",\$name); print "@"\$name"\\n"\$seq"\\n+\\n"\$qual}' \${fastq} | gzip ${args} > \${filename}
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
