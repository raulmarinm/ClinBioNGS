/*
 * FASTQ processing - Transfer UMI from a separate FASTQ to the new FASTQ header
 */

process FASTQ_UMI_TRANSFER {

	tag "${sample}"

	label 'umitransfer'
	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.fastqUmiTransfer_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.fastqUmiTransfer_publishMode,
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
		def suffix    = params.fastqUmiTransfer_suffix ?: ''
		def umiFile   = params.fastqUmiTransfer_umiFile
		def in1File   = params.fastqUmiTransfer_in1File
		def in2File   = params.fastqUmiTransfer_in2File
		def delim     = params.fastqUmiTransfer_delim
		def in2       = params.singleEnd ? "${in1File}" : "${in2File}"
		def out1      = params.singleEnd ? "${sample}_${type}${suffix}.fastq.gz" : "${sample}_${type}${suffix}.R1.fastq.gz"
		def out2      = params.singleEnd ? '/dev/null' : "${sample}_${type}${suffix}.R2.fastq.gz"
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		fastq1=\$(ls data/*.fastq* | grep ${in1File})
		fastq2=\$(ls data/*.fastq* | grep ${in2})
		umi=\$(ls data/*.fastq* | grep ${umiFile})
		
		umi-transfer external \\
			--gzip --correct_numbers --delim '${delim}' \\
			--in \${fastq1} \\
			--in2 \${fastq2} \\
			--umi \${umi} \\
			--out ${out1} \\
			--out2 ${out2} \\
			${args}
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
