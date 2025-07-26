/*
 * FASTQ processing - Merge FASTQ files from the same lane
 */

process FASTQ_MERGE_LANES {

	tag "${sample}"

	label 'minCpu'
	label 'minMem'

	publishDir(
		enabled: params.fastqMergeLanes_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.fastqMergeLanes_publishMode,
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
		def suffix    = params.fastqMergeLanes_suffix ?: ''
		def isUmi     = type.toUpperCase() == "RNA" ? params.umiRna: params.umiDna
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		if( params.singleEnd && isUmi && params.fastqUmiTransfer)
		"""
		zcat data/*1.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R1.fastq.gz
		zcat data/*2.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R2.fastq.gz
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
		else if( params.singleEnd )
		"""
		zcat data/*.fastq* | gzip ${args} > ${sample}_${type}${suffix}.fastq.gz
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

		else if( isUmi && params.fastqUmiTransfer )
		"""
		zcat data/*1.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R1.fastq.gz
		zcat data/*2.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R2.fastq.gz
		zcat data/*3.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R3.fastq.gz
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
		else
		"""
		zcat data/*1.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R1.fastq.gz
		zcat data/*2.fastq* | gzip ${args} > ${sample}_${type}${suffix}.R2.fastq.gz
		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
