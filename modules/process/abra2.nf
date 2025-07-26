/*
 * Local realignment - ABRA2
 */

process ABRA2 {

	tag "${sample}"

	label 'abra2'
	label 'highCpu'
	label 'lowMem'

	publishDir(
		enabled: params.abra2_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.abra2_publishMode,
		pattern: "${sample}_${type}*.bam"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path('input/*'), path('input/*')
		path('genome/*')
		path bed
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def single    = params.singleEnd ? '--single' : ''
		def amplicon  = params.amplicon ? '--sa --cons' : ''
		def suffix    = params.bamRealigned_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls input/*.bam)
		fasta=\$(ls genome/*genome.fa)
		abra2 --in \${bam} --ref \${fasta} --targets ${bed} --threads ${task.cpus} \\
			--out ${sample}_${type}${suffix}.bam ${single} ${amplicon} ${args}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
