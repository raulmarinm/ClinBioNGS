/*
 * FASTQ processing - preprocessing for FASTQ files using fastp
 */

process FASTP {

	tag "${sample}"

	label 'fastp'
	label 'lowCpu'
	label 'minMem'

	publishDir(
		enabled: params.fastp_publish,
		path: "${params.dataRunDir}/${sample}/FASTQ",
		mode: params.fastp_publishMode,
		pattern: "${sample}_${type}*.fastq*"
	)
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}",
		mode: 'copy',
		pattern: "${sample}_${type}_fastp.{html,json}",
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
		path "${sample}_${type}_fastp.{html,json}", emit: files
		path '*.{err,log,out,run,sh}'
	
	script:
		def inFastq    = params.singleEnd ? "--in1 \$(ls data/*.fastq*)" : 
			"--in1 \$(ls data/${sample}_${type}*1*.fastq*) --in2 \$(ls data/${sample}_${type}*2*.fastq*)"
		def isUmi      = type.toUpperCase() == "RNA" ? params.umiRna && params.fastp_umi: params.umiDna && params.fastp_umi
		def umi_loc    = params.fastp_umi_loc
		def umi_len    = params.fastp_umi_len
		def umi_prefix = params.fastp_umi_prefix
		def umi_skip   = params.fastp_umi_skip
		def umi_delim  = params.fastp_umi_delim
		def umiArgs    = isUmi ? "-U --umi_loc=${umi_loc} --umi_len=${umi_len} --umi_prefix=${umi_prefix} --umi_skip=${umi_skip} --umi_delim=${umi_delim}" : ""
		def suffix     = params.fastp_suffix ?: ''
		def out1       = params.singleEnd ? "${sample}_${type}${suffix}.fastq.gz" : "${sample}_${type}${suffix}.R1.fastq.gz"
		def out2       = params.singleEnd ? "" : "--out2 ${sample}_${type}${suffix}.R2.fastq.gz" 
		def args       = task.ext.args ?: ''
		def logPrefix  = "${step}_${type}_${sample}"

		"""
		fastp --dont_eval_duplication --thread ${task.cpus} \\
			${inFastq} \\
			--out1 ${out1} \\
			${out2} \\
			--json ${sample}_${type}_fastp.json \\
			--html ${sample}_${type}_fastp.html \\
			${umiArgs} \\
			${args}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
