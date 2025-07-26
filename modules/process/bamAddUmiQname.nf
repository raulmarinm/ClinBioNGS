/*
 * FASTQ generation - BAM add UMI to the Qname
 */

process BAM_ADD_UMI_QNAME {

	tag "${sample}"

	label 'samtools'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		enabled: params.bamAddUmiQname_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.bamAddUmiQname_publishMode,
		pattern: "${sample}_${type}*.bam"
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
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def umi1      = params.bamAddUmiQname_flagUmi1
		def umi2      = params.bamAddUmiQname_flagUmi2
		def umiDelim  = params.bamAddUmiQname_umiDelim
		def umiSep    = params.bamAddUmiQname_umiSep
		def suffix    = params.bamAddUmiQname_suffix ?: ''
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls data/*.bam)

		samtools view -h \${bam} \\
			| awk 'BEGIN {FS="\\t";OFS="\\t"} { if (/^@/) { print } else { for (i=11; i<=NF; ++i) { if (\$i ~ "^${umi1}:Z:|^${umi2}:Z:") { split(\$i, tag, ":Z:"); umi[tag[1]] = tag[2] } }; \$1=\$1"${umiDelim}"umi["${umi1}"]"${umiSep}"umi["${umi2}"]; print } }' \\
			| samtools view ${args} -bh -o ${sample}_${type}${suffix}.bam -

		

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""


}
