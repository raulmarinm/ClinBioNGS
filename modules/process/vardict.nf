/*
 * Small variant calling - VarDict
 */

process VARDICT {

	tag "${sample}_${chr}"

	label 'vardict'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai), val(chr)
		path('genome/*')
		path target
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.vcf"), emit: vcf
		path '*.{err,log,out,run,sh}'
	
	script:
		def targetFileExt = target.getExtension()
		def padding       = params.variantIntervalPadding ?: '0'
		def realign       = params.seqPlatform.toLowerCase() == 'iontorrent' ? '0' : '1'
		def minbq         = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.vardict_minbq_IonTorrent : params.vardict_minbq
		def suffix        = params.vcfCalling_suffix ?: ''
		def args          = task.ext.args ?: ''
		def args2         = task.ext.args2 ?: ''
		def logPrefix     = "${step}_DNA_${sample}_${chr}"

		"""
		fasta=\$(ls genome/*genome.fa)

		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^@/ {print; next}; {if (\$1 == "${chr}") { print }}' ${target} > tmp_${chr}.${targetFileExt}

		vardict-java \\
			-b ${bam} \\
			-G \${fasta} \\
			-N ${sample}_DNA \\
			-U \\
			-th ${task.cpus} \\
			-x ${padding} \\
			-c ${params.vardict_chrom} \\
			-S ${params.vardict_start} \\
			-E ${params.vardict_end} \\
			-g ${params.vardict_gene} \\
			-k ${params.vardict_localRealignment} \\
			-q ${minbq} \\
			-r ${params.vardict_minreads} \\
			-f ${params.vardict_minfreq} \\
			-X ${params.vardict_bpAfterIndel} \\
			${args} \\
			tmp_${chr}.${targetFileExt} | \\
			teststrandbias.R | \\
			var2vcf_valid.pl \\
				-N ${sample}_DNA \\
				-q ${minbq} \\
				-Q ${params.vardict_minmapq} \\
				-f ${params.vardict_minfreq} \\
				-d ${params.vardict_mindp} \\
				-v ${params.vardict_minad} \\
				${args2} > ${sample}_DNA_SmallVariant_vardict_${chr}${suffix}.vcf
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
