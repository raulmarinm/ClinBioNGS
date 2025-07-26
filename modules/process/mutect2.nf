/*
 * Small variant calling - Mutect2
 */

// Mutect2
process MUTECT2_CALL {

	tag "${sample}_${chr}"

	label 'gatk4'
	label 'minCpu'
	label 'lowMem'
	
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/Mutect2_Reports",
		mode: 'copy',
		pattern: '*.{stats,tar.gz}',
	)
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
		tuple val(sample), val(chr), path("*.vcf"), path("*.vcf.stats"), emit: vcf
		tuple val(sample), path("*f1r2.tar.gz"), emit: f1r2
		path '*.{err,log,out,run,sh}'
	
	script:
		def targetFileExt = target.getExtension()
		def args          = task.ext.args ?: ''
		def logPrefix     = "${step}_DNA_${sample}_${chr}"

		"""
		fasta=\$(ls genome/*.fa)

		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^@/ {print; next}; {if (\$1 == "${chr}") { print }}' ${target} > tmp_${chr}.${targetFileExt}

		gatk Mutect2 \\
			-I ${bam} \\
			-O tmp_${sample}_DNA_${chr}.vcf \\
			-R \${fasta} \\
			-L tmp_${chr}.${targetFileExt} \\
			--native-pair-hmm-threads ${task.cpus} \\
			--f1r2-tar-gz ${sample}_DNA_mutect2_${chr}_f1r2.tar.gz --f1r2-max-depth ${params.mutect2_F1r2MaxDepth} \\
			--max-reads-per-alignment-start ${params.mutect2_maxReadsPerAlignmentStart} \\
			--max-mnp-distance ${params.mutect2_maxMnpDistance} \\
			${args}
		
		mv tmp_${sample}_DNA_${chr}.vcf.stats ${sample}_DNA_mutect2_${chr}.vcf.stats

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

// LearnReadOrientationModel
process MUTECT2_MODEL {

	tag "${sample}"

	label 'gatk4'
	label 'minCpu'
	label 'medMem'
	
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/Mutect2_Reports",
		mode: 'copy',
		pattern: '*read-orientation-model.tar.gz',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(f1r2)
		val subworkflow
		val step

	output:
		tuple val(sample), path("*read-orientation-model.tar.gz"), emit: model
		path '*.{err,log,out,run,sh}'
	
	script:
		def input_list = f1r2.collect{ "-I ${it}"}.join(' ')
		def args       = task.ext.args ?: ''
		def logPrefix  = "${step}_DNA_${sample}"

		"""
		gatk LearnReadOrientationModel \\
			${input_list} \\
			-O ${sample}_DNA_mutect2_read-orientation-model.tar.gz \\
			--max-depth ${params.mutect2_F1r2MaxDepth} \\
			${args}
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

// FilterMutectCalls
process MUTECT2_FILTER {

	tag "${sample}_${chr}"

	label 'gatk4'
	label 'minCpu'
	label 'minMem'
	
	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/Mutect2_Reports",
		mode: 'copy',
		pattern: '*.tsv',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), val(chr), path(vcf), path(stats), path(model)
		path('genome/*')
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.vcf"), emit: vcf
		path '*.{tsv,err,log,out,run,sh}'
	
	script:
		def minbq     = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.mutect2_minbq_IonTorrent : params.mutect2_minbq
		def suffix    = params.vcfCalling_suffix ?: '' 
		def args      = task.ext.args ?: ''
		def logPrefix = "${step}_DNA_${sample}_${chr}"

		"""
		fasta=\$(ls genome/*.fa)

		gatk FilterMutectCalls \\
			-V ${vcf} \\
			-O ${sample}_DNA_SmallVariant_mutect2_${chr}${suffix}.vcf \\
			-R \${fasta} \\
			--stats ${stats} \\
			-ob-priors ${model} \\
			--min-median-base-quality ${minbq} \\
			--min-median-mapping-quality ${params.mutect2_minmapq} \\
			--min-allele-fraction ${params.mutect2_minfreq} \\
			--max-events-in-region ${params.mutect2_maxEventsInRegion} \\
			--max-alt-allele-count ${params.mutect2_maxAltAlleleCount} \\
			${args}
		
		mv ${sample}_DNA_SmallVariant_mutect2_${chr}${suffix}.vcf.filteringStats.tsv ${sample}_DNA_mutect2_${chr}.vcf.filteringStats.tsv
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
