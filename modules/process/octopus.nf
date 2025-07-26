/*
 * Small variant calling - Octopus
 */

process OCTOPUS {

	tag "${sample}_${chr}"

	label 'octopus'
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
		def targetFileExt    = target.getExtension()
		def callingModel     = params.octopus_callingModel ?: 'individual' 
		def errorModel       = params.octopus_errorModel ?: 'PCR' 
		def minbq            = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.octopus_minbq_IonTorrent : params.octopus_minbq
		def filterExpression = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.octopus_filterExpression_IonTorrent : params.octopus_filterExpression
		def suffix           = params.vcfCalling_suffix ?: '' 
		def args             = task.ext.args ?: ''
		def logPrefix        = "${step}_DNA_${sample}_${chr}"

		"""
		fasta=\$(ls genome/*genome.fa)

		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^@/ {print; next}; {if (\$1 == "${chr}") { print }}' ${target} > tmp_${chr}.${targetFileExt}

		octopus \\
			-R \${fasta} \\
			-I ${bam} \\
			-o tmp_${sample}_DNA_${chr}.vcf \\
			--regions-file tmp_${chr}.${targetFileExt} \\
			--threads ${task.cpus} \\
			--allow-octopus-duplicates \\
			--caller ${callingModel} \\
			--sequence-error-model ${errorModel} \\
			--max-reference-cache-memory ${params.octopus_maxRefCacheMem} \\
			--target-read-buffer-memory ${params.octopus_targetReadBufferMem} \\
			--downsample-above ${params.octopus_downsampleAbove} \\
			--downsample-target ${params.octopus_downsampleTarget} \\
			--min-pileup-base-quality ${minbq} \\
			--good-base-quality ${minbq} \\
			--min-mapping-quality ${params.octopus_minmapq} \\
			--filter-expression '${filterExpression}' \\
			--annotations ${params.octopus_annotations} \\
			${args}
		
		awk 'BEGIN {OFS="\\t"}; /^#/ {print; next}; {gsub(/:\\.,/,":0,",\$10); print}' tmp_${sample}_DNA_${chr}.vcf > ${sample}_DNA_SmallVariant_octopus_${chr}${suffix}.vcf
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
