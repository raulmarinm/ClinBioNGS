/*
 * Small variant calling - TVC
 */

process TVC {

	tag "${sample}_${chr}"

	label 'tvc'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/TVC_Reports",
		mode: 'copy',
		pattern: "tvc_${chr}_metrics.json",
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
		path json
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.vcf"), emit: vcf
		path "tvc_${chr}_metrics.json"
		path '*.{err,log,out,run,sh}'
	
	script:
		def targetFileExt = target.getExtension()
		def suffix        = params.vcfCalling_suffix ?: '' 
		def args          = task.ext.args ?: ''
		def args2         = task.ext.args2 ?: ''
		def logPrefix     = "${step}_DNA_${sample}_${chr}"

		"""
		fasta=\$(ls genome/*genome.fa)

		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^@/ {print; next}; {if (\$1 == "${chr}") { print }}' ${target} > tmp_${chr}.${targetFileExt}

		tvc \\
			--reference \${fasta} \\
			--input-bam ${bam} \\
			--target-file tmp_${chr}.${targetFileExt} \\
			--num-threads ${task.cpus} \\
			--parameters-file ${json} \\
			${args}
		
		tvcutils unify_vcf \\
			--reference-fasta \${fasta} \\
			--filter-by-target off \\
			--novel-tvc-vcf small_variants.vcf \\
			--novel-assembly-vcf indel_assembly.vcf \\
			--tvc-metrics tvc_metrics.json \\
			--output-vcf ${sample}_DNA_SmallVariant_tvc_${chr}${suffix}.vcf \\
			${args2}
		
		mv tvc_metrics.json tvc_${chr}_metrics.json
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
