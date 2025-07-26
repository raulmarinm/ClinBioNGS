/*
 * VCF - Concat files
 */

process VCF_CONCAT {

	tag "${sample}"

	label 'bcftools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "*.vcf"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(vcf)
		path('genome/*')
		val subworkflow
		val step
		val caller

	output:
		tuple val(sample), path("*.vcf"), emit: vcf
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix     = params.vcfCalling_suffix ?: ''
		def args       = task.ext.args ?: ''
		def args2      = task.ext.args2 ?: ''
		def logPrefix  = "${step}_DNA_${sample}"

		"""
		fai=\$(ls genome/*.fai)

		for f in ${vcf}
		do
			bcftools reheader --threads ${task.cpus} --fai \${fai} -o tmp_\${f} \${f}
		done

		vcf_list=\$(ls tmp_*.vcf | sort -V -k1,1 -k2,2)

		bcftools concat --threads ${task.cpus} ${args} -o ${sample}_DNA_SmallVariant_${caller}${suffix}.vcf \${vcf_list}

		rm tmp_*.vcf
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
