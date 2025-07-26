/*
 * Small variant calling - VCF processing with Bcftools
 *	- Add to the header contig names and their lengths
 *	- Keep PASS & TARGET variants
 */

process BCFTOOLS {

	tag "${sample}"

	label 'bcftools'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(vcf)
		path bed
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA.vcf"), emit: vcf
		path '*.{err,log,out,run,sh}'
	
	script:
		def args = task.ext.args ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		bcftools view -Oz --threads ${task.cpus} -o tmp_${sample}_DNA.vcf.gz ${vcf}

		bcftools index --threads ${task.cpus} tmp_${sample}_DNA.vcf.gz
		
		bcftools view --threads ${task.cpus} --regions-file ${bed} --apply-filters PASS -o ${sample}_DNA.vcf ${args} tmp_${sample}_DNA.vcf.gz

		rm tmp_*
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
