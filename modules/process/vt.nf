/*
 * Small variant calling - VCF normalization with Vt
 *	- Decompose multiallelic variants and biallelic block substitutions
 *	- Indel normalization and remove duplicates
 */

process VT {

	tag "${sample}"

	label 'vt'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.dataRunDir}/${sample}/VCF",
		mode: 'copy',
		pattern: "${sample}_DNA*.vcf"
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
		tuple val(sample), path("${sample}_DNA*.vcf"), emit: vcf
		path '*.{err,log,out,run,sh}'
	
	script:
		def argsDecompose = task.ext.argsDecompose ?: ''
		def argsView = task.ext.argsView ?: ''
		def argsDecomposeBlocksub = task.ext.argsDecomposeBlocksub ?: ''
		def argsNormalize = task.ext.argsNormalize ?: ''
		def argsUniq = task.ext.argsUniq ?: ''
		def decompose_blocksub = params.decomposeMnv ? "vt decompose_blocksub -a tmp_1_${sample}_DNA.vcf ${argsDecomposeBlocksub}" : ''
		def normalize = params.decomposeMnv ? "${decompose_blocksub} | vt normalize - -n -r \${fasta} -o tmp_2_${sample}_DNA.vcf ${argsNormalize}" :
			"vt normalize -n -r \${fasta} -o tmp_2_${sample}_DNA.vcf tmp_1_${sample}_DNA.vcf ${argsNormalize}"
		def suffix = params.vcfProcessing_suffix ?: ''
		def logPrefix = "${step}_DNA_${sample}"

		"""
		fasta=\$(ls genome/*genome.fa)

		vt decompose -s ${vcf} ${argsDecompose} | vt view - -f "ALT!='*'" -o tmp_1_${sample}_DNA.vcf ${argsView}

		${normalize}		
		
		vt uniq -o ${sample}_DNA_SmallVariant_${caller}${suffix}.vcf tmp_2_${sample}_DNA.vcf ${argsUniq}
		
		rm tmp_*
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
