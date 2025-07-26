/*
 * Prepare resources - CIViC
 */

process PREPARE_CIVIC_VCF{

	tag "${params.runName}"

	label 'gatk4'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path vcf
		path chain
		path('genome/*')
		val outName

	output:
		tuple path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf
		path '*.vcf'
		path '*.{err,log,out,run,sh}'

	script:
		def inName  = vcf.getSimpleName()
		def logPrefix = "PREPARE_CIVIC_VCF"

		"""
		fasta=\$(ls genome/*genome.fa)

		awk 'BEGIN {OFS = "\\t"}; /^#/ {print; next}; \$4 != \$5 {\$1="chr"\$1; print}' ${vcf} > tmp_1.vcf
		gatk LiftoverVcf -I tmp_1.vcf -O tmp_2.vcf -C ${chain} -R \${fasta} --REJECT ${inName}_rejected_liftover.vcf
		bgzip --stdout tmp_2.vcf > ${outName}
		tabix -p vcf ${outName}

		rm tmp_*


		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}


process PREPARE_CIVIC_EVIDENCE{

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path tumorNames
		path variant
		path molecular
		path evidence

	output:
		path '*_Oncogenic.txt', emit: oncogenic
		path '*_Clinical.txt', emit: clinical
		path '*.{err,log,out,run,sh,tsv}'

	script:
		def logPrefix = "PREPARE_CIVIC_EVIDENCE"

		"""
		prepareCivicEvidence.R \\
			--tumorNames=${tumorNames} \\
			--variantFile=${variant} \\
			--molecularFile=${molecular} \\
			--evidenceFile=${evidence} \\
			--date=${params.civicDate}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

