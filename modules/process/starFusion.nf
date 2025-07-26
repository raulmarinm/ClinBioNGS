/*
 * Fusion analysis - STAR-fusion
 */

process STAR_FUSION {

	tag "${sample}"

	label 'starFusion'
	label 'medCpu'
	label 'highMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/STAR-Fusion",
		mode: 'copy',
		pattern: '**.{tsv,wAnnot}*',
		saveAs: { fn -> "${sample}_RNA_${fn}" }
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.command.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai), path(junction)
		path(ctatLib)
		val subworkflow
		val step

	output:
		tuple val(sample), path('**.{tsv,wAnnot}*'), emit: results
		path '*.command.{err,log,out,run,sh}'
	
	script:
		def samtoolsFastqs    = params.singleEnd ? '-0 tmp.fastq.gz' : 
			'-1 tmp_R1.fastq.gz -2 tmp_R2.fastq.gz -s /dev/null -0 /dev/null'
		def fusionFastqs    = params.singleEnd ? '--left_fq tmp.fastq.gz' : 
			'--left_fq tmp_R1.fastq.gz --right_fq tmp_R2.fastq.gz'
		def args             = task.ext.args ?: ''
		def logPrefix        = "${step}_RNA_${sample}"

		"""
		samtools collate ${bam} -O | samtools fastq - ${samtoolsFastqs}

		STAR-Fusion --genome_lib_dir ${ctatLib} --chimeric_junction ${junction} ${fusionFastqs} \\
			--CPU ${task.cpus} --examine_coding_effect --FusionInspector validate --output_dir . ${args}
		
		rm -f *.fastq.gz

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
