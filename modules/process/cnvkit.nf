/*
 * CNA analysis - CNVkit: Genome-wide copy number from high-throughput sequencing
 */

process CNVKIT {

	tag "${sample}"

	label 'cnvkit'
	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/CNVkit_Reports",
		mode: 'copy',
		pattern: "${sample}_DNA*",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai)
		path baseline
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA.cnr"), emit: bins
		tuple val(sample), path("${sample}_DNA.cns"), path("${sample}_DNA_sex.txt"), emit: segments
		tuple val(sample), path("${sample}_DNA_sex.txt"), emit: sex
		path '*.{err,log,out,run,sh}'
	
	script:
		def ampliconArgs    = params.amplicon ? '--no-edge' : ''
		def smoothCbs       = params.cnvkit_smoothCbs ? '--smooth-cbs' : ''
		def dropLowCoverage = params.cnvkit_dropLowCoverage ? '--drop-low-coverage' : ''
		def argsCoverage    = task.ext.argsCoverage ?: ''
		def argsFix         = task.ext.argsFix ?: ''
		def argsSex         = task.ext.argsSex ?: ''
		def argsSegment     = task.ext.argsSegment ?: ''
		def logPrefix       = "${step}_DNA_${sample}"

		"""
		tail -n+2 ${baseline} | awk 'BEGIN {FS="\\t";OFS="\\t"}{if (\$4 != "Antitarget") { print \$1,\$2,\$3,\$4}}' > target.bed
		tail -n+2 ${baseline} | awk 'BEGIN {FS="\\t";OFS="\\t"}{if (\$4 == "Antitarget") { print \$1,\$2,\$3,\$4}}' > antitarget.bed

		cnvkit.py coverage \\
			${bam} \\
			target.bed \\
			--processes ${task.cpus} \\
			-o ${sample}_DNA_target.cnn \\
			--min-mapq ${params.cnvkit_minMapq} \\
			${argsCoverage}
		cnvkit.py coverage \\
			${bam} \\
			antitarget.bed \\
			--processes ${task.cpus} \\
			-o ${sample}_DNA_antitarget.cnn \\
			--min-mapq ${params.cnvkit_minMapq} \\
			${argsCoverage}

		cnvkit.py fix \\
			${sample}_DNA_target.cnn \\
			${sample}_DNA_antitarget.cnn \\
			${baseline} \\
			--sample-id ${sample}_DNA \\
			-o ${sample}_DNA.cnr \\
			${ampliconArgs} \\
			${argsFix}

		cnvkit.py sex \\
			${sample}_DNA.cnr \\
			-o ${sample}_DNA_sex.txt \\
			${argsSex}
		
		cnvkit.py segment \\
			${sample}_DNA.cnr \\
			--processes ${task.cpus} \\
			-o ${sample}_DNA.cns \\
			--method ${params.cnvkit_segmentMethod} \\
			--threshold ${params.cnvkit_segmentThreshold} \\
			${dropLowCoverage} \\
			${smoothCbs} \\
			${argsSegment}

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""


}
