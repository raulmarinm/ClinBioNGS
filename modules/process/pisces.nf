/*
 * Small variant calling - Pisces
 */

process PISCES {

	tag "${sample}_${chr}"

	label 'pisces'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${sample}/${subworkflow}/Pisces_Reports",
		mode: 'copy',
		pattern: "PiscesLogs_${chr}",
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), path(bam), path(bai), val(chr)
		path genomeDir
		path target
		val subworkflow
		val step

	output:
		tuple val(sample), path("${sample}_DNA*.vcf"), emit: vcf
		path "PiscesLogs_${chr}"
		path '*.{err,log,out,run,sh}'
	
	script:
		def targetFileExt = target.getExtension()
		def bamName       = bam.getBaseName()
		def minbq         = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.pisces_minbq_IonTorrent : params.pisces_minbq
		def suffix        = params.vcfCalling_suffix ?: '' 
		def args          = task.ext.args ?: ''
		def logPrefix     = "${step}_DNA_${sample}_${chr}"

		"""
		awk 'BEGIN {FS="\\t";OFS="\\t"}; /^@/ {print; next}; {if (\$1 == "${chr}") { print }}' ${target} > tmp_${chr}.${targetFileExt}

		Pisces \\
			-bam ${bam} \\
			-g ${genomeDir} \\
			-o . \\
			-i tmp_${chr}.${targetFileExt} \\
			--maxthreads ${task.cpus} \\
			--gvcf false \\
			--callmnvs ${params.pisces_callmnvs} \\
			--outputsbfiles ${params.pisces_outputsbfiles} \\
			--collapse ${params.pisces_collapse} \\
			--usestitchedxd ${params.pisces_usestitchedxd} \\
			--minbasecallquality ${minbq} \\
			--minmapquality ${params.pisces_minmapq} \\
			--minvariantqscore ${params.pisces_minqscore} \\
			--mindepth ${params.pisces_mindepth} \\
			--minimumfrequency ${params.pisces_minfreq} \\
			--variantqualityfilter ${params.pisces_minqscorefilter} \\
			--maxvariantqscore ${params.pisces_maxqscorefilter} \\
			--ssfilter ${params.pisces_ssfilter} \\
			--reportnocalls ${params.pisces_reportnocalls} \\
			${args}

		mv ${bamName}.vcf ${sample}_DNA_SmallVariant_pisces_${chr}${suffix}.vcf
		mv PiscesLogs PiscesLogs_${chr}
		
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
