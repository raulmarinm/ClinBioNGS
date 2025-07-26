/*
 * Prepare Ion Torrent files from the results folder
 */

process PREPARE_MANIFEST_IONTORRENT_DNA {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path inputDir

	output:
		path '*.bed', emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_MANIFEST_IONTORRENT_DNA"

		"""
		results=\$(find -L ${inputDir} -name "*.xz" | grep -v temp.tar.xz | head -n1)
		tar -xf \${results}

		prepareIonTorrentManifestDna.R --output=${params.seqPanel}_DNA_manifest.bed

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

process PREPARE_MANIFEST_IONTORRENT_RNA {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outManifestDir'

	input:
		path inputDir

	output:
		path '*.bed', emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_MANIFEST_IONTORRENT_RNA"

		"""
		results=\$(find -L ${inputDir} -name "*.xz" | grep -v temp.tar.xz | head -n1)
		tar -xf \${results}

		prepareIonTorrentManifestRna.R --output=${params.seqPanel}_RNA_manifest.bed --padding=${params.ionTorrentRnaPadding}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

process PREPARE_HOTSPOTS_IONTORRENT {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'
	label 'outSmallVariantDir'

	input:
		path inputDir
		path chain

	output:
		path '*.bed', emit: bed
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_HOTSPOTS_IONTORRENT"

		"""
		results=\$(find -L ${inputDir} -name "*.xz" | grep -v temp.tar.xz | head -n1)
		tar -xf \${results}

		prepareIonTorrentHotspots.R --output=${params.seqPanel}_DNA_${params.genomeVersionHg}.hotspots.bed --chain=${chain}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

process EXTRACT_UBAM_IONTORRENT {

	tag "${params.runName}"

	label 'minCpu'
	label 'lowMem'

	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path inputDir
		val subworkflow
		val step
		val type

	output:
		path '*rawlib.basecaller.bam', emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def logPrefix = "${step}_${type}"

		"""
		for dir in \$(find -L ${inputDir} -name "*.xz" | grep -v temp.tar.xz)
		do
			tar -xf \${dir}
		done

		all_bams=\$(find . -iname '*rawlib.basecaller.bam')
		for bam in \${all_bams}
		do
			ln -s \${bam} \$(basename \${bam})
		done

		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}


process PREPARE_UBAM_IONTORRENT {

	tag "${sample}"

	label 'minCpu'
	label 'lowMem'

	publishDir(
		enabled: params.ubamIonTorrent_publish,
		path: "${params.dataRunDir}/${sample}/BAM",
		mode: params.ubamIonTorrent_publishMode,
		pattern: "${sample}_${type}*.bam"
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${type}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		tuple val(sample), val(barcode)
		path('data/*')
		val subworkflow
		val step
		val type

	output:
		tuple val(sample), path("${sample}_${type}*.bam"), emit: bam
		path '*.{err,log,out,run,sh}'
	
	script:
		def suffix    = params.ubamIonTorrent_suffix ?: ''
		def logPrefix = "${step}_${type}_${sample}"

		"""
		bam=\$(ls data/*.bam | grep ${barcode})

		ln -s \${bam} ${sample}_${type}${suffix}.bam

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

