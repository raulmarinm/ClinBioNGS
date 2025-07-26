/*
 * Prepare sampleSheet
 * 	- Illumina: From 'SampleSheet.csv' file in startingDir
 * 	- IonTorrent: From 'Info.csv' file in startingDir
 */

process PREPARE_SAMPLESHEET {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: { params.sampleSheetDir },
		mode: 'copy',
		pattern: '*.csv'
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path inputDir

	output:
		path '*.csv', emit: csv
		path '*.{err,log,out,run,sh}'

	script:
		def run       = params.runName
		def logPrefix = "PREPARE_SAMPLESHEET"

		if( params.prepareIontorrentSamplesheet )
		"""
		for dir in \$(find -L ${inputDir} -name "*.xz" | grep -v temp.tar.xz)
		do
			tar -xf \${dir}
		done

		prepareIonTorrentSampleSheet.R --output=SampleSheet_${run}.csv --run=${run}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
		else
		"""
		cp ${inputDir}/SampleSheet.csv SampleSheet_${run}.csv

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}
