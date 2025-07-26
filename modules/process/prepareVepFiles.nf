/*
 * Prepare resources - VEP files
 */

// REVEL
process PREPARE_REVEL_FILE{

	tag "${params.runName}"

	label 'gatk4'
	label 'minCpu'
	label 'lowMem'
	label 'outAnnotDir'

	input:
		path zip
		val outName

	output:
		tuple path('*.tsv.gz'), path('*.tsv.gz.tbi'), emit: files
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "PREPARE_REVEL_FILE"

		"""
		unzip ${zip}
		cat revel_with_transcript_ids | tr "," "\\t" > tmp_tabbed_revel.tsv
		sed '1s/.*/#&/' tmp_tabbed_revel.tsv > tmp_new_tabbed_revel.tsv
		bgzip tmp_new_tabbed_revel.tsv
		zcat tmp_new_tabbed_revel.tsv.gz | head -n1 > tmp_h
		zgrep -h -v ^#chr tmp_new_tabbed_revel.tsv.gz | awk '\$3 != "." ' | sort -k1,1 -k3,3n - | cat tmp_h - | bgzip -c > ${outName}
		tabix -f -s 1 -b 3 -e 3 ${outName}

		rm revel_with_transcript_ids
		rm tmp_*

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

// Alphamissense
process INDEX_ALPHAMISSENSE_FILE{

	tag "${params.runName}"

	label 'gatk4'
	label 'minCpu'
	label 'minMem'
	label 'outAnnotDir'

	input:
		path input

	output:
		path "*.tbi", emit: file
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "INDEX_ALPHAMISSENSE_FILE"

		"""
		tabix -s 1 -b 2 -e 2 -f -S 1 ${input}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""
		
}

