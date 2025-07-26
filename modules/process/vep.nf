/*
 * Small variant annotation - VEP
 */

// VEP running
process VEP_RUN {

	tag "${params.runName}"

	label 'ensemblVep'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${subworkflow}",
		mode: 'copy',
		pattern: '*.{html,txt,vcf}',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path inVcf
		val outVcf
		path './*'
		tuple path(caddSnv), path(caddSnvTbi)
		tuple path(caddIndel), path(caddIndelTbi)
		tuple path(revelTsv), path(revelTbi)
		tuple path(alphaMissenseFileTsv), path(alphaMissenseFileTbi)
		tuple path(clinvarVcf), path(clinvarTbi)
		tuple path(civicVcf), path(civicTbi)
		val subworkflow
		val step

	output:
		path "${outVcf}", emit: vcf
		path '*.{html,txt,vcf}'
		path '*.{err,log,out,run,sh}'
	
	script:
		def args       = task.ext.args ?: ''
		def logPrefix  = "${step}_DNA"

		"""
		cat ${inVcf} | \\
			grep -v "^#" | cut -f1-5 |  sort -k1,1V -k2,2n | uniq | \\
			awk -F'\\t' -v OFS='\\t' '{print \$1,\$2,".",\$4,\$5,".",".","."}' > VEP_input.vcf
		
		vep \\
			--dir_cache . --offline --cache \\
			--assembly ${params.genomeVersionGrc} \\
			--fork ${task.cpus} \\
			-i VEP_input.vcf \\
			-o ${outVcf} \\
			--vcf -v --force_overwrite --pick \\
			--variant_class --total_length --numbers \\
			--sift b --polyphen b \\
			--hgvs --hgvsg --shift_hgvs ${params.vepShiftHgvs} \\
			--symbol --protein --ccds --biotype \\
			--canonical --mane --tsl --appris \\
			--af_gnomade --af_gnomadg \\
			--distance ${params.vepDistance} \\
			--check_existing \\
			--plugin NMD \\
			--plugin CADD,snv=${caddSnv},indels=${caddIndel} \\
			--plugin REVEL,${revelTsv} \\
			--plugin AlphaMissense,file=${alphaMissenseFileTsv} \\
			--custom ${clinvarVcf},ClinVar,vcf,exact,0,${params.vepClinvarFields} \\
			--custom ${civicVcf},CIViC,vcf,exact,0,${params.vepCivicFields} \\
			${args}
		 

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}

// VEP processing
process VEP_PROCESSING {

	tag "${params.runName}"

	label 'R'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: "${params.analysisDir}/${subworkflow}",
		mode: 'copy',
		pattern: '*.txt',
	)
	publishDir(
		path: { "${params.logsDir}/${subworkflow}/${step}" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		path vcf
		path hotspotsWhitelistBed
		path problematicRegionsBed
		path ctrRegionsBed
		path mmrGenes
		path cancerdrivers
		path cancerHotspotsResults
		path sopOncogenic
		path genieCounts
		path genieOncogenic
		path cgiOncogenic
		path civicOncogenic
		path panelRecurrentMutations
		val subworkflow
		val step

	output:
		path "*_DNA_SmallVariant_annot.txt", emit: results
		path '*.{err,log,out,run,sh}'
	
	script:
		def outFile = "${params.runName}_DNA_SmallVariant_annot.txt"
		def logPrefix  = "${step}_DNA"

		"""
		vepProcessing.R \\
			--vcf=${vcf} \\
			--outFile=${outFile} \\
			--hotspotsWhitelistBed=${hotspotsWhitelistBed} \\
			--problematicRegionsBed=${problematicRegionsBed} \\
			--ctrRegionsBed=${ctrRegionsBed} \\
			--mmrGenes=${mmrGenes} \\
			--cancerdrivers=${cancerdrivers} \\
			--cancerHotspotsResults=${cancerHotspotsResults} \\
			--sopOncogenic=${sopOncogenic} \\
			--genieCounts=${genieCounts} \\
			--genieOncogenic=${genieOncogenic} \\
			--cgiOncogenic=${cgiOncogenic} \\
			--civicOncogenic=${civicOncogenic} \\
			--panelRecurrentMutations=${panelRecurrentMutations} \\
			--caddCutoff=${params.caddCutoff} \\
			--revelCutoff=${params.revelCutoff} \\
			--alphamissenseCutoff=${params.alphamissenseCutoff} \\
			--somaticMaxpVAF=${params.variantSomaticMaxpVAF}
		
		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}


// Download VEP cache data
process DOWNLOAD_VEP_CACHE {

	tag "${params.runName}"

	label 'ensemblVep'
	label 'minCpu'
	label 'minMem'

	publishDir(
		path: params.annotationDir,
		mode: params.vepCache_publishMode,
		pattern: "${specie}/${version}_${genome}"
	)
	publishDir(
		path: { "${params.logsDir}/RESOURCES" },
		mode: 'copy',
		pattern: '*.{err,log,out,run,sh}'
	)

	input:
		val specie
		val genome
		val version

	output:
		path "${specie}", emit: dir
		path "${specie}/${version}_${genome}"
		path '*.{err,log,out,run,sh}'

	script:
		def logPrefix = "DOWNLOAD_VEP_CACHE"

		"""
		INSTALL.pl -c "." -a cf -s ${specie} -y ${genome}

		for file in .command.{err,log,out,run,sh}; do cp -a \${file} ${logPrefix}\${file}; done
		"""

}
