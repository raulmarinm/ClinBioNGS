/*
 * Small variant analysis - Calling
 */

include { BAM_CLIPPING                } from '../process/bamClipping'
include { BAM_REMOVE_TAGS             } from '../process/bamRemoveTags'
include { PISCES                      } from '../process/pisces'
include { 
	MUTECT2_CALL; MUTECT2_MODEL; 
	MUTECT2_FILTER                    } from '../process/mutect2'
include { VARDICT                     } from '../process/vardict'
include { OCTOPUS                     } from '../process/octopus'
include { TVC                         } from '../process/tvc'
include { 
	VCF_CONCAT as VCF_CONCAT_PISCES;
	VCF_CONCAT as VCF_CONCAT_MUTECT2;
	VCF_CONCAT as VCF_CONCAT_VARDICT;
	VCF_CONCAT as VCF_CONCAT_OCTOPUS;
	VCF_CONCAT as VCF_CONCAT_TVC;     } from '../process/vcfConcat'
include { 
	BCFTOOLS as BCFTOOLS_PISCES; 
	BCFTOOLS as BCFTOOLS_MUTECT2; 
	BCFTOOLS as BCFTOOLS_VARDICT; 
	BCFTOOLS as BCFTOOLS_OCTOPUS;
	BCFTOOLS as BCFTOOLS_TVC          } from '../process/bcftools'
include { 
	VT as VT_PISCES; 
	VT as VT_MUTECT2; 
	VT as VT_VARDICT; 
	VT as VT_OCTOPUS;
	VT as VT_TVC;                     } from '../process/vt'
include { VCF_CONSENSUS               } from '../process/vcfConsensus'


workflow VARIANT_CALLING {

	take:
		bam
		genomeFastaFiles
		genomeFastaDir
		targetBed
		targetBedExt
		targetIntervalListExt
		offTargetBed
		targetChr
		tvcJson
		chainFile
		panelHotspotsBed
		panelBlacklistBed
		vcfHeader
	
	main:

		chCallerNames = params.seqPlatform.toLowerCase() == 'iontorrent' ? params.variantCallers_IonTorrent.toLowerCase() : params.variantCallers.toLowerCase()
		chVcfCallers  = Channel.empty()
		chChromosomes = targetChr.splitCsv().flatten()

		if (params.startingDataType.toUpperCase() == 'VCF') {
			chVcfCallers = Channel.fromPath("${params.startingDataDir}/**_DNA*{${chCallerNames}}${params.vcfProcessing_suffix}.vcf")
				.map{ [it.simpleName.toString().split("_DNA")[0], it] }
				.groupTuple()
			
			// 1. VCF consensus
			VCF_CONSENSUS(chVcfCallers, chCallerNames, chainFile, panelHotspotsBed, panelBlacklistBed, vcfHeader, '06_VARIANT_CALLING', '01_VCF_CONSENSUS')

		} else {

			// 0. BAM clipping
			if (params.bamclipping) {
				bam = BAM_CLIPPING(bam, offTargetBed, '06_VARIANT_CALLING', '00_BAM_CLIPPING', 'DNA').bam
			}

			// 1. PISCES
			if (params.pisces && chCallerNames.split(',').contains('pisces')) {
				if (params.seqPlatform.toLowerCase() == 'iontorrent') {
					chPiscesBam = BAM_REMOVE_TAGS(bam, params.bamTags_IonTorrent, '06_VARIANT_CALLING', '00_PISCES_BAM_REMOVE_TAGS', 'DNA').bam
				} else { chPiscesBam = bam }
				PISCES(chPiscesBam.combine(chChromosomes), genomeFastaDir, targetIntervalListExt, '06_VARIANT_CALLING', '01_PISCES')
				VCF_CONCAT_PISCES(PISCES.out.vcf.groupTuple(), genomeFastaFiles, '06_VARIANT_CALLING', '02_VCF_CONCAT_PISCES', 'pisces')
				BCFTOOLS_PISCES(VCF_CONCAT_PISCES.out.vcf, targetBed, '06_VARIANT_CALLING', '03_BCFTOOLS_PISCES')
				VT_PISCES(BCFTOOLS_PISCES.out.vcf, genomeFastaFiles, '06_VARIANT_CALLING', '04_VT_PISCES', 'pisces')
				chVcfPisces = VT_PISCES.out.vcf
			} else { chVcfPisces = Channel.empty() }

			// 1. MUTECT2
			if (params.mutect2 && chCallerNames.split(',').contains('mutect2')) {
				MUTECT2_CALL(bam.combine(chChromosomes), genomeFastaFiles, targetIntervalListExt, '06_VARIANT_CALLING', '01_MUTECT2_CALL')
				MUTECT2_MODEL(MUTECT2_CALL.out.f1r2.groupTuple(), '06_VARIANT_CALLING', '01_MUTECT2_MODEL')
				MUTECT2_FILTER(MUTECT2_CALL.out.vcf.combine(MUTECT2_MODEL.out.model, by: 0), genomeFastaFiles, '06_VARIANT_CALLING', '01_MUTECT2_FILTER')
				VCF_CONCAT_MUTECT2(MUTECT2_FILTER.out.vcf.groupTuple(), genomeFastaFiles, '06_VARIANT_CALLING', '02_VCF_CONCAT_MUTECT2', 'mutect2')
				BCFTOOLS_MUTECT2(VCF_CONCAT_MUTECT2.out.vcf, targetBed, '06_VARIANT_CALLING', '03_BCFTOOLS_MUTECT2')
				VT_MUTECT2(BCFTOOLS_MUTECT2.out.vcf, genomeFastaFiles, '06_VARIANT_CALLING', '04_VT_MUTECT2', 'mutect2')
				chVcfMutect2 = VT_MUTECT2.out.vcf
			} else { chVcfMutect2 = Channel.empty() }

			// 1. VARDICT
			if (params.vardict && chCallerNames.split(',').contains('vardict')) {
				VARDICT(bam.combine(chChromosomes), genomeFastaFiles, targetBed, '06_VARIANT_CALLING', '01_VARDICT')
				VCF_CONCAT_VARDICT(VARDICT.out.vcf.groupTuple(), genomeFastaFiles, '06_VARIANT_CALLING', '02_VCF_CONCAT_VARDICT', 'vardict')
				BCFTOOLS_VARDICT(VCF_CONCAT_VARDICT.out.vcf, targetBed, '06_VARIANT_CALLING', '03_BCFTOOLS_VARDICT')
				VT_VARDICT(BCFTOOLS_VARDICT.out.vcf, genomeFastaFiles, '06_VARIANT_CALLING', '04_VT_VARDICT', 'vardict')
				chVcfVardict = VT_VARDICT.out.vcf
			} else { chVcfVardict = Channel.empty() }

			// 1. OCTOPUS
			if (params.octopus && chCallerNames.split(',').contains('octopus')) {
				OCTOPUS(bam.combine(chChromosomes), genomeFastaFiles, targetBedExt, '06_VARIANT_CALLING', '01_OCTOPUS')
				VCF_CONCAT_OCTOPUS(OCTOPUS.out.vcf.groupTuple(), genomeFastaFiles, '06_VARIANT_CALLING', '02_VCF_CONCAT_OCTOPUS', 'octopus')
				BCFTOOLS_OCTOPUS(VCF_CONCAT_OCTOPUS.out.vcf, targetBed, '06_VARIANT_CALLING', '03_BCFTOOLS_OCTOPUS')
				VT_OCTOPUS(BCFTOOLS_OCTOPUS.out.vcf, genomeFastaFiles, '06_VARIANT_CALLING', '04_VT_OCTOPUS', 'octopus')
				chVcfOctopus = VT_OCTOPUS.out.vcf
			} else { chVcfOctopus = Channel.empty() }

			// 1. TVC (Ion torrent data)
			if (params.tvc && chCallerNames.split(',').contains('tvc')) {
				TVC(bam.combine(chChromosomes), genomeFastaFiles, targetBedExt, tvcJson, '06_VARIANT_CALLING', '01_TVC')
				VCF_CONCAT_TVC(TVC.out.vcf.groupTuple(), genomeFastaFiles, '06_VARIANT_CALLING', '02_VCF_CONCAT_TVC', 'tvc')
				BCFTOOLS_TVC(VCF_CONCAT_TVC.out.vcf, targetBed, '06_VARIANT_CALLING', '03_BCFTOOLS_TVC')
				VT_TVC(BCFTOOLS_TVC.out.vcf, genomeFastaFiles, '06_VARIANT_CALLING', '04_VT_TVC', 'tvc')
				chVcfTvc = VT_TVC.out.vcf
			} else { chVcfTvc = Channel.empty() }

			chVcfCallers = chVcfCallers.concat(chVcfPisces, chVcfMutect2, chVcfVardict, chVcfOctopus, chVcfTvc).groupTuple()
			chEnabledCallerNames = [
				pisces  : params.pisces && chCallerNames.split(',').contains('pisces'),
				mutect2 : params.mutect2 && chCallerNames.split(',').contains('mutect2'),
				vardict : params.vardict && chCallerNames.split(',').contains('vardict'),
				octopus : params.octopus && chCallerNames.split(',').contains('octopus'),
				tvc     : params.tvc && chCallerNames.split(',').contains('tvc')
			].findAll{ k,v -> v }.keySet().join(',')

			// 2. VCF consensus
			VCF_CONSENSUS(chVcfCallers, chEnabledCallerNames, chainFile, panelHotspotsBed, panelBlacklistBed, vcfHeader, '06_VARIANT_CALLING', '05_VCF_CONSENSUS')
		}

	
	emit:
		vcf     = VCF_CONSENSUS.out.vcf
		results = VCF_CONSENSUS.out.results
		
}
