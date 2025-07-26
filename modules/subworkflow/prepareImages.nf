/*
 * Prepare singularity images workflow
 */

include { 
	SIF_DOWNLOAD as SIF_DOWNLOAD_ABRA2; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_BCFTOOLS; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_BEDTOOLS; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_BIOAWK; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_BWAMEM2; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_CNVKIT; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_CTATSPLICING; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_FASTP; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_FASTQC; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_GENCORE; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_MOSDEPTH; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_MSISENSORPRO; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_MULTIQC; 
	SIF_DOWNLOAD as SIF_DOWNLOAD_SAMTOOLS;
	SIF_DOWNLOAD as SIF_DOWNLOAD_STARFUSION;
	SIF_DOWNLOAD as SIF_DOWNLOAD_UCSCBIGBEDTOBED;
	SIF_DOWNLOAD as SIF_DOWNLOAD_UCSCLIFTOVER;
	SIF_DOWNLOAD as SIF_DOWNLOAD_UMITOOLS;
	SIF_DOWNLOAD as SIF_DOWNLOAD_UMITRANSFER;
	SIF_DOWNLOAD as SIF_DOWNLOAD_VARDICT;
	SIF_DOWNLOAD as SIF_DOWNLOAD_VT;

	SIF_BUILD as SIF_BUILD_AWSCLI; 
	SIF_BUILD as SIF_BUILD_BCLCONVERT;
	SIF_BUILD as SIF_BUILD_VEP;
	SIF_BUILD as SIF_BUILD_GATK;
	SIF_BUILD as SIF_BUILD_OCTOPUS;
	SIF_BUILD as SIF_BUILD_PISCES;
	SIF_BUILD as SIF_BUILD_R;
	SIF_BUILD as SIF_BUILD_TMAPTVC;
	SIF_BUILD as SIF_BUILD_XENGSORT;

} from '../process/prepareImages'


workflow PREPARE_TOOL_IMAGES {

	main:

		// Abra2
		if (!file(params.abra2Img).exists()) { SIF_DOWNLOAD_ABRA2(params.abra2Link, 'abra2', params.abra2Version) }

		// AWS CLI
		if (!file(params.awscliImg).exists()) { SIF_BUILD_AWSCLI(params.awscliLink, 'awscli', params.awscliVersion) }

		// Bcftools
		if (!file(params.bcftoolsImg).exists()) { SIF_DOWNLOAD_BCFTOOLS(params.bcftoolsLink, 'bcftools', params.bcftoolsVersion) }

		// BCL-Convert
		if (!file(params.bclconvertImg).exists()) { SIF_BUILD_BCLCONVERT(params.bclconvertLink, 'bclconvert', params.bclconvertVersion) }
		
		// Bedtools
		if (!file(params.bedtoolsImg).exists()) { SIF_DOWNLOAD_BEDTOOLS(params.bedtoolsLink, 'bedtools', params.bedtoolsVersion) }

		// Bioawk
		if (!file(params.bioawkImg).exists()) { SIF_DOWNLOAD_BIOAWK(params.bioawkLink, 'bioawk', params.bioawkVersion) }

		// BWA-mem2
		if (!file(params.bwamem2Img).exists()) { SIF_DOWNLOAD_BWAMEM2(params.bwamem2Link, 'bwa-mem2', params.bwamem2Version) }

		// CNVkit
		if (!file(params.cnvkitImg).exists()) { SIF_DOWNLOAD_CNVKIT(params.cnvkitLink, 'cnvkit', params.cnvkitVersion) }

		// CTAT-splicing
		if (!file(params.ctatSplicingImg).exists()) { SIF_DOWNLOAD_CTATSPLICING(params.ctatSplicingLink, 'ctat-splicing', params.ctatSplicingVersion) }

		// Ensembl-VEP
		if (!file(params.ensemblVepImg).exists()) { SIF_BUILD_VEP(params.ensemblVepLink, 'ensembl-vep', params.ensemblVepVersion) }

		// FastP
		if (!file(params.fastpImg).exists()) { SIF_DOWNLOAD_FASTP(params.fastpLink, 'fastp', params.fastpVersion) }

		// FastQC
		if (!file(params.fastqcImg).exists()) { SIF_DOWNLOAD_FASTQC(params.fastqcLink, 'fastqc', params.fastqcVersion) }

		// GATK4
		if (!file(params.gatk4Img).exists()) { SIF_BUILD_GATK(params.gatk4Link, 'gatk', params.gatk4Version) }

		// Gencore
		if (!file(params.gencoreImg).exists()) { SIF_DOWNLOAD_GENCORE(params.gencoreLink, 'gencore', params.gencoreVersion) }

		// Mosdepth
		if (!file(params.mosdepthImg).exists()) { SIF_DOWNLOAD_MOSDEPTH(params.mosdepthLink, 'mosdepth', params.mosdepthVersion) }

		// MSIsensor-pro
		if (!file(params.msisensorproImg).exists()) { SIF_DOWNLOAD_MSISENSORPRO(params.msisensorproLink, 'msisensor-pro', params.msisensorproVersion) }

		// MultiQC
		if (!file(params.multiqcImg).exists()) { SIF_DOWNLOAD_MULTIQC(params.multiqcLink, 'multiqc', params.multiqcVersion) }

		// Octopus
		if (!file(params.octopusImg).exists()) { SIF_BUILD_OCTOPUS(params.octopusLink, 'octopus', params.octopusVersion) }

		// Pisces
		if (!file(params.piscesImg).exists()) { SIF_BUILD_PISCES(params.piscesLink, 'pisces', params.piscesVersion) }

		// R
		if (!file(params.rImg).exists()) { SIF_BUILD_R(params.rLink, 'R', params.rVersion) }

		// Samtools
		if (!file(params.samtoolsImg).exists()) { SIF_DOWNLOAD_SAMTOOLS(params.samtoolsLink, 'samtools', params.samtoolsVersion) }

		// STAR-fusion
		if (!file(params.starFusionImg).exists()) { SIF_DOWNLOAD_STARFUSION(params.starFusionLink, 'star-fusion', params.starFusionVersion) }

		// TMAP-TVC
		if (!file(params.tmapTvcImg).exists()) { SIF_BUILD_TMAPTVC(params.tmapTvcLink, 'tmap-tvc', params.tmapTvcVersion) }

		// UCSC-bigBedToBed
		if (!file(params.ucscBigbedtobedImg).exists()) { SIF_DOWNLOAD_UCSCBIGBEDTOBED(params.ucscBigbedtobedLink, 'ucsc-bigbedtobed', params.ucscBigbedtobedVersion) }

		// UCSC-liftover
		if (!file(params.ucscLiftoverImg).exists()) { SIF_DOWNLOAD_UCSCLIFTOVER(params.ucscLiftoverLink, 'ucsc-liftover', params.ucscLiftoverVersion) }

		// UMI-tools
		if (!file(params.umitoolsImg).exists()) { SIF_DOWNLOAD_UMITOOLS(params.umitoolsLink, 'umitools', params.umitoolsVersion) }

		// UMI-transfer
		if (!file(params.umitransferImg).exists()) { SIF_DOWNLOAD_UMITRANSFER(params.umitransferLink, 'umi-transfer', params.umitransferVersion) }

		// VarDict
		if (!file(params.vardictImg).exists()) { SIF_DOWNLOAD_VARDICT(params.vardictLink, 'vardict-java', params.vardictVersion) }

		// Vt
		if (!file(params.vtImg).exists()) { SIF_DOWNLOAD_VT(params.vtLink, 'vt', params.vtVersion) }

		// Xengsort
		if (!file(params.xengsortImg).exists()) { SIF_BUILD_XENGSORT(params.xengsortLink, 'xengsort', params.xengsortVersion) }


}
