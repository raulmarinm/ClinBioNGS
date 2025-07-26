#!/usr/bin/env Rscript

# ============================================================================ #
# What: Collect run results and create a library
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(DBI)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
run <- args_list[["run"]]
dna_sample.file <- args_list[["dnaSampleFile"]]
rna_sample.file <- args_list[["rnaSampleFile"]]
vars_min_callers <- as.numeric(args_list[["variantMinCallers"]])
vars_min_ad <- as.numeric(args_list[["variantMinAD"]])
vars_min_dp <- as.numeric(args_list[["variantMinDP"]])
vars_min_vaf <- as.numeric(args_list[["variantMinVAF"]])
vars_somatic_max_pvaf <- as.numeric(args_list[["variantSomaticMaxpVAF"]])
vars_somatic_max_vaf <- as.numeric(args_list[["variantSomaticMaxVAF"]])
cna_min_target_bins <- as.numeric(args_list[["cnaMinTargetBins"]])
cna_CN_AMP_OK <- as.numeric(args_list[["cnaHighAmpCN"]])
cna_CN_DEL_OK <- as.numeric(args_list[["cnaHighDelCN"]])
cna_arm_call_threshold <- as.numeric(args_list[["cnaArmCallingThreshold"]])
fusion_calling_min_ad <- as.numeric(args_list[["fusionCallingMinAD"]])
fusion_calling_min_ffpm <- as.numeric(args_list[["fusionCallingMinFFPM"]])
fusion_flag_min_ad <- as.numeric(args_list[["fusionFlagMinAD"]])
fusion_flag_min_ad_nonfused <- as.numeric(args_list[["fusionFlagMinNonFused"]])
fusion_flag_min_dp <- as.numeric(args_list[["fusionFlagMinDP"]])
fusion_flag_min_vaf <- as.numeric(args_list[["fusionFlagMinVAF"]])
splicing_calling_min_ad <- as.numeric(args_list[["splicingCallingMinAD"]])
splicing_flag_min_ad <- as.numeric(args_list[["splicingFlagMinAD"]])
splicing_flag_min_dp <- as.numeric(args_list[["splicingFlagMinDP"]])
splicing_flag_min_vaf <- as.numeric(args_list[["splicingFlagMinVAF"]])
library.db <- paste0(run, "_", format(Sys.time(), "%Y%m%d"), ".db")


# Load samples ---------------------------------------------------------------

all_sample_info <- matrix(data = NA, nrow = 0, ncol = 1)
if(file.exists(dna_sample.file)){
	dna_sample_info <- read.table(dna_sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Library) %>% mutate(TYPE = "DNA")
	all_sample_info <- rbind(all_sample_info, dna_sample_info)
	rm(dna_sample_info)
}
if(file.exists(rna_sample.file)){
	rna_sample_info <- read.table(rna_sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Library) %>% mutate(TYPE = "RNA")
	all_sample_info <- rbind(all_sample_info, rna_sample_info)
	rm(rna_sample_info)
}

if(!nrow(all_sample_info)){
	library_conn <- dbConnect(RSQLite::SQLite(), library.db)
	dbWriteTable(library_conn, "Samples", as.data.frame(all_sample_info), overwrite = TRUE)
	dbDisconnect(library_conn)

	cat(paste0("\n", "0 samples for registering in ", run, " run", "\n"))
	
} else {

	all_sample_info <- all_sample_info %>% mutate(Sample = as.character(Sample))
	dna_samples <- (all_sample_info %>% filter(TYPE == "DNA"))$Sample
	rna_samples <- (all_sample_info %>% filter(TYPE == "RNA"))$Sample

	all_sample_info <- all_sample_info %>% distinct(Sample, .keep_all = TRUE) %>%
		mutate(
			ID = paste0(run, "_" , Sample), RUN = run, SAMPLE = Sample, SEX = Sex, AGE = Age, PURITY = Purity, 
			TUMOR_NAME = TumorName, TUMOR_DOID = TumorDOID, TUMOR_CODE = TumorCode,
			TOPNODE_NAME = TopNodeName, TOPNODE_DOID = TopNodeDOID, 
			BAMQC_DNA = FALSE, BAMQC_RNA = FALSE, SMALL_VARIANT = FALSE, CNA = FALSE, FUSION = FALSE, SPLICING = FALSE, 
			TMB = FALSE, MSI = FALSE) %>%
		select(
			ID, RUN, SAMPLE, SEX, AGE, PURITY, TUMOR_NAME, TUMOR_DOID, TUMOR_CODE, TOPNODE_NAME, TOPNODE_DOID,
			BAMQC_DNA, BAMQC_RNA, SMALL_VARIANT, CNA, FUSION, SPLICING, TMB, MSI)
	rownames(all_sample_info) <- all_sample_info$SAMPLE
	all_bamqc_dna <- matrix(data = NA, nrow = 0, ncol = 1)
	all_variants <- matrix(data = NA, nrow = 0, ncol = 1)
	all_variants_annot <- matrix(data = NA, nrow = 0, ncol = 1)
	all_cna_gene <- matrix(data = NA, nrow = 0, ncol = 1)
	all_cna_arm <- matrix(data = NA, nrow = 0, ncol = 1)
	all_tmb <- matrix(data = NA, nrow = 0, ncol = 1)
	all_msi <- matrix(data = NA, nrow = 0, ncol = 1)
	all_bamqc_rna <- matrix(data = NA, nrow = 0, ncol = 1)
	all_fusion <- matrix(data = NA, nrow = 0, ncol = 1)
	all_splicing <- matrix(data = NA, nrow = 0, ncol = 1)


	# Collect DNA results -------------------------------------------------------

	cat(paste0("\n", "Collecting DNA results from ", run, ":\n"))

	for (sample in dna_samples){

		sample_id <- paste0(run, "_" , sample)

		# Set result files
		bamqc_dna.file <- list.files(pattern = paste0(sample, "_DNA_BamQC.txt"))
		variants.file <- list.files(pattern = paste0(sample, "_DNA_SmallVariant_consensus.txt"))
		variants_annot.file <- list.files(pattern = paste0(sample, "_DNA_SmallVariant_annot.txt"))
		cna_gene.file <- list.files(pattern = paste0(sample, "_DNA_CNA_gene.txt"))
		cna_arm.file <- list.files(pattern = paste0(sample, "_DNA_CNA_arm.txt"))
		tmb.file <- list.files(pattern = paste0(sample, "_DNA_TMB_results.txt"))
		msi.file <- list.files(pattern = paste0(sample, "_DNA_MSI_results.txt"))

		# Load results
		if(length(bamqc_dna.file)){
			bamqc_dna <- read.table(bamqc_dna.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_bamqc_dna <- rbind(all_bamqc_dna, bamqc_dna)
			all_sample_info[sample, "BAMQC_DNA"] <- TRUE
			rm(bamqc_dna)
		}

		if(length(variants.file)){
			variants <- read.table(variants.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_variants <- rbind(all_variants, variants)
			all_sample_info[sample, "SMALL_VARIANT"] <- TRUE
			rm(variants)
		}

		if(length(variants_annot.file)){
			variants_annot <- read.table(variants_annot.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_variants_annot <- rbind(all_variants_annot, variants_annot)
			rm(variants_annot)
		}

		if(length(cna_gene.file)){
			cna_gene <- read.table(cna_gene.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_cna_gene <- rbind(all_cna_gene, cna_gene)
			all_sample_info[sample, "CNA"] <- TRUE
			rm(cna_gene)
		}

		if(length(cna_arm.file)){
			cna_arm <- read.table(cna_arm.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_cna_arm <- rbind(all_cna_arm, cna_arm)
			rm(cna_arm)
		}

		if(length(tmb.file)){
			tmb <- read.table(tmb.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_tmb <- rbind(all_tmb, tmb)
			all_sample_info[sample, "TMB"] <- TRUE
			rm(tmb)
		}

		if(length(msi.file)){
			msi <- read.table(msi.file, header = TRUE, sep = "\t") %>% 
				mutate(
					ID = sample_id, RUN = run, SAMPLE = sample, TotalSites = Total_Number_of_Sites,
					MsiSites = Number_of_Somatic_Sites, MSI_PCT = round(Number_of_Somatic_Sites * 100 / Total_Number_of_Sites, 1)) %>% 
				select(ID, RUN, SAMPLE, MSI_PCT, MsiSites, TotalSites)
			all_msi <- rbind(all_msi, msi)
			all_sample_info[sample, "MSI"] <- TRUE
			rm(msi)
		}

		rm(sample_id, bamqc_dna.file, variants.file, variants_annot.file, cna_gene.file, cna_arm.file, tmb.file, msi.file)
		
		cat(paste0("   - ", sample, "\n"))
	}
	rm(sample)


	# Collect RNA results -------------------------------------------------------

	cat(paste0("\n", "Collecting RNA results from ", run, ":\n"))

	for (sample in rna_samples){

		sample_id <- paste0(run, "_" , sample)

		# Set result files
		bamqc_rna.file <- list.files(pattern = paste0(sample, "_RNA_BamQC.txt"))
		fusion.file <- list.files(pattern = paste0(sample, "_RNA_Fusion_results.txt"))
		splicing.file <- list.files(pattern = paste0(sample, "_RNA_Splicing_results.txt"))

		# Load results
		if(length(bamqc_rna.file)){
			bamqc_rna <- read.table(bamqc_rna.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_bamqc_rna <- rbind(all_bamqc_rna, bamqc_rna)
			all_sample_info[sample, "BAMQC_RNA"] <- TRUE
			rm(bamqc_rna)
		}

		if(length(fusion.file)){
			fusion <- read.table(fusion.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_fusion <- rbind(all_fusion, fusion)
			all_sample_info[sample, "FUSION"] <- TRUE
			rm(fusion)
		}

		if(length(splicing.file)){
			splicing <- read.table(splicing.file, header = TRUE, sep = "\t") %>% mutate(ID = sample_id) %>% select(ID, everything())
			all_splicing <- rbind(all_splicing, splicing)
			all_sample_info[sample, "SPLICING"] <- TRUE
			rm(splicing)
		}

		rm(sample_id, bamqc_rna.file, fusion.file, splicing.file)

		cat(paste0("   - ", sample, "\n"))

	}
	rm(sample)


	# Save data ---------------------------------------------------------------

	# Write tables to the run DB
	library_conn <- dbConnect(RSQLite::SQLite(), library.db)

	dbWriteTable(library_conn, "Samples", all_sample_info, overwrite = TRUE)
	if(nrow(all_bamqc_dna)){dbWriteTable(library_conn, "BamQC_DNA", all_bamqc_dna, overwrite = TRUE)}
	if(nrow(all_bamqc_rna)){dbWriteTable(library_conn, "BamQC_RNA", all_bamqc_rna, overwrite = TRUE)}
	if(nrow(all_variants)){dbWriteTable(library_conn, "SmallVariant_calling", all_variants, overwrite = TRUE)}
	if(nrow(all_variants_annot)){dbWriteTable(library_conn, "SmallVariant_annot", all_variants_annot, overwrite = TRUE)}
	if(nrow(all_cna_gene)){dbWriteTable(library_conn, "CNA_gene", all_cna_gene, overwrite = TRUE)}
	if(nrow(all_cna_arm)){dbWriteTable(library_conn, "CNA_arm", all_cna_arm, overwrite = TRUE)}
	if(nrow(all_fusion)){dbWriteTable(library_conn, "Fusion", all_fusion, overwrite = TRUE)}
	if(nrow(all_splicing)){dbWriteTable(library_conn, "Splicing", all_splicing, overwrite = TRUE)}
	if(nrow(all_tmb)){dbWriteTable(library_conn, "TMB", all_tmb, overwrite = TRUE)}
	if(nrow(all_msi)){dbWriteTable(library_conn, "MSI", all_msi, overwrite = TRUE)}

	dbDisconnect(library_conn)

	cat(paste0("\n", "The run variant library is created!", "\n"))
	cat(paste0("    ", basename(library.db), "\n"))
	cat(paste0("    ", nrow(all_sample_info), " samples", "\n"))
	cat(paste0("      - ", length(dna_samples), " DNA", "\n"))
	cat(paste0("      - ", length(rna_samples), " RNA", "\n"))
	
	if(nrow(all_variants_annot)){
		ok_vars <- all_variants_annot%>%filter(NO_FLAGS)
		pass_vars <- all_variants_annot%>%filter(PASS)
		nocall_vars <- all_variants_annot%>%filter(!PASS)
		cat(paste0("    ", nrow(pass_vars), " PASS small variants", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(TYPE == "SNV")), " SNVs", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(TYPE != "SNV")), " InDels", "\n"))
		cat(paste0("      - ", nrow(ok_vars), " 'OK' variants (no flags)", "\n"))
		cat(paste0("        AMP/ASCO/CAP", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier IA")), " Tier IA", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier IB")), " Tier IB", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier IIC")), " Tier IIC", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier IID")), " Tier IID", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
		cat(paste0("        ClinGen/CGC/VICC", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(ONCOGENICITY %in% "Oncogenic")), " Oncogenic", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(ONCOGENICITY %in% "Likely Oncogenic")), " Likely Oncogenic", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(ONCOGENICITY %in% "VUS")), " VUS", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(ONCOGENICITY %in% "Likely Benign")), " Likely Benign", "\n"))
		cat(paste0("          - ", nrow(ok_vars%>%filter(ONCOGENICITY %in% "Benign")), " Benign", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(GERMLINE)), " 'Germline' (>", vars_somatic_max_pvaf*100, "% pVAF or >", vars_somatic_max_vaf*100, "% VAF) variants", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(!CTR_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT))), " 'OutCTR' variants", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(PROBLEMATIC_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT))), " 'ProblematicRegion' variants", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(RECURRENT & !CANCER_HOTSPOT)), " 'Recurrent' variants", "\n"))
		cat(paste0("    ", nrow(nocall_vars), " 'NoCall' small variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(TYPE == "SNV")), " SNVs", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(TYPE != "SNV")), " InDels", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowCallers)), " 'LowCallers' (<", vars_min_callers, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowAD)), " 'LowAD' (<", vars_min_ad, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowDP)), " 'LowDP' (<", vars_min_dp, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowVAF)), " 'LowVAF' (<", vars_min_vaf*100, "%) variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(PANEL_BLACKLIST & !PANEL_HOTSPOT)), " 'Blacklist' variants", "\n"))
	} else if(nrow(all_variants)) {
		pass_vars <- all_variants%>%filter(PASS)
		nocall_vars <- all_variants%>%filter(!PASS)
		cat(paste0("    ", nrow(pass_vars), " PASS small variants (NO annotation)", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(TYPE == "SNV")), " SNVs", "\n"))
		cat(paste0("      - ", nrow(pass_vars%>%filter(TYPE != "SNV")), " InDels", "\n"))
		cat(paste0("    ", nrow(nocall_vars), " 'NoCall' small variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(TYPE == "SNV")), " SNVs", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(TYPE != "SNV")), " InDels", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowCallers)), " 'LowCallers' (<", vars_min_callers, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowAD)), " 'LowAD' (<", vars_min_ad, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowDP)), " 'LowDP' (<", vars_min_dp, ") variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(LowVAF)), " 'LowVAF' (<", vars_min_vaf*100, "%) variants", "\n"))
		cat(paste0("      - ", nrow(nocall_vars%>%filter(PANEL_BLACKLIST & !PANEL_HOTSPOT)), " 'Blacklist' variants", "\n"))
	} else {
		cat(paste0("    ", nrow(all_variants), " small variants", "\n"))
	}
	
	cat(paste0("    ", nrow(all_cna_gene), " gene CNAs", "\n"))
	if(nrow(all_cna_gene)){
		ok_cna <- all_cna_gene%>%filter(NO_FLAGS)
		cat(paste0("      - ", nrow(ok_cna), " OK CNAs (no flags)", "\n"))
		cat(paste0("        - ", nrow(ok_cna%>%filter(CNA %in% "AMP")), " AMP (CN\u2265", cna_CN_AMP_OK, ")", "\n"))
		cat(paste0("        - ", nrow(ok_cna%>%filter(CNA %in% "DEL")), " DEL (CN\u2264", cna_CN_DEL_OK, ")", "\n"))
		cat(paste0("        AMP/ASCO/CAP", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IA")), " Tier IA", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IB")), " Tier IB", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IIC")), " Tier IIC", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IID")), " Tier IID", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
		cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
		cat(paste0("      - ", nrow(all_cna_gene%>%filter(LOW_AMP)), " 'LowAmp' CNAs (2<CN<", cna_CN_AMP_OK, ")", "\n"))
		cat(paste0("      - ", nrow(all_cna_gene%>%filter(LOW_DEL)), " 'LowDel' CNAs (", cna_CN_DEL_OK, "<CN<2)", "\n"))
		cat(paste0("      - ", nrow(all_cna_gene%>%filter(LOWBINS)), " 'LowBins' (<", cna_min_target_bins, ") CNAs", "\n"))
	}
	
	cat(paste0("    ", nrow(all_cna_arm), " arm CNAs", "\n"))
	if(nrow(all_cna_arm)){
		cat(paste0("      - ", nrow(all_cna_arm%>%filter(CNA %in% "AMP")), " AMP (copy ratio\u2265", cna_arm_call_threshold, ")", "\n"))
		cat(paste0("      - ", nrow(all_cna_arm%>%filter(CNA %in% "DEL")), " DEL (copy ratio\u2264", -cna_arm_call_threshold, ")", "\n"))
	}
	
	if(nrow(all_fusion)){
		ok_fusion <- all_fusion%>%filter(NO_FLAGS)
		pass_fusion <- all_fusion%>%filter(PASS)
		nocall_fusion <- all_fusion%>%filter(!PASS)
		cat(paste0("    ", nrow(pass_fusion), " PASS fusions", "\n"))
		cat(paste0("      - ", nrow(ok_fusion), " OK fusions (no flags)", "\n"))
		cat(paste0("        AMP/ASCO/CAP guidelines", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier IA")), " Tier IA", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier IB")), " Tier IB", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier IIC")), " Tier IIC", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier IID")), " Tier IID", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
		cat(paste0("          - ", nrow(ok_fusion%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
		cat(paste0("      - ", nrow(pass_fusion%>%filter(LowAD)), " 'LowAD' (<", fusion_flag_min_ad, ") fusions", "\n"))
		cat(paste0("      - ", nrow(pass_fusion%>%filter(LowNonFused)), " 'LowNonFused' (<", fusion_flag_min_ad_nonfused, ") fusions", "\n"))
		cat(paste0("      - ", nrow(pass_fusion%>%filter(LowDP)), " 'LowDP' (<", fusion_flag_min_dp, ") fusions", "\n"))
		cat(paste0("      - ", nrow(pass_fusion%>%filter(LowVAF)), " 'LowVAF' (<", fusion_flag_min_vaf*100, "%) fusions", "\n"))
		cat(paste0("      - ", nrow(pass_fusion%>%filter(UNKNOWN)), " 'Unknown' fusions", "\n"))
		cat(paste0("    ", nrow(nocall_fusion), " 'NoCall' fusions", "\n"))
		cat(paste0("      - ", nrow(nocall_fusion%>%filter(LowSupport)), " 'LowSupport' (AD<", fusion_calling_min_ad, " or FFPM<", fusion_calling_min_ffpm, ") fusions", "\n"))
		cat(paste0("      - ", nrow(nocall_fusion%>%filter(Normal)), " 'Normal' fusions", "\n"))
		cat(paste0("      - ", nrow(nocall_fusion%>%filter(RTartifact)), " 'RTartifact' fusions", "\n"))
		cat(paste0("      - ", nrow(nocall_fusion%>%filter(!InSilicoValid)), " 'NoInSilicoValid' fusions", "\n"))
	} else {
		cat(paste0("    ", nrow(all_fusion), " gene fusions", "\n"))
	}

	if(nrow(all_splicing)){
		ok_splicing <- all_splicing%>%filter(NO_FLAGS)
		pass_splicing <- all_splicing%>%filter(PASS)
		nocall_splicing <- all_splicing%>%filter(!PASS)
		cat(paste0("    ", nrow(pass_splicing), " PASS splicing variants", "\n"))
		cat(paste0("      - ", nrow(ok_splicing), " OK variants (no flags)", "\n"))
		cat(paste0("        AMP/ASCO/CAP guidelines", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier IA")), " Tier IA", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier IB")), " Tier IB", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier IIC")), " Tier IIC", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier IID")), " Tier IID", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
		cat(paste0("          - ", nrow(ok_splicing%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
		cat(paste0("      - ", nrow(pass_splicing%>%filter(!CancerEnriched & !WHITELIST_SPLICING)), " 'NoCancerEnriched' variants", "\n"))
		cat(paste0("      - ", nrow(pass_splicing%>%filter(LowAD)), " 'LowAD' (<", splicing_flag_min_ad, ") variants", "\n"))
		cat(paste0("      - ", nrow(pass_splicing%>%filter(LowDP)), " 'LowDP' (<", splicing_flag_min_dp, ") variants", "\n"))
		cat(paste0("      - ", nrow(pass_splicing%>%filter(LowVAF)), " 'LowVAF' (<", splicing_flag_min_vaf*100, "%) variants", "\n"))
		cat(paste0("    ", nrow(nocall_splicing), " 'NoCall' variants", "\n"))
		cat(paste0("      - ", nrow(nocall_splicing%>%filter(LowSupport)), " 'LowSupport' (AD<", splicing_calling_min_ad, ") variants", "\n"))
	} else {
		cat(paste0("    ", nrow(all_splicing), " splicing variants", "\n"))
	}

}
