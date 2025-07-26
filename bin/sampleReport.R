#!/usr/bin/env Rscript

# ============================================================================ #
# What: Collect sample results and create a report
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(knitr)))
suppressPackageStartupMessages(suppressWarnings(library(flexdashboard)))
suppressPackageStartupMessages(suppressWarnings(library(kableExtra)))
suppressPackageStartupMessages(suppressWarnings(library(DT)))
suppressPackageStartupMessages(suppressWarnings(library(crosstalk)))
suppressPackageStartupMessages(suppressWarnings(library(DBI)))
suppressPackageStartupMessages(suppressWarnings(library(scales)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
run <- args_list[["run"]]
sample <- args_list[["sample"]]
panel <- args_list[["panel"]]
is_iontorrent <- tolower(args_list[["platform"]]) %in% "iontorrent"
dna_sample.file <- args_list[["dnaSampleFile"]]
rna_sample.file <- args_list[["rnaSampleFile"]]
is_panel_hotspots <- file.exists(args_list[["panelHotspots"]])
civic.file <- args_list[["civicEvidence"]]
library.db <- args_list[["libraryDb"]]
vars_min_callers <- as.numeric(args_list[["variantMinCallers"]])
vars_min_ad <- as.numeric(args_list[["variantMinAD"]])
vars_min_dp <- as.numeric(args_list[["variantMinDP"]])
vars_min_vaf <- as.numeric(args_list[["variantMinVAF"]])
vars_somatic_max_pvaf <- as.numeric(args_list[["variantSomaticMaxpVAF"]])
vars_somatic_max_vaf <- as.numeric(args_list[["variantSomaticMaxVAF"]])
cna_min_target_bins <- as.numeric(args_list[["cnaMinTargetBins"]])
cna_CN_AMP_OK <- as.numeric(args_list[["cnaHighAmpCN"]])
cna_CN_DEL_OK <- as.numeric(args_list[["cnaHighDelCN"]])
cna_plot_bygene_all <- as.logical(args_list[["cnaPlotByGeneAll"]])
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
bamqc_cov_byexon_all_dna <- as.logical(args_list[["covByExonAllGenesDna"]])
bamqc_cov_byexon_all_rna <- as.logical(args_list[["covByExonAllGenesRna"]])
bamqc_plot_bygene_all_dna <- as.logical(args_list[["bamqcPlotByGeneAllDna"]])
bamqc_plot_bygene_all_rna <- as.logical(args_list[["bamqcPlotByGeneAllRna"]])

# Set QC arguments
qc_parameters <- c(
	"TOTAL_READS" = "totalReads", "ALIGNED_READS" = "alignedReads", "PCT_ALIGNED" = "alignedReadsPct",
	"ONTARGET_READS" = "ontargetReads", "PCT_ONTARGET" = "ontargetReadsPct",
	"HQ_ALIGNED_READS" = "hqAlignedReads", "PCT_HQ_ALIGNED" = "hqAlignedReadsPct",
	"MEDIAN_READ_LENGTH" = "medianReadLength", "MEDIAN_INSERT_SIZE" = "medianInsertSize",
	"UNIQUE_READS" = "uniqueReads", "PCT_DUP" = "duplicateReadsPct",
	"MEDIAN_COVERAGE" = "medianCoverage", "MEAN_COVERAGE" = "meanCoverage",
	"PCT_0.4X_MEAN" = "pct04xMeanCoverage", "PCT_100X" = "pct100xCoverage", "PCT_1000X" = "pct1000xCoverage"
)
dna_thresholds <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
for(param in qc_parameters){
	param_id <- paste0(param, "Dna")
	value <- args_list[[param_id]]
	aux_df <- data.frame(
		NAME = param, LOWER = ifelse(identical(value, param_id), NA, as.numeric(unlist(strsplit(value, "-"))[1])),
		UPPER = ifelse(identical(value, param_id), NA, as.numeric(unlist(strsplit(value, "-"))[2]))) %>%
		mutate(OUT = is.na(LOWER) | is.na(UPPER))
	dna_thresholds <- rbind(dna_thresholds, aux_df)
	rm(value, aux_df)
}
rm(param)
rna_thresholds <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
for(param in qc_parameters){
	param_id <- paste0(param, "Rna")
	value <- args_list[[paste0(param, "Rna")]]
	aux_df <- data.frame(
		NAME = param, LOWER = ifelse(identical(value, param_id), NA, as.numeric(unlist(strsplit(value, "-"))[1])),
		UPPER = ifelse(identical(value, param_id), NA, as.numeric(unlist(strsplit(value, "-"))[2]))) %>%
		mutate(OUT = is.na(LOWER) | is.na(UPPER))
	rna_thresholds <- rbind(rna_thresholds, aux_df)
	rm(value, aux_df)
}
rm(param)

# Set sample results
variants.file <- list.files(pattern = paste0(sample, "_DNA_SmallVariant_annot.txt"))
cna_gene.file <- list.files(pattern = paste0(sample, "_DNA_CNA_gene.txt"))
cna_arm.file <- list.files(pattern = paste0(sample, "_DNA_CNA_arm.txt"))
cna_sex.file <- list.files(pattern = paste0(sample, "_DNA_sex.txt"))
fusion.file <- list.files(pattern = paste0(sample, "_RNA_Fusion_allCandidates.txt"))
splicing.file <- list.files(pattern = paste0(sample, "_RNA_Splicing_allCandidates.txt"))
tmb.file <- list.files(pattern = paste0(sample, "_DNA_TMB_results.txt"))
tmb_trace.file <- list.files(pattern = paste0(sample, "_DNA_TMB_trace.txt"))
msi.file <- list.files(pattern = paste0(sample, "_DNA_MSI_results.txt"))
dna_metrics.file <- list.files(pattern = paste0(sample, "_DNA_BamQC.txt"))
dna_genecov.file <- list.files(pattern = paste0(sample, "_DNA_CoverageByGene.txt"))
dna_exoncov.file <- list.files(pattern = paste0(sample, "_DNA_CoverageByExon.txt"))
rna_metrics.file <- list.files(pattern = paste0(sample, "_RNA_BamQC.txt"))
rna_genecov.file <- list.files(pattern = paste0(sample, "_RNA_CoverageByGene.txt"))
rna_exoncov.file <- list.files(pattern = paste0(sample, "_RNA_CoverageByExon.txt"))

# Set other variables
is_vars <- as.logical(length(variants.file))
is_cna_gene <- as.logical(length(cna_gene.file))
is_cna_arm <- as.logical(length(cna_arm.file))
is_fusion <- as.logical(length(fusion.file))
is_splicing <- as.logical(length(splicing.file))
is_tmb <- as.logical(length(tmb.file))
is_msi <- as.logical(length(msi.file))
is_dna_metrics <- as.logical(length(dna_metrics.file))
is_rna_metrics <- as.logical(length(rna_metrics.file))
all_colors <- c(
	"primary" = "#2780E3", "ND" = "grey", "OK" = "mediumseagreen", "FAIL" = "lightcoral", "CLIN_RELEVANT" = "indigo",
	"Tier IA" = "indigo", "Tier IB" = "indigo", "Tier IIC" = "mediumpurple", "Tier IID" = "mediumpurple", 
	"Tier III" = "dodgerblue", "Tier IV" = "mediumseagreen",
	"SNV" = "blue", "INDEL" = "dodgerblue",
	"Oncogenic" = "red", "Likely Oncogenic" = "darkorange", "VUS" = "#FFDB6D", "Likely Benign" = "lightgreen", "Benign" = "mediumseagreen",
	"AMP" = "blue", "DEL" = "red", "NEUTRAL" = "mediumseagreen", "CONFLICT" = "#FFDB6D",
	"HighAmp" = "blue", "HighDel" = "red", "LowAmp" = "#00BFFF", "LowDel" = "lightcoral", "Neutral" = "mediumseagreen",
	"QC_PASS" = "mediumseagreen", "QC_WARN" = "darkorange", "QC_FAIL" = "red"
)
clinical_levels <- c("Tier IA", "Tier IB", "Tier IIC", "Tier IID", "Tier III", "Tier IV")
vars_chr_levels <- paste0("chr", c(1:22, "X", "Y"))
vars_onco_levels <- c("Oncogenic", "Likely Oncogenic", "VUS", "Likely Benign", "Benign")
vars_impact_levels <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
cna_chr_levels <- paste0("chr", c(1:22, "X"))
cna_CNA_levels <- c("AMP", "DEL")
cna_class_levels <- c("HighAmp", "HighDel", "LowAmp", "LowDel")
cna_arm_levels <- c("AMP", "DEL", "NEUTRAL", "CONFLICT")
fusion_type_levels <- c("INFRAME", "FRAMESHIFT")
qc_flags <- c("QC_PASS", "QC_WARN", "QC_FAIL")


# Load data ---------------------------------------------------------------

# Sample info
sample_info <- matrix(data = NA, nrow = 0, ncol = 1)
if(file.exists(dna_sample.file)){
	dna_samples <- read.table(dna_sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE)
	sample_info <- rbind(sample_info, dna_samples)
	rm(dna_samples)
}
if(file.exists(rna_sample.file)){
	rna_samples <- read.table(rna_sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE)
	sample_info <- rbind(sample_info, rna_samples)
	rm(rna_samples)
}
sample_info <- sample_info %>% mutate(Sample = as.character(Sample)) %>% distinct(Sample, .keep_all = TRUE) %>% filter(Sample %in% sample)
sex <- sample_info$Sex
if(sex == "."){
	sex_predict <- ifelse(length(cna_sex.file), read.table(cna_sex.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)$sex, ".")
}
age <- sample_info$Age
purity <- sample_info$Purity
whitelist <- sample_info$Whitelist
tumor <- sample_info$TumorCode
doid <- sample_info$TumorDOID

# Clinical data (CiVIC)
civic_data <- read.table(civic.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE) %>%
	mutate(LEVEL = EVIDENCE_LEVEL, RATING = EVIDENCE_RATING, TYPE = EVIDENCE_TYPE, TUMOR = TUMOR_NAME, SCORE = EVIDENCE_SCORE) %>%
	select(TYPE, ALTERATION, MOLECULAR_ID, LEVEL, RATING, EFFECT, DRUG, TUMOR, SCORE, ORIGIN)
civic_predictive <- civic_data %>% filter(TYPE %in% "Predictive") %>% select(-TYPE)
civic_prognostic <- civic_data %>% filter(TYPE %in% "Prognostic") %>% select(-c(TYPE, DRUG))
civic_diagnostic <- civic_data %>% filter(TYPE %in% "Diagnostic") %>% select(-c(TYPE, DRUG))

# Variant library
if(file.exists(library.db)){
	library_conn <- dbConnect(RSQLite::SQLite(), library.db)
	library_tables <- dbListTables(library_conn)
	
	if("Samples" %in% library_tables){
		samples_lib <- tbl(library_conn, "Samples") %>% collect()
		vars_lib_N <- nrow(samples_lib %>% filter(as.logical(SMALL_VARIANT)))
		cna_lib_N <- nrow(samples_lib %>% filter(as.logical(CNA)))
		fusion_lib_N <- nrow(samples_lib %>% filter(as.logical(FUSION)))
		splicing_lib_N <- nrow(samples_lib %>% filter(as.logical(SPLICING)))
	} else {
		vars_lib_N <- cna_lib_N <- fusion_lib_N <- splicing_lib_N <- 0
	}

	if("SmallVariant_calling" %in% library_tables){
		vars_lib <- tbl(library_conn, "SmallVariant_calling") %>% select(VAR, AC_SAMPLES, AF_SAMPLES) %>% collect()
	} else {
		vars_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
	}
	
	if("CNA_gene" %in% library_tables){
		cna_gene_lib <- tbl(library_conn, "CNA_gene") %>% select(VAR, AC_SAMPLES, AF_SAMPLES) %>% collect()
	} else {
		cna_gene_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
	}
	
	if("CNA_arm" %in% library_tables){
		cna_arm_lib <- tbl(library_conn, "CNA_arm") %>% select(VAR, AC_SAMPLES, AF_SAMPLES) %>% collect()
	} else {
		cna_arm_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
	}
	
	if("Fusion" %in% library_tables){
		fusion_lib <- tbl(library_conn, "Fusion") %>% select(VAR, AC_SAMPLES, AF_SAMPLES) %>% collect()
	} else {
		fusion_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
	}
	
	if("Splicing" %in% library_tables){
		splicing_lib <- tbl(library_conn, "Splicing") %>% select(VAR, AC_SAMPLES, AF_SAMPLES) %>% collect()
	} else {
		splicing_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
	}
	
	dbDisconnect(library_conn)
} else {
	vars_lib_N <- cna_lib_N <- fusion_lib_N <- splicing_lib_N <- 0
	vars_lib <- cna_gene_lib <- cna_arm_lib <- fusion_lib <- splicing_lib <- data.frame(VAR = NA, AC_SAMPLES = NA, AF_SAMPLES = NA)
}

# Small variants
#	- PASS or relevant variants are considered for reporting
#	- Tier-based reporting: Tier I-III variants (tier IV or benign/likely benign variants are NOT reported)
#	- Tier I, Tier II and Oncogenic variants always considered
# 	- OK variants (If there is a panel hotspot focus on that)
if(is_vars){
	all_variants <- read.table(variants.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

	if(nrow(all_variants)){
		
		all_variants <- all_variants %>% mutate(LABEL = gsub("\\*", "\\\\*", MUT_ID))
		variants <- all_variants %>% filter(PASS | CLIN_RELEVANT | ONCOGENICITY %in% "Oncogenic")

		nocall_vars <- all_variants %>% filter(!PASS)
		pass_vars <- all_variants %>% filter(PASS)
		ok_vars <- all_variants %>% filter(NO_FLAGS)
		main_vars <- ok_vars %>% filter(!CLIN_STATUS %in% "Tier IV" & ((!is_panel_hotspots) | (is_panel_hotspots & PANEL_HOTSPOT)))
		somatic_vars <- pass_vars %>% filter(!GERMLINE)
		germline_vars <- pass_vars %>% filter(GERMLINE)
		tier1_vars <- pass_vars %>% filter(CLIN_STATUS %in% c("Tier IA", "Tier IB"))
		tier2_vars <- pass_vars %>% filter(CLIN_STATUS %in% c("Tier IIC", "Tier IID"))
		tier3_vars <- pass_vars %>% filter(CLIN_STATUS %in% "Tier III")
		tier4_vars <- pass_vars %>% filter(CLIN_STATUS %in% "Tier IV")
		
		variants_plot <- ok_vars %>% filter(CLIN_RELEVANT | ONCOGENICITY %in% c("Oncogenic", "Likely Oncogenic"))

		variants_table <- variants %>% filter(!CLIN_STATUS %in% "Tier IV") %>%
			left_join(vars_lib, by = "VAR") %>%
			mutate(
				GENE = GENE_SYMBOL, FLAGS = gsub(";", ", ", FLAGS), CALLERS = gsub(";", ", ", CALLERS), 
				CONSEQUENCE = gsub("_", " ", CONSEQUENCE), EXISTING_VARIATION = gsub("&", ", ", EXISTING_VARIATION),
				ClinVar_CLNSIG = gsub("_", " ", gsub("/", " / ", ClinVar_CLNSIG)), 
				BIOTYPE = gsub("_", " ", BIOTYPE), 
				AC_SAMPLES = ifelse(is.na(AC_SAMPLES), paste0(0, "/", vars_lib_N), AC_SAMPLES),
				AF_SAMPLES = ifelse(is.na(AF_SAMPLES), 0, as.numeric(AF_SAMPLES)),
				N_REF = nchar(REF), N_ALT = nchar(ALT), INDEL_1BP = abs(N_REF - N_ALT) == 1,
				CLIN_STATUS = factor(CLIN_STATUS, levels = clinical_levels),
				ONCOGENICITY = factor(ONCOGENICITY, levels = vars_onco_levels),
				IMPACT = factor(IMPACT, levels = vars_impact_levels), CHROM = factor(CHROM, levels = vars_chr_levels),
				across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
			select(
				FLAGS, CLASS, GENE, MUTATION, AD_ALT, AF, CLIN_STATUS, ONCOGENICITY, FILTER, ONCO_POINTS, ONCO_CODES, 
				CONSEQUENCE, IMPACT, EXON, INTRON, DP, AD_REF, AC_SAMPLES, AF_SAMPLES, OVERLAP_CALLERS, MATCH_CALLERS, CALLERS,
				HGVSg, HGVSc_ENSEMBL, HGVSc_REFSEQ, HGVSp_ENSEMBL, dbSNP_ID, ClinVar_ID, EXISTING_VARIATION, gnomAD_MAX_AF, GENIE_CNT, 
				HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, ClinVar_CLNSIG, ClinVar_ORIGIN,
				CADD_TERM, CADD_SCORE, REVEL_TERM, REVEL_SCORE, AlphaMissense_TERM, AlphaMissense_SCORE, 
				GENE_ENSEMBL, GENE_HGNC, BIOTYPE, APPRIS, TSL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, HGVSc, cDNA_POS, CDS_POS, CCDS, 
				PROTEIN_ENSEMBL, HGVSp, HGVSp_SHORT, PROTEIN_POS, AA, CODONS, 
				
				CLIN_RELEVANT, NMD_ESCAPE, GENE_WHITELIST, PANEL_HOTSPOT, MMR_GENE, SOMATIC_WHITELIST, 
				ClinVar_PATHOGENIC, PREDICTOR_PATHOGENIC, PREDICTOR_BENIGN,
				ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL, CODING, NONCODING, INDEL_1BP,
				PASS, LowCallers, LowAD, LowDP, LowVAF, PANEL_BLACKLIST,
				NO_FLAGS, CANCER_HOTSPOT, GERMLINE, CTR_REGION, PROBLEMATIC_REGION, RECURRENT,

				MUT_ID, CHROM, START, END, REF, ALT, STRAND, VAR, VAR_HG19, VAR_MUTECT2, VAR_OCTOPUS, VAR_PISCES, VAR_VARDICT, VAR_TVC)
		
		somatic_table <- variants_table %>% filter(GERMLINE == "FALSE")
		germline_table <- variants_table %>% filter(GERMLINE == "TRUE")
	
	} else {
		variants <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
		is_vars <- FALSE
	}

} else {
	variants <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# CNAs
if(is_cna_gene){
	cna_all <- read.table(cna_gene.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

	if(nrow(cna_all)){
		cna_all <- cna_all %>%
			mutate(
				TIER_I = CLIN_STATUS %in% c("Tier IA", "Tier IB"), TIER_II = CLIN_STATUS %in% c("Tier IIC", "Tier IID"),
				TIER_III = CLIN_STATUS %in% "Tier III", TIER_IV = CLIN_STATUS %in% "Tier IV", LABEL = gsub("_", " ", ALTERATION))

		cna <- cna_all %>% filter(CNA %in% c("AMP", "DEL"))
		amp <- cna %>% filter(CNA %in% "AMP") %>% arrange(GENE)
		del <- cna %>% filter(CNA %in% "DEL") %>% arrange(GENE)

		cna_gene_table <- cna %>% left_join(cna_gene_lib, by = "VAR") %>%
			mutate(
				AC_SAMPLES = ifelse(is.na(AC_SAMPLES), paste0(0, "/", cna_lib_N), AC_SAMPLES),
				AF_SAMPLES = ifelse(is.na(AF_SAMPLES), 0, as.numeric(AF_SAMPLES)),
				FLAGS = gsub(";", ", ", FLAGS), CNA = factor(CNA, levels = cna_CNA_levels),
				CLIN_STATUS = factor(CLIN_STATUS, levels = clinical_levels),
				CLASS = factor(CLASS, levels = cna_class_levels), CHROM = factor(CHROM, levels = cna_chr_levels),
				across(c(CN, LOG2, DEPTH, WEIGHT), ~ as.numeric(.x)),
				PERC_COV_GENE_ALL = as.numeric(ifelse(PERC_COV_GENE_ALL != ".", PERC_COV_GENE_ALL, NA)),
				PERC_COV_GENE_EXON = as.numeric(ifelse(PERC_COV_GENE_EXON != ".", PERC_COV_GENE_EXON, NA)),
				across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
			select(
				FLAGS, CYTOBAND, GENE, CNA, CLASS, CLIN_STATUS, CN, LOG2, DEPTH, AC_SAMPLES, AF_SAMPLES, GENIE_CNT, GENIE_FREQ, 
				BINS_TOTAL, BINS_TARGET, WEIGHT, GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ,
				CLIN_RELEVANT, WHITELIST, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, NO_FLAGS, LOW_AMP, LOW_DEL, LOWBINS, NEUTRAL,
				ALTERATION, VAR, CHROM, START, END, STRAND, PERC_COV_GENE_ALL, PERC_COV_GENE_EXON)
	} else {
		cna <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
		is_cna_gene <- FALSE
	}
	
} else {
	cna <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

if(is_cna_arm){
	cna_arm_all <- read.table(cna_arm.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

	if(nrow(cna_arm_all)){
		cna_arm_all <- cna_arm_all %>% filter(!CNA %in% ".")
		cna_arm <- cna_arm_all %>% filter(CNA %in% c("AMP", "DEL"))

		cna_arm_table <- cna_arm_all %>% left_join(cna_arm_lib, by = "VAR") %>%
			mutate(
				AC_SAMPLES = ifelse(is.na(AC_SAMPLES), paste0(0, "/", cna_lib_N), AC_SAMPLES),
				AF_SAMPLES = ifelse(is.na(AF_SAMPLES), 0, as.numeric(AF_SAMPLES)),
				CHROM = factor(CHROM, levels = cna_chr_levels), CNA = factor(CNA, levels = cna_arm_levels), LOG2 = as.numeric(LOG2)) %>%
			select(CHR_ARM, CNA, LOG2, CHROM, SIDE, START, END, AC_SAMPLES, AF_SAMPLES, VAR)
	} else {
		cna_arm_all <- cna_arm <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
		is_cna_arm <- FALSE
	}
	
} else {
	cna_arm_all <- cna_arm <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# Gene fusions
#	- Whitelist fusions always considered for plotting
#	- Tier IV (unknown) fusions are NOT visualized
if(is_fusion){
	fusion <- read.table(fusion.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

	if(nrow(fusion)){
		fusion <- fusion %>%
			mutate(
				TIER_I = CLIN_STATUS %in% c("Tier IA", "Tier IB"), TIER_II = CLIN_STATUS %in% c("Tier IIC", "Tier IID"),
				TIER_III = CLIN_STATUS %in% "Tier III", TIER_IV = CLIN_STATUS %in% "Tier IV",
				LABEL = ifelse(FUSION_VARIANT != ".", FUSION_VARIANT, FUSION_SHORT))
		fusion_called <- fusion %>% filter(PASS) %>% mutate(VAR = as.character(VAR))
		fusion_plot <- fusion %>% filter(WHITELIST_FUSION | (PASS & !TIER_IV))
		fusion_table <- fusion_called %>% left_join(fusion_lib, by = "VAR") %>%
			mutate(
				FUSION = FUSION_SHORT, VARIANT = FUSION_VARIANT, TYPE = MODEL_FUSION_TYPE,
				AC_SAMPLES = ifelse(is.na(AC_SAMPLES), paste0(0, "/", fusion_lib_N), AC_SAMPLES),
				AF_SAMPLES = ifelse(is.na(AF_SAMPLES), 0, as.numeric(AF_SAMPLES)),
				FLAGS = gsub(";", ", ", FLAGS), RESOURCES = gsub(",", ", ", RESOURCES), 
				AF = as.numeric(AF), DP = as.numeric(DP), 
				COV_IN_A = as.numeric(ifelse(COV_IN_FUSION_A != ".", COV_IN_FUSION_A, NA)),
				COV_OUT_A = as.numeric(ifelse(COV_OUT_FUSION_A != ".", COV_OUT_FUSION_A, NA)),
				COV_IN_B = as.numeric(ifelse(COV_IN_FUSION_B != ".", COV_IN_FUSION_B, NA)),
				COV_OUT_B = as.numeric(ifelse(COV_OUT_FUSION_B != ".", COV_OUT_FUSION_B, NA)),
				MODEL_FUSION_TYPE = factor(MODEL_FUSION_TYPE, levels = fusion_type_levels),
				CLIN_STATUS = factor(CLIN_STATUS, levels = clinical_levels),
				across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
			select(
				FLAGS, FUSION, FUSION_NAME, VARIANT, AD, AF, CLIN_STATUS, CALLING, TYPE, DP, FFPM,
				COV_IN_A, COV_OUT_A, COV_IN_B, COV_OUT_B, AC_SAMPLES, AF_SAMPLES, RESOURCES, 
				MitelmanDB_COUNT, GENIE_CNT, LeftBreakDinuc, RightBreakDinuc, 
				JunctionReadCount, SpanningFragCount, AD_NonFused_A, AD_NonFused_B, 
				FUSION_RANGE_A, MODEL_FUSION_A, MODEL_CDS_A, MODEL_CDS_RANGE_A, 
				FUSION_RANGE_B, MODEL_FUSION_B, MODEL_CDS_B, MODEL_CDS_RANGE_B, 
				GENE_A, GENE_ENSEMBL_A, TRANSCRIPT_ENSEMBL_A, TRANSCRIPT_REFSEQ_A, EXON_A, INTRON_A, BASES_FROM_EXON_A,
				GENE_B, GENE_ENSEMBL_B, TRANSCRIPT_ENSEMBL_B, TRANSCRIPT_REFSEQ_B, EXON_B, INTRON_B, BASES_FROM_EXON_B, 

				CLIN_RELEVANT, WHITELIST_FUSION, WHITELIST_GENE, LargeAnchorSupport, RefSpliceSite, 
				ONCOGENE_A, TSG_A, DRIVER_A, CANONICAL_DRIVER_A, ONCOGENE_B, TSG_B, DRIVER_B, CANONICAL_DRIVER_B, 
				PASS, LowSupport, Normal, RTartifact, InSilicoValid, 
				NO_FLAGS, LowAD, LowNonFused, LowDP, LowVAF, UNKNOWN, MitelmanDB, GENIE,
				
				VAR, BREAKPOINT_A, CHROM_A, POS_A, STRAND_A, BREAKPOINT_B, CHROM_B, POS_B, STRAND_B, 
				VAR_HG19, BREAKPOINT_A_HG19, BREAKPOINT_B_HG19)
	} else {
		fusion <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
		is_fusion <- FALSE
	}
	
	
} else {
	fusion <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# Splicing variants
#	- Whitelist splicing always considered for plotting
if(is_splicing){
	splicing <- read.table(splicing.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

	if(nrow(splicing)){
		splicing <- splicing %>% rowwise() %>%
			mutate(
				TIER_I = CLIN_STATUS %in% c("Tier IA", "Tier IB"), TIER_II = CLIN_STATUS %in% c("Tier IIC", "Tier IID"),
				TIER_III = CLIN_STATUS %in% "Tier III", TIER_IV = CLIN_STATUS %in% "Tier IV",
				LABEL = ifelse(VAR_NAME != ".", VAR_NAME, VAR_GENE), VAR_GENE = gsub("_", " ", VAR_GENE)) %>% 
			as.data.frame()
		splicing_called <- splicing %>% filter(PASS)
		splicing_plot <- splicing_called %>% filter(WHITELIST_SPLICING | VAR_NAME != "." | (NO_FLAGS & CancerEnriched))
		splicing_table <- splicing_called %>% filter(WHITELIST_SPLICING | VAR_NAME != "." | CancerEnriched) %>%
			left_join(splicing_lib, by = "VAR") %>%
			mutate(
				VARIANT = VAR_NAME, REGION = REGION_AFFECTED,
				AC_SAMPLES = ifelse(is.na(AC_SAMPLES), paste0(0, "/", splicing_lib_N), AC_SAMPLES),
				AF_SAMPLES = ifelse(is.na(AF_SAMPLES), 0, as.numeric(AF_SAMPLES)),
				CALLING = gsub(";", ", ", CALLING), FLAGS = gsub(";", ", ", FLAGS),
				CHROM = factor(CHROM, levels = vars_chr_levels), MUTATED_SITE = MUTATION != ".",
				COV_EXONS_AFFECTED = as.numeric(ifelse(COV_EXONS_AFFECTED != ".", COV_EXONS_AFFECTED, NA)),
				COV_EXONS_FLANKING = as.numeric(ifelse(COV_EXONS_FLANKING != ".", COV_EXONS_FLANKING, NA)),
				CLIN_STATUS = factor(CLIN_STATUS, levels = clinical_levels),
				across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
			select(
				FLAGS, GENE, REGION, VARIANT, AD, AF, CLIN_STATUS, CALLING, MUTATION, 
				DP_MAX, COV_EXONS_AFFECTED, COV_EXONS_FLANKING, AC_SAMPLES, AF_SAMPLES, UNIQ_MAPPED_READS, 
				MULTI_MAPPED_READS, DP_START, DP_END, DP_MEAN, 
				GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, EXONS_AFFECTED, INTRONS_AFFECTED, 

				CLIN_RELEVANT, WHITELIST_SPLICING, WHITELIST_GENE, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, 
				MUTATED_SITE, PASS, LowSupport, NO_FLAGS, CancerEnriched, LowAD, LowDP, LowVAF,

				VAR_GENE, VAR, CHROM, START, END, STRAND, VAR_HG19, START_HG19, END_HG19)
	} else {
		splicing <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
		is_splicing <- FALSE
	}
	
} else {
	splicing <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# TMB result
if(is_tmb){
	tmb <- read.table(tmb.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
	tmb_trace <- read.table(tmb_trace.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% filter(IncludedInTmbNumerator) %>%
		mutate(across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
		select(
			VAR, VAR_HG19, GeneName, Mutation, VAF, TotalDepth, AltDepth, ClinStatus, Oncogenicity, 
			Consequence, dbSNP, MaxpVAF, IsGermline:IncludedInTmbNumerator,
			VariantType, VariantClass, Chromosome, Position, RefCall, AltCall)
} else {
	tmb <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# MSI result
if(is_msi){
	msi <- read.table(msi.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
		mutate(
			TotalSites = Total_Number_of_Sites, MsiSites = Number_of_Somatic_Sites, 
			MSI_PCT = round(Number_of_Somatic_Sites * 100 / Total_Number_of_Sites, 1)) %>% 
		select(MSI_PCT, MsiSites, TotalSites)
	if(is_vars){
		msi_table <- variants_table %>% filter(as.logical(MMR_GENE))
	} else {
		msi_table <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
	}
} else {
	msi <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
}

# DNA metrics
if(is_dna_metrics){
	dna_metrics <- read.table(dna_metrics.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% mutate(
		PCT_ALIGNED = PCT_ALIGNED_READS, PCT_ONTARGET = PCT_ONTARGET_ALIGNED_READS,
		PCT_HQ_ALIGNED = PCT_HQ_ALIGNED_READS, PCT_DUP = PCT_DUPLICATION)
	
	for(param in names(qc_parameters)){
		aux_df <- dna_thresholds %>% filter(NAME %in% qc_parameters[param])
		value <- dna_metrics[[param]]
		flag <- ifelse(aux_df$OUT, ".", ifelse(value > aux_df$UPPER, "QC_PASS", ifelse(value < aux_df$LOWER, "QC_FAIL", "QC_WARN")))
		dna_metrics[[paste0(qc_parameters[param], "_FLAG")]] <- flag
	}
	rm(param)
	
	dna_metrics <- dna_metrics %>% 
		mutate(
			across(matches("_READS$", ignore.case = FALSE),  ~ label_number(scale_cut = cut_short_scale())(.x)),
			PANEL_SIZE = label_number(accuracy = 0.01, suffix = "b", scale_cut = cut_short_scale())(PANEL_SIZE),
			across(matches("_COVERAGE$", ignore.case = FALSE), ~ format(.x, big.mark = ",")),
			across(matches("PCT_\\w+", ignore.case = FALSE), ~ paste0(.x, "%"))) %>%
		select(
			TOTAL_READS, ALIGNED_READS, PCT_ALIGNED, ONTARGET_READS, PCT_ONTARGET, HQ_ALIGNED_READS, PCT_HQ_ALIGNED, 
			MEDIAN_READ_LENGTH, MEDIAN_INSERT_SIZE, UNIQUE_READS, PCT_DUP, MEDIAN_COVERAGE, MEAN_COVERAGE,
			matches("PCT_\\d+", ignore.case = FALSE), MIN_COVERAGE, MAX_COVERAGE, PANEL_SIZE, contains("_FLAG")) %>%
		select(as.character(names(qc_parameters[!dna_thresholds$OUT])), everything())
	
	aux_metrics <- as.character(dna_metrics[1,])
	dna_qc <- ifelse("QC_FAIL" %in% aux_metrics, "QC_FAIL", ifelse("QC_WARN" %in% aux_metrics, "QC_WARN", "QC_PASS"))
	rm(aux_metrics)
	
	dna_genecov <- read.table(dna_genecov.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
		mutate(
			TARGET = DP_TARGET, CODING = DP_MANE_CODING, EXON = DP_MANE_EXON, WHOLE = DP_MANE_ALL, MIN = MIN_COVERAGE, 
			MAX = MAX_COVERAGE, across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
		arrange(WHITELIST) %>%
		select(
			GENE, TARGET, CODING, EXON, WHOLE, MIN, MAX, contains("PCT_"), PANEL_SIZE, 
			matches("_ENSEMBL|_REFSEQ"), WHITELIST, LowCoverage)
	dna_genecov_whitelist <- dna_genecov %>% filter(as.logical(WHITELIST))

	if(bamqc_plot_bygene_all_dna){dna_genes <- dna_genecov$GENE} else {dna_genes <- dna_genecov_whitelist$GENE}
	
	dna_exoncov <- read.table(dna_exoncov.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
		mutate(MIN = MIN_COVERAGE, MAX = MAX_COVERAGE) %>%
		select(GENE, EXON, DEPTH, MIN, MAX, LENGTH, contains("PCT_"), CHROM, START, END, contains("TRANSCRIPT_"), WHITELIST, LowCoverage) %>%
		mutate(across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
		filter(bamqc_cov_byexon_all_dna | as.logical(WHITELIST))
	
} else {
	dna_metrics <- dna_exoncov <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
	dna_qc <- "ND"
}

# RNA metrics
if(is_rna_metrics){
	rna_metrics <- read.table(rna_metrics.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% mutate(
		PCT_ALIGNED = PCT_ALIGNED_READS, PCT_ONTARGET = PCT_ONTARGET_ALIGNED_READS,
		PCT_HQ_ALIGNED = PCT_HQ_ALIGNED_READS, PCT_DUP = PCT_DUPLICATION)
	
	for(param in names(qc_parameters)){
		aux_df <- rna_thresholds %>% filter(NAME %in% qc_parameters[param])
		value <- rna_metrics[[param]]
		flag <- ifelse(aux_df$OUT, ".", ifelse(value > aux_df$UPPER, "QC_PASS", ifelse(value < aux_df$LOWER, "QC_FAIL", "QC_WARN")))
		rna_metrics[[paste0(qc_parameters[param], "_FLAG")]] <- flag
	}
	rm(param)
	
	rna_metrics <- rna_metrics %>% 
		mutate(
			across(matches("_READS$", ignore.case = FALSE),  ~ label_number(scale_cut = cut_short_scale())(.x)),
			PANEL_SIZE = label_number(accuracy = 0.01, suffix = "b", scale_cut = cut_short_scale())(PANEL_SIZE),
			across(matches("_COVERAGE$", ignore.case = FALSE), ~ format(.x, big.mark = ",")),
			across(matches("PCT_\\w+", ignore.case = FALSE), ~ paste0(.x, "%"))) %>%
		select(
			TOTAL_READS, ALIGNED_READS, PCT_ALIGNED, ONTARGET_READS, PCT_ONTARGET, HQ_ALIGNED_READS, PCT_HQ_ALIGNED, 
			MEDIAN_READ_LENGTH, MEDIAN_INSERT_SIZE, UNIQUE_READS, PCT_DUP, MEDIAN_COVERAGE, MEAN_COVERAGE,
			matches("PCT_\\d+", ignore.case = FALSE), MIN_COVERAGE, MAX_COVERAGE, PANEL_SIZE, contains("_FLAG")) %>%
		select(as.character(names(qc_parameters[!rna_thresholds$OUT])), everything())
	
	aux_metrics <- as.character(rna_metrics[1,])
	rna_qc <- ifelse("QC_FAIL" %in% aux_metrics, "QC_FAIL", ifelse("QC_WARN" %in% aux_metrics, "QC_WARN", "QC_PASS"))
	rm(aux_metrics)
	
	rna_genecov <- read.table(rna_genecov.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
		mutate(
			TARGET = DP_TARGET, CODING = DP_MANE_CODING, EXON = DP_MANE_EXON, WHOLE = DP_MANE_ALL, 
			MIN = MIN_COVERAGE, MAX = MAX_COVERAGE, across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
		arrange(WHITELIST) %>%
		select(
			GENE, TARGET, CODING, EXON, WHOLE, MIN, MAX, contains("PCT_"), PANEL_SIZE, 
			matches("_ENSEMBL|_REFSEQ"), WHITELIST, LowCoverage)
	rna_genecov_whitelist <- rna_genecov %>% filter(as.logical(WHITELIST))

	if(bamqc_plot_bygene_all_rna){rna_genes <- rna_genecov$GENE} else {rna_genes <- rna_genecov_whitelist$GENE}
	
	rna_exoncov <- read.table(rna_exoncov.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
		mutate(MIN = MIN_COVERAGE, MAX = MAX_COVERAGE) %>%
		select(GENE, EXON, DEPTH, MIN, MAX, LENGTH, contains("PCT_"), CHROM, START, END, contains("TRANSCRIPT_"), WHITELIST, LowCoverage) %>%
		mutate(across(where(is.logical), ~ factor(.x, levels = c(TRUE, FALSE)))) %>%
		filter(bamqc_cov_byexon_all_dna | as.logical(WHITELIST))

} else {
	rna_metrics <- rna_exoncov <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 1))
	rna_qc <- "ND"
}


# Make report ---------------------------------------------------------------

suppressWarnings(
	rmarkdown::render(input = "sampleReport.Rmd", 
	output_file = paste0(run, "_", sample, "_FinalReport_", format(Sys.time(), "%Y%m%d"), ".html")))

