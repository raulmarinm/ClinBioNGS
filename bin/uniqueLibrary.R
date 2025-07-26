#!/usr/bin/env Rscript

# ============================================================================ #
# What: Collect runs results and create a library of unique variants
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
unique_library_name <- args_list[["libraryName"]]
run.db <- args_list[["runLibrary"]]
results.dir <- args_list[["resultsDir"]]

# Set input unique library
unique.db <- list.files("input", unique_library_name, full.names = TRUE)

# Set other variables
new_unique.db <- paste0(unique_library_name, "_", format(Sys.time(), "%Y%m%d"), ".db")
run_library_path <- paste0(results.dir, "/", basename(run.db))
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
chr_arm_levels <- unlist(lapply(c(1:22, "X", "Y"), function(x) paste0(x, c("p", "q"))))


# Load data ---------------------------------------------------------------

# Load run library DB
if(file.exists(run.db)){
	run_conn <- dbConnect(RSQLite::SQLite(), run.db)
	run_samples <- as.data.frame(tbl(run_conn, "Samples"))
	dbDisconnect(run_conn)
} else {
	run_samples <- matrix(data = NA, nrow = 0, ncol = 1)
}

# Load unique library DB
if(length(unique.db)){
	unique_conn <- dbConnect(RSQLite::SQLite(), unique.db)
	runs_lib <- as.data.frame(tbl(unique_conn, "Runs"))
	dbDisconnect(unique_conn)
} else {
	runs_lib <- data.frame(RUN = NA, DB = NA) %>% filter(!is.na(RUN))
}


# Data processing ---------------------------------------------------------

# Registrate or update the RUN library into the UNIQUE one (in case of having samples for library)
if(nrow(run_samples)){
	if(run %in% runs_lib$RUN){
		runs_lib <- runs_lib %>% mutate(DB = ifelse(RUN %in% run, run_library_path, DB))
	} else {
		runs_lib <- rbind(runs_lib, data.frame(RUN = run, DB = run_library_path))
	}
} else {
	# Remove the run if it is registered in the UNIQUE library
	if(run %in% runs_lib$RUN){runs_lib <- runs_lib %>% filter(!RUN %in% run)}
}

# Collect the unique results from each registered run library
all_sample_info <- matrix(data = NA, nrow = 0, ncol = 1)
all_bamqc_dna <- matrix(data = NA, nrow = 0, ncol = 1)
all_bamqc_rna <- matrix(data = NA, nrow = 0, ncol = 1)
all_variants <- matrix(data = NA, nrow = 0, ncol = 1)
all_variants_annot <- matrix(data = NA, nrow = 0, ncol = 1)
all_cna_gene <- matrix(data = NA, nrow = 0, ncol = 1)
all_cna_arm <- matrix(data = NA, nrow = 0, ncol = 1)
all_fusion <- matrix(data = NA, nrow = 0, ncol = 1)
all_splicing <- matrix(data = NA, nrow = 0, ncol = 1)
all_tmb <- matrix(data = NA, nrow = 0, ncol = 1)
all_msi <- matrix(data = NA, nrow = 0, ncol = 1)

cat(paste0("\n", "Collecting results from:", "\n"))

for (run_name in runs_lib$RUN){

	aux_lib <- runs_lib %>% filter(RUN == run_name)

	# Load the RUN library and collect the results
	run_conn <- dbConnect(RSQLite::SQLite(), aux_lib$DB)
	library_tables <- dbListTables(run_conn)

	all_sample_info <- rbind(all_sample_info, as.data.frame(tbl(run_conn, "Samples")))
	if("BamQC_DNA" %in% library_tables){all_bamqc_dna <- rbind(all_bamqc_dna, as.data.frame(tbl(run_conn, "BamQC_DNA")))}
	if("BamQC_RNA" %in% library_tables){all_bamqc_rna <- rbind(all_bamqc_rna, as.data.frame(tbl(run_conn, "BamQC_RNA")))}
	if("TMB" %in% library_tables){all_tmb <- rbind(all_tmb, as.data.frame(tbl(run_conn, "TMB")))}
	if("MSI" %in% library_tables){all_msi <- rbind(all_msi, as.data.frame(tbl(run_conn, "MSI")))}		

	if("SmallVariant_calling" %in% library_tables){
		variants <- tbl(run_conn, "SmallVariant_calling") %>% filter(PASS) %>%
			select(CHROM, START, END, REF, ALT, VAR, TYPE, START_HG19, END_HG19, VAR_HG19) %>%
			collect() %>% group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_variants <- rbind(all_variants, variants) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(variants)
	}	
	
	if("SmallVariant_annot" %in% library_tables){
		variants_annot <- tbl(run_conn, "SmallVariant_annot") %>% filter(PASS) %>% 
			select(
				VAR, GENE_SYMBOL, MUTATION, CLIN_RELEVANT, CLIN_STATUS, ONCOGENICITY, ONCO_POINTS, ONCO_CODES, 
				CONSEQUENCE, IMPACT, TYPE, CLASS, EXON, INTRON, 
				HGVSg, HGVSc, HGVSp, HGVSp_SHORT, HGVSc_ENSEMBL, HGVSc_REFSEQ, HGVSp_ENSEMBL, 
				dbSNP_ID, ClinVar_ID, EXISTING_VARIATION, gnomAD_MAX_AF,
				NMD_ESCAPE, GENIE_CNT, SOMATIC_WHITELIST, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, BENIGN,
				GENE_WHITELIST, MMR_GENE, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, GENE_ENSEMBL, GENE_HGNC, BIOTYPE, CANONICAL, 
				MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, TSL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, 
				cDNA_POS, CDS_POS, CCDS, PROTEIN_ENSEMBL, PROTEIN_POS, AA, CODONS, CODING, NONCODING, 
				contains("gnomAD"),
				POLYPHEN_TERM, POLYPHEN_SCORE, SIFT_TERM, SIFT_SCORE, REVEL_TERM, REVEL_SCORE, AlphaMissense_TERM, AlphaMissense_SCORE, 
				CADD_TERM, CADD_SCORE, PREDICTOR_PATHOGENIC, PREDICTOR_BENIGN,
				contains("ClinVar_"),
				CANCER_HOTSPOT, CTR_REGION, PROBLEMATIC_REGION, 
				NULL_VARIANT, STOPLOSS_VARIANT, INFRAME_INDEL, MISSENSE_VARIANT, SILENT_VARIANT,
				ONCOGENIC_SOP_MUT, ONCOGENIC_SOP_POS, ONCOGENIC_VALID_VAR, ONCOGENIC_VALID_MUT, SBVS1:OVS1,
				MUT_ID, CIViC_ALTERATION, CIViC_VARIANT_ID, CIViC_MOLECULAR_ID, CIViC_VARIANT_GROUP,
				CHROM, START, END, REF, ALT, STRAND, VAR_HG19) %>%
			collect() %>% 
			mutate(CLIN_STATUS = ifelse(CLIN_RELEVANT, "Tier I/II", CLIN_STATUS)) %>%
			group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_variants_annot <- rbind(all_variants_annot, variants_annot) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(variants_annot)
	}
	
	if("CNA_gene" %in% library_tables){
		cna_gene <- tbl(run_conn, "CNA_gene") %>% filter(!NEUTRAL) %>%
			select(
				ALTERATION, VAR, CLIN_RELEVANT, CLIN_STATUS, CYTOBAND, GENE, CNA, CLASS, WHITELIST, GENIE_CNT, GENIE_FREQ, 
				ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ,
				NO_FLAGS, LOW_AMP, LOW_DEL, LOWBINS, CHROM, START, END, STRAND, 
				CIViC_ALTERATION, CIViC_VARIANT_ID, CIViC_MOLECULAR_ID) %>% 
			collect() %>% 
			mutate(CLIN_STATUS = ifelse(CLIN_RELEVANT, "Tier I/II", CLIN_STATUS)) %>%
			group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_cna_gene <- rbind(all_cna_gene, cna_gene) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(cna_gene)
	}
	
	if("CNA_arm" %in% library_tables){
		cna_arm <- tbl(run_conn, "CNA_arm") %>% filter(CNA %in% c("AMP", "DEL")) %>%
			select(VAR, CHROM, START, END, CHR_ARM, SIDE, CNA) %>% collect() %>%
			group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_cna_arm <- rbind(all_cna_arm, cna_arm) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(cna_arm)
	}
	
	if("Fusion" %in% library_tables){
		fusion <- tbl(run_conn, "Fusion") %>% filter(NO_FLAGS) %>%
			select(
				VAR, FUSION_SHORT, FUSION_NAME, FUSION_VARIANT, CLIN_RELEVANT, CLIN_STATUS, WHITELIST_FUSION, 
				RESOURCES, MitelmanDB_COUNT, GENIE_CNT, WHITELIST_GENE, FUSION_RANGE_A, FUSION_RANGE_B, 
				ONCOGENE_A, TSG_A, DRIVER_A, CANONICAL_DRIVER_A, GENE_A, GENE_ENSEMBL_A, 
				TRANSCRIPT_ENSEMBL_A, TRANSCRIPT_REFSEQ_A, EXON_A, INTRON_A, BASES_FROM_EXON_A, 
				ONCOGENE_B, TSG_B, DRIVER_B, CANONICAL_DRIVER_B, GENE_B, GENE_ENSEMBL_B, 
				TRANSCRIPT_ENSEMBL_B, TRANSCRIPT_REFSEQ_B, EXON_B, INTRON_B, BASES_FROM_EXON_B, 
				BREAKPOINT_A, CHROM_A, POS_A, STRAND_A, BREAKPOINT_B, CHROM_B, POS_B, STRAND_B, 
				VAR_HG19, BREAKPOINT_A_HG19, BREAKPOINT_B_HG19, UNKNOWN, MitelmanDB, GENIE, 
				CIViC_ALTERATION, CIViC_VARIANT_ID, CIViC_MOLECULAR_ID, CIViC_VARIANT_GROUP) %>% 
			collect() %>% 
			mutate(CLIN_STATUS = ifelse(CLIN_RELEVANT, "Tier I/II", CLIN_STATUS)) %>%
			group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_fusion <- rbind(all_fusion, fusion) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(fusion)
	}
	
	if("Splicing" %in% library_tables){
		splicing <- tbl(run_conn, "Splicing") %>% filter(NO_FLAGS) %>%
			select(
				VAR, VAR_NAME, VAR_GENE, CLIN_RELEVANT, CLIN_STATUS, WHITELIST_SPLICING, WHITELIST_GENE, GENE, REGION_AFFECTED, 
				ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, 
				EXONS_AFFECTED, INTRONS_AFFECTED, CHROM, START, END, STRAND, VAR_HG19, START_HG19, END_HG19, 
				CancerEnriched, CIViC_ALTERATION, CIViC_VARIANT_ID, CIViC_MOLECULAR_ID) %>% 
			collect() %>% 
			mutate(CLIN_STATUS = ifelse(CLIN_RELEVANT, "Tier I/II", CLIN_STATUS)) %>%
			group_by(VAR) %>% mutate(AC_SAMPLES = n()) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		all_splicing <- rbind(all_splicing, splicing) %>% group_by(VAR) %>% 
			mutate(AC_SAMPLES = sum(AC_SAMPLES)) %>% distinct(VAR, .keep_all = TRUE) %>% as.data.frame()
		rm(splicing)
	}

	dbDisconnect(run_conn)

	rm(aux_lib, run_conn, library_tables)

	cat(paste0("   - ", run_name, "\n"))
}
rm(run_name)


# Save data ---------------------------------------------------------------

# Load the new DB and write tables
# Calculate de sample frequencies of each alteration
new_unique_conn <- dbConnect(RSQLite::SQLite(), new_unique.db)
dbWriteTable(new_unique_conn, "Runs", runs_lib, overwrite = TRUE)
if(nrow(all_sample_info)){dbWriteTable(new_unique_conn, "Samples", all_sample_info, overwrite = TRUE)}
if(nrow(all_bamqc_dna)){dbWriteTable(new_unique_conn, "BamQC_DNA", all_bamqc_dna, overwrite = TRUE)}
if(nrow(all_bamqc_rna)){dbWriteTable(new_unique_conn, "BamQC_RNA", all_bamqc_rna, overwrite = TRUE)}
if(nrow(all_variants)){
	variants_N <- nrow(all_sample_info %>% filter(as.logical(SMALL_VARIANT)))
	all_variants <- all_variants %>% 
		mutate(
			AF_SAMPLES = round(AC_SAMPLES / variants_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", variants_N), 
			CHROM = factor(CHROM, levels = chr_levels)) %>%
		arrange(CHROM, START, END, REF, ALT)
	dbWriteTable(new_unique_conn, "SmallVariant_calling", all_variants, overwrite = TRUE)
}
if(nrow(all_variants_annot)){
	variants_N <- nrow(all_sample_info %>% filter(as.logical(SMALL_VARIANT)))
	all_variants_annot <- all_variants_annot %>% 
		mutate(
			AF_SAMPLES = round(AC_SAMPLES / variants_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", variants_N), 
			CHROM = factor(CHROM, levels = chr_levels)) %>%
		arrange(CHROM, START, END, REF, ALT)
	dbWriteTable(new_unique_conn, "SmallVariant_annot", all_variants_annot, overwrite = TRUE)
}
if(nrow(all_cna_gene)){
	cna_N <- nrow(all_sample_info %>% filter(as.logical(CNA)))
	all_cna_gene <- all_cna_gene %>% 
		mutate(AF_SAMPLES = round(AC_SAMPLES / cna_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", cna_N)) %>% 
		arrange(GENE, CNA)
	dbWriteTable(new_unique_conn, "CNA_gene", all_cna_gene, overwrite = TRUE)
}
if(nrow(all_cna_arm)){
	cna_N <- nrow(all_sample_info %>% filter(as.logical(CNA)))
	all_cna_arm <- all_cna_arm %>% 
		mutate(
			AF_SAMPLES = round(AC_SAMPLES / cna_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", cna_N), 
			CHR_ARM = factor(CHR_ARM, levels = chr_arm_levels)) %>% arrange(CHR_ARM, CNA)
	dbWriteTable(new_unique_conn, "CNA_arm", all_cna_arm, overwrite = TRUE)
}
if(nrow(all_fusion)){
	fusion_N <- nrow(all_sample_info %>% filter(as.logical(FUSION)))
	all_fusion <- all_fusion %>% 
		mutate(AF_SAMPLES = round(AC_SAMPLES / fusion_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", fusion_N)) %>% 
		arrange(GENE_A, POS_A, GENE_B, POS_B)
	dbWriteTable(new_unique_conn, "Fusion", all_fusion, overwrite = TRUE)
}
if(nrow(all_splicing)){
	splicing_N <- nrow(all_sample_info %>% filter(as.logical(SPLICING)))
	all_splicing <- all_splicing %>% 
		mutate(AF_SAMPLES = round(AC_SAMPLES / splicing_N, 3), AC_SAMPLES = paste0(AC_SAMPLES, "/", splicing_N)) %>% 
		arrange(GENE, START, END)
	dbWriteTable(new_unique_conn, "Splicing", all_splicing, overwrite = TRUE)
}
if(nrow(all_tmb)){dbWriteTable(new_unique_conn, "TMB", all_tmb, overwrite = TRUE)}
if(nrow(all_msi)){dbWriteTable(new_unique_conn, "MSI", all_msi, overwrite = TRUE)}
dbDisconnect(new_unique_conn)

cat(paste0("\n", "The global variant library is created! (unique variants)", "\n"))
cat(paste0("    ", basename(new_unique.db), "\n"))
cat(paste0("    ", nrow(runs_lib), " runs", "\n"))
cat(paste0("    ", nrow(all_sample_info), " samples", "\n"))
cat(paste0("      - ", nrow(all_bamqc_dna), " DNA", "\n"))
cat(paste0("      - ", nrow(all_bamqc_rna), " RNA", "\n"))
if(nrow(all_variants_annot)){
	cat(paste0("    ", nrow(all_variants_annot), " small variants", "\n"))
	cat(paste0("      - ", nrow(all_variants_annot%>%filter(TYPE == "SNV")), " SNVs", "\n"))
	cat(paste0("      - ", nrow(all_variants_annot%>%filter(TYPE != "SNV")), " InDels", "\n"))
	cat(paste0("      AMP/ASCO/CAP", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(as.logical(CLIN_RELEVANT))), " Tier I/II", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
	cat(paste0("      ClinGen/CGC/VICC", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(ONCOGENICITY %in% "Oncogenic")), " Oncogenic", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(ONCOGENICITY %in% "Likely Oncogenic")), " Likely Oncogenic", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(ONCOGENICITY %in% "VUS")), " VUS", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(ONCOGENICITY %in% "Likely Benign")), " Likely Benign", "\n"))
	cat(paste0("        - ", nrow(all_variants_annot%>%filter(ONCOGENICITY %in% "Benign")), " Benign", "\n"))
} else if(nrow(all_variants)) {
	cat(paste0("    ", nrow(all_variants), " small variants (NO annotation)", "\n"))
	cat(paste0("      - ", nrow(all_variants%>%filter(TYPE == "SNV")), " SNVs", "\n"))
	cat(paste0("      - ", nrow(all_variants%>%filter(TYPE != "SNV")), " InDels", "\n"))
} else {
	cat(paste0("    ", nrow(all_variants), " small variants", "\n"))
}
cat(paste0("    ", nrow(all_cna_gene), " gene CNAs", "\n"))
if(nrow(all_cna_gene)){
	ok_cna <- all_cna_gene%>%filter(as.logical(NO_FLAGS))
	cat(paste0("      - ", nrow(ok_cna), " OK CNAs (no flags)", "\n"))
	cat(paste0("        - ", nrow(ok_cna%>%filter(CNA %in% "AMP")), " AMP", "\n"))
	cat(paste0("        - ", nrow(ok_cna%>%filter(CNA %in% "DEL")), " DEL", "\n"))
	cat(paste0("        AMP/ASCO/CAP", "\n"))
	cat(paste0("          - ", nrow(ok_cna%>%filter(as.logical(CLIN_RELEVANT))), " Tier I/II", "\n"))
	cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
	cat(paste0("          - ", nrow(ok_cna%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
	cat(paste0("      - ", nrow(all_cna_gene%>%filter(as.logical(LOW_AMP))), " 'LowAmp' CNAs", "\n"))
	cat(paste0("      - ", nrow(all_cna_gene%>%filter(as.logical(LOW_DEL))), " 'LowDel' CNAs", "\n"))
	cat(paste0("      - ", nrow(all_cna_gene%>%filter(as.logical(LOWBINS))), " 'LowBins' CNAs", "\n"))
}
cat(paste0("    ", nrow(all_cna_arm), " arm CNAs", "\n"))
if(nrow(all_cna_arm)){
	cat(paste0("      - ", nrow(all_cna_arm%>%filter(CNA %in% "AMP")), " AMP", "\n"))
	cat(paste0("      - ", nrow(all_cna_arm%>%filter(CNA %in% "DEL")), " DEL", "\n"))
}
cat(paste0("    ", nrow(all_fusion), " gene fusions", "\n"))
if(nrow(all_fusion)){
	cat(paste0("      AMP/ASCO/CAP", "\n"))
	cat(paste0("        - ", nrow(all_fusion%>%filter(as.logical(CLIN_RELEVANT))), " Tier I/II", "\n"))
	cat(paste0("        - ", nrow(all_fusion%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
	cat(paste0("        - ", nrow(all_fusion%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
}
cat(paste0("    ", nrow(all_splicing), " splicing variants", "\n"))
if(nrow(all_splicing)){
	cat(paste0("      AMP/ASCO/CAP", "\n"))
	cat(paste0("        - ", nrow(all_splicing%>%filter(as.logical(CLIN_RELEVANT))), " Tier I/II", "\n"))
	cat(paste0("        - ", nrow(all_splicing%>%filter(CLIN_STATUS %in% "Tier III")), " Tier III", "\n"))
	cat(paste0("        - ", nrow(all_splicing%>%filter(CLIN_STATUS %in% "Tier IV")), " Tier IV", "\n"))
}
