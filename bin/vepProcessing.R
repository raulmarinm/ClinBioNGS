#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of vep results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(VariantAnnotation)))
suppressWarnings(suppressPackageStartupMessages(library(DBI)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
vep.file <- args_list[["vcf"]]
out.file <- args_list[["outFile"]]
hotspots_whitelist.bed <- args_list[["hotspotsWhitelistBed"]]
problematic_regions.bed <- args_list[["problematicRegionsBed"]]
ctr_regions.bed <- args_list[["ctrRegionsBed"]]
mmr.file <- args_list[["mmrGenes"]]
cancerdrivers.file <- args_list[["cancerdrivers"]]
hotspots_results.file <- args_list[["cancerHotspotsResults"]]
sop_oncogenic.file <- args_list[["sopOncogenic"]]
genie_counts.file <- args_list[["genieCounts"]]
genie_oncogenic.file <- args_list[["genieOncogenic"]]
cgi_oncogenic.file <- args_list[["cgiOncogenic"]]
civic_oncogenic.file <- args_list[["civicOncogenic"]]
panel_recurrent.file <- args_list[["panelRecurrentMutations"]]
cadd_cutoff <- as.numeric(args_list[["caddCutoff"]])
revel_cutoff <- as.numeric(args_list[["revelCutoff"]])
alphamissense_cutoff <- as.numeric(args_list[["alphamissenseCutoff"]])
somatic_max_pvaf <- as.numeric(args_list[["somaticMaxpVAF"]])

# Set other variables
aa_names <- c(
	"Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Cys" = "C", "Glu" = "E", "Gln" = "Q", 
	"Gly" = "G", "His" = "H", "Ile" = "I", "Leu" = "L", "Lys" = "K", "Met" = "M", "Phe" = "F", 
	"Pro" = "P", "Ser" = "S", "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V", "Ter" = "*"
)
specific_biomarkers <- c(
	"EGFR_Exon_19_inframe_insertion" = "EGFR_Exon_19_Deletion", 
	"EGFR_Exon_19_inframe_deletion" = "EGFR_Exon_19_Deletion", 
	"EGFR_Exon_19_protein_altering_variant" = "EGFR_Exon_19_Deletion", 
	"EGFR_Exon_20_inframe_insertion" = "EGFR_Exon_20_Insertion", 
	"EGFR_Exon_20_inframe_deletion" = "EGFR_Exon_20_Insertion", 
	"EGFR_Exon_20_protein_altering_variant" = "EGFR_Exon_20_Insertion",
	"ERBB2_Exon_20_inframe_insertion" = "ERBB2_Exon_20_Insertion", 
	"ERBB2_Exon_20_inframe_deletion" = "ERBB2_Exon_20_Insertion", 
	"ERBB2_Exon_20_protein_altering_variant" = "ERBB2_Exon_20_Insertion"
)
bottlenecked_pops <- c("AMI", "ASJ", "FIN", "MID", "OTH", "REMAINING")
clinvar_orig <- c("1" = "Germline", "2" = "Somatic", "3" = "Both")
onco_points <- c(
	"OVS1" = 8,
	"OS1" = 4, "OS2" = 4, "OS3" = 4,
	"OM2" = 2, "OM3" = 2, "OM4" = 2,
	"OP1" = 1, "OP3" = 1, "OP4" = 1,
	"SBP1" = -1, "SBP2" = -1,
	"SBS1" = -4,
	"SBVS1" = -8
)
oncogenic_sop_mut <- c()
oncogenic_sop_pos <- c()
oncogenic_valid_mut <- c()


# Load data ---------------------------------------------------------------

# VEP VCF
vep.vcf <- readVcf(vep.file)
vep_vars <- rownames(vep.vcf)
vep_fields <- unlist(strsplit(info(header(vep.vcf))["CSQ", "Description"], split = "|", fixed = TRUE))
vep_fields[1] <- "Allele"
vep_info <- info(vep.vcf)[, "CSQ"]
names(vep_info) <- vep_vars
vep_ranges <- rowRanges(vep.vcf)

# Genomic info
hotspots_whitelist_ranges <- makeGRangesFromDataFrame(fread(
	hotspots_whitelist.bed, header = FALSE, col.names = c("seqnames", "start", "end"),
	sep = "\t", quote = "", fill = FALSE, nThread = 1), starts.in.df.are.0based = TRUE)
problematic_regions_ranges <- makeGRangesFromDataFrame(fread(
	problematic_regions.bed, header = FALSE, col.names = c("seqnames", "start", "end"), 
	sep = "\t", quote = "", fill = FALSE, nThread = 1), starts.in.df.are.0based = TRUE)
ctr_regions_ranges <- makeGRangesFromDataFrame(fread(
	ctr_regions.bed, header = FALSE, col.names = c("seqnames", "start", "end"), 
	sep = "\t", quote = "", fill = FALSE, nThread = 1), starts.in.df.are.0based = TRUE)

# MMR genes
if(file.exists(mmr.file)){mmr_genes <- read.table(mmr.file, header = TRUE, sep = ";")[, 1]} else {mmr_genes <- c()}

# Network of Cancer Genes (NCG) data
cancerdrivers <- read.table(cancerdrivers.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

# Cancer Hotspots results
hotspots_results <- read.table(hotspots_results.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
	select(MUT_ID, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT)

# GENIE mutation counts
if(file.exists(genie_counts.file)){
	genie_counts <- read.table(genie_counts.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
} else {
	genie_counts <- data.frame(MUT_ID = NA, GENIE_CNT = NA)
}

# GENIE oncogenic mutations with AA change
if(file.exists(genie_oncogenic.file)){
	genie_oncogenic <- read.table(genie_oncogenic.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
		filter(ONCOGENICITY %in% "Oncogenic" & !is.na(AA_POS)) %>%
		mutate(AA_POS_ID = paste0(GENE, "_", AA_POS))
	oncogenic_sop_mut <- na.omit(unique(c(oncogenic_sop_mut, genie_oncogenic$MUT_ID)))
	oncogenic_sop_pos <- na.omit(unique(c(oncogenic_sop_pos, genie_oncogenic$AA_POS_ID)))
	rm(genie_oncogenic)
}

# ClinGen/CGC/VICC SOP oncogenic mutations
if(file.exists(sop_oncogenic.file)){
	sop_mutations <- read.table(sop_oncogenic.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
		mutate(AA_POS_ID = paste0(GENE, "_", AA_POS))
	oncogenic_sop_mut <- na.omit(unique(c(oncogenic_sop_mut, sop_mutations$MUT_ID)))
	oncogenic_sop_pos <- na.omit(unique(c(oncogenic_sop_pos, sop_mutations$AA_POS_ID)))
	rm(sop_mutations)
}

# CGI: Catalog of Validated Oncogenic Mutations
cgi_oncogenic <- read.table(cgi_oncogenic.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
oncogenic_valid_mut <- na.omit(unique(c(oncogenic_valid_mut, cgi_oncogenic$MUT_ID)))

# CIViC: oncogenic evidences
civic_oncogenic <- read.table(civic_oncogenic.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
oncogenic_valid_mut <- na.omit(unique(c(oncogenic_valid_mut, civic_oncogenic$ALTERATION)))

# Panel recurrent mutations
if(file.exists(panel_recurrent.file)){
	panel_recurrent <- read.table(panel_recurrent.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)$VAR
} else {panel_recurrent <- c()}


# Data processing ---------------------------------------------------------

# Collect and prepare VEP annotation
#	- Add VAR (chr:pos_ref/alt)
#	- Create a MUT_ID (GEN_MUTATION)
#	- Annotate variant impact: VEP, ClinVar, predictors
#	- Annotate if variants are located in PROBLEMATIC and known/panel HOTSPOT regions
#	- Calculate max gnomAD AF of non-bottleneck populations
#	- Annotate the gene role (DRIVER/ONCOGENE/TSG)
#	- Annotate Cancer Hotspots and GENIE mutation counts
#	- Annotate oncogenic variants
#	- Annotate if variants are recurrently reported in this panel
vep_annot <- unlist(vep_info)
vep_vars <- names(vep_annot)
vars_annot <- as.data.frame(matrix(unlist(vep_annot), nrow = length(vep_annot), byrow = TRUE)) %>%
	separate(col = 1, into = vep_fields, sep = "\\|", remove = TRUE) %>% 
	separate(EXON, c("EXON_NUM", NA), sep = "/", fill = "left", remove = FALSE) %>%
	separate(Consequence, "CONSEQUENCE", sep = "&", extra = "drop") %>%
	mutate(
		VAR = vep_vars, isGENE = SYMBOL != "", isHGVSc = HGVSc != "", isHGVSp = HGVSp != "", CLASS = VARIANT_CLASS,
		CANONICAL = CANONICAL %in% "YES", TRANSCRIPT_REFSEQ = paste0(MANE_SELECT, MANE_PLUS_CLINICAL), 
		MANE_SELECT = MANE_SELECT != "", MANE_PLUS_CLINICAL = MANE_PLUS_CLINICAL != "", STRAND = ifelse(STRAND == "1", "+", "-"), 
		GENE_SYMBOL = SYMBOL, GENE_ENSEMBL = Gene, TRANSCRIPT_ENSEMBL = Feature, PROTEIN_ENSEMBL = ENSP, 
		cDNA_POS = cDNA_position, CDS_POS = CDS_position, PROTEIN_POS = Protein_position, AA = Amino_acids, CODONS = Codons,
		CODING = CDS_POS != "", NONCODING = !CODING,
		EXISTING_VARIATION = Existing_variation, NMD_ESCAPE = NMD != "", 
		MMR_GENE = GENE_SYMBOL %in% mmr_genes,
		CIViC_VARIANT_ID = CIViC, CIViC_GENE = CIViC_GN, CIViC_MUT = CIViC_VT,

		GENE_HGNC = gsub(".*:", "", HGNC_ID),
		HGVSc_ENSEMBL = HGVSc, HGVSc = gsub(".*:", "", HGVSc),
		HGVSc_REFSEQ = ifelse(isHGVSc & (MANE_SELECT | MANE_PLUS_CLINICAL), paste0(TRANSCRIPT_REFSEQ, ":", HGVSc), NA),
		HGVSp = gsub("%3D", "=", HGVSp), HGVSp_ENSEMBL = HGVSp, HGVSp = gsub(".*:", "", HGVSp),
		HGVSp_SHORT = str_replace_all(gsub("p\\.", "", HGVSp), pattern = aa_names),
		AA_POS = as.integer(gsub("\\D", "", gsub("\\*.*", "", gsub("_.*", "", HGVSp_SHORT)))),
		MUTATION = ifelse(isHGVSp, HGVSp_SHORT, ifelse(isHGVSc, HGVSc, HGVSg)),
		MUT_ID = ifelse(isGENE & (isHGVSp | isHGVSc), paste0(GENE_SYMBOL, "_", MUTATION), ifelse(
				(!isGENE) & (isHGVSp | isHGVSc), paste0(GENE_ENSEMBL, "_", MUTATION), MUTATION)),
		AA_POS_ID = ifelse(is.na(AA_POS), NA, paste0(GENE_SYMBOL, "_", AA_POS)),
		OTHER_ID = str_replace_all(paste0(GENE_SYMBOL, "_Exon_", EXON_NUM, "_", CONSEQUENCE), pattern = specific_biomarkers),
		
		STOPLOSS_VARIANT = CONSEQUENCE %in% "stop_lost",
		NULL_VARIANT = IMPACT %in% "HIGH" & !STOPLOSS_VARIANT,
		INFRAME_INDEL = CONSEQUENCE %in% c("inframe_deletion", "inframe_insertion"),
		MISSENSE_VARIANT = CONSEQUENCE %in% "missense_variant",
		SILENT_VARIANT = IMPACT %in% c("LOW", "MODIFIER"),
				
		SIFT_TERM = gsub("\\(.*\\)", "", SIFT), SIFT_SCORE = as.numeric(gsub(".*\\(|\\)", "", SIFT)),
		POLYPHEN_TERM = gsub("\\(.*\\)", "", PolyPhen), POLYPHEN_SCORE = as.numeric(gsub(".*\\(|\\)", "", PolyPhen)),
		CADD_SCORE = as.numeric(CADD_PHRED), CADD_TERM = ifelse(CADD_SCORE >= cadd_cutoff, "likely_pathogenic", "likely_benign"),
		REVEL_SCORE = as.numeric(REVEL), REVEL_TERM = ifelse(REVEL_SCORE >= revel_cutoff, "likely_pathogenic", "likely_benign"),
		AlphaMissense_SCORE = as.numeric(am_pathogenicity),
		AlphaMissense_TERM = ifelse(AlphaMissense_SCORE >= alphamissense_cutoff, "likely_pathogenic", "likely_benign"),
		
		ClinVar_ID = ClinVar,
		ClinVar_ORIGIN = ifelse(ClinVar_ORIGIN == "", "", ifelse(ClinVar_ORIGIN %in% names(clinvar_orig), clinvar_orig[ClinVar_ORIGIN], "Other")),
		ClinVar_PATHOGENIC = grepl("Pathogenic|Likely_pathogenic|drug_response", ClinVar_CLNSIG),

		ONCOGENIC_SOP_MUT = MUT_ID %in% oncogenic_sop_mut, ONCOGENIC_SOP_POS = AA_POS_ID %in% oncogenic_sop_pos,
		ONCOGENIC_VALID_MUT = MUT_ID %in% oncogenic_valid_mut,
		ONCOGENIC_VALID_VAR = HGVSg %in% cgi_oncogenic$HGVSg | CIViC_VARIANT_ID %in% civic_oncogenic$VARIANT_ID,

		RECURRENT = VAR %in% panel_recurrent,
		SOMATIC_WHITELIST = overlapsAny(vep_ranges, hotspots_whitelist_ranges),
		PROBLEMATIC_REGION = overlapsAny(vep_ranges, problematic_regions_ranges),
		CTR_REGION = overlapsAny(vep_ranges, ctr_regions_ranges),

		across(everything(), ~ ifelse(.x == "", NA, .x))) %>% 
	left_join(cancerdrivers, by = c("GENE_SYMBOL" = "GENE")) %>%
	left_join(hotspots_results, by = "MUT_ID") %>%
	left_join(genie_counts, by = "MUT_ID") %>%
	left_join(genie_counts, by = c("VAR" = "MUT_ID")) %>%
	mutate(
		Existing_variation = ifelse(is.na(Existing_variation), "", Existing_variation),
		across(c(ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER), ~ .x & !is.na(.x)),
		across(contains(c("gnomAD", "CNT", "EVIDENCES")), ~ ifelse(is.na(.x), 0, as.numeric(.x))),
		CANCER_HOTSPOT = SOMATIC_WHITELIST | HOTSPOT_MUT_CNT >= 1) %>%
	rowwise() %>%
	mutate(
		GENIE_CNT = sum(GENIE_CNT.x, GENIE_CNT.y),
		IMPACT_PREDICTOR = paste(na.omit(unique(c(CADD_TERM, REVEL_TERM, AlphaMissense_TERM))), collapse = ""),
		dbSNP_ID = ifelse(Existing_variation != "", paste0("", grep("rs", unlist(strsplit(Existing_variation, "&")), value = TRUE)), ""),
		gnomAD_MAX_AF = max(across(contains("gnomAD") & !contains(bottlenecked_pops)))) %>% 
	as.data.frame() %>%
	mutate(PREDICTOR_PATHOGENIC = IMPACT_PREDICTOR %in% "likely_pathogenic", PREDICTOR_BENIGN = IMPACT_PREDICTOR %in% "likely_benign")


# Oncogenicity classification (ClinGen/CGC/VICC SOP)
#	- Oncogenic: >= 10 points
#	- Likely Oncogenic: 6-9 points
#	- VUS: 0-5 points
#	- Likely Benign: -1-(-6) points
#	- Benign: <= -7 points
vars_annot <- vars_annot %>%
	mutate(
		OVS1 = NULL_VARIANT & TSG & TSG_EVIDENCES > 1,
		OS1 = ONCOGENIC_SOP_MUT,
		OS2 = (ONCOGENIC_VALID_MUT & !OS1) | ONCOGENIC_VALID_VAR,
		OS3 = HOTSPOT_MUT_CNT >= 10 & HOTSPOT_POS_CNT >= 50 & !OS1, 
		OM2 = (INFRAME_INDEL & (ONCOGENE|TSG)) | (STOPLOSS_VARIANT & TSG),
		OM4 = MISSENSE_VARIANT & ONCOGENIC_SOP_POS & !(OS1 | OS3),
		OM3 = HOTSPOT_MUT_CNT >= 10 & HOTSPOT_POS_CNT < 50 & !OM4,
		OP1 = PREDICTOR_PATHOGENIC, 
		OP3 = HOTSPOT_MUT_CNT >= 1 & HOTSPOT_MUT_CNT < 10,
		OP4 = gnomAD_MAX_AF <= somatic_max_pvaf,
		SBP1 = PREDICTOR_BENIGN,
		SBP2 = SILENT_VARIANT & !PREDICTOR_PATHOGENIC,
		SBS1 = gnomAD_MAX_AF > 0.01 & gnomAD_MAX_AF <= 0.05,
		SBVS1 = gnomAD_MAX_AF > 0.05) %>%
	rowwise() %>%
	mutate(
		ONCO_POINTS = sum(onco_points[c(OVS1, OS1, OS2, OS3, OM2, OM3, OM4, OP1, OP3, OP4, SBP1, SBP2, SBS1, SBVS1)]),
		ONCO_CODES = paste(collapse = ", ", names(onco_points)[
			c(OVS1, OS1, OS2, OS3, OM2, OM3, OM4, OP1, OP3, OP4, SBP1, SBP2, SBS1, SBVS1)])) %>% 
	as.data.table() %>%
	mutate(
		ONCOGENICITY = ifelse(is.na(ONCO_POINTS), ".", ifelse(
			ONCO_POINTS >= 10, "Oncogenic", ifelse(
				ONCO_POINTS >= 6, "Likely Oncogenic", ifelse(
					ONCO_POINTS >= 0, "VUS", ifelse(
						ONCO_POINTS >= -6, "Likely Benign", ifelse(
							ONCO_POINTS <= -7, "Benign", ".")))))),
		BENIGN = ONCOGENICITY %in% c("Benign", "Likely Benign"))

# Create final table
vars_annot <- vars_annot %>%
	select(
		CLASS, VAR, MUT_ID, CANCER_HOTSPOT, PROBLEMATIC_REGION, CTR_REGION, RECURRENT, ONCOGENICITY, 
		ONCO_POINTS, ONCO_CODES, CIViC_VARIANT_ID, ClinVar_ID, ClinVar_CLNSIG, 
		GENE_SYMBOL, MUTATION, OTHER_ID, CONSEQUENCE, IMPACT, NMD_ESCAPE, MMR_GENE, 
		SOMATIC_WHITELIST, GENIE_CNT, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, BENIGN,
		ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, gnomAD_MAX_AF, dbSNP_ID, EXISTING_VARIATION,
		HGVSg, HGVSc, HGVSp, HGVSp_SHORT, HGVSc_ENSEMBL, HGVSc_REFSEQ, HGVSp_ENSEMBL,
		
		GENE_ENSEMBL, GENE_HGNC, BIOTYPE, STRAND, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, PROTEIN_ENSEMBL, CCDS,
		EXON, INTRON, cDNA_POS, CDS_POS, PROTEIN_POS, AA, AA_POS, CODONS, CODING, NONCODING, 
		CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, TSL, ONCOGENE_EVIDENCES, TSG_EVIDENCES,

		contains("gnomAD"),

		POLYPHEN_TERM, POLYPHEN_SCORE, SIFT_TERM, SIFT_SCORE, REVEL_TERM, REVEL_SCORE,
		AlphaMissense_TERM, AlphaMissense_SCORE, CADD_SCORE, CADD_TERM, PREDICTOR_PATHOGENIC, PREDICTOR_BENIGN,
		ClinVar_CLNREVSTAT, ClinVar_CLNDN, ClinVar_ORIGIN, ClinVar_PATHOGENIC, 
		CIViC_GENE, CIViC_MUT,

		NULL_VARIANT, STOPLOSS_VARIANT, INFRAME_INDEL, MISSENSE_VARIANT, SILENT_VARIANT,
		ONCOGENIC_SOP_MUT, ONCOGENIC_SOP_POS, ONCOGENIC_VALID_VAR, ONCOGENIC_VALID_MUT,
		
		SBVS1:OVS1
	)

vars_annot[is.na(vars_annot) | vars_annot == "NA" | vars_annot == ""] <- "."



# Save data ---------------------------------------------------------------

write.table(vars_annot, file = out.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

