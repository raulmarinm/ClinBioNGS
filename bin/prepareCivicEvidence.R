#!/usr/bin/env Rscript

# ============================================================================ #
# What: Prepare the CIViC evidence file
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
tumor_names.file <- args_list[["tumorNames"]]
var.file <- args_list[["variantFile"]]
molecular.file <- args_list[["molecularFile"]]
evidence.file <- args_list[["evidenceFile"]]
civic_date <- args_list[["date"]]

# Set other variables
specific_biomarkers <- c(
	"v::" = ".::", "\\?::" = ".::", "::v" = "::.", "::\\?" = "::.",
	"MET_Exon_14_Skipping_Mutation" = "METx14del", "EGFR_VIII" = "EGFRvIII,EGFRvIIIb", "AR_V7_EXPRESSION" = "AR-V7")
evidence_levels <- c("A", "B", "C", "D")
evidence_rating <- c(3:5)
evidence_types <- c("Predictive", "Prognostic", "Diagnostic")
prognostic_effect <- c("Better Outcome", "Poor Outcome")
diagnostic_effect <- c("Positive")
predictive_effect <- c("Sensitivity/Response", "Reduced Sensitivity", "Resistance", "Adverse Response")



# Load data ---------------------------------------------------------------

# Tumor names
tumor_names <- read.table(tumor_names.file, header = TRUE, sep = ";", quote = "\"", na.strings = c("NA", ".", ""), fill = TRUE) %>%
	mutate(TUMOR_NAME = tolower(TUMOR_NAME))
tumor_doid <- tumor_names %>% filter(!is.na(TUMOR_DOID)) %>% distinct(TUMOR_DOID, .keep_all = TRUE)

# CIViC data
var_data <- read.table(var.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
molecular_data <- read.table(molecular.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
evidence_data <- read.table(evidence.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)



# Data processing ---------------------------------------------------------

# CNAs
var_cna <- var_data %>% filter(tolower(variant) %in% c("amplification", "deletion", "loss")) %>%
	mutate(
		VARIANT_ID = as.character(variant_id), VARIANT_GROUP = variant_groups, VARIANT_TYPE = variant_types, GENE = feature_name,
		ALTERATION = paste0(GENE, "_", ifelse(tolower(variant) %in% "amplification", "AMP", "DEL"))) %>%
	select(VARIANT_ID, VARIANT_TYPE, ALTERATION, VARIANT_GROUP)

# Fusions
var_fus <- var_data %>% filter(tolower(feature_type) %in% "fusion") %>% 
	mutate(
		VARIANT_ID = as.character(variant_id), VARIANT_GROUP = variant_groups, 
		VARIANT_TYPE = variant_types, GENE = feature_name, ALTERATION = feature_name) %>%
	select(VARIANT_ID, VARIANT_TYPE, ALTERATION, VARIANT_GROUP)

# Join all variants
all_var_data <- var_data %>% filter(!as.character(variant_id) %in% c(var_cna$VARIANT_ID, var_fus$VARIANT_ID)) %>% 
	mutate(
		VARIANT_ID = as.character(variant_id), VARIANT_GROUP = variant_groups, VARIANT_TYPE = variant_types, GENE = feature_name,
		ALTERATION = paste0(GENE, "_", gsub(" ", "_", gsub("\\(.*", "", variant)))) %>%
	select(VARIANT_ID, VARIANT_TYPE, ALTERATION, VARIANT_GROUP) %>%
	rbind(var_cna, var_fus) %>%
	mutate(ALTERATION = str_replace_all(ALTERATION, pattern = specific_biomarkers))

# Molecular profiles
molecular_data <- molecular_data %>%
	separate_rows(variant_ids, sep = ", ") %>%
	mutate(VARIANT_ID = variant_ids, MOLECULAR_ID = molecular_profile_id, EVIDENCE_SCORE = evidence_score) %>%
	left_join(all_var_data, by = "VARIANT_ID") %>% 
	select(MOLECULAR_ID, VARIANT_ID, VARIANT_TYPE, ALTERATION, VARIANT_GROUP, EVIDENCE_SCORE)

# Evidences
#	- Levels: A-D; Rating: 3:5 stars; Direction: Supports
#	- Discard evidences with NO tumor names and germline origin
#	- Oncogenic: small variants with oncogenic effect
#	- Clinical: 
#		- Type: Predictive, Prognostic, Diagnostic
#		- A drug is only associated with one ALTERATION-TUMOR-EFFECT-LEVEL entry (up-to-date data)
evidence <- evidence_data %>% 
	filter(
		evidence_level %in% evidence_levels & rating %in% evidence_rating & evidence_direction %in% "Supports" & 
		disease != "" & !variant_origin %in% c("Common Germline", "Rare Germline")) %>% 
	mutate(disease = tolower(disease)) %>%
	left_join(tumor_doid, by = c("doid" = "TUMOR_DOID")) %>%
	left_join(tumor_names, by = c("disease" = "TUMOR_NAME")) %>%
	mutate(
		MOLECULAR_ID = molecular_profile_id, EVIDENCE_TYPE = evidence_type, EVIDENCE_RATING = rating, 
		ORIGIN = variant_origin, EFFECT = significance, DRUG = therapies, 
		TUMOR_DOID = ifelse(!is.na(TUMOR_NAME), doid, TUMOR_DOID),
		TUMOR_NAME = ifelse(!is.na(TUMOR_NAME), TUMOR_NAME, disease),
		TUMOR_CODE = ifelse(!is.na(TUMOR_CODE.x), TUMOR_CODE.x, TUMOR_CODE.y),
		TOPNODE_NAME = ifelse(!is.na(TOPNODE_NAME.x), TOPNODE_NAME.x, TOPNODE_NAME.y),
		TOPNODE_DOID = ifelse(!is.na(TOPNODE_DOID.x), TOPNODE_DOID.x, TOPNODE_DOID.y),
		EVIDENCE_LEVEL = factor(evidence_level, levels = evidence_levels)) %>% 
	left_join(molecular_data, by = "MOLECULAR_ID", relationship = "many-to-many") %>%
	separate_rows(ALTERATION, sep = ",") %>% as.data.frame() %>%
	filter(!is.na(ALTERATION)) %>%
	select(
		ALTERATION, EVIDENCE_LEVEL, EVIDENCE_RATING, EVIDENCE_TYPE, EFFECT, DRUG, 
		TUMOR_NAME, TUMOR_DOID, TUMOR_CODE, TOPNODE_NAME, TOPNODE_DOID, MOLECULAR_ID, 
		VARIANT_ID, VARIANT_TYPE, VARIANT_GROUP, EVIDENCE_SCORE, ORIGIN) %>%
	arrange(EVIDENCE_LEVEL, desc(EVIDENCE_RATING), desc(EVIDENCE_SCORE), ALTERATION)

oncogenic <- evidence %>% filter(EFFECT %in% "Oncogenicity") %>% distinct(MOLECULAR_ID, .keep_all = TRUE) %>%
	select(
		ALTERATION, EVIDENCE_LEVEL, EVIDENCE_RATING, EVIDENCE_TYPE, EFFECT, MOLECULAR_ID, 
		VARIANT_ID, VARIANT_TYPE, VARIANT_GROUP, EVIDENCE_SCORE, ORIGIN)

clinical_prognostic <- evidence %>% filter(EVIDENCE_TYPE %in% "Prognostic" & EFFECT %in% prognostic_effect)
clinical_diagnostic <- evidence %>% filter(EVIDENCE_TYPE %in% "Diagnostic" & EFFECT %in% diagnostic_effect)
clinical_predictive <- evidence %>% filter(EVIDENCE_TYPE %in% "Predictive" & EFFECT %in% predictive_effect) %>%
	separate_rows(DRUG, sep = ",") %>% as.data.frame() %>%
	distinct(MOLECULAR_ID, EFFECT, TUMOR_DOID, TOPNODE_DOID, DRUG, .keep_all = TRUE) %>%
	group_by(MOLECULAR_ID, EVIDENCE_LEVEL, EFFECT, TUMOR_DOID, TOPNODE_DOID) %>%
	mutate(DRUG = paste(DRUG, collapse = ", ")) %>% as.data.frame() %>% 
	distinct(MOLECULAR_ID, EVIDENCE_LEVEL, EFFECT, TUMOR_DOID, TOPNODE_DOID, .keep_all = TRUE)
clinical <- rbind(clinical_predictive, clinical_prognostic, clinical_diagnostic) %>%
	arrange(EVIDENCE_LEVEL, desc(EVIDENCE_RATING), desc(EVIDENCE_SCORE))



# Save data ---------------------------------------------------------------

write.table(oncogenic, file = paste0("CIViC_", civic_date, "_", "Oncogenic.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(clinical, file = paste0("CIViC_", civic_date, "_", "Clinical.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
