#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of Network of Cancer Genes (NCG) file
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(HGNChelper)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
input.file <- args_list[["in"]]
output.file <- args_list[["out"]]


# Load data ---------------------------------------------------------------

cancer_drivers <- read.table(input.file, header = TRUE, sep = "\t", na.strings = c("NA", ""))

# Data processing ---------------------------------------------------------

# Update gene symbols
genes_updated <- suppressWarnings(checkGeneSymbols(cancer_drivers$symbol))
cancer_drivers$gene <- genes_updated$Suggested.Symbol

# Prepare list of gene role status
#	- Set ONCOGENE/TSG from NCG status
#	- Set number of evidences: 1-3
gene_role <- cancer_drivers %>%
	mutate(
		GENE = ifelse(!is.na(gene), gene, symbol), GENE_ENTREZ = entrez,
		ONCOGENE = as.logical(NCG_oncogene), TSG = as.logical(NCG_tsg),
		DRIVER = TRUE, CANONICAL_DRIVER = type %in% "Canonical Cancer Driver",
		cgc_oncogene = grepl("oncogene", cgc_annotation),
		cgc_tsg = grepl("TSG", cgc_annotation),
		vogelstein_oncogene = vogelstein_annotation %in% "Oncogene",
		vogelstein_tsg = vogelstein_annotation %in% "TSG",
		saito_oncogene = saito_annotation %in% "Oncogene",
		saito_tsg = saito_annotation %in% "TSG",
		across(c(ONCOGENE, TSG), ~ .x & !is.na(.x))) %>%
	group_by(GENE) %>% 
	mutate(
		ONCOGENE = any(ONCOGENE), TSG = any(TSG),
		CANONICAL_DRIVER = any(CANONICAL_DRIVER),
		cgc_oncogene = any(cgc_oncogene), cgc_tsg = any(cgc_tsg), 
		vogelstein_oncogene = any(vogelstein_oncogene), vogelstein_tsg = any(vogelstein_tsg),
		saito_oncogene = any(saito_oncogene), saito_tsg = any(saito_tsg)) %>%
	ungroup() %>%
	distinct(GENE, .keep_all = TRUE) %>%
	rowwise() %>%
	mutate(
		ONCOGENE_EVIDENCES = sum(cgc_oncogene, vogelstein_oncogene, saito_oncogene),
		TSG_EVIDENCES = sum(cgc_tsg, vogelstein_tsg, saito_tsg)) %>%
	as.data.frame() %>%
	select(GENE, GENE_ENTREZ, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, ONCOGENE_EVIDENCES, TSG_EVIDENCES) %>%
	arrange(GENE)


# Save data ---------------------------------------------------------------

write.table(gene_role, file = output.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
