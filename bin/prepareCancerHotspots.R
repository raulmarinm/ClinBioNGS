#!/usr/bin/env Rscript

# ============================================================================ #
# What: Prepare the results of Cancer Hotspot resource
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(HGNChelper)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
excel.file <- args_list[["excel"]]
out.file <- args_list[["outFile"]]


# Load data ---------------------------------------------------------------

hotspots_snv <- read_xls(excel.file, na = c("NA", ""), sheet = "SNV-hotspots")
hotspots_indel <- read_xls(excel.file, na = c("NA", ""), sheet = "INDEL-hotspots")


# Data processing ---------------------------------------------------------

# Update gene symbols
snv_genes_updated <- suppressWarnings(checkGeneSymbols(hotspots_snv$Hugo_Symbol))
hotspots_snv$gene <- snv_genes_updated$Suggested.Symbol
indel_genes_updated <- suppressWarnings(checkGeneSymbols(hotspots_indel$Hugo_Symbol))
hotspots_indel$gene <- indel_genes_updated$Suggested.Symbol

# Obtain AA changes and counts
snv_metrics <- hotspots_snv %>%
	separate(Variant_Amino_Acid, c("HOTSPOT_ALT", "HOTSPOT_MUT_CNT"), sep = ":") %>%
	mutate(
		HOTSPOT_GENE = ifelse(!is.na(gene), gene, Hugo_Symbol), HOTSPOT_POS = Amino_Acid_Position, 
		HOTSPOT_POS_CNT = Mutation_Count, HOTSPOT_MUT = ifelse(ref %in% "splice", HOTSPOT_POS, paste0(ref, HOTSPOT_POS, HOTSPOT_ALT))) %>%
	select(HOTSPOT_GENE, HOTSPOT_MUT, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, HOTSPOT_POS, HOTSPOT_ALT)
indel_metrics <- hotspots_indel %>%
	separate(Variant_Amino_Acid, c("HOTSPOT_ALT", "HOTSPOT_MUT_CNT"), sep = ":") %>%
	mutate(
		HOTSPOT_GENE = ifelse(!is.na(gene), gene, Hugo_Symbol), HOTSPOT_POS = Amino_Acid_Position, 
		HOTSPOT_POS_CNT = Mutation_Count, HOTSPOT_MUT = HOTSPOT_ALT) %>%
	select(HOTSPOT_GENE, HOTSPOT_MUT, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, HOTSPOT_POS, HOTSPOT_ALT)


hotspots_metrics <- rbind(snv_metrics, indel_metrics) %>% arrange(HOTSPOT_GENE) %>%
	mutate(MUT_ID = paste0(HOTSPOT_GENE, "_", HOTSPOT_MUT)) %>% select(MUT_ID, everything())


# Save data ---------------------------------------------------------------

write.table(hotspots_metrics, file = out.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

