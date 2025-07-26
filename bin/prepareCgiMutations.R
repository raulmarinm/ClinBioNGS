#!/usr/bin/env Rscript

# ============================================================================ #
# What: Prepare the CGI mutations file
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(HGNChelper)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))


# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
chain.file <- args_list[["chain"]]
input.file <- args_list[["in"]]
output.file <- args_list[["out"]]


# Load data ---------------------------------------------------------------

cgi_data <- read.table(input.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)


# Data processing ---------------------------------------------------------

# Discard OncoKB resource and germline context
# Separate multiple gDNAs into multiple rows (sep="__")
mutations <- cgi_data %>% filter(source != "OncoKB" & context %in% "somatic") %>%
	separate_rows(gdna, sep = "__") %>%
	mutate(HGVSg_HG19 = gdna) %>%
	separate(gdna, c("CHROM", "gdna"), sep = ":") %>%
	mutate(
		GENE = gene, MUTATION = gsub("p\\.", "", protein), 
		AA_POS = as.integer(gsub("\\D", "", gsub("\\*.*", "", gsub("_.*", "", MUTATION)))),
		TRANSCRIPT = transcript, gdna = gsub("g\\.", "", gdna)) %>%
	select(GENE, MUTATION, AA_POS, TRANSCRIPT, HGVSg_HG19, CHROM, gdna) %>%
	as.data.frame()

# Update gene symbols
genes_updated <- suppressWarnings(checkGeneSymbols(mutations$GENE))
mutations$gene <- genes_updated$Suggested.Symbol
mutations <- mutations %>% mutate(GENE = ifelse(!is.na(gene), gene, GENE)) %>% select(-gene)

# Separate SNVs from InDels and extract genomic regions gDNA
snvs <- mutations %>% filter(grepl(">", mutations$gdna)) %>%
	mutate(POS = as.integer(gsub("\\D", "", gdna)), ALT = gsub("\\d", "", gdna))
indels <- mutations %>% filter(!grepl(">", mutations$gdna)) %>%
	mutate(ALT = gsub("\\d+_\\d+|\\d+", "", gdna), POS = gsub("del\\D+|dup\\D+|ins\\D+", "", gdna)) %>%
	separate(POS, c("START", "END"), sep = "_", fill = "right") %>%
	mutate(aux_end = ifelse(is.na(END), START, END))

# Convert hg19 to hg38 coordinates
snvs_hg38 <- as.data.frame(unlist(
	liftOver(GRanges(paste0(snvs$CHROM, ":", snvs$POS), HGVSg_HG19 = snvs$HGVSg_HG19), import.chain(chain.file)))) %>%
	select(HGVSg_HG19, start)
indels_hg38_start <- as.data.frame(unlist(
	liftOver(GRanges(paste0(indels$CHROM, ":", indels$START), HGVSg_HG19 = indels$HGVSg_HG19), import.chain(chain.file)))) %>%
	select(HGVSg_HG19, start)
indels_hg38_end <- as.data.frame(unlist(
	liftOver(GRanges(paste0(indels$CHROM, ":", indels$aux_end), HGVSg_HG19 = indels$HGVSg_HG19), import.chain(chain.file)))) %>%
	select(HGVSg_HG19, end)

# Add hg38 HGVSg
snvs <- snvs %>% left_join(snvs_hg38, by = "HGVSg_HG19") %>%
	mutate(HGVSg = ifelse(is.na(start), NA, paste0(CHROM, ":g.", start, ALT))) %>%
	select(GENE, MUTATION, AA_POS, TRANSCRIPT, HGVSg, HGVSg_HG19)
indels <- indels %>% left_join(indels_hg38_start, by = "HGVSg_HG19") %>% left_join(indels_hg38_end, by = "HGVSg_HG19") %>%
	mutate(
		HGVSg = ifelse(is.na(start) | is.na(end), NA, ifelse(
			is.na(END), paste0(CHROM, ":g.", start, ALT), paste0(CHROM, ":g.", start, "_", end, ALT)))) %>%
	select(GENE, MUTATION, AA_POS, TRANSCRIPT, HGVSg, HGVSg_HG19)

# Join snvs and indels
mutations_ok <- rbind(snvs, indels) %>% arrange(GENE, MUTATION) %>% 
	mutate(MUT_ID = paste0(GENE, "_", MUTATION)) %>% select(MUT_ID, everything())


# Save data ---------------------------------------------------------------

write.table(mutations_ok, file = output.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
