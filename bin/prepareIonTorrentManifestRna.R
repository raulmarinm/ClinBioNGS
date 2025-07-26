#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate RNA manifest from IonTorrent results
# Last modification: 2023/03/12 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
output.file <- args_list[["output"]]
padding <- round(as.numeric(args_list[["padding"]]))

# Set other variables
chr_levels <- paste0("chr", c(1:22, "X", "Y"))

# Set '<...>_Gen.amplicon.cov.xls' amplicon file
amplicon.file <- list.files(pattern = "Gen.amplicon.cov.xls$", full.names = TRUE, recursive = TRUE)[1]

# Set 'Fusion.tsv' file
fusion.file <- list.files(pattern = "Fusion.tsv$", full.names = TRUE, recursive = TRUE)[1]


# Manifest generation ---------------------------------------------------

# Load DNA amplicons and create a BED file
#  - Split the attributes column (seb by ";") and keep only the gene name (sep by "GENE_ID=")
dna_amplicon_bed <- read.table(amplicon.file, header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
	separate(attributes, c("GENE"), sep = ";", extra = "drop") %>% separate(GENE, c(NA, "GENE"), sep = "GENE_ID=") %>%
	mutate(START = contig_srt-1, END = contig_end, CHROM = factor(contig_id, levels = chr_levels)) %>% 
	select(CHROM, START, END, GENE)

# Load the fusion file, extract the regions and genes, and create a BED file
#  - Add padding to each position
rna_non_fusions <- read.table(fusion.file, header = TRUE, sep = "\t", quote = "", comment.char = "") %>% 
	filter(Type != "Fusion") %>% 
	separate(Locus, c("LocusA", "LocusB"), sep = "-", fill = "right")
rna_fusions <- read.table(fusion.file, header = TRUE, sep = "\t", quote = "", comment.char = "") %>% 
	filter(Type == "Fusion") %>%
	separate(Genes..Exons., c("GeneA", "GeneB"), sep = " - ", fill = "right") %>%
	separate(Locus, c("LocusA", "LocusB"), sep = "-", fill = "right") %>% rowwise() %>%
	mutate(LOCUS = ifelse(
		grepl(paste0(Oncomine.Driver.Gene, "\\("), GeneA), LocusA, ifelse(
			grepl(paste0(Oncomine.Driver.Gene, "\\("), GeneB), LocusB, NA))) %>% 
	as.data.frame()

rna_amplicon_bed <- data.frame(
	LOCUS = c(rna_non_fusions$LocusA, rna_non_fusions$LocusB, rna_fusions$LOCUS), 
	GENE = c(rna_non_fusions$Oncomine.Driver.Gene, rna_non_fusions$Oncomine.Driver.Gene, rna_fusions$Oncomine.Driver.Gene)) %>%
	filter(!is.na(LOCUS)) %>% mutate(LOCUS = gsub(" ", "", LOCUS)) %>% distinct(LOCUS, .keep_all = TRUE) %>%
	separate(LOCUS, c("CHROM", "POS"), sep = ":") %>%
	mutate(START = as.integer(POS) - 1 - padding, END = as.integer(POS) + padding, CHROM = factor(CHROM, levels = chr_levels)) %>% 
	select(CHROM, START, END, GENE)

# Merge DNA and RNA regions for RNA manifest
rna_amplicon_bed <- rbind(dna_amplicon_bed, rna_amplicon_bed) %>% arrange(CHROM, START, END)


# Save data ---------------------------------------------------------------

write.table(rna_amplicon_bed, file = output.file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
