#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of target genes
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
options(dplyr.summarise.inform = FALSE)
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(VariantAnnotation)))
suppressWarnings(suppressPackageStartupMessages(library(HGNChelper)))
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
panel <- args_list[["panel"]]
mane.file <- args_list[["mane"]]
cytoband.bed <- args_list[["cytoband"]]
target.bed <- args_list[["target"]]
type <- args_list[["type"]]

# Set other parameters
chr_names <- paste0("chr", c(1:22, "X", "Y"))


# Load data ---------------------------------------------------------------

# MANE genes
mane_genes <- fread(mane.file, header = TRUE, sep = "\t", quote = "", nThread = 1)

# UCSC CytoBand BED file
cytoband <- fread(cytoband.bed, header = FALSE, sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% 
	select(V1:V4) %>% filter(V1 %in% chr_names) %>% mutate(LABEL = sub("chr", "", paste0(V1, V4))) %>% as.data.frame()
colnames(cytoband) <- c("CHROM", "START", "END", "CYTOBAND", "LABEL")
cytoband_ranges <- makeGRangesFromDataFrame(cytoband, keep.extra.columns = TRUE, ignore.strand = TRUE, starts.in.df.are.0based = TRUE)

# Target BED file
target <- read.table(target.bed, header = FALSE, sep = "\t", col.names = c("CHROM", "START", "END", "GENE")) %>% 
	filter(!(is.na(GENE) | GENE %in% c("", ".", "iIndel", "iSNP", "chrY", "Y")))
target_ranges <- makeGRangesFromDataFrame(target, ignore.strand = TRUE, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)

# Target genes
genes <- sort(unique(target$GENE))


# Data processing ---------------------------------------------------------

# Update gene symbols
genes_updated <- suppressWarnings(checkGeneSymbols(genes)); colnames(genes_updated) <- c("OLD_GENE", "APPROVED", "NEW_GENE")

# Add CytoBand info
target$CYTOBAND <- unlist(lapply(seq_len(nrow(target)), function(x){
	paste(cytoband[overlapsAny(cytoband_ranges, target_ranges[x]), "LABEL"], collapse = "-")}))

# Create a table by gene
genes_cytoband <- target %>% group_by(GENE) %>% 
	summarise(CYTOBAND = paste(sort(unique(unlist(strsplit(CYTOBAND, "-")))), collapse = "-")) %>% as.data.frame()

# Add the Cytoband column to the gene file
genes_updated <- left_join(genes_updated, genes_cytoband, by = c("OLD_GENE" = "GENE"))

# Add the MANE info to genes
genes_updated <- left_join(genes_updated, mane_genes, by = c("NEW_GENE" = "GENE_SYMBOL"))


# Save data ---------------------------------------------------------------

write.table(genes_updated, file = paste0(panel, "_", type, "_genes.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

