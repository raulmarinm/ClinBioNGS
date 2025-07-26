#!/usr/bin/env Rscript

# ============================================================================ #
# What: Annotate genes to the BED file
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
options(dplyr.summarise.inform = FALSE)
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(VariantAnnotation)))
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
target.bed <- args_list[["bed"]]
genes.file <- args_list[["genes"]]
out.file <- args_list[["outFile"]]


# Load data ---------------------------------------------------------------

# Target BED file
target <- read.table(target.bed, header = FALSE, sep = "\t")[, 1:3]
target_ranges <- makeGRangesFromDataFrame(
	target, seqnames.field = "V1", start.field = "V2", end.field = "V3", 
	ignore.strand = TRUE, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)

# Genes
genes <- fread(genes.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% as.data.frame()
genes_ranges <- makeGRangesFromDataFrame(genes, ignore.strand = TRUE)


# Data processing ---------------------------------------------------------

# Annotate gene
target$GENE <- unlist(lapply(seq_len(nrow(target)), function(x){
	paste(genes[overlapsAny(genes_ranges, target_ranges[x]), "GENE_SYMBOL"], collapse = "_")}))
target[target == ""] <- "."


# Save data ---------------------------------------------------------------

write.table(target, file = out.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

