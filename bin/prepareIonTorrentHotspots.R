#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate DNA hotspots BED from IonTorrent results
# Last modification: 2023/03/01 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
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
chain.file <- args_list[["chain"]]

# Set other variables
chr_levels <- paste0("chr", c(1:22, "X", "Y"))

# Set 'Snvindel.tsv' file
snvindel.file <- list.files(pattern = "Snvindel.tsv$", full.names = TRUE, recursive = TRUE)[1]


# Manifest generation ---------------------------------------------------

# Load Snvindel file and extract the HOTSPOT regions 
hotspots <- read.table(snvindel.file, header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
	filter(Variant.ID != ".") %>% separate(Locus, c("CHROM", "POS")) %>% rowwise() %>%
	mutate(START = as.integer(POS), END = as.integer(POS) + nchar(Ref) - 1, GENE = Gene, CHROM = factor(CHROM, levels = chr_levels)) %>%
	as.data.frame() %>% arrange(CHROM, START, END) %>% distinct(CHROM, START, END, .keep_all = TRUE) %>% select(CHROM, START, END, GENE)
hotspots_ranges <- makeGRangesFromDataFrame(hotspots, ignore.strand = TRUE, keep.extra.columns = TRUE)

# Convert hg19 to hg38 coordinates and create a BED file
hotspots_hg38_bed <- as.data.frame(unlist(liftOver(hotspots_ranges, import.chain(chain.file)))) %>%
	mutate(start = start-1) %>% select(seqnames, start, end, GENE)


# Save data ---------------------------------------------------------------

write.table(hotspots_hg38_bed, file = output.file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
