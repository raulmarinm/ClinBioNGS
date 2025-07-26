#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate DNA manifest from IonTorrent results
# Last modification: 2023/03/01 by RM
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

# Set other variables
chr_levels <- paste0("chr", c(1:22, "X", "Y"))

# Set '<...>_Gen.amplicon.cov.xls' amplicon file
amplicon.file <- list.files(pattern = "Gen.amplicon.cov.xls$", full.names = TRUE, recursive = TRUE)[1]


# Manifest generation ---------------------------------------------------

# Load amplicons and create a BED file
#  - Split the attributes column (seb by ";") and keep only the gene name (sep by "GENE_ID=")
amplicon_bed <- read.table(amplicon.file, header = TRUE, sep = "\t", quote = "", comment.char = "") %>%
	separate(attributes, c("GENE"), sep = ";", extra = "drop") %>% separate(GENE, c(NA, "GENE"), sep = "GENE_ID=") %>%
	mutate(START = contig_srt-1, END = contig_end, CHROM = factor(contig_id, levels = chr_levels)) %>% 
	select(CHROM, START, END, GENE) %>% arrange(CHROM, START, END)


# Save data ---------------------------------------------------------------

write.table(amplicon_bed, file = output.file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
