#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate arm genomic coordinates from the UCSC cytoband file
# Last modification: 2024/01/09 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
options(dplyr.summarise.inform = FALSE)
.libPaths("opt/R/lib/R/library")

# Load libraries
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
cytoband.bed <- args_list[["cytoband"]]
output.file <- args_list[["fileName"]]

# Set other parameters
chr_names <- paste0("chr", c(1:22, "X", "Y"))


# Load data ---------------------------------------------------------------

# UCSC CytoBand BED file
cytoband <- fread(cytoband.bed, header = FALSE, sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% 
	select(V1:V4) %>% filter(V1 %in% chr_names)

# Data processing ---------------------------------------------------------

# Generate arm genomic coordinates
arm_coord <- cytoband %>% separate(V4, c("side", NA), sep = 1) %>% mutate(arm = paste0(sub("chr", "", paste0(V1, side)))) %>% 
	group_by(arm) %>% summarise(chrom = unique(V1), start = min(V2), end = max(V3), side = unique(side)) %>%
	select(chrom, arm, side, start, end) %>%
	mutate(chrom = factor(chrom, levels = chr_names)) %>% arrange(chrom)

# Save data ---------------------------------------------------------------

write.table(arm_coord, file = output.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
