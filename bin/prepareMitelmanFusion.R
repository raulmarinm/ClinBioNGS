#!/usr/bin/env Rscript

# ============================================================================ #
# What: Prepare the MitelmanDB fusion file
# Last modification: 2024/01/30 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
options(dplyr.summarise.inform = FALSE)
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
input.file <- args_list[["in"]]
output.file <- args_list[["out"]]


# Load data ---------------------------------------------------------------

mitelman <- read.csv(input.file, header = TRUE)


# Data processing ---------------------------------------------------------

# Keep fusions and count them
mitelman <- mitelman %>% 
	filter(grepl("::", GeneShort)) %>%
	mutate(FUSION = gsub("\\+", "", GeneShort)) %>%
    group_by(FUSION) %>% summarise(MitelmanDB_COUNT = n()) %>% arrange(desc(MitelmanDB_COUNT))


# Save data ---------------------------------------------------------------

write.table(mitelman, file = output.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

