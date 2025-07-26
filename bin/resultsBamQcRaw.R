#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of raw BAM QC
# Last modification: 2024/03/21 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

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
run <- args_list[["run"]]
sample <- args_list[["sample"]]
type <- args_list[["type"]]
isPaired <- !as.logical(args_list[["singleEnd"]])

# Set other variables
sample_id <- paste0(sample, "_", type)

# Set metrics files
alignmentSummary.file <- list.files(pattern = "AlignmentSummaryMetrics.txt")
insertSize.file <- list.files(pattern = "InsertSizeMetrics.txt")
onTargetReads.file <- list.files(pattern = "OnTargetReads.txt")


# Load data ---------------------------------------------------------------

onTargetReads <- read.table(onTargetReads.file, header = TRUE, sep = "\t") %>% mutate(SAMPLE = sample)
if (isPaired) {
	alignmentSummary <- read.table(alignmentSummary.file, header = TRUE, sep = "\t", comment.char = "#", nrows = 3) %>% 
		filter(CATEGORY == "PAIR") %>% mutate(SAMPLE = sample)
	insertSize <- read.table(insertSize.file, header = TRUE, sep = "\t", comment.char = "#", nrows = 1) %>% 
		mutate(SAMPLE = sample) %>% select(SAMPLE, MEDIAN_INSERT_SIZE)
} else {
	alignmentSummary <- read.table(alignmentSummary.file, header = TRUE, sep = "\t", comment.char = "#", nrows = 1) %>% 
		mutate(SAMPLE = sample)
	insertSize <- data.frame(SAMPLE = sample, MEDIAN_INSERT_SIZE = ".")
}


# Data processing ---------------------------------------------------------

bamqc_metrics <- onTargetReads %>%
	left_join(alignmentSummary, by = "SAMPLE") %>% left_join(insertSize, by = "SAMPLE") %>%
	mutate(
		RUN = run, TOTAL_READS = PF_READS, ALIGNED_READS = PF_READS_ALIGNED, 
		PCT_ALIGNED_READS = round(ALIGNED_READS*100/TOTAL_READS, 1), 
		PCT_ONTARGET_ALIGNED_READS = round(ONTARGET_READS*100/ALIGNED_READS, 1),
		HQ_ALIGNED_READS = PF_HQ_ALIGNED_READS, PCT_HQ_ALIGNED_READS = round(HQ_ALIGNED_READS*100/TOTAL_READS, 1)) %>% 
	select(
		RUN, SAMPLE, TOTAL_READS, ALIGNED_READS, PCT_ALIGNED_READS, ONTARGET_READS, PCT_ONTARGET_ALIGNED_READS, 
		HQ_ALIGNED_READS, PCT_HQ_ALIGNED_READS, MEDIAN_READ_LENGTH, MEDIAN_INSERT_SIZE)


# Save data ---------------------------------------------------------------

write.table(bamqc_metrics, file = paste0(sample, "_", type, "_RawBamQC.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)  

