#!/usr/bin/env Rscript

# ============================================================================ #
# What: Arm-level Somatic Copy-number Events in Targeted Sequencing (ASCETS)
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(VariantAnnotation)))
suppressPackageStartupMessages(suppressWarnings(library(DBI)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

# Load arguments
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments (--arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args_df <- as.data.frame(do.call("rbind", parseArgs(args)))
args_list <- as.list(as.character(args_df$V2))
names(args_list) <- args_df$V1

# Set arguments
run <- args_list[["run"]]
sample_name <- args_list[["sample"]]
type <- args_list[["type"]]
sample.file <- args_list[["samplesFile"]]
resources.file <- args_list[["resources"]]
segments.file <- args_list[["segments"]]
sex.file <- args_list[["sex"]]
arm.file <- args_list[["armCoordinates"]]
min_boc <- args_list[["minBoc"]]
call_threshold <- args_list[["callingThreshold"]]
fraction_threshold <- args_list[["fractionThreshold"]]

# Set other variables
sample_id <- paste0(sample_name, "_", type)


# Load data ---------------------------------------------------------------

# Sample info
sample_info <- read.table(sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Sample == sample_name)
if(sample_info$Sex != "."){
	is_male <- sample_info$Sex %in% "Male"
} else {
	sex_df <- read.table(sex.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
	is_male <- sex_df$sex %in% "Male"
}

# Ascets resources
source(resources.file)

# Segments (remove antitarget, iSNP/iIndel and chrY segments)
segments <- read.table(segments.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
	filter(!(gene %in% c("-", "iSNP", "iIndel") | chromosome %in% "chrY")) %>% rowwise() %>% 
	mutate(sample = sample_name, start = start+1, log2 = ifelse(is_male & chromosome %in% "chrX", log2+1, log2)) %>%
	select(sample, chromosome, start, end, probes, log2)

# Arm genomic coordinates (remove chrY)
arm_coord <- read.table(arm.file, header = TRUE) %>% filter(chrom != "chrY")


# Data processing ---------------------------------------------------------

ascets_output <- ascets(
	cna = segments, cytoband = arm_coord, name = sample_id, min_boc = min_boc, 
	threshold = call_threshold, alteration_threshold = fraction_threshold
)
aux_cn <- as.data.frame(t(ascets_output$calls %>% select(-sample))) %>% mutate(CNA = ifelse(V1 == "NC", "CONFLICT", V1)) %>% select(CNA)
aux_cn$arm <- rownames(aux_cn)
aux_log2 <- as.data.frame(t(ascets_output$weight_ave %>% select(-sample)))
aux_log2$arm <- rownames(aux_log2)
aux_log2 <- aux_log2 %>% mutate(LOG2 = ifelse(is_male & arm %in% c("Xp", "Xq"), V1-1, V1)) %>% select(LOG2, arm)

arm_results <- arm_coord %>% left_join(aux_cn, by = "arm") %>% left_join(aux_log2, by = "arm") %>% 
	mutate(
		RUN = run, SAMPLE = sample_name, CHROM = chrom, START = start+1, END = end, 
		CHR_ARM = arm, SIDE = side, VAR = paste0(CHR_ARM, "_", CNA)) %>% 
	select(RUN, SAMPLE, VAR, CHROM, START, END, CHR_ARM, SIDE, CNA, LOG2)
arm_results[is.na(arm_results) | arm_results == "NA" | arm_results == "NaN" | arm_results == ""] <- "."


# Save data ---------------------------------------------------------------

suppressWarnings(dir.create("Ascets_Reports"))
write_outputs_to_file(ascets_output, location = "Ascets_Reports/")

write.table(arm_results, file = paste0(sample_id, "_CNA_arm.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

