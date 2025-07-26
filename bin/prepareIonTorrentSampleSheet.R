#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate sampleSheet from IonTorrent results
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
run <- args_list[["run"]]
output.file <- args_list[["output"]]

# Set 'Info.csv' files
info.files <- list.files(pattern = "Info.csv$", full.names = TRUE, recursive = TRUE)


# SampleSheet generation ---------------------------------------------------

samplesheet <- data.frame(V1 = c("[Header]", "RunName", "", "[Data]"), V2 = c("", run, "", "" ), V3 = rep("", 4))
samplesheet <- rbind(samplesheet, c("Sample_ID", "Barcode", "Sample_Type"))

for (file in info.files) {

	# Load the input csv file and remove blank rows
	aux_info <- read.csv(file, header = FALSE) %>% filter(V1 != "")

	# Parse sample info
	aux_list <- as.list(as.character(aux_info$V2))
	names(aux_list) <- aux_info$V1

	# Set sample name
	sample_name <- aux_list[["Sample Name"]]

	# Set DNA/RNA barcodes
	dna_start <- grep("^DNA$", aux_list)
	rna_start <- grep("^RNA$", aux_list)
	if(dna_start < rna_start){
		dna_barcode <- aux_list[dna_start:rna_start][["Barcode Id"]]
		rna_barcode <- aux_list[rna_start:length(aux_list)][["Barcode Id"]]
	} else {
		dna_barcode <- aux_list[dna_start:length(aux_list)][["Barcode Id"]]
		rna_barcode <- aux_list[rna_start:dna_start][["Barcode Id"]]
	}

	# Add sample info
	samplesheet <- rbind(samplesheet, c(paste0(sample_name, "_DNA"), dna_barcode, "DNA"))
	samplesheet <- rbind(samplesheet, c(paste0(sample_name, "_RNA"), rna_barcode, "RNA"))
	
	rm(aux_info, aux_list, sample_name, dna_start, rna_start, dna_barcode, rna_barcode)
}
rm(file)


# Save data ---------------------------------------------------------------

write.table(samplesheet, file = output.file, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
