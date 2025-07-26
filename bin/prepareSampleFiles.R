#!/usr/bin/env Rscript

# ============================================================================ #
# What: Generate sample files from the sample sheet
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
sample_sheet.file <- args_list[["sampleSheet"]]
sample_info.file <- args_list[["sampleInfo"]]
tumor_names.file <- args_list[["tumorNames"]]
panel <- args_list[["panel"]]
run_name <- args_list[["run"]]
type <- args_list[["type"]]
isIontorrentSamplesheet <- as.logical(args_list[["isIontorrentSamplesheet"]])
isBcl <- as.logical(args_list[["isBcl"]])
isUmi <- as.logical(args_list[["isUmi"]])

# Other variables
reads_header <- c("[Reads]")
data_header <- c("[Data]", "[data]", "[BCLConvert_Data]")
settings_header <- c("[Settings]", "[settings]", "[BCLConvert_Settings]")


# Load data ---------------------------------------------------------------

# Tumor names
tumor_names <- read.table(tumor_names.file, header = TRUE, sep = ";", quote = "\"", na.strings = c("NA", ".", ""), fill = TRUE, colClasses = "character")
tumor_doid <- tumor_names %>% filter(!is.na(TUMOR_DOID)) %>% distinct(TUMOR_DOID, .keep_all = TRUE)
tumor_codes <- tumor_names %>% filter(is.na(TUMOR_DOID)) %>% distinct(TUMOR_CODE, .keep_all = TRUE) %>% select(-TUMOR_DOID)

# Sample info
sample_info <- read.table(sample_info.file, header = TRUE, sep = ";", quote = "\"", na.strings = c("NA", ".", ""), fill = TRUE, colClasses = "character") %>% 
    filter(Run %in% run_name) %>% 
    left_join(tumor_doid, by = c("DOID" = "TUMOR_DOID")) %>%
    left_join(tumor_names, by = c("Tumor" = "TUMOR_NAME")) %>%
    left_join(tumor_codes, by = c("Tumor" = "TUMOR_CODE")) %>%
    mutate(
        Sample = as.character(Sample),
        Sex = ifelse(tolower(Sex) %in% c("m", "male"), "Male", ifelse(tolower(Sex) %in% c("f", "female"), "Female", Sex)), 
        TumorName = ifelse(!is.na(TUMOR_NAME.x), TUMOR_NAME.x, ifelse(!is.na(TUMOR_NAME.y), TUMOR_NAME.y, Tumor)),
        TumorDOID = ifelse(!is.na(DOID), DOID, TUMOR_DOID),
        TumorCode = ifelse(!is.na(TUMOR_CODE.x), TUMOR_CODE.x, ifelse(!is.na(TUMOR_CODE.y), TUMOR_CODE.y, Tumor)),
        TopNodeName = ifelse(!is.na(TOPNODE_NAME.x), TOPNODE_NAME.x, ifelse(!is.na(TOPNODE_NAME.y), TOPNODE_NAME.y, TOPNODE_NAME)),
        TopNodeDOID = ifelse(!is.na(TOPNODE_DOID.x), TOPNODE_DOID.x, ifelse(!is.na(TOPNODE_DOID.y), TOPNODE_DOID.y, TOPNODE_DOID))) %>% 
    select(Sample, Sex, Age, Purity, Whitelist, TumorName, TumorDOID, TumorCode, TopNodeName, TopNodeDOID, LibraryDna, LibraryRna)
sample_info[is.na(sample_info)] <- "."

# Load the sample sheet and Remove blank columns
sampleSheet <- read.csv(sample_sheet.file, header = FALSE)
sampleSheet <- sampleSheet[, colSums(!is.na(sampleSheet)) > 0 ]
sampleSheet[nrow(sampleSheet)+1, ] <- ""


# Data processing ---------------------------------------------------------

# Get header info
header_pos <- grep(TRUE, sampleSheet[,1] %in% "[Header]")[1]
if (is.na(header_pos)) {
    header_info <- sampleSheet[0, 1:2]
} else {
    header_end <- header_pos + grep("^$|^\\[.*\\]$", sampleSheet[header_pos+1:nrow(sampleSheet), 1])[1] - 1
    header_info <- sampleSheet[header_pos:header_end, 1:2]
}

# Get data info
data_pos <- grep(TRUE, sampleSheet[,1] %in% data_header)[1]
data_end <- data_pos + grep("^$|^\\[.*\\]$", sampleSheet[data_pos+1:nrow(sampleSheet), 1])[1] - 1
data_info <- sampleSheet[data_pos:data_end, ]
sample_table <- data_info[-1, ]
sample_table <- as.data.frame(sample_table[, colSums(sample_table != "") > 0])
colnames(sample_table) <- sample_table[1,]
sample_table <- sample_table %>% filter(Sample_ID != "")

# Split "_[D/R]NA" from sample ID (If no sample type is specified, DNA by default)
if(! "Sample_Type" %in% colnames(sample_table)){sample_table$Sample_Type <- str_extract(sample_table$Sample_ID, "[DR]NA$")}
sample_table <- sample_table %>% mutate(
    Sample_ID = ifelse(is.na(Sample_Type) & Sample_ID != "Sample_ID", paste0(Sample_ID, "_DNA"), Sample_ID),
    Sample_Type = ifelse(is.na(Sample_Type), "DNA", Sample_Type),
    Sample = gsub(".[DR]NA$", "", sample_table$Sample_ID))

# Split Sample table by sample type
dna_table <- sample_table[-1, ] %>% filter(toupper(Sample_Type) %in% "DNA")
rna_table <- sample_table[-1, ] %>% filter(toupper(Sample_Type) %in% "RNA")

# Create a specific sample sheet for BCL convert (BCL input data type)
if(isBcl){

    aux_sampleSheet <- header_info

    # Add reads info
    if(any(reads_header %in% sampleSheet[,1])){
        reads_pos <- grep(TRUE, sampleSheet[,1] %in% reads_header)[1]
        reads_end <- reads_pos + grep("^$|^\\[.*\\]$", sampleSheet[reads_pos+1:nrow(sampleSheet), 1])[1] - 1
        reads_info <- sampleSheet[reads_pos:reads_end, 1:2]

        aux_sampleSheet <- rbind(aux_sampleSheet, reads_info)
    }

    # Add settings info
    if(any(settings_header %in% sampleSheet[,1])){
        settings_pos <- grep(TRUE, sampleSheet[,1] %in% settings_header)[1]
        settings_end <- settings_pos + grep("^$|^\\[.*\\]$", sampleSheet[settings_pos+1:nrow(sampleSheet), 1])[1] - 1
        settings_info <- sampleSheet[settings_pos:settings_end, 1:2]

        aux_sampleSheet <- rbind(aux_sampleSheet, settings_info)
    }
    
    # Add data info
    aux_data <- sample_table[, colnames(sample_table) %in% c("Sample_ID", "index", "index2", "Index", "Index2")]
    colnames(aux_data) <- paste0("V", seq_len(ncol(aux_data)))

    sampleSheet <- bind_rows(aux_sampleSheet, data_info[1, 1:2], aux_data)
    sampleSheet[is.na(sampleSheet)] <- ""

}

### DNA processing ###
if(type %in% "DNA"){
    # Create specific DNA files and save them in the Run directory
    dna_info <- dna_table %>% select(Sample) %>% left_join(sample_info, by = "Sample") %>% 
        mutate(Library = !tolower(LibraryDna) %in% c("no", "false")) %>% select(-LibraryDna, -LibraryRna)
    sampleSheet_dna <- sampleSheet %>% filter(!(V1 %in% rna_table$Sample_ID))

    # Ion Torrent samplesheet
    if(isIontorrentSamplesheet){dna_info$Barcode <- dna_table$Barcode}
    
    # TSO500 samplesheet
    if(panel %in% "tso500" & isBcl){
        if(!isUmi){
            OverrideCycles <- (settings_info %>% filter(V1 %in% "OverrideCycles"))[,2]
            sampleSheet_dna[,2] <- sub(OverrideCycles, "Y101;I8;I8;Y101", sampleSheet_dna[,2])
        }
    }

    write.table(sampleSheet_dna, file = "SampleSheet_DNA.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

    dna_info[dna_info == "" | is.na(dna_info)] <- "."
    write.table(dna_info, file = "SamplesDNA.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}



### RNA processing ###
if(type %in% "RNA"){
    # Create specific RNA files and save them in the Run directory
    rna_info <- rna_table %>% select(Sample) %>% left_join(sample_info, by = "Sample") %>% 
        mutate(Library = !tolower(LibraryRna) %in% c("no", "false")) %>% select(-LibraryDna, -LibraryRna)
    sampleSheet_rna <- sampleSheet %>% filter(!(V1 %in% dna_table$Sample_ID))

    # Ion Torrent samplesheet
    if(isIontorrentSamplesheet){rna_info$Barcode <- rna_table$Barcode}

    # TSO500 samplesheet
    if(panel %in% "tso500" & isBcl){
        if(!isUmi){
            OverrideCycles <- (settings_info %>% filter(V1 %in% "OverrideCycles"))[,2]
            sampleSheet_rna[,2] <- sub(OverrideCycles, "Y101;I8;I8;Y101", sampleSheet_rna[,2])
        }
    }
        
    write.table(sampleSheet_rna, file = "SampleSheet_RNA.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

    rna_info[rna_info == "" | is.na(rna_info)] <- "."
    write.table(rna_info, file = "SamplesRNA.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}


