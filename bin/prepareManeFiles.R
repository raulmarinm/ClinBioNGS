#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of MANE files
# Last modification: 2024/01/29 by RM
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
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
mane.gtf <- args_list[["gtf"]]
mane_summary.file <- args_list[["summary"]]
genes.file <- args_list[["genes"]]
exons.file <- args_list[["exons"]]
introns.file <- args_list[["introns"]]
coding.file <- args_list[["coding"]]

# Set other parameters
chr_names <- paste0("chr", c(1:22, "X", "Y")); names(chr_names) <- paste0("NC_0000", c(paste0("0", 1:9), 10:24))


# Load data ---------------------------------------------------------------

# Summary of MANE info (only MANE SELECT)
mane_summary <- fread(mane_summary.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% 
    filter(MANE_status == "MANE Select" & grepl("NC_", GRCh38_chr))

# Genomic data
genomic_data <- as.data.frame(rtracklayer::import(mane.gtf)) %>% filter(gene_id %in% mane_summary$Ensembl_Gene) %>% 
    mutate(GENE_SYMBOL = gene_name, EXON = as.integer(exon_number))
genomic_genes <- unique(genomic_data$GENE_SYMBOL); names(genomic_genes) <- unique(genomic_data$gene_id)


# Data processing ---------------------------------------------------------

# GENES 
gene_data <- mane_summary %>%
    separate(GRCh38_chr, c("GRCh38_chr", NA), sep = "\\.") %>%
    separate(`#NCBI_GeneID`, c(NA, "GENE_NCBI"), sep = ":") %>%
    separate(HGNC_ID, c(NA, "GENE_HGNC"), sep = ":", fill = "right") %>%
    mutate(
        GENE_SYMBOL = genomic_genes[Ensembl_Gene], GENE_ENSEMBL = Ensembl_Gene, GENE_NAME = name, TRANSCRIPT_REFSEQ = RefSeq_nuc, TRANSCRIPT_ENSEMBL = Ensembl_nuc,
        PROTEIN_REFSEQ = RefSeq_prot, PROTEIN_ENSEMBL = Ensembl_prot, CHROM = chr_names[GRCh38_chr], START = chr_start, END = chr_end, 
        STRAND = chr_strand) %>%
    select(
        GENE_SYMBOL, CHROM, START, END, STRAND, GENE_NCBI, GENE_ENSEMBL, GENE_HGNC, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL, PROTEIN_REFSEQ, 
        PROTEIN_ENSEMBL, GENE_NAME) %>%
    arrange(GENE_SYMBOL)

# EXONS
exon_data <- genomic_data %>% filter(type == "exon" & transcript_id %in% gene_data$TRANSCRIPT_ENSEMBL) %>% 
	left_join(gene_data, by = "GENE_SYMBOL") %>% 
    mutate(CHROM = seqnames, START = start, END = end, STRAND = strand) %>% arrange(GENE_SYMBOL, EXON) %>% 
    select(GENE_SYMBOL, EXON, CHROM, START, END, STRAND, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL)

# INTRONS
intron_data <- matrix(data = NA, nrow = 0, ncol = 1)
for(gene_name in unique(exon_data$GENE_SYMBOL)){
    aux_df <- exon_data %>% filter(GENE_SYMBOL == gene_name)
    rownames(aux_df) <- aux_df$EXON
    introns_n <- max(aux_df$EXON)-1
    aux_chr <- unique(aux_df$CHROM)
    aux_strand <- unique(aux_df$STRAND)
    aux_transcript_refseq <- unique(aux_df$TRANSCRIPT_REFSEQ)
    aux_transcript_ensembl <- unique(aux_df$TRANSCRIPT_ENSEMBL)
    if(!introns_n){next}
    for(n in seq_len(introns_n)){
        aux_intron <- data.frame(
            GENE_SYMBOL = gene_name, INTRON = n, CHROM = aux_chr, 
            START = ifelse(aux_strand == "-", aux_df[n+1, "END"]+1, aux_df[n, "END"]+1),
            END = ifelse(aux_strand == "-", aux_df[n, "START"]-1, aux_df[n+1, "START"]-1),
            STRAND = aux_strand, TRANSCRIPT_REFSEQ = aux_transcript_refseq, TRANSCRIPT_ENSEMBL = aux_transcript_ensembl)
        intron_data <- rbind(intron_data, aux_intron)
        rm(aux_intron)
    }
    rm(n, aux_df, introns_n, aux_chr, aux_strand, aux_transcript_refseq, aux_transcript_ensembl)
}
rm(gene_name)


# CODING
coding_data <- genomic_data %>% filter(type == "CDS" & transcript_id %in% gene_data$TRANSCRIPT_ENSEMBL) %>% 
	left_join(gene_data, by = "GENE_SYMBOL") %>% 
    mutate(CHROM = seqnames, START = start, END = end, STRAND = strand) %>% arrange(GENE_SYMBOL, EXON) %>% 
    select(GENE_SYMBOL, EXON, CHROM, START, END, STRAND, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL)


# Save data ---------------------------------------------------------------

write.table(gene_data, file = genes.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(exon_data, file = exons.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(intron_data, file = introns.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(coding_data, file = coding.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
