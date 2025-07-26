#!/usr/bin/env Rscript

# ============================================================================ #
# What: Calculate TMB from small variant annotation
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(VariantAnnotation)))
suppressPackageStartupMessages(suppressWarnings(library(openxlsx)))
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
sample <- args_list[["sample"]]
type <- args_list[["type"]]
vars.file <- args_list[["variants"]]
perBaseCoverage.bed <- args_list[["perBaseCovBed"]]
target.bed <- args_list[["targetBed"]]
mane_coding.file <- args_list[["maneCoding"]]
problematic_regions.bed <- args_list[["problematicRegionsBed"]]
min_ad <- as.numeric(args_list[["minAD"]])
min_dp <- as.numeric(args_list[["minDP"]])
min_vaf <- as.numeric(args_list[["minVAF"]])
somatic_max_pvaf <- as.numeric(args_list[["somaticMaxpVAF"]])
somatic_max_vaf <- as.numeric(args_list[["somaticMaxVAF"]])

# Set other variables
chr_names <- paste0("chr", c(1:22, "X", "Y"))
sample_id <- paste0(sample, "_", type)


# Load data ---------------------------------------------------------------

# Variants (exclude LowCallers)
variants <- read.table(vars.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% filter(!LowCallers)
variants_ranges <- makeGRangesFromDataFrame(variants)

# Per-base coverage
perBaseCov_ranges <- makeGRangesFromDataFrame(fread(
	perBaseCoverage.bed, header = FALSE, col.names = c("CHROM", "START", "END", "DEPTH"),
	sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% filter(CHROM %in% chr_names),
	starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)

# Target regions
target_ranges <- makeGRangesFromDataFrame(
	read.table(target.bed, header = FALSE, sep = "\t", col.names = c("CHROM", "START", "END", "GENE")), starts.in.df.are.0based = TRUE)
target_perbase_ranges <- makeGRangesFromDataFrame(
	target_ranges %>% as.data.table() %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
	mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.table(), keep.extra.columns = TRUE)

# MANE coding regions
mane_coding_ranges <- makeGRangesFromDataFrame(fread(mane_coding.file, header = TRUE, sep = "\t", quote = "", nThread = 1))

# Problematic regions
problematic_ranges <- makeGRangesFromDataFrame(fread(
	problematic_regions.bed, header = FALSE, col.names = c("seqnames", "start", "end"), 
	sep = "\t", quote = "", fill = FALSE, nThread = 1), starts.in.df.are.0based = TRUE)


# Data processing ---------------------------------------------------------

# Keep effective target regions (Mb) with coverage >= minDP
#	- Coding (from MANE regions)
#	- Confidence regions (not found in problematic genomic regions)
effective_ranges <- subsetByOverlaps(subsetByOverlaps(target_perbase_ranges, mane_coding_ranges), problematic_ranges, invert = TRUE)
bases_info <- mergeByOverlaps(effective_ranges, perBaseCov_ranges) %>% as.data.table() %>% filter(DEPTH >= min_dp)
bases_ranges <- GRanges(paste(bases_info$effective_ranges.seqnames, bases_info$effective_ranges.start, sep = ":"))
region_size <- nrow(bases_info) / 1e6

# Create a TMB trace table
tmb_trace <- variants %>% 
	mutate(
		Position = START, RefCall = REF, AltCall = ALT, VAF = AF, TotalDepth = DP, AltDepth = AD_ALT, 
		VariantType = TYPE, VariantClass = CLASS, Consequence = CONSEQUENCE, GeneName = GENE_SYMBOL, Mutation = MUTATION, 
		ClinStatus = CLIN_STATUS, Oncogenicity = ONCOGENICITY, dbSNP = dbSNP_ID, MaxpVAF = gnomAD_MAX_AF, 
		IsGermline = (!(CANCER_HOTSPOT | PANEL_HOTSPOT)) & (gnomAD_MAX_AF > somatic_max_pvaf | AF > somatic_max_vaf), 
		LowConfidenceRegion = PROBLEMATIC_REGION, CancerHotspot = CANCER_HOTSPOT | PANEL_HOTSPOT, Recurrent = RECURRENT, CodingVariant = CODING, 
		ClinRelevant = CLIN_RELEVANT, IsOncogenic = ONCOGENICITY %in% c("Oncogenic", "Likely Oncogenic"),
		Nonsynonymous = (IMPACT == "MODERATE" | IMPACT == "HIGH"), 
		Chromosome = factor(variants$CHROM, levels = chr_names), 
		InEffectiveRegions = overlapsAny(variants_ranges, bases_ranges)) %>% 
	select(
		RUN, SAMPLE, Chromosome, Position, RefCall, AltCall, VAR, VAR_HG19, VAF, TotalDepth, AltDepth, VariantType, VariantClass, Consequence, 
		GeneName, Mutation, ClinStatus, Oncogenicity, dbSNP, MaxpVAF, IsGermline, LowConfidenceRegion, CancerHotspot, Recurrent, CodingVariant, 
		ClinRelevant, IsOncogenic, InEffectiveRegions, Nonsynonymous) %>% 
	mutate(
		IncludedInTmbNumerator = (!IsGermline) & (!LowConfidenceRegion) & VAF>=min_vaf & TotalDepth>=min_dp & AltDepth>=min_ad & 
		(!CancerHotspot) & CodingVariant & (!Recurrent) & (!ClinRelevant) & (!IsOncogenic) & InEffectiveRegions) %>%
	arrange(desc(IncludedInTmbNumerator), Chromosome, Position)

# Create a TMB report
vars_count <- sum(tmb_trace$IncludedInTmbNumerator)
nonsynvars_count <- sum(tmb_trace$IncludedInTmbNumerator & tmb_trace$Nonsynonymous)
results <- data.frame(
	RUN = run, SAMPLE = sample, TmbPerMb = round(vars_count / region_size, 1), 
	NonsynonymousTmbPerMb = round(nonsynvars_count / region_size, 1), EligibleVars = vars_count, 
	EligibleNonsynonymousVars = nonsynvars_count, EffectiveRegionsSizeMb = region_size)


# Save data ---------------------------------------------------------------

write.table(tmb_trace, file = paste0(sample_id, "_TMB_trace.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(results, file = paste0(sample_id, "_TMB_results.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.xlsx(
	list("Results" = results, "Trace" = tmb_trace), asTable = TRUE, overwrite = TRUE,
	file = paste0(run, "_", sample_id, "_TMB_", format(Sys.time(), "%Y%m%d"), ".xlsx"))


