#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of fusion results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(cowplot)))
suppressPackageStartupMessages(suppressWarnings(library(ggplotify)))
suppressPackageStartupMessages(suppressWarnings(library(karyoploteR)))
suppressPackageStartupMessages(suppressWarnings(library(GenomicFeatures)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
suppressPackageStartupMessages(suppressWarnings(library(VariantAnnotation)))
suppressPackageStartupMessages(suppressWarnings(library(DBI)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(openxlsx)))
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
sample.file <- args_list[["samplesFile"]]
target.file <- args_list[["targetBed"]]
genes.file <- args_list[["targetGenes"]]
genes_whitelist.file <- args_list[["whitelistGenes"]]
fusion_whitelist.file <- args_list[["whitelistFusions"]]
mane.gtf <- args_list[["maneGtf"]]
mane_gene.file <- args_list[["maneGene"]]
mane_exon.file <- args_list[["maneExon"]]
mane_intron.file <- args_list[["maneIntron"]]
cancerdrivers.file <- args_list[["cancerdrivers"]]
mitelman.file <- args_list[["mitelmanFusion"]]
genie.file <- args_list[["genieFusion"]]
civic.file <- args_list[["civicClinical"]]
chain.file <- args_list[["chain"]]
fasta.file <- args_list[["fasta"]]
fasta_name <- args_list[["fastaName"]]
headerVcf.file <- args_list[["vcfHeader"]]
calling_min_junction <- as.numeric(args_list[["callingMinJunctionReads"]])
calling_min_spanning <- as.numeric(args_list[["callingMinSpanningReads"]])
calling_min_ad <- as.numeric(args_list[["callingMinAD"]])
calling_min_ffpm <- as.numeric(args_list[["callingMinFFPM"]])
flag_min_ad <- as.numeric(args_list[["flagMinAD"]])
flag_min_ad_nonfused <- as.numeric(args_list[["flagMinADNonFused"]])
flag_min_dp <- as.numeric(args_list[["flagMinDP"]])
flag_min_vaf <- as.numeric(args_list[["flagMinVAF"]])

# Set fusion results
fusion_all.file <- list.files(pattern = "wAnnot$", full.names = TRUE, recursive = TRUE)
fusion_normal.file <- list.files(pattern = "wAnnot.annot_filter.fail$", full.names = TRUE, recursive = TRUE)
fusion_RTartifact.file <- list.files(pattern = "RTartifact.filtered$", full.names = TRUE, recursive = TRUE)
fusion_coding_effect.file <- list.files(pattern = "coding_effect.tsv", full.names = TRUE)
fusion_validated.file <- list.files(pattern = "FusionInspector.fusions.abridged.tsv$", full.names = TRUE, recursive = TRUE)

# Set other files
perBaseCoverage.file <- list.files(pattern = "per-base.bed.gz$")

# Set other variables
sample_id <- paste0(sample, "_", type)
chr_names <- paste0("chr", c(1:22, "X", "Y"))
calling_status <- c("PASS", "LowSupport")
fusion_flags <- c("OK", "NoCall", "Normal", "RTartifact", "NoInSilicoValid", "LowAD", "LowNonFused", "LowDP", "LowVAF", "Unknown")
resources_names <- c("MitelmanDB", "GENIE")
results_excel <- createWorkbook()
results_excel.file <- paste0(run, "_", sample_id, "_Fusion_", format(Sys.time(), "%Y%m%d"), ".xlsx")


# Load data ---------------------------------------------------------------

# Sample info
sample_info <- read.table(sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Sample == sample)
tumor_doid <- sample_info[1, "TumorDOID"]
tumor_topnode_doid <- sample_info[1, "TopNodeDOID"]
tumor_whitelist <- sample_info[1, "Whitelist"]

# FASTA file
fasta <- FaFile(fasta.file)

# Target regions
target_bed <- read.table(target.file, header = FALSE, sep = "\t", col.names = c("CHROM", "START", "END", "GENE"))
target_ranges <- makeGRangesFromDataFrame(target_bed, ignore.strand = TRUE, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)

# Genes
target_genes <- read.table(genes.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
	mutate(GENE = ifelse(!is.na(NEW_GENE), NEW_GENE, OLD_GENE))
target_genes_names <- target_genes$GENE
target_genes <- target_genes %>% select(GENE, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL)
genes_whitelist <- read.table(genes_whitelist.file, header = TRUE, sep = ";", fill = TRUE, na.strings = c("NA", ".", ""))
genes_whitelist <- as.character(na.omit(unique(genes_whitelist[[tumor_whitelist]])))
genes_whitelist <- genes_whitelist[genes_whitelist %in% target_genes_names]

# Known fusion variants from whitelist
fusion_whitelist <- read.table(fusion_whitelist.file, header = TRUE, sep = ";", quote = "") %>% select(VAR, FUSION_VARIANT)

# MANE data
mane_TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(mane.gtf, format = "gtf")))
mane_gene <- fread(mane_gene.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% select(-GENE_NAME) %>% as.data.frame()
mane_gene_names <- mane_gene$GENE_ENSEMBL; names(mane_gene_names) <- mane_gene$GENE_SYMBOL
mane_gene_ranges <- makeGRangesFromDataFrame(mane_gene, keep.extra.columns = TRUE)
mane_exon <- fread(mane_exon.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% as.data.frame()
mane_exon_ranges <- makeGRangesFromDataFrame(mane_exon, keep.extra.columns = TRUE)
mane_intron <- fread(mane_intron.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% as.data.frame()
mane_intron_ranges <- makeGRangesFromDataFrame(mane_intron, keep.extra.columns = TRUE)

# Per-base coverage
if(length(perBaseCoverage.file)){
	perBase_ranges <- makeGRangesFromDataFrame(fread(
		perBaseCoverage.file, header = FALSE, col.names = c("CHROM", "START", "END", "DEPTH"),
		sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% filter(CHROM %in% chr_names), 
		starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
	perBase_ranges$y1 <- perBase_ranges$DEPTH
}

# Network of Cancer Genes (NCG) data
cancerdrivers <- read.table(cancerdrivers.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

# Fusion DBs
mitelman <- read.table(mitelman.file, header = TRUE, sep = "\t", quote = "")
genie <- read.table(genie.file, header = TRUE, sep = "\t", quote = "") %>% distinct(FUSION, GENIE_CNT)

# CIViC clinical data
civic_data <- read.table(civic.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)
colnames(civic_data) <- paste0("CIViC_", colnames(civic_data))

# STAR-Fusion: all fusion candidates
fusion_all <- read.table(fusion_all.file, header = TRUE, sep = "\t", comment.char = "") %>% select(-JunctionReads, -SpanningFrags) %>%
	separate(LeftGene, c("GENE_A", "GENE_ENSEMBL_A"), "\\^", extra = "drop") %>%
	separate(RightGene, c("GENE_B", "GENE_ENSEMBL_B"), "\\^", extra = "drop") %>%
	separate(LeftBreakpoint, c("CHROM_A", "POS_A", "STRAND_A"), ":", extra = "drop", remove = FALSE) %>%
	separate(RightBreakpoint, c("CHROM_B", "POS_B", "STRAND_B"), ":", extra = "drop", remove = FALSE) %>%
	mutate(
		ID = paste0(LeftBreakpoint, "::", RightBreakpoint), BREAKPOINT_A = paste0(CHROM_A, ":", POS_A), 
		BREAKPOINT_B = paste0(CHROM_B, ":", POS_B), FUSION_RANGE_A = BREAKPOINT_A, FUSION_RANGE_B = BREAKPOINT_B, 
		VAR = paste0(BREAKPOINT_A, "::", BREAKPOINT_B))
colnames(fusion_all)[1] <- "FUSION_SHORT"

# STAR-Fusion: normal fusions
fusion_normal_ids <- (read.table(fusion_normal.file, header = TRUE, sep = "\t", comment.char = "") %>% 
	mutate(ID = paste0(LeftBreakpoint, "::", RightBreakpoint)))$ID

# STAR-Fusion: RT artifacts
fusion_RTartifacts <- (read.table(fusion_RTartifact.file, header = TRUE, sep = "\t", comment.char = "") %>% 
	mutate(ID = paste0(LeftBreakpoint, "::", RightBreakpoint)))$ID

# STAR-Fusion: coding effect
fusion_coding_effect <- read.table(fusion_coding_effect.file, header = TRUE, sep = "\t", comment.char = "") %>% 
	filter(FUSION_MODEL != ".") %>%
	mutate(
		ID = paste0(LeftBreakpoint, "::", RightBreakpoint), MODEL_CDS_A = CDS_LEFT_ID, MODEL_CDS_RANGE_A = CDS_LEFT_RANGE, 
		MODEL_CDS_B = CDS_RIGHT_ID, MODEL_CDS_RANGE_B = CDS_RIGHT_RANGE, MODEL_FUSION_TYPE = PROT_FUSION_TYPE,
		FUSION_MODEL = gsub("chr\\d+|chr[XY]|\\|[+-]\\||\\[[012]\\]", "", FUSION_MODEL)) %>%
	separate(FUSION_MODEL, c("MODEL_FUSION_A", "MODEL_FUSION_B"), "<==>", extra = "drop") %>%
	select(ID, MODEL_CDS_A, MODEL_CDS_RANGE_A, MODEL_CDS_B, MODEL_CDS_RANGE_B, MODEL_FUSION_TYPE, MODEL_FUSION_A, MODEL_FUSION_B)

# STAR-Fusion: in-silico validated fusions
if(length(fusion_validated.file)){
	fusion_validated <- read.table(fusion_validated.file, header = TRUE, sep = "\t", comment.char = "") %>%
		mutate(ID = paste0(LeftBreakpoint, "::", RightBreakpoint)) %>%
		select(ID, NumCounterFusionLeft, NumCounterFusionRight)
} else {
	fusion_validated <- data.frame(ID = NA, NumCounterFusionLeft = NA, NumCounterFusionRight = NA)
}


# Data processing ---------------------------------------------------------

# Remove atypical chromosomes
# Remove off-target genes (almost 1 target name)
# Annotate if found in fusion DBs: Mitelman, GENIE
# Annotate known fusion variants from the whitelist
# AD (JunctionReadCount + SpanningFragCount)
# AD_NonFused_{Left,Right}
# DP: AD + AD_NonFused_{A,B}
# AF: AD / DP
# CALLING: PASS or LowSupport
# Add fusion coding effect
# FLAGS: OK or NoCall/Normal/RTartifact/NoInSilicoValid/LowAD/LowNonFused/LowDP/LowVAF/UNKNOWN
fusion_target <- fusion_all %>%
	filter((CHROM_A %in% chr_names & CHROM_B %in% chr_names) & (GENE_A %in% target_genes_names | GENE_B %in% target_genes_names)) %>%
	mutate(FUSION_SHORT = sub("--", "::", FUSION_SHORT)) %>% 
	left_join(mitelman, c("FUSION_SHORT" = "FUSION")) %>% 
	left_join(genie, c("FUSION_SHORT" = "FUSION")) %>% 
	left_join(fusion_whitelist, "VAR") %>%
	mutate(
		RUN = run, SAMPLE = sample, WHITELIST_GENE = GENE_A %in% genes_whitelist | GENE_B %in% genes_whitelist,
		RefSpliceSite = SpliceType %in% "ONLY_REF_SPLICE", LargeAnchorSupport = LargeAnchorSupport %in% "YES_LDAS",
		Normal = ID %in% fusion_normal_ids, RTartifact = ID %in% fusion_RTartifacts, InSilicoValid = ID %in% fusion_validated$ID,
		across(c("MitelmanDB_COUNT", "GENIE_CNT"), ~ ifelse(is.na(.x), 0, as.numeric(.x))),
		MitelmanDB = MitelmanDB_COUNT > 0, GENIE = GENIE_CNT > 0, WHITELIST_FUSION = !is.na(FUSION_VARIANT)) %>%
	left_join(mane_gene, by = c("GENE_A" = "GENE_SYMBOL")) %>% left_join(mane_gene, by = c("GENE_B" = "GENE_SYMBOL")) %>%
	left_join(cancerdrivers, by = c("GENE_A" = "GENE")) %>% left_join(cancerdrivers, by = c("GENE_B" = "GENE")) %>%
	mutate(
		TRANSCRIPT_ENSEMBL_A = TRANSCRIPT_ENSEMBL.x, TRANSCRIPT_REFSEQ_A = TRANSCRIPT_REFSEQ.x,
		ONCOGENE_A = ONCOGENE.x, TSG_A = TSG.x, DRIVER_A = DRIVER.x, CANONICAL_DRIVER_A = CANONICAL_DRIVER.x,
		TRANSCRIPT_ENSEMBL_B = TRANSCRIPT_ENSEMBL.y, TRANSCRIPT_REFSEQ_B = TRANSCRIPT_REFSEQ.y,
		ONCOGENE_B = ONCOGENE.y, TSG_B = TSG.y, DRIVER_B = DRIVER.y, CANONICAL_DRIVER_B = CANONICAL_DRIVER.y,
		across(contains(c("ONCOGENE_", "TSG_", "DRIVER_", "CANONICAL_DRIVER_")), ~ .x & !is.na(.x))) %>%
	left_join(fusion_coding_effect, by = "ID") %>% rowwise() %>%
	mutate(
		MODEL_FUSION_A = ifelse(!is.na(MODEL_FUSION_TYPE),
			paste0(CHROM_A, ":", min(unlist(strsplit(MODEL_FUSION_A, "-|\\|"))), "-", max(unlist(strsplit(MODEL_FUSION_A, "-|\\|")))), NA),
		MODEL_FUSION_B = ifelse(!is.na(MODEL_FUSION_TYPE),
			paste0(CHROM_B, ":", min(unlist(strsplit(MODEL_FUSION_B, "-|\\|"))), "-", max(unlist(strsplit(MODEL_FUSION_B, "-|\\|")))), NA)) %>% 
	ungroup() %>% left_join(fusion_validated, by = "ID") %>% rowwise() %>%
	mutate(
		AD = JunctionReadCount + SpanningFragCount,
		AD_NonFused_A = NumCounterFusionLeft, AD_NonFused_B = NumCounterFusionRight,
		DP = AD + AD_NonFused_A + AD_NonFused_B, AF = round(AD / DP, 4),
		LowSupport = JunctionReadCount < calling_min_junction | SpanningFragCount < calling_min_spanning | 
			AD < calling_min_ad | FFPM < calling_min_ffpm | (SpanningFragCount == 0 & !LargeAnchorSupport),
		PASS = !LowSupport,
		CALLING = paste(calling_status[c(PASS, LowSupport)], collapse = ";"),
		LowAD = AD < flag_min_ad, LowNonFused = ifelse(is.na(AD_NonFused_A), FALSE, (AD_NonFused_A + AD_NonFused_B) < flag_min_ad_nonfused),
        LowDP = ifelse(is.na(DP), FALSE, DP < flag_min_dp), LowVAF = ifelse(is.na(AF), FALSE, AF < flag_min_vaf), 
		UNKNOWN = !(MitelmanDB | GENIE), RESOURCES = paste(resources_names[c(MitelmanDB, GENIE)], collapse = ","),
		NO_FLAGS = PASS & InSilicoValid & !(Normal | RTartifact | LowAD | LowNonFused | LowDP | LowVAF | UNKNOWN),
		FLAGS = paste(fusion_flags[c(NO_FLAGS, !PASS, Normal, RTartifact, !InSilicoValid, LowAD, LowNonFused, LowDP, LowVAF, UNKNOWN)], collapse = ";"),
		EXON_A = NA, INTRON_A = NA, BASES_FROM_EXON_A = NA, COV_IN_FUSION_A = NA, COV_OUT_FUSION_A = NA, 
		EXON_B = NA, INTRON_B = NA, BASES_FROM_EXON_B = NA, COV_IN_FUSION_B = NA, COV_OUT_FUSION_B = NA, 
		VAR_HG19 = NA, BREAKPOINT_A_HG19 = NA, BREAKPOINT_B_HG19 = NA) %>% 
	as.data.frame() %>% arrange(desc(AD)) %>% distinct(VAR, .keep_all = TRUE)

if(nrow(fusion_target)){

	# Annotate the EXON/INTRON of GENE A and B
	fusion_ranges_A <- GRanges(fusion_target$LeftBreakpoint, VAR = fusion_target$VAR, GENE = fusion_target$GENE_A)
	fusion_exons_A <- suppressWarnings(as.data.frame(mergeByOverlaps(fusion_ranges_A, mane_exon_ranges)) %>%
		filter(GENE %in% GENE_SYMBOL) %>% group_by(VAR) %>% summarise(exon_A = EXON))
	fusion_introns_A <- suppressWarnings(as.data.frame(mergeByOverlaps(fusion_ranges_A, mane_intron_ranges)) %>% 
		filter(GENE %in% GENE_SYMBOL) %>% group_by(VAR) %>% summarise(intron_A = INTRON))
	fusion_target <- left_join(left_join(fusion_target, fusion_exons_A, by = "VAR"), fusion_introns_A, by = "VAR") %>%
		mutate(EXON_A = exon_A, INTRON_A = intron_A)
	
	fusion_ranges_B <- GRanges(fusion_target$RightBreakpoint, VAR = fusion_target$VAR, GENE = fusion_target$GENE_B)
	fusion_exons_B <- suppressWarnings(as.data.frame(mergeByOverlaps(fusion_ranges_B, mane_exon_ranges)) %>% 
		filter(GENE %in% GENE_SYMBOL) %>% group_by(VAR) %>% summarise(exon_B = EXON))
	fusion_introns_B <- suppressWarnings(as.data.frame(mergeByOverlaps(fusion_ranges_B, mane_intron_ranges)) %>% 
		filter(GENE %in% GENE_SYMBOL) %>% group_by(VAR) %>% summarise(intron_B = INTRON))
	fusion_target <- left_join(left_join(fusion_target, fusion_exons_B, by = "VAR"), fusion_introns_B, by = "VAR") %>%
		mutate(EXON_B = exon_B, INTRON_B = intron_B) %>% as.data.frame()
	
	# Annotate previous/next exon of intronic breakpoints (depeding of A or B gene)
	# Calculate number of bases from breakpoint to exon
	# Annotate the genomic region of each part of the fusion
	for(n in seq_len(nrow(fusion_target))){
		aux_fusion <- fusion_target[n, ]

		aux_strand_A <- aux_fusion$STRAND_A
		aux_chrom_A <- aux_fusion$CHROM_A
		aux_pos_A <- as.integer(aux_fusion$POS_A)
		aux_gene_A <- aux_fusion$GENE_A
		aux_mane_A <- mane_exon %>% filter(GENE_SYMBOL %in% aux_gene_A)
		aux_intron_A <- aux_fusion$INTRON_A
		aux_exon_A <- ifelse(is.na(aux_fusion$EXON_A) & !is.na(aux_intron_A), aux_intron_A, aux_fusion$EXON_A)
		aux_mane_exon_A <- aux_mane_A %>% filter(EXON %in% aux_exon_A)
		aux_bases_exon_A <- ifelse(is.na(aux_exon_A) & is.na(aux_intron_A), NA, ifelse(
			aux_strand_A %in% "-", aux_mane_exon_A$START - aux_pos_A, aux_pos_A - aux_mane_exon_A$END))
		aux_bases_exon_A <- ifelse(is.na(aux_bases_exon_A) | aux_bases_exon_A == 0, "", ifelse(
			aux_bases_exon_A > 0, paste0("ins", aux_bases_exon_A), sub("-", "del", aux_bases_exon_A)))
		aux_range_A <- ifelse(is.na(aux_exon_A), aux_fusion$FUSION_RANGE_A, ifelse(
			aux_strand_A %in% "-", paste0(aux_chrom_A, ":", aux_pos_A, "-", aux_mane_A[1, "END"]),
			paste0(aux_chrom_A, ":", aux_mane_A[1, "START"], "-", aux_pos_A)))
		if(eval(parse(text = unlist(strsplit(aux_range_A, ":"))[2])) > 0){aux_range_A <- aux_fusion$FUSION_RANGE_A}

		aux_strand_B <- aux_fusion$STRAND_B
		aux_chrom_B <- aux_fusion$CHROM_B
		aux_pos_B <- as.integer(aux_fusion$POS_B)
		aux_gene_B <- aux_fusion$GENE_B
		aux_mane_B <- mane_exon %>% filter(GENE_SYMBOL %in% aux_gene_B)
		aux_intron_B <- aux_fusion$INTRON_B
		aux_exon_B <- ifelse(is.na(aux_fusion$EXON_B) & !is.na(aux_intron_B), aux_intron_B+1, aux_fusion$EXON_B)
		aux_mane_exon_B <- mane_exon %>% filter(GENE_SYMBOL %in% aux_gene_B & EXON %in% aux_exon_B)
		aux_bases_exon_B <- ifelse(is.na(aux_exon_B) & is.na(aux_intron_B), NA, ifelse(
			aux_strand_B %in% "-", aux_pos_B - aux_mane_exon_B$END, aux_mane_exon_B$START - aux_pos_B))
		aux_bases_exon_B <- ifelse(is.na(aux_bases_exon_B) | aux_bases_exon_B == 0, "", ifelse(
			aux_bases_exon_B > 0, paste0("ins", aux_bases_exon_B), sub("-", "del", aux_bases_exon_B)))
		aux_range_B <- ifelse(is.na(aux_exon_B), aux_fusion$FUSION_RANGE_B, ifelse(
			aux_strand_B %in% "-", paste0(aux_chrom_B, ":", aux_mane_B[nrow(aux_mane_B), "START"], "-", aux_pos_B),
			paste0(aux_chrom_B, ":", aux_pos_B, "-", aux_mane_B[nrow(aux_mane_B), "END"])))
		if(eval(parse(text = unlist(strsplit(aux_range_B, ":"))[2])) > 0){aux_range_B <- aux_fusion$FUSION_RANGE_B}
					
		fusion_target[n, ] <- aux_fusion %>% mutate(
			EXON_A = aux_exon_A, FUSION_RANGE_A = aux_range_A, BASES_FROM_EXON_A = aux_bases_exon_A, 
			EXON_B = aux_exon_B, FUSION_RANGE_B = aux_range_B, BASES_FROM_EXON_B = aux_bases_exon_B)
		rm(
			aux_fusion, aux_strand_A, aux_chrom_A, aux_pos_A, aux_gene_A, aux_mane_A, aux_intron_A, aux_exon_A, 
			aux_mane_exon_A, aux_bases_exon_A, aux_range_A, aux_strand_B, aux_chrom_B, aux_pos_B, aux_gene_B, aux_mane_B, 
			aux_intron_B, aux_exon_B, aux_mane_exon_B, aux_bases_exon_B, aux_range_B)
	}
	rm(n)
	
	# Annotate the exon coverage before/after the breakpoint
	if(length(perBaseCoverage.file)){
		mane_exon_perbase_ranges <- makeGRangesFromDataFrame(
			mane_exon_ranges %>% as.data.table() %>% 
				filter(TRANSCRIPT_ENSEMBL %in% c(fusion_target$TRANSCRIPT_ENSEMBL_A, fusion_target$TRANSCRIPT_ENSEMBL_B)) %>% 
				mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% mutate(pos = 1:n(), end = start+pos-1, start=end) %>% 
				as.data.table() %>% select(-ID, -pos), keep.extra.columns = TRUE)
		exon_coverage <- makeGRangesFromDataFrame(
			suppressWarnings(mergeByOverlaps(mane_exon_perbase_ranges, perBase_ranges)) %>% as.data.table(), 
			seqnames.field = "mane_exon_perbase_ranges.seqnames", start.field = "mane_exon_perbase_ranges.start",
			end.field = "mane_exon_perbase_ranges.end", ignore.strand = TRUE, keep.extra.columns = TRUE)

		for(n in seq_len(nrow(fusion_target))){

			# Gene A
			aux_range_A <- fusion_target[n, "FUSION_RANGE_A"]
			aux_transcript_A <- fusion_target[n, "TRANSCRIPT_ENSEMBL_A"]
			aux_exon_A <- fusion_target[n, "EXON_A"]
			aux_coverage_A <- exon_coverage[exon_coverage$TRANSCRIPT_ENSEMBL %in% aux_transcript_A]
			if(!is.na(aux_exon_A)){
				CovInExons_A <- suppressWarnings((subsetByOverlaps(aux_coverage_A, GRanges(aux_range_A)) %>% as.data.table() %>% 
					summarise(COV = round(mean(DEPTH))))$COV)
				CovOutExons_A <- suppressWarnings((subsetByOverlaps(aux_coverage_A, GRanges(aux_range_A), invert = TRUE) %>% as.data.table() %>% 
					summarise(COV = round(mean(DEPTH))))$COV)
				fusion_target[n, "COV_IN_FUSION_A"] <- CovInExons_A
				fusion_target[n, "COV_OUT_FUSION_A"] <- CovOutExons_A

				rm(CovInExons_A, CovOutExons_A)
			}
			rm(aux_range_A, aux_transcript_A, aux_exon_A, aux_coverage_A)

			# Gene B
			aux_range_B <- fusion_target[n, "FUSION_RANGE_B"]
			aux_transcript_B <- fusion_target[n, "TRANSCRIPT_ENSEMBL_B"]
			aux_exon_B <- fusion_target[n, "EXON_B"]
			aux_coverage_B <- exon_coverage[exon_coverage$TRANSCRIPT_ENSEMBL %in% aux_transcript_B]
			if(!is.na(aux_exon_B)){
				CovInExons_B <- suppressWarnings((subsetByOverlaps(aux_coverage_B, GRanges(aux_range_B)) %>% as.data.table() %>% 
					summarise(COV = round(mean(DEPTH))))$COV)
				CovOutExons_B <- suppressWarnings((subsetByOverlaps(aux_coverage_B, GRanges(aux_range_B), invert = TRUE) %>% as.data.table() %>% 
					summarise(COV = round(mean(DEPTH))))$COV)
				fusion_target[n, "COV_IN_FUSION_B"] <- CovInExons_B
				fusion_target[n, "COV_OUT_FUSION_B"] <- CovOutExons_B

				rm(CovInExons_B, CovOutExons_B)
			}
			rm(aux_range_B, aux_transcript_B, aux_exon_B, aux_coverage_B)
		}
		rm(n, mane_exon_perbase_ranges, exon_coverage)
	}

	# Convert Hg38 to Hg19 coordinates
	fusion_hg19_A <- as.data.frame(unlist(liftOver(fusion_ranges_A, import.chain(chain.file)))) %>% distinct(VAR, .keep_all = TRUE)
	if(nrow(fusion_hg19_A)){
		fusion_hg19_A <- fusion_hg19_A %>% mutate(breakpoint_A_hg19 = paste0(seqnames, ":", start)) %>% select(VAR, breakpoint_A_hg19)
	} else {
		fusion_hg19_A <- fusion_hg19_A %>% mutate(VAR = "", breakpoint_A_hg19 = "") %>% select(VAR, breakpoint_A_hg19)
	}
	fusion_hg19_B <- as.data.frame(unlist(liftOver(fusion_ranges_B, import.chain(chain.file)))) %>% distinct(VAR, .keep_all = TRUE)
	if(nrow(fusion_hg19_B)){
		fusion_hg19_B <- fusion_hg19_B %>% mutate(breakpoint_B_hg19 = paste0(seqnames, ":", start)) %>% select(VAR, breakpoint_B_hg19)
	} else {
		fusion_hg19_B <- fusion_hg19_B %>% mutate(VAR = "", breakpoint_B_hg19 = "") %>% select(VAR, breakpoint_B_hg19)
	}
	fusion_target <- left_join(fusion_target, left_join(fusion_hg19_A, fusion_hg19_B, by = "VAR"), by = "VAR") %>% mutate(
		BREAKPOINT_A_HG19 = breakpoint_A_hg19, BREAKPOINT_B_HG19 = breakpoint_B_hg19, 
		VAR_HG19 = paste0(BREAKPOINT_A_HG19, "::", BREAKPOINT_B_HG19))
}

# Create final output
# Sorting: WHITELIST_FUSION > NO_FLAGS > PASS > GENE WHITELIST > GENE > AD
fusion_results <- fusion_target %>%
	mutate(EXON_A = ifelse(is.na(EXON_A), "", EXON_A), EXON_B = ifelse(is.na(EXON_B), "", EXON_B)) %>% rowwise() %>% 
	mutate(FUSION_NAME = paste0(
		GENE_A, "::", GENE_B, " ", if(isTRUE(!is.na(FUSION_VARIANT))){paste0(FUSION_VARIANT, " ")}, "(", unlist(strsplit(GENE_A, ""))[1], 
		EXON_A, BASES_FROM_EXON_A, "::", BASES_FROM_EXON_B, unlist(strsplit(GENE_B, ""))[1], EXON_B, ")")) %>% as.data.frame() %>% 
	arrange(desc(WHITELIST_FUSION), desc(NO_FLAGS), desc(PASS), desc(WHITELIST_GENE), desc(AD)) %>%
	select(
		RUN, SAMPLE, VAR, FUSION_SHORT, FUSION_NAME, FUSION_VARIANT, AD, AF, CALLING, FLAGS, DP, FFPM, 
		COV_IN_FUSION_A, COV_OUT_FUSION_A, COV_IN_FUSION_B, COV_OUT_FUSION_B, WHITELIST_FUSION, 
		RESOURCES, MitelmanDB_COUNT, GENIE_CNT, WHITELIST_GENE, MODEL_FUSION_TYPE, 
		LargeAnchorSupport, RefSpliceSite, LeftBreakDinuc, RightBreakDinuc, 
		JunctionReadCount, SpanningFragCount, AD_NonFused_A, AD_NonFused_B, 
		FUSION_RANGE_A, MODEL_FUSION_A, MODEL_CDS_A, MODEL_CDS_RANGE_A, 
		FUSION_RANGE_B, MODEL_FUSION_B, MODEL_CDS_B, MODEL_CDS_RANGE_B, 
		ONCOGENE_A, TSG_A, DRIVER_A, CANONICAL_DRIVER_A, GENE_A, GENE_ENSEMBL_A, 
		TRANSCRIPT_ENSEMBL_A, TRANSCRIPT_REFSEQ_A, EXON_A, INTRON_A, BASES_FROM_EXON_A,
		ONCOGENE_B, TSG_B, DRIVER_B, CANONICAL_DRIVER_B, GENE_B, GENE_ENSEMBL_B, 
		TRANSCRIPT_ENSEMBL_B, TRANSCRIPT_REFSEQ_B, EXON_B, INTRON_B, BASES_FROM_EXON_B,  
		BREAKPOINT_A, CHROM_A, POS_A, STRAND_A, BREAKPOINT_B, CHROM_B, POS_B, STRAND_B, 
		VAR_HG19, BREAKPOINT_A_HG19, BREAKPOINT_B_HG19,
		PASS, LowSupport, Normal, RTartifact, InSilicoValid, 
		NO_FLAGS, LowAD, LowNonFused, LowDP, LowVAF, UNKNOWN, MitelmanDB, GENIE)

# Clinical classification (AMP/ASCO/CAP guidelines)
#	- Predictive, prognostic & diagnostic evidences in CIViC resource
# 	- Tier IA: CIViC evidence A (guidelines) in the SAME tumor
# 	- Tier IB: CIViC evidence B (trials) in the SAME tumor
# 	- Tier IIC: CIViC evidence A/B in OTHER tumor OR CIViC evidence C
# 	- Tier IID: CIViC evidence D
# 	- Tier III: Observed in cancer resources (GENIE, MitelmanDB)
# 	- Tier IV: Not observed in cancer resources
fusion_civic <- fusion_results %>% 
	mutate(FUSION_A = paste0(GENE_A, "::."), FUSION_B = paste0(".::", GENE_B)) %>%
	filter(FUSION_SHORT %in% civic_data$CIViC_ALTERATION | FUSION_A %in% civic_data$CIViC_ALTERATION | FUSION_B %in% civic_data$CIViC_ALTERATION)
civic_annot <- civic_data %>% mutate(VAR = ".", CLIN_STATUS = ".") %>% filter(CIViC_ALTERATION == ".")
for(n in seq_len(nrow(fusion_civic))){
	aux_var <- fusion_civic[n, ]
	aux_civic <- civic_data %>% filter(CIViC_ALTERATION %in% c(aux_var$FUSION, aux_var$FUSION_A, aux_var$FUSION_B)) %>%
		mutate(
			ID = aux_var$ID, VAR = aux_var$VAR,
			isTumor = CIViC_TUMOR_DOID %in% tumor_doid | CIViC_TOPNODE_DOID %in% tumor_topnode_doid | CIViC_TUMOR_NAME %in% "cancer",
			CLIN_STATUS = ifelse(
				isTumor & CIViC_EVIDENCE_LEVEL %in% "A", "Tier IA", ifelse(
					isTumor & CIViC_EVIDENCE_LEVEL %in% "B", "Tier IB", ifelse(
						(CIViC_EVIDENCE_LEVEL %in% c("A", "B") & !isTumor) | CIViC_EVIDENCE_LEVEL %in% "C", "Tier IIC", ifelse(
							CIViC_EVIDENCE_LEVEL %in% "D", "Tier IID", NA))))) %>%
		select(-isTumor) %>%
		arrange(CLIN_STATUS) %>%
		head(1)
	civic_annot <- rbind(civic_annot, aux_civic)
	rm(aux_var, aux_civic)
}
rm(n)
fusion_results <- fusion_results %>% left_join(civic_annot, by = "VAR") %>%
	mutate(
		CLIN_STATUS = ifelse(!is.na(CLIN_STATUS), CLIN_STATUS, ifelse(!UNKNOWN, "Tier III", "Tier IV")),
		CLIN_RELEVANT = CLIN_STATUS %in% c("Tier IA", "Tier IB", "Tier IIC", "Tier IID")) %>%
	select(RUN:FLAGS, CLIN_RELEVANT, CLIN_STATUS, everything()) %>%
	arrange(CLIN_STATUS)
fusion_results[is.na(fusion_results) | fusion_results == "NaN" | fusion_results == ""] <- "."

# Generate a VCF file
aux_header <- readLines(headerVcf.file)
aux_header[grep("##fileDate", aux_header)] <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
aux_header[grep("##reference", aux_header)] <- paste0("##reference=file:", fasta_name)
if(nrow(fusion_results)){
	aux_vars <- fusion_results %>% 
		mutate(
			POS_A = as.integer(POS_A), POS_B = as.integer(POS_B), SVTYPE = "BND", 
			CHROM_A = factor(CHROM_A, levels = chr_names), CHROM_B= factor(CHROM_B, levels = chr_names)) %>% 
		arrange(CHROM_A, POS_A, CHROM_B, POS_B)
	aux_ranges <- GRanges(aux_vars$BREAKPOINT_A)
	aux_ref <- as.data.frame(getSeq(fasta, aux_ranges))$x
	aux_info <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c(
			"VAR", "RUN", "SAMPLE", "SVTYPE", "FUSION_SHORT", "FUSION_NAME", "FUSION_VARIANT", "AD", "AF", "CLIN_STATUS", 
			"DP", "FFPM", "COV_IN_FUSION_A", "COV_OUT_FUSION_A", "COV_IN_FUSION_B", "COV_OUT_FUSION_B",
			"RESOURCES", "MitelmanDB_COUNT", "GENIE_CNT", "MODEL_FUSION_TYPE", 
			"LeftBreakDinuc", "RightBreakDinuc", "JunctionReadCount", "SpanningFragCount", "AD_NonFused_A", "AD_NonFused_B",
			"FUSION_RANGE_A", "MODEL_FUSION_A", "MODEL_CDS_A", "MODEL_CDS_RANGE_A", 
			"FUSION_RANGE_B", "MODEL_FUSION_B", "MODEL_CDS_B", "MODEL_CDS_RANGE_B", 
			"GENE_A", "GENE_ENSEMBL_A", "TRANSCRIPT_ENSEMBL_A", "TRANSCRIPT_REFSEQ_A", "EXON_A", "BASES_FROM_EXON_A",
			"GENE_B", "GENE_ENSEMBL_B", "TRANSCRIPT_ENSEMBL_B", "TRANSCRIPT_REFSEQ_B", "EXON_B", "BASES_FROM_EXON_B",  
			"BREAKPOINT_A", "CHROM_A", "POS_A", "STRAND_A", "BREAKPOINT_B", "CHROM_B", "POS_B", "STRAND_B", 
			"VAR_HG19", "BREAKPOINT_A_HG19", "BREAKPOINT_B_HG19",
			"CIViC_ALTERATION", "CIViC_VARIANT_ID", "CIViC_EVIDENCE_LEVEL", "CIViC_EVIDENCE_RATING", "CIViC_EVIDENCE_TYPE", 
			"CIViC_EFFECT", "CIViC_DRUG", "CIViC_TUMOR_NAME", "CIViC_TUMOR_DOID", "CIViC_TUMOR_CODE", "CIViC_TOPNODE_NAME",
			"CIViC_TOPNODE_DOID", "CIViC_MOLECULAR_ID", "CIViC_VARIANT_GROUP", "CIViC_EVIDENCE_SCORE", "CIViC_ORIGIN"),
		function(x) paste0(x, "=", aux_vars[n, x]) )), collapse = ";")  }))
	aux_flags <- paste0(aux_vars$FLAGS, ";", unlist(lapply(seq_len(nrow(aux_vars)), function(n){paste(
		unlist(lapply(c(
			"LargeAnchorSupport", "RefSpliceSite", "WHITELIST_GENE", 
			"ONCOGENE_A", "TSG_A", "DRIVER_A", "CANONICAL_DRIVER_A", 
			"ONCOGENE_B", "TSG_B", "DRIVER_B", "CANONICAL_DRIVER_B", 
			"MitelmanDB", "GENIE", "WHITELIST_FUSION", "CLIN_RELEVANT"), 
		function(x) ifelse(aux_vars[n, x], paste0(x, ";"), "")  )), collapse = "")  })))
	aux_info <- paste0(aux_info, ";", aux_flags)
	aux_format <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("AD", "AF", "DP"), function(x) paste0(aux_vars[n, x]) )), collapse = ":")  }))
	aux_df <- aux_vars %>% mutate(
		CHROM_A = as.character(CHROM_A), POS = POS_A, ID = ".", REF = aux_ref, QUAL = ".", FILTER = CALLING, INFO = aux_info,
		FORMAT = paste(c("AD", "AF", "DP"), collapse = ":"), SAMPLE = aux_format,
		ALT = ifelse(STRAND_B %in% "-", paste0(REF, "]", BREAKPOINT_B, "]"), paste0(REF, "[", BREAKPOINT_B, "["))) %>%
		select(CHROM_A, POS_A, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)
	vcf_text <- c(aux_header,
		paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"),
		unlist(lapply(seq_len(nrow(aux_df)), function(x) paste(aux_df[x, ], collapse = "\t"))) )
	
	rm(aux_header, aux_vars, aux_ranges, aux_ref, aux_info, aux_flags, aux_format, aux_df)

} else {
	vcf_text <- c(aux_header, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"))
}


# Save data ---------------------------------------------------------------

write.table(
	fusion_results, file = paste0(sample_id, "_Fusion_allCandidates.txt"), 
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(
	fusion_results %>% filter(PASS), file = paste0(sample_id, "_Fusion_results.txt"), 
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

addWorksheet(results_excel, "OK")
writeDataTable(results_excel, "OK", fusion_results %>% filter(NO_FLAGS))
addWorksheet(results_excel, "PASS")
writeDataTable(results_excel, "PASS", fusion_results %>% filter(PASS))
addWorksheet(results_excel, "FAIL")
writeDataTable(results_excel, "FAIL", fusion_results %>% filter(!PASS))
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

writeLines(con = paste0(sample_id, "_Fusion.vcf"), text = vcf_text)


# Data visualization ---------------------------------------------------------------

# Visualization of gene coverage: WHITELIST_FUSION | PASS variants
plot_vars <- fusion_results %>% filter(WHITELIST_FUSION | PASS)

if(nrow(plot_vars)) {

	addWorksheet(results_excel, "Plots")

	plot.params <- getDefaultPlotParams(plot.type = 3)
	plot.params$data2height <- 0
	plot.params$topmargin <- 40

	excel_startRow <- 1

	for(n in seq_len(nrow(plot_vars))){
		
		aux_fusion_vars <- plot_vars[n, ]
		fusion_id <- paste0(
			gsub("::", "-", aux_fusion_vars$FUSION_SHORT), "_", 
			gsub(":", "_", aux_fusion_vars$BREAKPOINT_A), "-", gsub(":", "_", aux_fusion_vars$BREAKPOINT_B))
		fusion_name <- aux_fusion_vars$FUSION_NAME
		fusion_AD <- aux_fusion_vars$AD
		fusion_DP <- aux_fusion_vars$DP
		fusion_AF <- ifelse(aux_fusion_vars$AF != ".", round(as.numeric(aux_fusion_vars$AF), 3), aux_fusion_vars$AF) 

		# Gene A
		transcript_A <- aux_fusion_vars$TRANSCRIPT_ENSEMBL_A
		if(transcript_A != "."){
			gene_A <- aux_fusion_vars$GENE_A
			aux_gene_A <- mane_gene_names[gene_A]
			gene_ranges_A <- mane_gene_ranges[mane_gene_ranges$TRANSCRIPT_ENSEMBL %in% transcript_A]
			exon_ranges_A <- mane_exon_ranges[mane_exon_ranges$TRANSCRIPT_ENSEMBL %in% transcript_A & mane_exon_ranges$EXON %in% aux_fusion_vars$EXON_A]
			width_A <- width(gene_ranges_A)
			zoom_A <- gene_ranges_A + width_A*0.01

			title_A <- paste0(ifelse(transcript_A != ".", paste0(gene_A, " (", transcript_A, ")"), paste0(aux_fusion_vars$GENE_A)))

			p1 <- as.ggplot(function(){
				kp <- plotKaryotype("hg38", cex = 1.5, plot.params = plot.params, zoom = zoom_A)
				kpAddBaseNumbers(kp, cex = 0.7, minor.ticks = TRUE, tick.len = 4, minor.tick.len = 1, add.units = TRUE, digits = 3,
					tick.dist=round(width(zoom_A)/5, -2), minor.tick.dist=round(width(zoom_A)/40, -2))
				title(paste0(fusion_name, "   AD=", fusion_AD, " DP=", fusion_DP, " AF=", fusion_AF, ifelse(
					aux_fusion_vars$NO_FLAGS, "", paste0("  *", aux_fusion_vars$FLAGS, "*"))), cex = 0.9, line = 0.5)
				title(title_A, cex = 0.8, adj = 0, line = -0.5)
				genes.data_A <- makeGenesDataFromTxDb(mane_TxDb, kp)
				is_gene <- names(genes.data_A$genes) %in% aux_gene_A
				genes.data_A$genes <- genes.data_A$genes[is_gene]
				genes.data_A$transcripts <- genes.data_A$transcripts[is_gene]
				is_transcript <- genes.data_A[["transcripts"]][[aux_gene_A]]$tx_name %in% transcript_A
				genes.data_A[["transcripts"]][[aux_gene_A]] <- genes.data_A[["transcripts"]][[aux_gene_A]][is_transcript]
				genes.data_A[["coding.exons"]] <- genes.data_A[["coding.exons"]][is_transcript]
				genes.data_A[["non.coding.exons"]] <- genes.data_A[["non.coding.exons"]][is_transcript]
				if(length(exon_ranges_A)){
					kpPlotMarkers(
						kp, exon_ranges_A, labels = exon_ranges_A$EXON, lwd = 0.4, cex = 0.7, text.orientation = "horizontal",
						data.panel = 1, r0 = 0.85, r1 = 0.97, line.color = "black", label.color = "black")
				}
				suppressWarnings(kpPlotGenes(
					kp, genes.data_A, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.8, r1 = 0.95, add.strand.marks=TRUE, 
					mark.height = 0.2, coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50", 
					non.coding.exons.border.col="grey50"))
				kpRect(kp, GRanges(aux_fusion_vars$FUSION_RANGE_A), y0 = 0, y1 = 1, border = "red", data.panel = 1, r0 = 0.79, r1 = 0.91, lwd = 1)
				if(length(perBaseCoverage.file)){
					perBase_ranges_A <- subsetByOverlaps(perBase_ranges, zoom_A)
					ymax_A <- plyr::round_any(max(perBase_ranges_A$y1), 100, f = ceiling)
					ymax_sep_A <- plyr::round_any(ymax_A/5, 100, f = ceiling)
					kpAxis(kp, ymin = 0, ymax = ymax_A, cex = 0.65, data.panel = 1, r0 = 0.03, r1 = 0.7, tick.pos = seq(0, ymax_A, ymax_sep_A), tick.len = width_A/800)
					kpAddLabels(kp, labels = "Coverage", cex = 0.8, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.03, r1 = 0.7)
					kpBars(kp, perBase_ranges_A, ymin = 0, ymax = ymax_A, data.panel = 1, r0 = 0.03, r1 = 0.7, col = "black", border = "black")
				}
				kpAbline(kp, v = as.integer(aux_fusion_vars$POS_A), data.panel = 1, r0 = 0, r1 = 0.9, col = "red", lwd = 1, lty = 2)
				if(overlapsAny(zoom_A, target_ranges)){
					kpPlotRegions(kp, target_ranges, col = "black", data.panel = 1, r0 = -0.05, r1 = 0)
				}
			})

			rm(gene_A, aux_gene_A, gene_ranges_A, exon_ranges_A, width_A, zoom_A, title_A)

		} else {p1 <- NULL}
		
		
		# Gene B
		transcript_B <- aux_fusion_vars$TRANSCRIPT_ENSEMBL_B
		if(transcript_B != "."){
			gene_B <- aux_fusion_vars$GENE_B
			aux_gene_B <- mane_gene_names[gene_B]
			gene_ranges_B <- mane_gene_ranges[mane_gene_ranges$TRANSCRIPT_ENSEMBL %in% transcript_B]
			exon_ranges_B <- mane_exon_ranges[mane_exon_ranges$TRANSCRIPT_ENSEMBL %in% transcript_B & mane_exon_ranges$EXON %in% aux_fusion_vars$EXON_B]
			width_B <- width(gene_ranges_B)
			zoom_B <- gene_ranges_B + width_B*0.01

			title_B <- paste0(ifelse(transcript_B != ".", paste0(gene_B, " (", transcript_B, ")"), paste0(aux_fusion_vars$GENE_B)))

			p2 <- as.ggplot(function(){
				kp <- plotKaryotype("hg38", cex = 1.5, plot.params = plot.params, zoom = zoom_B)
				kpAddBaseNumbers(kp, cex = 0.7, minor.ticks = TRUE, tick.len = 4, minor.tick.len = 1, add.units = TRUE, digits = 3,
					tick.dist=round(width(zoom_B)/5, -2), minor.tick.dist=round(width(zoom_B)/40, -2))
				if(transcript_A == "."){title(paste0(fusion_name, "   AD=", fusion_AD, " DP=", fusion_DP, " AF=", fusion_AF, ifelse(
					aux_fusion_vars$NO_FLAGS, "", paste0("  *", aux_fusion_vars$FLAGS, "*"))), cex = 0.9, line = 0.5)}
				title(title_B, cex = 0.8, adj = 0, line = -0.5)
				genes.data_B <- makeGenesDataFromTxDb(mane_TxDb, kp)
				is_gene <- names(genes.data_B$genes) %in% aux_gene_B
				genes.data_B$genes <- genes.data_B$genes[is_gene]
				genes.data_B$transcripts <- genes.data_B$transcripts[is_gene]
				is_transcript <- genes.data_B[["transcripts"]][[aux_gene_B]]$tx_name %in% transcript_B
				genes.data_B[["transcripts"]][[aux_gene_B]] <- genes.data_B[["transcripts"]][[aux_gene_B]][is_transcript]
				genes.data_B[["coding.exons"]] <- genes.data_B[["coding.exons"]][is_transcript]
				genes.data_B[["non.coding.exons"]] <- genes.data_B[["non.coding.exons"]][is_transcript]
				if(length(exon_ranges_B)){
					kpPlotMarkers(
						kp, exon_ranges_B, labels = exon_ranges_B$EXON, lwd = 0.4, cex = 0.7, text.orientation = "horizontal",
						data.panel = 1, r0 = 0.85, r1 = 0.97, line.color = "black", label.color = "black")
				}
				suppressWarnings(kpPlotGenes(
					kp, genes.data_B, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.8, r1 = 0.95, add.strand.marks=TRUE,
					mark.height = 0.2, coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50",
					non.coding.exons.border.col="grey50"))
				kpRect(kp, GRanges(aux_fusion_vars$FUSION_RANGE_B), y0 = 0, y1 = 1, border = "red", data.panel = 1, r0 = 0.79, r1 = 0.91, lwd = 1)
				if(length(perBaseCoverage.file)){
					perBase_ranges_B <- subsetByOverlaps(perBase_ranges, zoom_B)
					ymax_B <- plyr::round_any(max(perBase_ranges_B$y1), 100, f = ceiling)
					ymax_sep_B <- plyr::round_any(ymax_B/5, 100, f = ceiling)
					kpAxis(kp, ymin = 0, ymax = ymax_B, cex = 0.65, data.panel = 1, r0 = 0.03, r1 = 0.7, tick.pos = seq(0, ymax_B, ymax_sep_B), tick.len = width_B/800)
					kpAddLabels(kp, labels = "Coverage", cex = 0.8, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.03, r1 = 0.7)
					kpBars(kp, perBase_ranges_B, ymin = 0, ymax = ymax_B, data.panel = 1, r0 = 0.03, r1 = 0.7, col = "black", border = "black")
				}
				kpAbline(kp, v = as.integer(aux_fusion_vars$POS_B), data.panel = 1, r0 = 0, r1 = 0.9, col = "red", lwd = 1, lty = 2)
				if(overlapsAny(zoom_B, target_ranges)){
					kpPlotRegions(kp, target_ranges, col = "black", data.panel = 1, r0 = -0.05, r1 = 0)
				}
			})

			rm(gene_B, aux_gene_B, gene_ranges_B, exon_ranges_B, width_B, zoom_B, title_B)

		} else {p2 <- NULL}

		# Combine plots
		png(paste0(sample_id, "_Fusion_", fusion_id, ".png"), width = 3000, height = 1200, res = 200)
		try(print(plot_grid(p1, p2, ncol = 1)))
		dev.off()

		# Add the plot to the excel of results
		try(insertImage(
				results_excel, "Plots", paste0(sample_id, "_Fusion_", fusion_id, ".png"), 
				units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow))
		saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)
		
		excel_startRow <- excel_startRow + 41

		rm(aux_fusion_vars, fusion_id, fusion_name, fusion_AD, fusion_DP, fusion_AF, transcript_A, p1, transcript_B, p2)  
			
	}
	rm(n, excel_startRow)
}


