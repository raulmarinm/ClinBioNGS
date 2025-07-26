#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of CNA results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(CopyNumberPlots))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(suppressWarnings(library(VariantAnnotation)))
suppressPackageStartupMessages(suppressWarnings(library(DBI)))
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
sample.file <- args_list[["samplesFile"]]
target.file <- args_list[["targetBed"]]
fasta.file <- args_list[["fasta"]]
fasta_name <- args_list[["fastaName"]]
arm.file <- args_list[["armCoordinates"]]
genes.file <- args_list[["targetGenes"]]
whitelist.file <- args_list[["whitelistGenes"]]
mane.gtf <- args_list[["maneGtf"]]
mane_gene.file <- args_list[["maneGene"]]
mane_exon.file <- args_list[["maneExon"]]
cancerdrivers.file <- args_list[["cancerdrivers"]]
genie.file <- args_list[["genie"]]
civic.file <- args_list[["civicClinical"]]
headerVcf.file <- args_list[["vcfHeader"]]
drop_low_coverage <- as.logical(args_list[["dropLowCoverage"]])
calling <- as.numeric(unlist(strsplit(args_list[["callingThresholds"]], ",")))
min_target_bins <- as.numeric(args_list[["minBins"]])
CN_AMP_OK <- as.numeric(args_list[["highAmpCN"]])
CN_DEL_OK <- as.numeric(args_list[["highDelCN"]])
plot_global_label_all_genes <- as.logical(args_list[["plotGlobalLabelAllGenes"]])
plot_bygene_all <- as.logical(args_list[["plotByGeneAll"]])

# Set other files
perBaseCoverage.file <- list.files(pattern = "per-base.bed.gz$")
cnvkit_bins.file <- list.files(pattern = "\\.cnr$")
cnvkit_sex.file <- list.files(pattern = "sex.txt")
ascets_results.file <- list.files(pattern = "CNA_arm.txt")

# Set other variables
sample_id <- paste0(sample, "_", type)
chr_levels <- paste0("chr", c(1:22, "X"))
CN_thresholds <- data.frame(
	LOWER = unique(c(-Inf, calling, log2((length(calling):99 + .5) / 2))),
	UPPER = unique(c(calling, log2( (length(calling):99 + .5) / 2), Inf)),
	CN = 0:100
) %>% mutate(
	CNA = ifelse(CN > 2, "AMP", ifelse(CN < 2, "DEL", "NEUTRAL")),
	CLASS = ifelse(
		CN >= CN_AMP_OK, "HighAmp", ifelse(
			CN > 2, "LowAmp", ifelse(
				CN == 2, "Neutral", ifelse(
					CN > CN_DEL_OK, "LowDel", ifelse(
						CN <= CN_DEL_OK, "HighDel", NA))))),
	OK = CN >= CN_AMP_OK | CN <= CN_DEL_OK)
log2_AMP <- (CN_thresholds %>% filter(CN == 3))$LOWER
log2_AMP_OK <- min((CN_thresholds %>% filter(CN >= CN_AMP_OK))$LOWER)
log2_DEL <- (CN_thresholds %>% filter(CN == 1))$UPPER
log2_DEL_OK <- max((CN_thresholds %>% filter(CN <= CN_DEL_OK))$UPPER)
CN_flags <- c("OK", "LowAmp", "LowDel", "LowBins", "Neutral")
CN_levels <- c("HighAmp", "HighDel", "LowAmp", "LowDel", "Neutral")
CN_colors <- c(
	"HighAmp" = "blue", "LowAmp" = "#00BFFF", "LowDel" = "lightcoral", "HighDel" = "red", 
	"Neutral" = "grey30", "AMP/ASCO/CAP Tier I/II" = "orange")
CN_legend <- c(
	paste0("HighAmp (CN\U2265", CN_AMP_OK, ")"), paste0("LowAmp (2<CN<", CN_AMP_OK, ")"), 
	paste0("LowDel (", CN_DEL_OK, "<CN<2)"), paste0("HighDel (CN\U2264", CN_DEL_OK, ")"), 
	"Neutral (CN=2)", "AMP/ASCO/CAP Tier I/II")
arm_colors <- c("AMP" = "steelblue", "DEL" = "brown", "CONFLICT" = "#FFDB6D")
results_excel <- createWorkbook()
results_excel.file <- paste0(run, "_", sample_id, "_CNA_", format(Sys.time(), "%Y%m%d"), ".xlsx")


# Load data ---------------------------------------------------------------

# Sample info
sample_info <- read.table(sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Sample == sample)
tumor_doid <- sample_info[1, "TumorDOID"]
tumor_topnode_doid <- sample_info[1, "TopNodeDOID"]
tumor_whitelist <- sample_info[1, "Whitelist"]
if(sample_info$Sex != "."){
	is_male <- sample_info$Sex %in% "Male"
} else {
	sex_df <- read.table(cnvkit_sex.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)
	is_male <- sex_df$sex %in% "Male"
}

# FASTA file
fasta <- FaFile(fasta.file)

# Arm genomic coordinates (remove chrY)
arm_coord <- read.table(arm.file, header = TRUE) %>% filter(chrom != "chrY")
arm_coord_autosomal_ranges <- GRanges(filter(arm_coord, !chrom %in% c("chrX", "chrY")))
arm_coord_X_ranges <- GRanges(filter(arm_coord, chrom %in% "chrX"))

# Target regions
target_bed <- read.table(target.file, header = FALSE, sep = "\t", col.names = c("CHROM", "START", "END", "GENE"))
target_ranges <- makeGRangesFromDataFrame(target_bed, ignore.strand = TRUE, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)

# Genes
panel_genes <- read.table(genes.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
	mutate(GENE = ifelse(!is.na(NEW_GENE), NEW_GENE, OLD_GENE))
panel_genes_names <- panel_genes$GENE; names(panel_genes_names) <- panel_genes$OLD_GENE
panel_genes <- panel_genes %>% select(GENE, CYTOBAND, STRAND, GENE_ENSEMBL, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL)
genes_whitelist <- read.table(whitelist.file, header = TRUE, sep = ";", fill = TRUE, na.strings = c("NA", ".", ""))
genes_whitelist <- as.character(na.omit(unique(genes_whitelist[[tumor_whitelist]])))
genes_whitelist <- genes_whitelist[genes_whitelist %in% panel_genes_names]

# MANE data (only target genes)
mane_TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(mane.gtf, format = "gtf")))
mane_gene <- fread(mane_gene.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% 
	filter(GENE_SYMBOL %in% panel_genes_names) %>% select(-GENE_NAME)
mane_gene_names <- mane_gene$GENE_ENSEMBL; names(mane_gene_names) <- mane_gene$GENE_SYMBOL
mane_gene_ranges <- makeGRangesFromDataFrame(mane_gene, keep.extra.columns = TRUE)
mane_gene_width <- width(mane_gene_ranges); names(mane_gene_width) <- mane_gene_ranges$GENE_SYMBOL
mane_perbase_ranges <- makeGRangesFromDataFrame(
	mane_gene_ranges %>% as.data.table() %>% select(seqnames:GENE_SYMBOL) %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
	mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.table() %>% select(-ID, -pos), keep.extra.columns = TRUE)
mane_exon <- fread(mane_exon.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% panel_genes_names)
mane_exon_ranges <- makeGRangesFromDataFrame(mane_exon, keep.extra.columns = TRUE)
mane_exon_width <- as.data.table(mane_exon_ranges) %>% group_by(GENE_SYMBOL) %>% summarise(WIDTH = sum(width)) %>% as.data.frame()
rownames(mane_exon_width) <- mane_exon_width$GENE_SYMBOL
mane_perbase_exon <- subsetByOverlaps(mane_perbase_ranges, mane_exon_ranges)

# Network of Cancer Genes (NCG) data
cancerdrivers <- read.table(cancerdrivers.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

# Load GENIE CNAs
if(file.exists(genie.file)){
	genie_cna <- read.table(genie.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% select(-SAMPLES)
} else {
	genie_cna <- data.frame(GENE = NA, CNA = NA, GENIE_CNT = NA, GENIE_FREQ = NA)
}

# CIViC clinical data
civic_data <- read.table(civic.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)
colnames(civic_data) <- paste0("CIViC_", colnames(civic_data))

# Per-base coverage
if(length(perBaseCoverage.file)){
	perBase_ranges <- makeGRangesFromDataFrame(fread(
		perBaseCoverage.file, header = FALSE, col.names = c("CHROM", "START", "END", "DEPTH"),
		sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% filter(CHROM %in% chr_levels), 
		starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
	perBase_ranges$y1 <- perBase_ranges$DEPTH
}

# CNVkit bins (remove bins of iSNP, iIndel, and chrY)
bins <- read.table(cnvkit_bins.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
	filter(!gene %in% c("iSNP", "iIndel", "chrY")) %>%
	mutate(start = start+1, y0 = log2, y1 = log2, BIN_ID = paste0(chromosome, ":", start, "-", end))

# Ascets results
if(length(ascets_results.file)){
	arm_results <- read.table(ascets_results.file, header = TRUE, sep = "\t")
	arm_results_ranges <- makeGRangesFromDataFrame(arm_results %>% filter(CNA != "."), keep.extra.columns = TRUE)
}


# Data processing ---------------------------------------------------------

# Set the updated gene name
bins <- bins %>% mutate(gene = ifelse(gene == "Antitarget", gene, panel_genes_names[gene]))

# If drop_low_coverage = TRUE, discard bins with log2<-15
if(drop_low_coverage) {bins <- bins %>% filter(log2 >= -15)}
bins_ranges <- makeGRangesFromDataFrame(bins, keep.extra.columns = TRUE, ignore.strand = TRUE)

# Keep only on-target bins
bins_target <- bins %>% filter(gene != "Antitarget") %>% group_by(gene) %>% summarise(
	chromosome = unique(chromosome), start = min(start), end = max(end), 
	DEPTH = weighted.mean(depth, weight), BINS_TARGET = n())
bins_target_ranges <- makeGRangesFromDataFrame(bins_target, ignore.strand = TRUE, keep.extra.columns = TRUE)

# Processing of bins by gene (considering intronic off-target bins)
#	- Set the updated gene name
#	- Calculate the weighted mean of log2 ratio by GENE 
#	- Calculate the number of bins in each gene (more bins more confidence)
#	- Calculate the CN (considering the sample sex)
bins_genes <- as.data.frame(mergeByOverlaps(bins_ranges, bins_target_ranges)) %>%
	mutate(gene = ifelse(gene == "Antitarget", bins_target_ranges.gene, gene)) %>%
	filter(gene == bins_target_ranges.gene)
bins_genes_ranges <- makeGRangesFromDataFrame(
	bins_genes, start.field = "bins_ranges.start", end.field = "bins_ranges.end", seqnames.field = "bins_ranges.seqnames", 
	ignore.strand = TRUE, keep.extra.columns = TRUE)
genes_ratio <- bins_genes %>% group_by(gene) %>% 
	summarise(
		LOG2 = weighted.mean(log2, weight), DEPTH = unique(DEPTH), WEIGHT = sum(weight),
		BINS_TOTAL = n(), BINS_TARGET = unique(BINS_TARGET), CHROM = unique(bins_ranges.seqnames), 
		START = min(bins_ranges.start), END = max(bins_ranges.end), .groups = "drop") %>% rowwise() %>% 
	mutate(
		aux_log2 = ifelse(is_male & CHROM == "chrX", LOG2+1, LOG2),
		CN = (CN_thresholds%>%filter(data.table::between(aux_log2, LOWER, UPPER)))$CN[1],
		CNA = (CN_thresholds%>%filter(data.table::between(aux_log2, LOWER, UPPER)))$CNA[1],
		CLASS = (CN_thresholds%>%filter(data.table::between(aux_log2, LOWER, UPPER)))$CLASS[1],
		NO_FLAGS = (CN_thresholds%>%filter(data.table::between(aux_log2, LOWER, UPPER)))$OK[1]) %>% as.data.frame()

# Calculate the % of covered gene (MANE) by target bins: ALL and EXON
aux_genes <- subsetByOverlaps(mane_perbase_ranges, bins_ranges[bins_ranges$gene != "Antitarget"]) %>% as.data.table() %>% 
	group_by(GENE_SYMBOL) %>% summarise(BIN_WIDTH = sum(width)) %>% rowwise() %>% 
	mutate(TOTAL_WIDTH = mane_gene_width[GENE_SYMBOL], PERC_COV_GENE_ALL = round(BIN_WIDTH*100/TOTAL_WIDTH,1)) %>% 
	select(GENE_SYMBOL, PERC_COV_GENE_ALL)
genes_ratio <- left_join(genes_ratio, aux_genes, by = c("gene" = "GENE_SYMBOL"))
rm(aux_genes)
aux_genes <- subsetByOverlaps(mane_perbase_exon, bins_ranges[bins_ranges$gene != "Antitarget"]) %>% as.data.table() %>% 
	group_by(GENE_SYMBOL) %>% summarise(BIN_WIDTH = sum(width)) %>% rowwise() %>% 
	mutate(TOTAL_WIDTH = mane_exon_width[GENE_SYMBOL, "WIDTH"], PERC_COV_GENE_EXON = round(BIN_WIDTH*100/TOTAL_WIDTH,1)) %>% 
	select(GENE_SYMBOL, PERC_COV_GENE_EXON)
genes_ratio <- left_join(genes_ratio, aux_genes, by = c("gene" = "GENE_SYMBOL"))
rm(aux_genes)

# Create final output
#	- Exclude genes out of the analysis
#	- LOWBINS: <=3 bins (by default)
#	- FLAGS: OK vs Neutral/LowAmp/LowDel/LowBins
cna_results <- left_join(panel_genes, genes_ratio, by = c("GENE" = "gene")) %>% 
	filter(!is.na(LOG2)) %>% 
	left_join(cancerdrivers, by = "GENE") %>%
	mutate(across(c(ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER), ~ .x & !is.na(.x))) %>%
	rowwise() %>%
	mutate(
		RUN = run, SAMPLE = sample, VAR = paste0(GENE, "_", CLASS), ALTERATION = paste0(GENE, "_", CNA), 
		LOG2 = round(LOG2, 4), DEPTH = round(DEPTH), WEIGHT = round(WEIGHT, 4), BINS_TOTAL = ifelse(!is.na(BINS_TOTAL), BINS_TOTAL, 0), 
		BINS_TARGET = ifelse(!is.na(BINS_TARGET), BINS_TARGET, 0), LOWBINS = BINS_TARGET < min_target_bins,
		LOW_AMP = CLASS %in% "LowAmp", LOW_DEL = CLASS %in% "LowDel", NEUTRAL = CLASS %in% "Neutral", NO_FLAGS = NO_FLAGS & !LOWBINS, 
		FLAGS = paste(CN_flags[c(NO_FLAGS, LOW_AMP, LOW_DEL, LOWBINS, NEUTRAL)], collapse = ";"),
		WHITELIST = GENE %in% genes_whitelist, CLASS = factor(CLASS, levels = CN_levels)) %>% as.data.frame() %>%
	left_join(genie_cna, by = c("GENE", "CNA")) %>% 
	mutate(across(c("GENIE_CNT", "GENIE_FREQ"), ~ ifelse(is.na(.x), 0, as.numeric(.x)))) %>%
	arrange(desc(NO_FLAGS), CLASS, LOWBINS, desc(WHITELIST), GENE) %>% 
	mutate(across(where(is.factor), as.character)) %>%
	select(
		RUN, SAMPLE, ALTERATION, VAR, FLAGS, CYTOBAND, GENE, CN, CNA, CLASS, 
		LOG2, DEPTH, WHITELIST, GENIE_CNT, GENIE_FREQ, BINS_TOTAL, BINS_TARGET, 
		WEIGHT, PERC_COV_GENE_ALL, PERC_COV_GENE_EXON, ONCOGENE, TSG, DRIVER, 
		CANONICAL_DRIVER, GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ,
		CHROM, START, END, STRAND, NO_FLAGS, LOW_AMP, LOW_DEL, LOWBINS, NEUTRAL)

# Clinical classification (AMP/ASCO/CAP guidelines)
#	- Predictive, prognostic & diagnostic evidences in CIViC resource
# 	- Tier IA: CIViC evidence A (guidelines) in the SAME tumor
# 	- Tier IB: CIViC evidence B (trials) in the SAME tumor
# 	- Tier IIC: CIViC evidence A/B in OTHER tumor OR CIViC evidence C
# 	- Tier IID: CIViC evidence D
# 	- Tier III: Observed in pan-cancer resource (GENIE) >= 0.1% Freq
# 	- Tier IV: Not observed in pan-cancer resource (GENIE) or < 0.1% Freq
cna_civic <- cna_results %>% filter(ALTERATION %in% civic_data$CIViC_ALTERATION)
civic_annot <- civic_data %>% mutate(ALTERATION = ".", CLIN_STATUS = ".") %>% filter(CIViC_ALTERATION == ".")
for(n in seq_len(nrow(cna_civic))){
	aux_var <- cna_civic[n, ]
	aux_civic <- civic_data %>% filter(CIViC_ALTERATION %in% aux_var$ALTERATION) %>%
		mutate(
			ALTERATION = aux_var$ALTERATION,
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
cna_results <- cna_results %>% left_join(civic_annot, by = "ALTERATION") %>%
	mutate(
		CLIN_STATUS = ifelse(!is.na(CLIN_STATUS), CLIN_STATUS, ifelse(GENIE_FREQ >= 0.001, "Tier III", "Tier IV")),
		CLIN_RELEVANT = CLIN_STATUS %in% c("Tier IA", "Tier IB", "Tier IIC", "Tier IID")) %>%
	select(RUN:FLAGS, CLIN_RELEVANT, CLIN_STATUS, everything()) %>%
	arrange(CLIN_STATUS)
cna_results[is.na(cna_results) | cna_results == "NA" | cna_results == "NaN" | cna_results == ""] <- "."

# Kepp OK (no flags) CNAs
cna_ok <- cna_results %>% filter(NO_FLAGS)

# Generate a VCF file
aux_header <- readLines(headerVcf.file)
aux_header[grep("##fileDate", aux_header)] <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
aux_header[grep("##reference", aux_header)] <- paste0("##reference=file:", fasta_name)
if(nrow(cna_results %>% filter(!NEUTRAL))){
	aux_vars <- cna_results %>% filter(!NEUTRAL) %>% mutate(
		START = as.integer(START), END = as.integer(END),
		SVTYPE = ifelse(CNA == "AMP", "DUP", ifelse(CNA == "DEL", "DEL", ".")), 
		ALT = ifelse(SVTYPE == ".", SVTYPE, paste0("<", SVTYPE, ">")))
	aux_vars$CHROM <- factor(aux_vars$CHROM, levels = chr_levels)
	aux_vars <- aux_vars %>% arrange(CHROM, START, END)
	aux_ranges <- makeGRangesFromDataFrame(aux_vars, start.field = "START", end.field = "START", seqnames.field = "CHROM")
	aux_ref <- as.data.frame(getSeq(fasta, aux_ranges))$x
	aux_info <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c(
			"ALTERATION", "VAR", "RUN", "SAMPLE", "CLIN_STATUS", "SVTYPE", "CYTOBAND", "GENE", "CN", "CNA", "CLASS", 
			"LOG2", "DEPTH", "GENIE_CNT", "GENIE_FREQ", "BINS_TOTAL", "BINS_TARGET", "WEIGHT", "PERC_COV_GENE_ALL", "PERC_COV_GENE_EXON", 
			"GENE_ENSEMBL", "TRANSCRIPT_ENSEMBL", "TRANSCRIPT_REFSEQ", "CHROM", "START", "END", "STRAND", 
			"CIViC_ALTERATION", "CIViC_VARIANT_ID", "CIViC_EVIDENCE_LEVEL", "CIViC_EVIDENCE_RATING", "CIViC_EVIDENCE_TYPE", 
			"CIViC_EFFECT", "CIViC_DRUG", "CIViC_TUMOR_NAME", "CIViC_TUMOR_DOID", "CIViC_TUMOR_CODE", "CIViC_TOPNODE_NAME",
			"CIViC_TOPNODE_DOID", "CIViC_MOLECULAR_ID", "CIViC_VARIANT_GROUP", "CIViC_EVIDENCE_SCORE", "CIViC_ORIGIN"),
			function(x) paste0(x, "=", aux_vars[n, x]) )), collapse = ";")  }))
	aux_flags <- paste0(aux_vars$FLAGS, ";", unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("WHITELIST", "ONCOGENE", "TSG", "DRIVER", "CANONICAL_DRIVER", "CLIN_RELEVANT"), 
			function(x) ifelse(aux_vars[n, x], paste0(x, ";"), "")  )), collapse = "")  })))
	aux_info <- paste0(aux_info, ";", aux_flags)
	aux_format <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("LOG2", "CN"), function(x) paste0(aux_vars[n, x]) )), collapse = ":")  }))
	aux_df <- aux_vars %>% mutate(
		CHROM = as.character(CHROM), POS = START, ID = ".", REF = aux_ref, QUAL = ".", FILTER = "PASS", INFO = aux_info,
		FORMAT = paste(c("LOG2", "CN"), collapse = ":"), SAMPLE = aux_format) %>%
		select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)
	vcf_text <- c(aux_header,
		paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"),
		unlist(lapply(seq_len(nrow(aux_df)), function(x) paste(aux_df[x, ], collapse = "\t"))))

	rm(aux_header, aux_vars, aux_ranges, aux_ref, aux_info, aux_flags, aux_format, aux_df)
} else {
	vcf_text <- c(aux_header, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"))
}



# Save data ---------------------------------------------------------------

write.table(cna_results, file = paste0(sample_id, "_CNA_gene.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	
addWorksheet(results_excel, "OK")
writeDataTable(results_excel, "OK", cna_ok)
addWorksheet(results_excel, "Gene")
writeDataTable(results_excel, "Gene", cna_results)
if(length(ascets_results.file)){
	addWorksheet(results_excel, "Arm")
	writeDataTable(results_excel, "Arm", arm_results)
}
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

writeLines(con = paste0(sample_id, "_CNA.vcf"), text = vcf_text)


# Data visualization ---------------------------------------------------------------

# Set ranges
#	- Exclude LOWBINS genes (excepting clinical relevant CNAs)
#   - Set whitelist to clinical CNAs to show them in the plot
cna_plots <- cna_results %>% filter(CLIN_RELEVANT | !LOWBINS) %>% mutate(
	LOG2 = as.numeric(LOG2), CLASS = ifelse(CLIN_RELEVANT, "AMP/ASCO/CAP Tier I/II", CLASS),
	WHITELIST = WHITELIST | CLIN_RELEVANT | plot_global_label_all_genes, PLOT = WHITELIST | CLIN_RELEVANT | plot_bygene_all)
genes_ranges <- makeGRangesFromDataFrame(cna_plots[rev(seq_len(nrow(cna_plots))),], keep.extra.columns = TRUE, ignore.strand = TRUE)

# Set colors and shapes
bin_colors <- unlist(lapply(genes_ranges$CLASS, function(x) CN_colors[x]))
whitelist_colors <- unlist(lapply(genes_ranges[genes_ranges$WHITELIST]$CLASS, function(x) CN_colors[x]))


# Genome visualization
png(paste0(sample_id, "_CNA.png"), width = 3000, height = 1600, res = 200)

plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$data2height <- 0
plot.params$topmargin <- 5
plot.params$bottommargin <- 90

kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 1.4, plot.params = plot.params, chromosomes = paste0("chr", c(1:22, "X")))
kpAddChromosomeNames(kp, srt = 70, cex = 1.3, yoffset = -235, xoffset = -0.003)
kpAddChromosomeSeparators(kp, lwd = 0.3, col = "#666666")

ymin <- suppressWarnings(floor(min(genes_ranges$LOG2)))
ymin <- suppressWarnings(floor(ifelse(is_male, ifelse(ymin>log2_DEL-1, log2_DEL-1, ymin), ifelse(ymin>log2_DEL, log2_DEL, ymin))))
ymax <- suppressWarnings(ceiling(max(genes_ranges$LOG2))); ymax <- ifelse(ymax<log2_AMP, log2_AMP, ymax)
kpRect(
	kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
	data.panel = 1, r0 = 0, r1 = 0.9)
kpSegments(
	kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
	data.panel = 1, r0 = 0, r1 = 0.9)
kpSegments(
	kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
	data.panel = 1, r0 = 0, r1 = 0.9)
if(ymin<=log2_DEL_OK){kpSegments(
	kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
	data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
if(ymax>=log2_AMP_OK){kpSegments(
	kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
	data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
if(is_male){
	kpRect(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_AMP-1, border = NA, col = "grey97", 
		data.panel = 1, r0 = 0, r1 = 0.9)
	kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_DEL-1, col = CN_colors["LowDel"], 
		data.panel = 1, r0 = 0, r1 = 0.9)
	kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP-1, y1 = log2_AMP-1, col = CN_colors["LowAmp"], 
		data.panel = 1, r0 = 0, r1 = 0.9)
	if(ymin<=log2_DEL_OK-1){kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK-1, y1 = log2_DEL_OK-1, col = CN_colors["HighDel"], 
		data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
	if(ymax>=log2_AMP_OK-1){kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK-1, y1 = log2_AMP_OK-1, col = CN_colors["HighAmp"], 
		data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
} else {
	kpRect(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
		data.panel = 1, r0 = 0, r1 = 0.9)
	kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
		data.panel = 1, r0 = 0, r1 = 0.9)
	kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
		data.panel = 1, r0 = 0, r1 = 0.9)
	if(ymin<=log2_DEL_OK){kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
		data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
	if(ymax>=log2_AMP_OK){kpSegments(
		kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
		data.panel = 1, r0 = 0, r1 = 0.9, lwd = 0.5)}
}
if(any(genes_ranges$WHITELIST)){kpPlotMarkers(
	kp, genes_ranges[genes_ranges$WHITELIST], labels = genes_ranges[genes_ranges$WHITELIST]$GENE, 
	y = 1, cex = 0.8, lwd = 0.5, adjust.label.position = TRUE, label.margin = 0.01, label.dist = 0.001,
	data.panel = 1, r0 = -0.05, r1 = 0.9, ignore.chromosome.ends = TRUE, marker.parts = c(0.95, 0.045, 0.005),
	line.color = whitelist_colors, label.color = whitelist_colors)
}
kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 0.4), cex = 0.8, data.panel = 1, r0 = 0, r1 = 0.9, tick.len = 10e6)
plotLRR(
	kp, genes_ranges, lrr.column = "LOG2", labels = "Copy ratio (log2)", ymin = ymin, ymax = ymax, label.margin = 0.035,
	line.at.0 = TRUE, line.at.0.col = "black", data.panel = 1, out.of.range = "points", out.of.range.col = "white",
	points.cex = 0.6, points.col = "black", points.pch = 1, label.cex = 1.4, add.axis = FALSE, r0 = 0, r1 = 0.9)
kpPoints(
	kp, genes_ranges, y = genes_ranges$LOG2, ymin = ymin, ymax = ymax, data.panel = 1, r0 = 0, r1 = 0.9, 
	cex = 0.6, pch = 16, col = bin_colors)	

if(length(ascets_results.file)){
	if(any(arm_results_ranges$CNA %in% names(arm_colors))){
		suppressWarnings(kpPlotRegions(
			kp, arm_results_ranges[arm_results_ranges$CNA == "AMP"], col = arm_colors["AMP"], data.panel = 1, r0 = -0.01, r1 = -0.04))
		suppressWarnings(kpPlotRegions(
			kp, arm_results_ranges[arm_results_ranges$CNA == "DEL"], col = arm_colors["DEL"], data.panel = 1, r0 = -0.01, r1 = -0.04))
		suppressWarnings(kpPlotRegions(
			kp, arm_results_ranges[arm_results_ranges$CNA == "CONFLICT"], col = arm_colors["CONFLICT"], data.panel = 1, r0 = -0.01, r1 = -0.04))
		
		legend("bottomleft", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
		legend("bottomright", legend = names(arm_colors), ncol = 3, fill = arm_colors, border = arm_colors, cex = 0.9, pt.cex = 1.2, title = "Arm CNAs")
	} else {
		legend("bottom", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
	}
} else {
	legend("bottom", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
}

rm(kp)

dev.off()

addWorksheet(results_excel, "Plot_Global")
insertImage(results_excel, "Plot_Global", paste0(sample_id, "_CNA.png"), units = "px", width = 3000, height = 1600, dpi = 200)
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)


# By chromosome visualization
addWorksheet(results_excel, "Plot_ByChr")

plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$data2height <- 0
plot.params$topmargin <- 30
plot.params$bottommargin <- 100

excel_startRow <- 1

for(chr in paste0("chr", c(1:22, "X"))){

	is_gene <- as.logical(seqnames(genes_ranges) == chr)
	aux_bins <- genes_ranges[is_gene]
	aux_colors <- bin_colors[is_gene]

	if(!length(aux_bins)){next}

	png(paste0(sample_id, "_CNA_ByChr_", chr, ".png"), width = 3000, height = 1600, res = 200)
	
	kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 1.5, chromosomes = chr, plot.params = plot.params)
	kpAddBaseNumbers(kp, cex = 1, minor.ticks = TRUE, tick.len = 5, minor.tick.len = 2, add.units = TRUE)
	kpAddChromosomeNames(kp, cex = 1.5, yoffset = -225, xoffset = -0.47)	
	kpAddMainTitle(kp, main = chr, cex = 2)

	kpRect(
		kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
		data.panel = 1, r0 = 0.25, r1 = 1)
	kpSegments(
		kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
		data.panel = 1, r0 = 0.25, r1 = 1)
	kpSegments(
		kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
		data.panel = 1, r0 = 0.25, r1 = 1)
	if(ymin<=log2_DEL_OK){kpSegments(
		kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
		data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
	if(ymax>=log2_AMP_OK){kpSegments(
		kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
		data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
	if(is_male & chr == "chrX"){
		kpRect(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_AMP-1, border = NA, col = "grey97", 
			data.panel = 1, r0 = 0.25, r1 = 1)
		kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_DEL-1, col = CN_colors["LowDel"], 
			data.panel = 1, r0 = 0.25, r1 = 1)
		kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP-1, y1 = log2_AMP-1, col = CN_colors["LowAmp"], 
			data.panel = 1, r0 = 0.25, r1 = 1)
		if(ymin<=log2_DEL_OK-1){kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK-1, y1 = log2_DEL_OK-1, col = CN_colors["HighDel"], 
			data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
		if(ymax>=log2_AMP_OK-1){kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK-1, y1 = log2_AMP_OK-1, col = CN_colors["HighAmp"], 
			data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
	} else {
		kpRect(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
			data.panel = 1, r0 = 0.25, r1 = 1)
		kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
			data.panel = 1, r0 = 0.25, r1 = 1)
		kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
			data.panel = 1, r0 = 0.25, r1 = 1)
		if(ymin<=log2_DEL_OK){kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
			data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
		if(ymax>=log2_AMP_OK){kpSegments(
			kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
			data.panel = 1, r0 = 0.25, r1 = 1, lwd = 0.5)}
	}

	suppressWarnings(kpPlotMarkers(
		kp, aux_bins, labels = aux_bins$GENE, y = 1, cex = 0.9, adjust.label.position = TRUE, ignore.chromosome.ends = TRUE,
		data.panel = 1, r0 = -0.05, r1 = 0.04, label.margin = 0.01, label.dist = 0.001, marker.parts = c(0.95, 0.045, 0.005),
		line.color = aux_colors, label.color = aux_colors))
	
	kpAxis(kp, ymin = ymin, ymax = ymax, tick.pos = seq(ymin, ymax, 0.4), cex = 0.7, data.panel = 1, r0 = 0.25, r1 = 1, tick.len = 50e4)
	plotLRR(
		kp, aux_bins, lrr.column = "LOG2", labels = "Copy ratio (log2)", ymin = ymin, ymax = ymax, label.margin = 0.035,
		line.at.0 = TRUE, line.at.0.col = "black", data.panel = 1, out.of.range = "points", out.of.range.col = "white",
		points.cex = 1, points.col = aux_colors, points.pch = 16, label.cex = 1.2, add.axis = FALSE, r0 = 0.25, r1 = 1)
	kpPoints(
		kp, aux_bins, y = aux_bins$LOG2, ymin = ymin, ymax = ymax, data.panel = 1, r0 = 0.25, r1 = 1, 
		cex = 1, pch = 1, col = "black")
	
	if(length(ascets_results.file)){
		if(any(arm_results_ranges[seqnames(arm_results_ranges) == chr]$CNA %in% names(arm_colors))){
			suppressWarnings(kpPlotRegions(
				kp, arm_results_ranges[arm_results_ranges$CNA == "AMP"], col = arm_colors["AMP"], data.panel = 1, r0 = -0.01, r1 = -0.047))
			suppressWarnings(kpPlotRegions(
				kp, arm_results_ranges[arm_results_ranges$CNA == "DEL"], col = arm_colors["DEL"], data.panel = 1, r0 = -0.01, r1 = -0.047))
			suppressWarnings(kpPlotRegions(
				kp, arm_results_ranges[arm_results_ranges$CNA == "CONFLICT"], col = arm_colors["CONFLICT"], data.panel = 1, r0 = -0.01, r1 = -0.047))
	
			legend("bottomleft", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
			legend("bottomright", legend = names(arm_colors), ncol = 3, fill = arm_colors, border = arm_colors, cex = 0.9, pt.cex = 1.2, title = "Arm CNAs")
		} else {
			legend("bottom", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
		}
	} else {
		legend("bottom", legend = CN_legend, col = CN_colors, ncol = 3, cex = 0.9, title = "Gene CNAs", pch = 16, pt.cex = 1.2)
	}
	
	dev.off()

	insertImage(
		results_excel, "Plot_ByChr", paste0(sample_id, "_CNA_ByChr_", chr, ".png"), 
		units = "px", width = 3000, height = 1600, dpi = 200, startRow = excel_startRow)
	saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

	excel_startRow <- excel_startRow + 35

	rm(kp, aux_bins, aux_colors)

}
rm(chr, excel_startRow)


# By gene (from Whitelist)
if(any(cna_plots$PLOT)) {
	
	addWorksheet(results_excel, "Plot_ByGene")

	# Setup karyoplot parameters
	plot.params <- getDefaultPlotParams(plot.type = 3)
	plot.params$data2height <- 0
	plot.params$topmargin <- 30
	plot.params$leftmargin <- 0.06

	excel_startRow <- 1

	for(gene in cna_plots[cna_plots$PLOT, "GENE"]){

		aux_bins <- bins_genes_ranges[bins_genes_ranges$gene == gene]
		aux_gene_cna <- genes_ranges[genes_ranges$GENE == gene]; aux_gene_cna$y0 <- aux_gene_cna$LOG2; aux_gene_cna$y1 <- aux_gene_cna$LOG2
		aux_colors <- bin_colors[genes_ranges$GENE == gene]

		if(!length(aux_bins)){next}

		png(paste0(sample_id, "_CNA_ByGene_", gene, ".png"), width = 3000, height = 1200, res = 200)
		
		aux_chr <- unique(seqnames(aux_bins))
		aux_gene_ranges <- mane_gene_ranges[mane_gene_ranges$GENE_SYMBOL == gene]
		aux_exon_ranges <- mane_exon_ranges[mane_exon_ranges$GENE_SYMBOL == gene]
		aux_gene <- mane_gene_names[unique(aux_exon_ranges$GENE_SYMBOL)]
		transcript <- unique(aux_exon_ranges$TRANSCRIPT_ENSEMBL)

		aux_ranges <- c(bins_target_ranges[bins_target_ranges$gene == gene], aux_gene_ranges)
		aux_region <- toGRanges(unique(seqnames(aux_ranges)), min(start(aux_ranges)), max(end(aux_ranges)))
		width <- width(aux_region)
		zoom <- aux_region + width*0.05
		aux_perBase_ranges <- subsetByOverlaps(perBase_ranges, zoom)

		kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 2, plot.params = plot.params, zoom = zoom)
		kpAddBaseNumbers(kp, cex = 1, minor.ticks = TRUE, tick.len = 5, minor.tick.len = 2, add.units = TRUE, digits = 3,
			tick.dist=round(width(zoom)/5, -2), minor.tick.dist=round(width(zoom)/50, -2))
		kpAddChromosomeNames(kp, cex = 1.5, yoffset = -225, xoffset = -0.47)
		kpAddMainTitle(kp, main = paste0(gene, " (", transcript, ")"), cex = 2)

		ymin <- suppressWarnings(floor(min(aux_bins$log2)))
		ymin <- suppressWarnings(floor(ifelse(is_male & aux_chr == "chrX", ifelse(ymin>log2_DEL-1, log2_DEL-1, ymin), ifelse(ymin>log2_DEL, log2_DEL, ymin))))
		ymax <- suppressWarnings(ceiling(max(aux_bins$log2)))
		ymax <- suppressWarnings(ceiling(ifelse(is_male & aux_chr == "chrX", ifelse(ymax<log2_AMP-1, log2_AMP-1, ymax), ifelse(ymax<log2_AMP, log2_AMP, ymax))))
		kpRect(
			kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
			data.panel = 1, r0 = 0.60, r1 = 0.95)
		kpSegments(
			kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
			data.panel = 1, r0 = 0.60, r1 = 0.95)
		kpSegments(
			kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
			data.panel = 1, r0 = 0.60, r1 = 0.95)
		if(ymin<=log2_DEL_OK){kpSegments(
			kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
			data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
		if(ymax>=log2_AMP_OK){kpSegments(
			kp, arm_coord_autosomal_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
			data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
		if(is_male & aux_chr == "chrX"){
			kpRect(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_AMP-1, border = NA, col = "grey97", 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL-1, y1 = log2_DEL-1, col = CN_colors["LowDel"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP-1, y1 = log2_AMP-1, col = CN_colors["LowAmp"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			if(ymin<=log2_DEL_OK-1){kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK-1, y1 = log2_DEL_OK-1, col = CN_colors["HighDel"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
			if(ymax>=log2_AMP_OK-1){kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK-1, y1 = log2_AMP_OK-1, col = CN_colors["HighAmp"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
		} else {
			kpRect(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_AMP, border = NA, col = "grey97", 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL, y1 = log2_DEL, col = CN_colors["LowDel"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP, y1 = log2_AMP, col = CN_colors["LowAmp"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95)
			if(ymin<=log2_DEL_OK){kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_DEL_OK, y1 = log2_DEL_OK, col = CN_colors["HighDel"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
			if(ymax>=log2_AMP_OK){kpSegments(
				kp, arm_coord_X_ranges, ymin = ymin, ymax = ymax, y0 = log2_AMP_OK, y1 = log2_AMP_OK, col = CN_colors["HighAmp"], 
				data.panel = 1, r0 = 0.60, r1 = 0.95, lwd = 0.5)}
		}
		kpAxis(kp, ymin = ymin, ymax = ymax, numticks = 5, cex = 0.7, data.panel = 1, r0 = 0.60, r1 = 0.95, tick.len = width/200)

		genes.data <- makeGenesDataFromTxDb(mane_TxDb, kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
		is_gene <- names(genes.data$genes) == aux_gene
		genes.data$genes <- genes.data$genes[is_gene]
		genes.data$transcripts <- genes.data$transcripts[is_gene]
		is_transcript <- genes.data[["transcripts"]][[aux_gene]]$tx_name == transcript
		genes.data[["transcripts"]][[aux_gene]] <- genes.data[["transcripts"]][[aux_gene]][is_transcript]
		genes.data[["coding.exons"]] <- genes.data[["coding.exons"]][is_transcript]
		genes.data[["non.coding.exons"]] <- genes.data[["non.coding.exons"]][is_transcript]

		kpPlotMarkers(
			kp, aux_exon_ranges, labels = aux_exon_ranges$EXON, adjust.label.position = TRUE, label.dist = 0.0003, lwd = 0.4, cex = 0.7,
			text.orientation = "horizontal", marker.parts = c(0.95, 0.045, 0.005), data.panel = 1, r0 = 0, r1 = 1.3, 
			line.color = "grey75", label.color = "black", ignore.chromosome.ends = TRUE)
		suppressWarnings(kpPlotGenes(
			kp, genes.data, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.55, r1 = 0.60, add.strand.marks=TRUE, mark.height = 0.25,
			coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50", non.coding.exons.border.col="grey50"))
		
		plotLRR(kp, bins_ranges[NULL], lrr.column = "log2", ymin = ymin, ymax = ymax, data.panel = 1, r0 = 0.60, r1 = 0.95, label.margin = 0.035,
				labels = "Copy ratio (log2)", label.cex = 1.2, line.at.0 = TRUE, line.at.0.col = "black", add.axis = FALSE)
		kpSegments(kp, aux_gene_cna, ymin = ymin, ymax = ymax, col = aux_colors, lwd = 6, data.panel = 1, r0 = 0.60, r1 = 0.95)
		kpSegments(kp, aux_bins, ymin = ymin, ymax = ymax, col = "grey70", lwd = 5, data.panel = 1, r0 = 0.60, r1 = 0.95)
		
		ymax <- plyr::round_any(max(aux_perBase_ranges$y1), 100, f = ceiling)
		ymax_sep <- plyr::round_any(ymax/5, 100, f = ceiling)
		kpAxis(kp, ymin = 0, ymax = ymax, cex = 0.9, data.panel = 1, r0 = 0.05, r1 = 0.5, tick.pos = seq(0, ymax, ymax_sep), tick.len = width/200)
		kpAddLabels(kp, labels = "Coverage", cex = 1.2, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.05, r1 = 0.5)
		kpAbline(kp, h = 0, data.panel = 1, r0 = 0.05, r1 = 0.5, col = "black", ymin = 0, ymax = ymax)
		kpBars(kp, aux_perBase_ranges, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.05, r1 = 0.5, col = "black", border = "black")

		kpAddChromosomeNames(kp, cex = 1, yoffset = -210, xoffset = -0.47, chr.names = "Targets")
		kpPlotRegions(kp, target_ranges, col = "black", data.panel = 1, r0 = 0, r1 = 0.04)

		dev.off()

		insertImage(
			results_excel, "Plot_ByGene", paste0(sample_id, "_CNA_ByGene_", gene, ".png"), 
			units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow)
		saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)
		
		excel_startRow <- excel_startRow + 31

	rm(
		aux_bins, aux_chr, aux_gene_cna, aux_colors, aux_gene_ranges, aux_exon_ranges, aux_gene, transcript, aux_ranges, aux_region, width, 
		zoom, aux_perBase_ranges, kp, ymin, ymax, genes.data, is_gene, is_transcript, ymax_sep)
	
	}
	rm(gene, excel_startRow)	
}

