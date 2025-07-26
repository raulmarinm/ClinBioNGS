#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of BAM QC results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(GenomicFeatures)))
suppressWarnings(suppressPackageStartupMessages(library(karyoploteR)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(openxlsx)))
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
sample.file <- args_list[["samplesFile"]]
target.file <- args_list[["targetBed"]]
genes.file <- args_list[["targetGenes"]]
whitelist.file <- args_list[["whitelistGenes"]]
mane.gtf <- args_list[["maneGtf"]]
mane_gene.file <- args_list[["maneGene"]]
mane_exon.file <- args_list[["maneExon"]]
mane_coding.file <- args_list[["maneCoding"]]
pct_X_coverage <- as.numeric(unlist(strsplit(args_list[["pctXCoverage"]], ",")))
low_coverage <- as.numeric(args_list[["lowCoverage"]])
plot_mark_chr_coverage <- as.logical(args_list[["plotMarkChrCoverage"]])
plot_global_label_all_genes <- as.logical(args_list[["plotGlobalLabelAllGenes"]])
plot_bygene_all <- as.logical(args_list[["plotByGeneAll"]])

# Set other variables
sample_id <- paste0(sample, "_", type)
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
whitelist_colors_1 <- c("FALSE" = "grey", "TRUE" = "orange")
whitelist_colors_2 <- c("FALSE" = "black", "TRUE" = "orange")
lowGenes_colors <- c("FALSE" = "black", "TRUE" = "red")
results_excel <- createWorkbook()
results_excel.file <- paste0(run, "_", sample_id, "_Metrics_", format(Sys.time(), "%Y%m%d"), ".xlsx")

# Set metrics files
alignmentSummary.file <- list.files(pattern = "AlignmentSummaryMetrics.txt")
insertSize.file <- list.files(pattern = "InsertSizeMetrics.txt")
onTargetReads.file <- list.files(pattern = "OnTargetReads.txt")
perBaseCoverage.file <- list.files(pattern = "per-base.bed.gz$")
raw_bamqc.file <- list.files(pattern = "RawBamQC.txt")


# Load data ---------------------------------------------------------------

# Sample info
sample_info <- read.table(sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Sample == sample)
tumor_whitelist <- sample_info[1, "Whitelist"]

# Genes
panel_genes <- read.table(genes.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
	mutate(NEW_GENE = ifelse(!is.na(NEW_GENE), NEW_GENE, OLD_GENE))
panel_genes_names <- panel_genes$NEW_GENE; names(panel_genes_names) <- panel_genes$OLD_GENE
panel_genes <- panel_genes %>% mutate(GENE = NEW_GENE) %>% select(GENE, GENE_ENSEMBL, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL)
genes_whitelist <- read.table(whitelist.file, header = TRUE, sep = ";", fill = TRUE, na.strings = c("NA", ".", ""))
genes_whitelist <- as.character(na.omit(unique(genes_whitelist[[tumor_whitelist]])))
genes_whitelist <- genes_whitelist[genes_whitelist %in% panel_genes_names]

# Target regions
target_bed <- read.table(target.file, header = FALSE, sep = "\t", col.names = c("CHROM", "START", "END", "GENE")) %>%
	filter(GENE %in% names(panel_genes_names)) %>%
	mutate(GENE = panel_genes_names[GENE], LOCUS = paste0(CHROM, ":", START, "-", END))
target_ranges <- makeGRangesFromDataFrame(target_bed, ignore.strand = TRUE, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
target_perbase_ranges <- makeGRangesFromDataFrame(
	target_ranges %>% as.data.table() %>% select(seqnames:GENE, LOCUS) %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
	mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.table() %>% select(-ID, -pos), keep.extra.columns = TRUE)

# MANE data (only target genes)
mane_TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(mane.gtf, format = "gtf")))
mane_gene <- fread(mane_gene.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% 
	filter(GENE_SYMBOL %in% panel_genes_names) %>% select(-GENE_NAME)
mane_gene_ranges <- makeGRangesFromDataFrame(mane_gene, keep.extra.columns = TRUE)
mane_perbase_ranges <- makeGRangesFromDataFrame(
	mane_gene_ranges %>% as.data.table() %>% select(seqnames:GENE_SYMBOL) %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
	mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.table() %>% select(-ID, -pos), keep.extra.columns = TRUE)
mane_exon <- fread(mane_exon.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% 
	filter(GENE_SYMBOL %in% panel_genes_names) %>% arrange(CHROM, START, END)
mane_exon_ranges <- makeGRangesFromDataFrame(mane_exon, keep.extra.columns = TRUE)
mane_coding <- fread(mane_coding.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% 
	filter(GENE_SYMBOL %in% panel_genes_names) %>% arrange(CHROM, START, END)
mane_coding_ranges <- makeGRangesFromDataFrame(mane_coding, keep.extra.columns = TRUE)

# Per-base coverage
perBase_ranges <- makeGRangesFromDataFrame(fread(
	perBaseCoverage.file, header = FALSE, col.names = c("CHROM", "START", "END", "DEPTH"),
	sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% filter(CHROM %in% chr_levels),
	starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
perBase_ranges$y1 <- perBase_ranges$DEPTH
target_perbase_coverage <- mergeByOverlaps(target_perbase_ranges, perBase_ranges) %>% as.data.table() %>% mutate(
	CHROM = target_perbase_ranges.seqnames, START = target_perbase_ranges.start, END = target_perbase_ranges.end,
	LENGTH = target_perbase_ranges.width)
target_genes <- as.data.frame(target_ranges) %>% group_by(GENE) %>% 
	summarise(CHROM = unique(seqnames), START = mean(c(start, end)), END = START, STRAND = unique(strand)) %>%
	filter(GENE %in% unique(target_perbase_coverage$GENE))

# BAM QC metrics
if (length(raw_bamqc.file)) {
	raw_bamqc <- read.table(raw_bamqc.file, header = TRUE, sep = "\t") %>% mutate(RUN = as.character(RUN), SAMPLE = as.character(SAMPLE))
}
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

# BAM coverage metrics
coverage_metrics <- target_perbase_coverage %>% summarise(
	MIN_COVERAGE = min(DEPTH), MAX_COVERAGE = max(DEPTH), MEDIAN_COVERAGE = round(median(DEPTH)), MEAN_COVERAGE = round(mean(DEPTH)), 
	PANEL_SIZE = sum(LENGTH), PCT_0.4X_MEAN = round((sum(LENGTH[DEPTH >= 0.4*MEAN_COVERAGE])*100/PANEL_SIZE)))
for(cov in pct_X_coverage){
	aux_df <- target_perbase_coverage %>%
		summarise(PCT_X = sum(LENGTH[DEPTH >= cov]), LENGTH = sum(LENGTH)) %>% mutate(PCT_X = round(PCT_X*100/LENGTH))
	coverage_metrics[[paste0("PCT_", cov, "X")]] <- aux_df$PCT_X
	rm(aux_df)
}
rm(cov)
coverage_metrics <- coverage_metrics %>% mutate(SAMPLE = sample) %>%
	select(SAMPLE, MEDIAN_COVERAGE, MEAN_COVERAGE, contains("PCT_"), MIN_COVERAGE, MAX_COVERAGE, PANEL_SIZE)



# Data processing ---------------------------------------------------------

# Join BAM QC metrics
bamqc_metrics <- onTargetReads %>% left_join(alignmentSummary, by = "SAMPLE") %>% left_join(insertSize, by = "SAMPLE") %>% mutate(RUN = run)
if (length(raw_bamqc.file)) {
	bamqc_metrics <- bamqc_metrics %>%
		mutate(
			UNIQUE_READS = PF_READS, UNIQUE_ALIGNED_READS = PF_READS_ALIGNED,
			UNIQUE_ONTARGET_READS = ONTARGET_READS) %>%
		select(-TOTAL_READS, -ONTARGET_READS, -MEDIAN_READ_LENGTH, -MEDIAN_INSERT_SIZE) %>%
		left_join(raw_bamqc, by = c("RUN", "SAMPLE"))		
} else {
	bamqc_metrics <- bamqc_metrics %>% mutate(
		TOTAL_READS = PF_READS, ALIGNED_READS = PF_READS_ALIGNED, PCT_ALIGNED_READS = round(ALIGNED_READS*100/TOTAL_READS, 1), 
		PCT_ONTARGET_ALIGNED_READS = round(ONTARGET_READS*100/ALIGNED_READS, 1),
		HQ_ALIGNED_READS = PF_HQ_ALIGNED_READS, PCT_HQ_ALIGNED_READS = round(HQ_ALIGNED_READS*100/TOTAL_READS, 1),
		UNIQUE_READS = TOTAL_READS, UNIQUE_ALIGNED_READS = ALIGNED_READS, UNIQUE_ONTARGET_READS = ONTARGET_READS)
}
bamqc_metrics <- bamqc_metrics %>% left_join(coverage_metrics, by = "SAMPLE") %>%
	mutate(
		PCT_DUPLICATION = round(100 - (UNIQUE_READS*100/TOTAL_READS), 1),
		PCT_UNIQUE_ALIGNED_READS = round(UNIQUE_ALIGNED_READS*100/UNIQUE_READS, 1),
		PCT_UNIQUE_ONTARGET_ALIGNED_READS = round(UNIQUE_ONTARGET_READS*100/UNIQUE_ALIGNED_READS, 1)) %>%
	select(
		RUN, SAMPLE, TOTAL_READS, ALIGNED_READS, PCT_ALIGNED_READS, ONTARGET_READS, PCT_ONTARGET_ALIGNED_READS,
		HQ_ALIGNED_READS, PCT_HQ_ALIGNED_READS, MEDIAN_READ_LENGTH, MEDIAN_INSERT_SIZE, UNIQUE_READS, PCT_DUPLICATION, 
		UNIQUE_ALIGNED_READS, PCT_UNIQUE_ALIGNED_READS, UNIQUE_ONTARGET_READS, PCT_UNIQUE_ONTARGET_ALIGNED_READS,
		colnames(coverage_metrics)[-1])

# Coverage by TARGET region
target_coverage <- target_perbase_coverage %>% group_by(LOCUS) %>% summarise(
	CHROM = unique(CHROM), START = min(START), END = max(END), GENE = unique(GENE), LENGTH = sum(LENGTH), 
	MIN_COVERAGE = min(DEPTH), MAX_COVERAGE = max(DEPTH), DEPTH = round(mean(DEPTH)))
for(cov in pct_X_coverage){
	aux_df <- target_perbase_coverage %>% group_by(LOCUS) %>%
		summarise(PCT_X = sum(LENGTH[DEPTH >= cov]), LENGTH = sum(LENGTH)) %>% mutate(PCT_X = round(PCT_X*100/LENGTH)) %>%
		ungroup() %>% select(LOCUS, PCT_X)
	colnames(aux_df) <- c("LOCUS", paste0("PCT_", cov, "X"))
	target_coverage <- left_join(target_coverage, aux_df, by = "LOCUS")
	rm(aux_df)
}
rm(cov)
target_coverage <- target_coverage %>% arrange(CHROM, START, END) %>%
	select(CHROM, START, END, GENE, DEPTH, LENGTH, contains("PCT_"), MIN_COVERAGE, MAX_COVERAGE)

# Coverage by CHROMOSOME (TARGET regions)
chrom_coverage <- target_perbase_coverage %>% group_by(CHROM) %>% summarise(
	START = min(START), END = max(END), MIN_COVERAGE = min(DEPTH), 
	MAX_COVERAGE = max(DEPTH), DEPTH = round(mean(DEPTH)), PANEL_SIZE = sum(LENGTH))
for(cov in pct_X_coverage){
	aux_df <- target_perbase_coverage %>% group_by(CHROM) %>%
		summarise(PCT_X = sum(LENGTH[DEPTH >= cov]), LENGTH = sum(LENGTH)) %>% mutate(PCT_X = round(PCT_X*100/LENGTH)) %>%
		ungroup() %>% select(CHROM, PCT_X)
	colnames(aux_df) <- c("CHROM", paste0("PCT_", cov, "X"))
	chrom_coverage <- left_join(chrom_coverage, aux_df, by = "CHROM")
	rm(aux_df)
}
rm(cov)
chrom_coverage <- chrom_coverage %>% select(CHROM, DEPTH, PANEL_SIZE, contains("PCT_"), MIN_COVERAGE, MAX_COVERAGE, START, END)

# Coverage by GENE: 4 levels
#   - TARGET: TSO500 regions
#   - MANE_CODING: CODING (CDS)
#   - MANE_EXON: EXONS
#   - MANE_ALL: EXONS + INTRONS
gene_DP_target <- target_perbase_coverage %>% group_by(GENE) %>% 
	summarise(MIN_COVERAGE = min(DEPTH), MAX_COVERAGE = max(DEPTH), DP_TARGET = round(mean(DEPTH)), PANEL_SIZE = sum(LENGTH)) %>% 
	ungroup()
gene_DP_mane_coding <- mergeByOverlaps(subsetByOverlaps(mane_perbase_ranges, mane_coding_ranges), perBase_ranges) %>% as.data.table() %>% 
	group_by(GENE_SYMBOL) %>% summarise(DP_MANE_CODING = round(mean(DEPTH))) %>% ungroup()
gene_DP_mane_exon <- mergeByOverlaps(subsetByOverlaps(mane_perbase_ranges, mane_exon_ranges), perBase_ranges) %>% as.data.table() %>% 
	group_by(GENE_SYMBOL) %>% summarise(DP_MANE_EXON = round(mean(DEPTH))) %>% ungroup()
gene_DP_mane <- mergeByOverlaps(mane_perbase_ranges, perBase_ranges) %>% as.data.table() %>% 
	group_by(GENE_SYMBOL) %>% summarise(DP_MANE_ALL = round(mean(DEPTH))) %>% ungroup()
gene_DP_thresholds <- panel_genes %>% select(GENE)
for(cov in pct_X_coverage){
	aux_df <- target_perbase_coverage %>% group_by(GENE) %>% 
		summarise(PCT_X = sum(LENGTH[DEPTH >= cov]), LENGTH = sum(LENGTH)) %>% mutate(PCT_X = round(PCT_X*100/LENGTH)) %>%
		ungroup() %>% select(GENE, PCT_X)
	colnames(aux_df) <- c("GENE", paste0("PCT_", cov, "X"))
	gene_DP_thresholds <- left_join(gene_DP_thresholds, aux_df, by = "GENE")
	rm(aux_df)
}
rm(cov)
gene_coverage <- left_join(left_join(by = c("GENE" = "GENE_SYMBOL"),
	gene_DP_target, left_join(gene_DP_mane_coding, left_join(gene_DP_mane_exon, left_join(
		gene_DP_mane, gene_DP_thresholds, by = c("GENE_SYMBOL" = "GENE")), by = "GENE_SYMBOL"), by = "GENE_SYMBOL")), panel_genes, by = "GENE") %>% 
	mutate(WHITELIST = GENE %in% genes_whitelist, LowCoverage = DP_TARGET < low_coverage) %>%
	select(GENE, contains("DP_"), MIN_COVERAGE, MAX_COVERAGE, PANEL_SIZE, contains("PCT_"), everything())

# Coverage by exon (from MANE transcript)
exon_metrics_ranges <- makeGRangesFromDataFrame(mergeByOverlaps(mane_perbase_ranges, mane_exon_ranges) %>% as.data.table() %>% 
	mutate(
		seqnames = mane_perbase_ranges.seqnames, start = mane_perbase_ranges.start, end = mane_perbase_ranges.end,
		strand = mane_perbase_ranges.strand, GENE = GENE_SYMBOL) %>%
	select(seqnames, start, end, strand, GENE, EXON, TRANSCRIPT_REFSEQ, TRANSCRIPT_ENSEMBL), keep.extra.columns = TRUE)
exon_metrics <- mergeByOverlaps(exon_metrics_ranges, perBase_ranges) %>% as.data.table() %>%
	mutate(
		CHROM = exon_metrics_ranges.seqnames, START = exon_metrics_ranges.start, END = exon_metrics_ranges.end, LENGTH = exon_metrics_ranges.width) %>% 
	group_by(GENE, EXON) %>%
	summarise(
		MIN_COVERAGE = min(DEPTH), MAX_COVERAGE = max(DEPTH), DEPTH = round(mean(DEPTH)), CHROM = unique(CHROM), START = min(START), END = max(END), 
		LENGTH = sum(LENGTH), TRANSCRIPT_REFSEQ = unique(TRANSCRIPT_REFSEQ), TRANSCRIPT_ENSEMBL = unique(TRANSCRIPT_ENSEMBL)) %>% ungroup()
for(cov in pct_X_coverage){
	aux_df <- mergeByOverlaps(exon_metrics_ranges, perBase_ranges) %>% as.data.table() %>% mutate(LENGTH = exon_metrics_ranges.width) %>%
		group_by(GENE, EXON) %>% summarise(PCT_X = sum(LENGTH[DEPTH >= cov]), LENGTH = sum(LENGTH)) %>% 
		mutate(PCT_X = round(PCT_X*100/LENGTH)) %>% ungroup() %>% select(GENE, EXON, PCT_X)
	colnames(aux_df) <- c("GENE", "EXON", paste0("PCT_", cov, "X"))
	exon_metrics <- left_join(exon_metrics, aux_df, by = c("GENE", "EXON"))
	rm(aux_df)
}
rm(cov)
exon_metrics <- exon_metrics %>% 
	mutate(WHITELIST = GENE %in% genes_whitelist, LowCoverage = DEPTH < low_coverage) %>%
	select(GENE, EXON, DEPTH, contains("PCT_"), everything())



# Save data ---------------------------------------------------------------

# BAM QC metrics
write.table(bamqc_metrics, file = paste0(sample_id, "_BamQC.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
addWorksheet(results_excel, "BamQC")
writeDataTable(results_excel, "BamQC", bamqc_metrics)

# Coverage by target region
write.table(target_coverage, file = paste0(sample_id, "_CoverageByTarget.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
addWorksheet(results_excel, "CoverageByTarget")
writeDataTable(results_excel, "CoverageByTarget", target_coverage)

# Coverage by exon
write.table(exon_metrics, file = paste0(sample_id, "_CoverageByExon.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
addWorksheet(results_excel, "CoverageByExon")
writeDataTable(results_excel, "CoverageByExon", exon_metrics)

# Coverage by gene
write.table(gene_coverage, file = paste0(sample_id, "_CoverageByGene.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
addWorksheet(results_excel, "CoverageByGene")
writeDataTable(results_excel, "CoverageByGene", gene_coverage)

# Coverage by chromosome (TARGET regions)
write.table(chrom_coverage, file = paste0(sample_id, "_CoverageByChr.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
addWorksheet(results_excel, "CoverageByChr")
writeDataTable(results_excel, "CoverageByChr", chrom_coverage)

saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)


# Data visualization ---------------------------------------------------------------

# Set ranges
chr_ranges <- makeGRangesFromDataFrame(chrom_coverage, keep.extra.columns = TRUE)
chr_ranges$y1 <- chr_ranges$y0 <- chr_ranges$DEPTH
gene_ranges <- makeGRangesFromDataFrame(
	(gene_DP_target %>% mutate(seqnames = target_genes$CHROM, start = target_genes$START, end = target_genes$END)),
	keep.extra.columns = TRUE, ignore.strand = TRUE)
gene_ranges$y1 <- gene_ranges$y <- gene_ranges$DP_TARGET
gene_ranges$Whitelist <- gene_ranges$GENE %in% genes_whitelist
gene_ranges$DP_LOW <- gene_ranges$DP_TARGET < low_coverage
low_genes <- gene_ranges[gene_ranges$DP_LOW]
low_genes$y1 <- low_genes$y <- low_genes$DP_TARGET


# Genome visualization

png(paste0(sample_id, "_GlobalCoverage.png"), width = 3000, height = 1200, res = 200)

plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$data2height <- 0
plot.params$topmargin <- 5
plot.params$bottommargin <- 20
plot.params$leftmargin <- 0.06

kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 1.4, chromosomes = chr_levels, plot.params = plot.params)
kpAddChromosomeNames(kp, srt = 70, cex = 1.3, yoffset = -235, xoffset = -0.003)
kpAddChromosomeSeparators(kp, lwd = 0.3, col = "#666666")

ymax <- plyr::round_any(max(gene_ranges$y), 100, f = ceiling)
ymax_sep <- plyr::round_any(ymax/10, 100, f = ceiling)
kpAxis(kp, ymin = 0, ymax = ymax, cex = 0.9, data.panel = 1, r0 = 0.1, r1 = 0.85, tick.pos = seq(0, ymax, ymax_sep), tick.len = 10e6)
kpAddLabels(kp, labels = "Coverage", cex = 1.2, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.1, r1 = 0.85)
kpAbline(kp, h = 0, data.panel = 1, r0 = 0.1, r1 = 0.85, col = "black", ymin = 0, ymax = ymax)
kpAbline(kp, h = low_coverage, data.panel = 1, r0 = 0.1, r1 = 0.85, col = "red", ymin = 0, ymax = ymax)
if(plot_mark_chr_coverage){kpSegments(kp, chr_ranges, ymin = 0, ymax = ymax, lwd = 2, data.panel = 1, r0 = 0.1, r1 = 0.85, col = "black")}
if(plot_global_label_all_genes){
	aux_colors <- unlist(lapply(as.character(gene_ranges$DP_LOW), function(x) lowGenes_colors[x]))
	kpPoints(kp, gene_ranges, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.1, r1 = 0.85, cex = 0.8, col = aux_colors, pch = 16)
	kpPlotMarkers(
		kp, gene_ranges, labels = gene_ranges$GENE, y = 1, cex = 0.6, adjust.label.position = TRUE,
		data.panel = 1, r0 = -0.05, r1 = -0.015, label.margin = 0.01, label.dist = 0.001,
		line.color = aux_colors, label.color = aux_colors, ignore.chromosome.ends = TRUE)
} else {
	kpPoints(
		kp, gene_ranges[!gene_ranges$Whitelist], ymin = 0, ymax = ymax, 
		data.panel = 1, r0 = 0.1, r1 = 0.85, cex = 0.8, col = whitelist_colors_1["FALSE"], pch = 16)
	if(length(genes_whitelist)){
		kpPoints(
			kp, gene_ranges[gene_ranges$Whitelist], ymin = 0, ymax = ymax, 
			data.panel = 1, r0 = 0.1, r1 = 0.85, cex = 0.8, col = whitelist_colors_1["TRUE"], pch = 16)
		kpPlotMarkers(
			kp, gene_ranges[gene_ranges$Whitelist], labels = gene_ranges[gene_ranges$Whitelist]$GENE, 
			y = 1, cex = 0.8, lwd = 0.5, adjust.label.position = TRUE, label.margin = 0.01, label.dist = 0.001,
			data.panel = 1, r0 = -0.05, r1 = 0.9, ignore.chromosome.ends = TRUE, marker.parts = c(0.95, 0.045, 0.005),
			line.color = whitelist_colors_1["TRUE"], label.color = "black")
	}
	if(length(low_genes)){
		kpPoints(kp, low_genes, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.1, r1 = 0.85, cex = 0.8, col = "red", pch = 16)
		kpPlotMarkers(
			kp, low_genes, labels = low_genes$GENE, y = 1, cex = 0.6, adjust.label.position = TRUE,
			data.panel = 1, r0 = -0.05, r1 = -0.015, label.margin = 0.01, label.dist = 0.001,
			line.color = "red", label.color = "red", ignore.chromosome.ends = TRUE)
	}
}

dev.off()

addWorksheet(results_excel, "Plot_GlobalCoverage")
insertImage(results_excel, "Plot_GlobalCoverage", paste0(sample_id, "_GlobalCoverage.png"), units = "px", width = 3000, height = 1200, dpi = 200)
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

rm(kp, ymax, ymax_sep)


# By chromosome visualization

addWorksheet(results_excel, "Plot_CoverageByChr")

plot.params <- getDefaultPlotParams(plot.type = 3)
plot.params$data2height <- 0
plot.params$leftmargin <- 0.06

excel_startRow <- 1

for(chr in paste0("chr", c(1:22, "X"))){

	aux_ranges <- gene_ranges[seqnames(gene_ranges) == chr]
	aux_low_genes <- low_genes[seqnames(low_genes) == chr]

	if(!length(aux_ranges)){next}

	png(paste0(sample_id, "_CoverageByChr_", chr, ".png"), width = 3000, height = 1200, res = 200)
		
	kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 1.5, chromosomes = chr, plot.params = plot.params)
	kpAddBaseNumbers(kp, cex = 1, minor.ticks = TRUE, tick.len = 5, minor.tick.len = 2, add.units = TRUE)
	kpAddMainTitle(kp, main = chr, cex = 2)

	ymax <- plyr::round_any(max(aux_ranges$y), 100, f = ceiling)
	ymax_sep <- plyr::round_any(ymax/10, 100, f = ceiling)
	kpAxis(kp, ymin = 0, ymax = ymax, cex = 0.9, data.panel = 1, r0 = 0.2, r1 = 1, tick.pos = seq(0, ymax, ymax_sep), tick.len = 50e4)
	kpAddLabels(kp, labels = "Coverage", cex = 1.2, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.2, r1 = 1)
	kpAbline(kp, h = 0, data.panel = 1, r0 = 0.2, r1 = 1, col = "black", ymin = 0, ymax = ymax)
	kpAbline(kp, h = low_coverage, data.panel = 1, r0 = 0.2, r1 = 1, col = "red", ymin = 0, ymax = ymax)
	if(plot_mark_chr_coverage){
      kpAbline(kp, h = chr_ranges[seqnames(chr_ranges) == chr]$y0, data.panel = 1, r0 = 0.2, r1 = 1, col = "grey", ymin = 0, ymax = ymax)
    }
	kpBars(
		kp, aux_ranges[!aux_ranges$Whitelist], ymin = 0, ymax = ymax, 
		data.panel = 1, r0 = 0.2, r1 = 1, cex = 0.8, col = whitelist_colors_2["FALSE"], border = whitelist_colors_2["FALSE"])
	kpPoints(
		kp, aux_ranges[!aux_ranges$Whitelist], ymin = 0, ymax = ymax, 
		data.panel = 1, r0 = 0.2, r1 = 1, cex = 1, col = whitelist_colors_2["FALSE"], pch = 16)

	if(length(aux_low_genes)){
		kpBars(
			kp, aux_low_genes[!aux_low_genes$Whitelist], ymin = 0, ymax = ymax, 
			data.panel = 1, r0 = 0.2, r1 = 1, cex = 0.8, col = "red", border = "red")
		kpPoints(
			kp, aux_low_genes[!aux_low_genes$Whitelist], ymin = 0, ymax = ymax, 
			data.panel = 1, r0 = 0.2, r1 = 1, cex = 1, col = "red", pch = 16) 
	}

	if(sum(aux_ranges$Whitelist)){
		kpBars(
			kp, aux_ranges[aux_ranges$Whitelist], ymin = 0, ymax = ymax, 
			data.panel = 1, r0 = 0.2, r1 = 1, cex = 0.8, col = whitelist_colors_2["TRUE"], border = whitelist_colors_2["TRUE"])
		kpPoints(
			kp, aux_ranges[aux_ranges$Whitelist], ymin = 0, ymax = ymax, 
			data.panel = 1, r0 = 0.2, r1 = 1, cex = 1, col = whitelist_colors_2["TRUE"], pch = 16) 
	}
		
	aux_col <- unlist(lapply(seq_len(length(aux_ranges)), function(x) ifelse(
		aux_ranges[x]$Whitelist, whitelist_colors_2["TRUE"], ifelse(
		aux_ranges[x]$GENE %in% aux_low_genes$GENE, "red", whitelist_colors_2["FALSE"]))))
	kpPlotMarkers(
		kp, aux_ranges, labels = aux_ranges$GENE, y = 1, cex = 0.9, adjust.label.position = TRUE,
		data.panel = 1, r0 = -0.05, r1 = 0, label.margin = 0.01, label.dist = 0.001,
		line.color = aux_col, label.color = aux_col, ignore.chromosome.ends = TRUE)

	dev.off()

	insertImage(
		results_excel, "Plot_CoverageByChr", paste0(sample_id, "_CoverageByChr_", chr, ".png"), 
		units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow)
	saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

	excel_startRow <- excel_startRow + 30

	rm(kp, aux_ranges, aux_low_genes, ymax, aux_col)

}
rm(chr, excel_startRow)


# By gene (from Whitelist)

if(plot_bygene_all){genes_whitelist <- panel_genes_names}
if(length(genes_whitelist)) {

	addWorksheet(results_excel, "Plot_CoverageByGene")

	plot.params <- getDefaultPlotParams(plot.type = 3)
	plot.params$data2height <- 0
	plot.params$leftmargin <- 0.06

	excel_startRow <- 1

	for(gene in genes_whitelist){
		
		aux_ranges <- mane_gene_ranges[mane_gene_ranges$GENE_SYMBOL == gene]

		if(!length(aux_ranges)){next}
		
		aux_exon_ranges <- mane_exon_ranges[mane_exon_ranges$GENE_SYMBOL == gene]
		aux_gene <- aux_ranges$GENE_ENSEMBL
		transcript <- aux_ranges$TRANSCRIPT_ENSEMBL
		width <- width(aux_ranges)
		zoom <- aux_ranges + width*0.05
		aux_perBase_ranges <- subsetByOverlaps(perBase_ranges, zoom)

		if(!length(aux_perBase_ranges)){next}

		png(paste0(sample_id, "_CoverageByGene_", gene, ".png"), width = 3000, height = 1200, res = 200)

		kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 2, plot.params = plot.params, zoom = zoom)
		kpAddBaseNumbers(
			kp, cex = 1, minor.ticks = TRUE, tick.len = 5, minor.tick.len = 2, add.units = TRUE, digits = 3,
			tick.dist=round(width(zoom)/5, -2), minor.tick.dist=round(width(zoom)/50, -2))
		kpAddChromosomeNames(kp, cex = 1.5, yoffset = -225, xoffset = -0.47)
		kpAddMainTitle(kp, main = paste0(gene, " (", transcript, ")"), cex = 1.8)

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
			kp, genes.data, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.85, r1 = 0.95, add.strand.marks=TRUE, mark.height = 0.25,
			coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50", non.coding.exons.border.col="grey50"))
		
		ymax <- plyr::round_any(max(aux_perBase_ranges$y1), 100, f = ceiling)
		ymax_sep <- plyr::round_any(ymax/10, 100, f = ceiling)
		kpAxis(kp, ymin = 0, ymax = ymax, cex = 0.9, data.panel = 1, r0 = 0.05, r1 = 0.8, tick.pos = seq(0, ymax, ymax_sep), tick.len = width/200)
		kpAddLabels(kp, labels = "Coverage", cex = 1.2, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.05, r1 = 0.8)
		kpAbline(kp, h = 0, data.panel = 1, r0 = 0.05, r1 = 0.8, col = "black", ymin = 0, ymax = ymax)
		kpBars(kp, aux_perBase_ranges, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.05, r1 = 0.8, col = "black", border = "black")

		kpAddChromosomeNames(kp, cex = 1, yoffset = -210, xoffset = -0.47, chr.names = "Targets")
		kpPlotRegions(kp, target_ranges, col = "black", data.panel = 1, r0 = 0, r1 = 0.05)

		# Close PNG object
		dev.off()

		insertImage(
			results_excel, "Plot_CoverageByGene", paste0(sample_id, "_CoverageByGene_", gene, ".png"), 
			units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow)
		saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)
		
		excel_startRow <- excel_startRow + 30
		
		rm(aux_ranges, aux_exon_ranges, aux_perBase_ranges, aux_gene, transcript, width, zoom, kp, genes.data, is_transcript, ymax, ymax_sep)
	}
	rm(gene, excel_startRow)

}

