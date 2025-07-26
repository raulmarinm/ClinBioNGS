#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of splicing results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(GenomicFeatures)))
suppressPackageStartupMessages(suppressWarnings(library(karyoploteR)))
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
splicing_whitelist.file <- args_list[["whitelistSplicing"]]
mane.gtf <- args_list[["maneGtf"]]
mane_gene.file <- args_list[["maneGene"]]
mane_exon.file <- args_list[["maneExon"]]
mane_intron.file <- args_list[["maneIntron"]]
cancerdrivers.file <- args_list[["cancerdrivers"]]
civic.file <- args_list[["civicClinical"]]
chain.file <- args_list[["chain"]]
fasta.file <- args_list[["fasta"]]
fasta_name <- args_list[["fastaName"]]
headerVcf.file <- args_list[["vcfHeader"]]
calling_min_ad <- as.numeric(args_list[["callingMinAD"]])
flag_min_ad <- as.numeric(args_list[["flagMinAD"]])
flag_min_dp <- as.numeric(args_list[["flagMinDP"]])
flag_min_vaf <- as.numeric(args_list[["flagMinVAF"]])

# Set other variables
sample_id <- paste0(sample, "_", type)
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
calling_status <- c("PASS", "LowSupport")
splicing_flags <- c("OK", "NoCall", "NoCancerEnriched", "LowAD", "LowDP", "LowVAF")
results_excel <- createWorkbook()
results_excel.file <- paste0(run, "_", sample_id, "_Splicing_", format(Sys.time(), "%Y%m%d"), ".xlsx")

# Set splicing results
splicing_all.file <- list.files(pattern = paste0(sample_id, ".introns$"))
splicing_cancer.file <- list.files(pattern = paste0(sample_id, ".cancer.introns.prelim$"))

# Set other files
perBaseCoverage.file <- list.files(pattern = "per-base.bed.gz$")
mutation.file <- list.files(pattern = "SmallVariant_annot.txt$")


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

# Known splicing variants from whitelist
splicing_whitelist <- read.table(splicing_whitelist.file, header = TRUE, sep = ";", quote = "") %>% select(VAR, VAR_NAME)

# MANE data
mane_TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(mane.gtf, format = "gtf")))
mane_gene <- fread(mane_gene.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% target_genes_names)
mane_gene_names <- mane_gene$GENE_ENSEMBL; names(mane_gene_names) <- mane_gene$GENE_SYMBOL
mane_exon <- fread(mane_exon.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% target_genes_names)
mane_exon_ranges <- makeGRangesFromDataFrame(mane_exon, keep.extra.columns = TRUE)
mane_intron <- fread(mane_intron.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% target_genes_names)
mane_intron_ranges <- makeGRangesFromDataFrame(mane_intron, keep.extra.columns = TRUE)

# Per-base coverage
if(length(perBaseCoverage.file)){
	perBase_ranges <- makeGRangesFromDataFrame(fread(
		perBaseCoverage.file, header = FALSE, col.names = c("CHROM", "START", "END", "DEPTH"),
		sep = "\t", quote = "", fill = FALSE, nThread = 1) %>% filter(CHROM %in% chr_levels), 
		starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
	perBase_ranges$y1 <- perBase_ranges$DEPTH
}

# Network of Cancer Genes (NCG) data
cancerdrivers <- read.table(cancerdrivers.file, header = TRUE, sep = "\t", quote = "", fill = FALSE)

# CIViC clinical data
civic_data <- read.table(civic.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)
colnames(civic_data) <- paste0("CIViC_", colnames(civic_data))

# CTAT-Splicing: all splicing variants
splicing_all <- read.table(splicing_all.file, header = TRUE, sep = "\t") %>%
	mutate(VAR = as.character(intron)) %>% 
	separate(genes, c("GENE_A", "GENE_B"), ",", fill = "right", extra = "drop") %>%
	separate(GENE_A, c("GENE_A", "GENE_A_ID"), "\\^", extra = "drop") %>%
	separate(GENE_B, c("GENE_B", "GENE_B_ID"), "\\^", extra = "drop") %>% select(-intron) %>%
	mutate(
		GENE = as.character(ifelse(is.na(GENE_B) | GENE_A %in% target_genes_names, GENE_A, GENE_B)),
		GENE_ENSEMBL = ifelse(is.na(GENE_B) | GENE_A %in% target_genes_names, GENE_A_ID, GENE_B_ID))

# CTAT-Splicing: cancer splicing variants
splicing_cancer <- read.table(splicing_cancer.file, header = TRUE, sep = "\t") %>%
	mutate(TCGA = TCGA_sample_counts, GTEx = GTEx_sample_counts, VAR_NAME = as.character(variant_name), VAR = as.character(intron)) %>%
	select(VAR, VAR_NAME, TCGA, GTEx)
splicing_cancer_ids <- splicing_cancer$VAR

# Mutations
if(length(mutation.file)){
	mutation_info <- read.table(mutation.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
		filter(GENE_SYMBOL %in% target_genes_names & CONSEQUENCE %in% c("splice_donor_variant", "splice_acceptor_variant")) %>%
		mutate(VAR_MUT = VAR) %>% select(CHROM, START, END, VAR_MUT)
} else {
	mutation_info <- data.frame()
}


# Data processing ---------------------------------------------------------

# Remove off-target genes
# Annotate known relevant splice variants from the whitelist
# AD: uniq_mapped + multi_mapped
# DP: AD of junctions with the same START or END
# AF: AD/DP
# LowSupport: low AD
# CancerEnriched: TRUE or FALSE
# CALLING: LowSupport or PASS
# FLAGS: OK or NoCancerEnriched/LowAD/LowDP/LowVAF
# Add TCGA and GTEX counts
splicing_target <- splicing_all %>% filter(GENE %in% target_genes_names) %>%
	left_join(splicing_cancer, by = "VAR") %>%
	left_join(splicing_whitelist, by = "VAR") %>%
	separate(VAR, c("CHROM", "START", "END"), remove = FALSE) %>% rowwise() %>%
	mutate(
		VAR_NAME = ifelse(is.na(VAR_NAME.x), VAR_NAME.y, VAR_NAME.x), 
		WHITELIST_SPLICING = !is.na(VAR_NAME),
		AD = uniq_mapped + multi_mapped, ID_START = ifelse(strand == "-", paste0(CHROM, ":", END), paste0(CHROM, ":", START)), 
		ID_END = ifelse(strand == "-", paste0(CHROM, ":", START), paste0(CHROM, ":", END)),
		START = as.numeric(START), END = as.numeric(END), UNIQ_MAPPED_READS = uniq_mapped,
		MULTI_MAPPED_READS = multi_mapped, WHITELIST_GENE = GENE %in% genes_whitelist, STRAND = strand, MUTATION = NA) %>% 
	group_by(ID_START) %>% mutate(DP_START = sum(AD)) %>% group_by(ID_END) %>% mutate(DP_END = sum(AD)) %>% ungroup() %>% rowwise() %>%
	mutate(
		DP_MAX = max(c(DP_START, DP_END)), DP_MEAN = mean(c(DP_START, DP_END)), AF = round(AD/DP_MAX, 4),
		LowSupport = AD < calling_min_ad, PASS = !LowSupport, CALLING = paste(calling_status[c(PASS, LowSupport)], collapse = ";"),
		CancerEnriched = VAR %in% splicing_cancer_ids, LowAD = AD < flag_min_ad, LowDP = DP_MAX < flag_min_dp, LowVAF = AF < flag_min_vaf, 
		NO_FLAGS = PASS & (CancerEnriched | WHITELIST_SPLICING) & !(LowAD | LowDP | LowVAF),
		FLAGS = paste(splicing_flags[c(NO_FLAGS, !PASS, (!CancerEnriched & !WHITELIST_SPLICING), LowAD, LowDP, LowVAF)], collapse = ";")) %>%
	as.data.frame()
splicing_ranges <- GRanges(splicing_target$VAR, VAR = splicing_target$VAR, GENE = splicing_target$GENE)


if(nrow(splicing_target)){
	# Annotate the affected EXONS/INTRONS
	splicing_exons <- as.data.frame(mergeByOverlaps(splicing_ranges, mane_exon_ranges)) %>% filter(GENE == GENE_SYMBOL) %>% 
		group_by(VAR) %>% 
		summarise(
			EXONS_AFFECTED = paste(unique(c(min(EXON), max(EXON))), collapse = "-"),
			transcript_exon = unique(TRANSCRIPT_ENSEMBL), COV_EXONS_AFFECTED = NA, COV_EXONS_FLANKING = NA, .groups = "drop")
	splicing_introns <- as.data.frame(mergeByOverlaps(splicing_ranges, mane_intron_ranges)) %>% filter(GENE == GENE_SYMBOL) %>% 
		group_by(VAR) %>% 
		summarise(
			INTRONS_AFFECTED = paste(unique(c(min(INTRON), max(INTRON))), collapse = "-"),
			transcript_intron = unique(TRANSCRIPT_ENSEMBL), .groups = "drop")
	splicing_target <- left_join(left_join(splicing_target, splicing_exons, by = "VAR"), splicing_introns, by = "VAR") %>%
		mutate(
			REGION_AFFECTED = ifelse(!is.na(EXONS_AFFECTED), paste0("Exon", EXONS_AFFECTED), ifelse(
				!is.na(INTRONS_AFFECTED), paste0("Intron", INTRONS_AFFECTED), ".")))
	
	# Calculate coverage of AFFECTED and FLANKING exons
	if(length(perBaseCoverage.file)){
		mane_exon_perbase_ranges <- makeGRangesFromDataFrame(
			mane_exon_ranges %>% as.data.table() %>% filter(GENE_SYMBOL %in% unique(splicing_target$GENE)) %>% 
				mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% mutate(pos = 1:n(), end = start+pos-1, start=end) %>% 
				as.data.table() %>% select(-ID, -pos), keep.extra.columns = TRUE)

		exon_coverage <- mergeByOverlaps(mane_exon_perbase_ranges, perBase_ranges) %>% as.data.table() %>% select(GENE_SYMBOL, EXON, DEPTH)

		for(n in seq_len(nrow(splicing_target))){
			aux_gene <- splicing_target[n, "GENE"]
			aux_exons <- eval(parse(text = sub("-", ":", splicing_target[n, "EXONS_AFFECTED"])))
			aux_introns <- eval(parse(text = sub("-", ":", splicing_target[n, "INTRONS_AFFECTED"])))
			if(all(is.na(aux_exons))){aux_exons_flanking <- c(aux_introns, max(aux_introns)+1)
			} else {aux_exons_flanking <- c(min(aux_exons)-1, max(aux_exons)+1)}

			aux_coverage <- exon_coverage %>% filter(GENE_SYMBOL == aux_gene)
			aux_depth <- aux_coverage %>% filter(EXON %in% aux_exons) %>% summarise(COV_EXONS_AFFECTED = round(mean(DEPTH)))
			aux_depth_flanking <- aux_coverage %>% filter(EXON %in% aux_exons_flanking) %>% summarise(COV_EXONS_FLANKING = round(mean(DEPTH)))
			splicing_target[n, "COV_EXONS_AFFECTED"] <- aux_depth$COV_EXONS_AFFECTED
			splicing_target[n, "COV_EXONS_FLANKING"] <- aux_depth_flanking$COV_EXONS_FLANKING
			
			rm(aux_gene, aux_exons, aux_exons_flanking, aux_coverage, aux_depth, aux_depth_flanking)
		}
		rm(n, mane_exon_perbase_ranges, exon_coverage)
	}

	# Annotate if there is an associated mutation (ONLY splice donor/acceptor)
	if(nrow(mutation_info)){
		mutation_ranges <- makeGRangesFromDataFrame(mutation_info, ignore.strand = TRUE, keep.extra.columns = TRUE)
		aux_mutation <- as.data.frame(mergeByOverlaps(splicing_ranges, mutation_ranges)) %>% select(VAR, VAR_MUT)
		splicing_target <- left_join(splicing_target, aux_mutation, by = "VAR") %>% mutate(MUTATION = VAR_MUT)
	}

	# Convert Hg38 to Hg19 coordinates
	splicing_start_ranges <- GRanges(paste(splicing_target$CHROM, splicing_target$START, sep = ":"), VAR = splicing_target$VAR)
	splicing_end_ranges <- GRanges(paste(splicing_target$CHROM, splicing_target$END, sep = ":"), VAR = splicing_target$VAR)
	splicing_hg19_start <- as.data.frame(unlist(liftOver(splicing_start_ranges, import.chain(chain.file)))) %>% 
		distinct(VAR, .keep_all = TRUE) %>% mutate(START_HG19 = start)
	splicing_hg19_end <- as.data.frame(unlist(liftOver(splicing_end_ranges, import.chain(chain.file)))) %>% 
		distinct(VAR, .keep_all = TRUE) %>% mutate(END_HG19 = start)
	splicing_target <- left_join(splicing_target, left_join(splicing_hg19_start, splicing_hg19_end, by = "VAR"), by = "VAR") %>%
		mutate(VAR_HG19 = paste0(CHROM, ":", START_HG19, "-", END_HG19))

} else {
	splicing_target <- splicing_target %>% mutate(
		EXONS_AFFECTED = NA, COV_EXONS_AFFECTED = NA, COV_EXONS_FLANKING = NA, INTRONS_AFFECTED = NA, REGION_AFFECTED = NA,
		START_HG19 = NA, END_HG19 = NA, VAR_HG19 = NA)
}

# Create final output
splicing_target$CHROM <- factor(splicing_target$CHROM, levels = chr_levels)
splicing_results <- splicing_target %>% 
	left_join(target_genes, by = "GENE") %>%
	left_join(cancerdrivers, by = "GENE") %>%
	mutate(
		RUN = run, SAMPLE = sample, VAR_GENE = paste0(GENE, "_", REGION_AFFECTED),
		across(c(ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER), ~ .x & !is.na(.x))) %>%
	arrange(desc(WHITELIST_SPLICING), VAR_NAME, desc(NO_FLAGS), desc(PASS), desc(CancerEnriched), desc(WHITELIST_GENE), CHROM, START, END) %>%
	select(
		RUN, SAMPLE, VAR, VAR_NAME, VAR_GENE, AD, AF, FLAGS, MUTATION, DP_MAX, 
		COV_EXONS_AFFECTED, COV_EXONS_FLANKING, WHITELIST_SPLICING, WHITELIST_GENE, 
		DP_START, DP_END, DP_MEAN, UNIQ_MAPPED_READS, MULTI_MAPPED_READS, 
		GENE, REGION_AFFECTED, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, GENE_ENSEMBL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, 
		EXONS_AFFECTED, INTRONS_AFFECTED, CHROM, START, END, STRAND, VAR_HG19, START_HG19, END_HG19, 
		CALLING, PASS, LowSupport, NO_FLAGS, CancerEnriched, LowAD, LowDP, LowVAF, TCGA, GTEx) %>%
	mutate(across(where(is.factor), as.character))

# Clinical classification (AMP/ASCO/CAP guidelines)
#	- Predictive, prognostic & diagnostic evidences in CIViC resource
# 	- Tier IA: CIViC evidence A (guidelines) in the SAME tumor
# 	- Tier IB: CIViC evidence B (trials) in the SAME tumor
# 	- Tier IIC: CIViC evidence A/B in OTHER tumor OR CIViC evidence C
# 	- Tier IID: CIViC evidence D
# 	- Tier III: Tumor-enriched (TCGA vs GTEX)
# 	- Tier IV: Non tumor-enriched
splicing_civic <- splicing_results %>% filter(!is.na(VAR_NAME) & VAR_NAME %in% civic_data$CIViC_ALTERATION)
civic_annot <- civic_data %>% mutate(VAR = ".", CLIN_STATUS = ".") %>% filter(CIViC_ALTERATION == ".")
for(n in seq_len(nrow(splicing_civic))){
	aux_var <- splicing_civic[n, ]
	aux_civic <- civic_data %>% filter(CIViC_ALTERATION %in% aux_var$VAR_NAME) %>%
		mutate(
			VAR = aux_var$VAR,
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
splicing_results <- splicing_results %>% left_join(civic_annot, by = "VAR") %>%
	mutate(
		CLIN_STATUS = ifelse(!is.na(CLIN_STATUS), CLIN_STATUS, ifelse(WHITELIST_SPLICING | CancerEnriched, "Tier III", "Tier IV")),
		CLIN_RELEVANT = CLIN_STATUS %in% c("Tier IA", "Tier IB", "Tier IIC", "Tier IID")) %>%
	select(RUN:FLAGS, CLIN_RELEVANT, CLIN_STATUS, everything()) %>%
	arrange(CLIN_STATUS)
splicing_results[is.na(splicing_results) | splicing_results == "NaN"] <- "."

# Generate a VCF file
aux_header <- readLines(headerVcf.file)
aux_header[grep("##fileDate", aux_header)] <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
aux_header[grep("##reference", aux_header)] <- paste0("##reference=file:", fasta_name)
if(nrow(splicing_results)){
	aux_vars <- splicing_results %>% 
		mutate(SVTYPE = "DEL", ALT = "<DEL>", DP = DP_MAX, CHROM = factor(CHROM, levels = chr_levels)) %>%
		arrange(CHROM, START, END)
	aux_ranges <- makeGRangesFromDataFrame(aux_vars, start.field = "START", end.field = "START", seqnames.field = "CHROM")
	aux_ref <- as.data.frame(getSeq(fasta, aux_ranges))$x
	aux_info <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){paste(unlist(lapply(c(
		"VAR", "RUN", "SAMPLE", "SVTYPE", "VAR_NAME", "VAR_GENE", "AD", "AF", "CLIN_STATUS", "MUTATION", "DP", 
		"COV_EXONS_AFFECTED", "COV_EXONS_FLANKING", "DP_START", "DP_END", "DP_MEAN", "UNIQ_MAPPED_READS", "MULTI_MAPPED_READS", 
		"GENE", "REGION_AFFECTED", "GENE_ENSEMBL", "TRANSCRIPT_ENSEMBL", "TRANSCRIPT_REFSEQ", 
		"EXONS_AFFECTED", "INTRONS_AFFECTED", "CHROM", "START", "END", "STRAND", "VAR_HG19", "START_HG19", "END_HG19", 
		"TCGA", "GTEx", "CIViC_ALTERATION", "CIViC_VARIANT_ID", "CIViC_EVIDENCE_LEVEL", "CIViC_EVIDENCE_RATING", "CIViC_EVIDENCE_TYPE", 
		"CIViC_EFFECT", "CIViC_DRUG", "CIViC_TUMOR_NAME", "CIViC_TUMOR_DOID", "CIViC_TUMOR_CODE", "CIViC_TOPNODE_NAME",
		"CIViC_TOPNODE_DOID", "CIViC_MOLECULAR_ID", "CIViC_VARIANT_GROUP", "CIViC_EVIDENCE_SCORE", "CIViC_ORIGIN"),
		function(x) paste0(x, "=", aux_vars[n, x]) )), collapse = ";")  }))
	aux_flags <- paste0(aux_vars$FLAGS, ";", unlist(lapply(seq_len(nrow(aux_vars)), function(n){paste(
		unlist(lapply(c(
			"WHITELIST_GENE", "ONCOGENE", "TSG", "DRIVER", "CANONICAL_DRIVER", "CancerEnriched", "WHITELIST_SPLICING", "CLIN_RELEVANT"), 
			function(x) ifelse(aux_vars[n, x], paste0(x, ";"), "")  )), 
		collapse = "")  })))
	aux_info <- paste0(aux_info, ";", aux_flags)
	aux_format <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("AD", "DP", "AF"), function(x) paste0(aux_vars[n, x]) )), collapse = ":")  }))
	aux_df <- aux_vars %>% mutate(
		CHROM = CHROM, POS = START, ID = ".", REF = aux_ref, QUAL = ".", FILTER = CALLING, INFO = aux_info,
		FORMAT = paste(c("AD", "DP", "AF"), collapse = ":"), SAMPLE = aux_format) %>%
		select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)
	vcf_text <- c(aux_header,
		paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"),
		unlist(lapply(seq_len(nrow(aux_df)), function(x) paste(aux_df[x, ], collapse = "\t"))) )
	
	rm(aux_header, aux_vars, aux_ranges, aux_ref, aux_info, aux_flags, aux_format, aux_df)
} else {
	vcf_text <- c(aux_header, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"))
}


# Save data ---------------------------------------------------------------

write.table(
	splicing_results, file = paste0(sample_id, "_Splicing_allCandidates.txt"), 
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(
	splicing_results %>% filter(PASS), file = paste0(sample_id, "_Splicing_results.txt"), 
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

addWorksheet(results_excel, "OK")
writeDataTable(results_excel, "OK", splicing_results %>% filter(NO_FLAGS))
addWorksheet(results_excel, "PASS")
writeDataTable(results_excel, "PASS", splicing_results %>% filter(PASS))
addWorksheet(results_excel, "FAIL")
writeDataTable(results_excel, "FAIL", splicing_results %>% filter(!PASS))
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

writeLines(con = paste0(sample_id, "_Splicing.vcf"), text = vcf_text)


# Data visualization ---------------------------------------------------------------

# Visualization of gene coverage: WHITELIST_SPLICING | (PASS & CancerEnriched) variants
plot_vars <- splicing_results %>% filter(WHITELIST_SPLICING | VAR_NAME != "." | (PASS & CancerEnriched))

if(nrow(plot_vars)) {
	
	addWorksheet(results_excel, "Plots")

	plot.params <- getDefaultPlotParams(plot.type = 3)
	plot.params$data2height <- 0
	plot.params$leftmargin <- 0.06

	excel_startRow <- 1        

	for(n in seq_len(nrow(plot_vars))){

		aux_splicing_vars <- plot_vars[n, ] %>% mutate(across(everything(), ~ ifelse(.x != ".", .x, NA)))
		aux_var <- aux_splicing_vars$VAR
		aux_var_name <- aux_splicing_vars$VAR_NAME
		aux_gene <- aux_splicing_vars$GENE
		aux_gene_id <- mane_gene_names[aux_gene]
		aux_region_affected <- aux_splicing_vars$REGION_AFFECTED
		aux_af <- aux_splicing_vars$AF
		aux_splicing_id <- sub(":", "-", aux_var)

		aux_transcript <- aux_splicing_vars$TRANSCRIPT_ENSEMBL
		aux_splicing_ranges <- GRanges(aux_var)
		aux_splicing_ranges_1 <- toGRanges(aux_splicing_vars$CHROM, aux_splicing_vars$START)
		aux_splicing_ranges_2 <- toGRanges(aux_splicing_vars$CHROM, aux_splicing_vars$END)
		
		aux_introns <- eval(parse(text = sub("-", ":", aux_splicing_vars$INTRONS_AFFECTED)))
		aux_exons <- eval(parse(text = sub("-", ":", aux_splicing_vars$EXONS_AFFECTED)))
		if(all(is.na(aux_exons))){
			aux_exons <- c(aux_introns, max(aux_introns)+1)
		} else {
			aux_exons <- (min(aux_exons)-1):(max(aux_exons)+1)
		}

		other_cancer_vars <- splicing_results %>% filter(
			GENE %in% aux_gene & VAR != aux_var & 
			((WHITELIST_SPLICING | VAR_NAME != ".") | (CancerEnriched & PASS & (START == aux_splicing_vars$START | END == aux_splicing_vars$END))))
		other_cancer_ranges <- toGRanges(other_cancer_vars$CHROM, other_cancer_vars$START, other_cancer_vars$END)
		other_cancer_vars_1 <- toGRanges(other_cancer_vars$CHROM, other_cancer_vars$START)
		other_cancer_vars_2 <- toGRanges(other_cancer_vars$CHROM, other_cancer_vars$END)

		if(!is.na(aux_transcript)){
			aux_normal_vars <- splicing_results %>% filter(
				GENE %in% aux_gene & PASS & (!VAR %in% c(aux_var, other_cancer_vars$VAR)) &
				(INTRONS_AFFECTED %in% aux_introns | START == aux_splicing_vars$START | END == aux_splicing_vars$END))
		} else {
			aux_normal_vars <- splicing_target %>% filter(GENE %in% aux_gene & VAR != aux_var & !WHITELIST_SPLICING & PASS)
		}
		aux_normal_ranges <- toGRanges(aux_normal_vars$CHROM, aux_normal_vars$START, aux_normal_vars$END)
		aux_normal_ranges_1 <- toGRanges(aux_normal_vars$CHROM, aux_normal_vars$START)
		aux_normal_ranges_2 <- toGRanges(aux_normal_vars$CHROM, aux_normal_vars$END)

		aux_mane_ranges <- mane_exon_ranges[mane_exon_ranges$GENE_SYMBOL == aux_gene & mane_exon_ranges$EXON %in% aux_exons]
		aux_region <- c(aux_splicing_ranges, aux_mane_ranges)
		zoom <- toGRanges(unique(seqnames(aux_region)), min(start(aux_region)), max(end(aux_region)))+200
		
		title_label <- paste0(
			ifelse(!is.na(aux_var_name), aux_var_name, paste0(aux_gene, " ", ifelse(!is.na(aux_region_affected), aux_region_affected, ""))),
			" ", ifelse(!is.na(aux_transcript), paste0("(", aux_transcript, ") "), ""), "AF=", aux_af,
			ifelse(aux_splicing_vars$NO_FLAGS, "", paste0(" *", aux_splicing_vars$FLAGS, "*"))
		)
		
		plot_name <- paste0(
			sample_id, "_Splicing_", aux_gene, "_", 
			ifelse(!is.na(aux_region_affected), paste0(aux_region_affected, "_"), ""),
			aux_splicing_id, ".png")

		png(plot_name, width = 3000, height = 1200, res = 200)

		kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 2, plot.params = plot.params, zoom = zoom)
		kpAddBaseNumbers(kp, cex = 1, minor.ticks = TRUE, tick.len = 5, minor.tick.len = 2, add.units = TRUE, digits = 3,
			tick.dist=round(width(zoom)/5, -1), minor.tick.dist=round(width(zoom)/50, -1))
		kpAddChromosomeNames(kp, cex = 1.5, yoffset = -225, xoffset = -0.47)
		kpAddMainTitle(kp, cex = 1.5, main = title_label)

		if(!is.na(aux_transcript)){
			genes.data <- makeGenesDataFromTxDb(mane_TxDb, kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
			is_gene <- names(genes.data$genes) == aux_gene_id
			genes.data$genes <- genes.data$genes[is_gene]
			genes.data$transcripts <- genes.data$transcripts[is_gene]
			is_transcript <- genes.data[["transcripts"]][[aux_gene_id]]$tx_name == aux_transcript
			genes.data[["transcripts"]][[aux_gene_id]] <- genes.data[["transcripts"]][[aux_gene_id]][is_transcript]
			genes.data[["coding.exons"]] <- genes.data[["coding.exons"]][is_transcript]
			genes.data[["non.coding.exons"]] <- genes.data[["non.coding.exons"]][is_transcript]

			suppressWarnings(kpPlotGenes(
				kp, genes.data, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.76, r1 = 0.85, add.strand.marks=TRUE, mark.height = 0.25, 
				coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50", non.coding.exons.border.col="grey50"))
			if(length(aux_mane_ranges)){
				kpText(kp, aux_mane_ranges, y = -0.15, labels = aux_mane_ranges$EXON, data.panel = 1, r0 = 0.75, r1 = 0.85, cex = 0.8, col = "black")
			}
			

			rm(genes.data, is_gene, is_transcript)
		}
		
		kpPlotLinks(kp, aux_normal_ranges_1, aux_normal_ranges_2, data.panel = 1, r0 = 0.8, r1 = 1, arch.height = 0.35, border = "blue")
		kpPlotLinks(kp, other_cancer_vars_1, other_cancer_vars_2, data.panel = 1, r0 = 0.8, r1 = 1, arch.height = 0.6, border = "lightsalmon")
		kpPlotLinks(kp, aux_splicing_ranges_1, aux_splicing_ranges_2, data.panel = 1, r0 = 0.8, r1 = 1, arch.height = 1, col = "red", border = "red")
		if(length(aux_normal_ranges)){suppressWarnings(kpPlotMarkers(
			kp, aux_normal_ranges, labels = aux_normal_vars$AD, y = 0.3, cex = 0.7, text.orientation = "vertical",
			adjust.label.position = TRUE, label.margin = 0.01, label.dist = 0.001, marker.parts = c(0, 0, 0),
			data.panel = 1, r0 = 0.8, r1 = 1, ignore.chromosome.ends = TRUE, label.color = "blue"))}
		
		if(length(other_cancer_ranges)){suppressWarnings(kpPlotMarkers(
			kp, other_cancer_ranges, labels = other_cancer_vars$AD, y = 0.5, cex = 0.7, text.orientation = "vertical",
			adjust.label.position = TRUE, label.margin = 0.01, label.dist = 0.001, marker.parts = c(0, 0, 0),
			data.panel = 1, r0 = 0.8, r1 = 1, ignore.chromosome.ends = TRUE, label.color = "lightsalmon"))}
		suppressWarnings(kpPlotMarkers(
			kp, aux_splicing_ranges, labels = aux_splicing_vars$AD, y = 0.8, cex = 0.7, text.orientation = "horizontal",
			adjust.label.position = TRUE, label.margin = 0.01, label.dist = 0.001, marker.parts = c(0, 0, 0),
			data.panel = 1, r0 = 0.8, r1 = 1, ignore.chromosome.ends = TRUE, label.color = "red"))

		if(length(perBaseCoverage.file)){
			aux_perBase_ranges <- subsetByOverlaps(perBase_ranges, zoom)
			ymax <- plyr::round_any(max(aux_perBase_ranges$y1), 100, f = ceiling)
			ymax_sep <- plyr::round_any(ymax/10, 100, f = ceiling)
			kpAxis(kp, ymin = 0, ymax = ymax, cex = 0.9, data.panel = 1, r0 = 0.03, r1 = 0.7, tick.pos = seq(0, ymax, ymax_sep), tick.len = width(zoom)/200)
			kpAddLabels(kp, labels = "Coverage", cex = 1.2, label.margin = 0.05, srt = 90, pos = 1, data.panel = 1, r0 = 0.03, r1 = 0.7)
			kpAbline(kp, h = 0, data.panel = 1, r0 = 0.03, r1 = 0.7, col = "black", ymin = 0, ymax = ymax)
			kpBars(kp, aux_perBase_ranges, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.03, r1 = 0.7, col = "black", border = "black")

			rm(aux_perBase_ranges, ymax, ymax_sep)
		}
		
		kpAbline(kp, v = aux_splicing_vars$START, data.panel = 1, r0 = 0, r1 = 1, col = "red", lwd = 0.8, lty = 2)
		kpAbline(kp, v = aux_splicing_vars$END, data.panel = 1, r0 = 0, r1 = 1, col = "red", lwd = 0.8, lty = 2)
		kpAddChromosomeNames(kp, cex = 1, yoffset = -210, xoffset = -0.47, chr.names = "Targets")
		kpPlotRegions(kp, target_ranges, col = "black", data.panel = 1, r0 = -0.05, r1 = 0.02)

		dev.off()

		insertImage(results_excel, "Plots", file = plot_name, units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow)
		saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)

		excel_startRow <- excel_startRow + 35

		rm(
			aux_splicing_vars, aux_var, aux_var_name, aux_gene, aux_gene_id, aux_region_affected, aux_splicing_id, aux_transcript, 
			aux_splicing_ranges, aux_splicing_ranges_1, aux_splicing_ranges_2, aux_introns, aux_exons, other_cancer_vars, 
			other_cancer_ranges, other_cancer_vars_1, other_cancer_vars_2, aux_normal_vars, aux_mane_ranges, aux_region, zoom, 
			title_label, plot_name, kp)
	}
	rm(n, excel_startRow)
	
}


