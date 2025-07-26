#!/usr/bin/env Rscript

# ============================================================================ #
# What: Processing of Small variant results
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")

# Load libraries
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(karyoploteR))
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
sample.file <- args_list[["samplesFile"]]
consensus.file <- args_list[["consensus"]]
annotation.file <- args_list[["annotation"]]
genes.file <- args_list[["targetGenes"]]
whitelist.file <- args_list[["whitelistGenes"]]
mane.gtf <- args_list[["maneGtf"]]
mane_gene.file <- args_list[["maneGene"]]
mane_exon.file <- args_list[["maneExon"]]
genie.file <- args_list[["genieOncogenic"]]
civic.file <- args_list[["civicClinical"]]
is_panel_hotspots <- file.exists(args_list[["panelHotspots"]])
headerVcf.file <- args_list[["vcfHeader"]]
fasta_name <- args_list[["fastaName"]]
min_vaf <- as.numeric(args_list[["minVAF"]])
somatic_max_pvaf <- as.numeric(args_list[["somaticMaxpVAF"]])
somatic_max_vaf <- as.numeric(args_list[["somaticMaxVAF"]])

# Set other variables
sample_id <- paste0(sample, "_", type)
var_flags <- c("OK", "NoCall", "Germline", "OutCTR", "ProblematicRegion", "Recurrent")
oncogenicity_levels <- c("Oncogenic", "Likely Oncogenic", "VUS", "Likely Benign", "Benign")
oncogenicity_colors <- c(
	"Oncogenic" = "red", "Likely Oncogenic" = "darkorange", "VUS" = "#FFDB6D", "Likely Benign" = "lightgreen", "Benign" = "mediumseagreen")
results_excel <- createWorkbook()
results_excel.file <- paste0(run, "_", sample_id, "_SmallVariant_", format(Sys.time(), "%Y%m%d"), ".xlsx")


# Load data ---------------------------------------------------------------

# Load sample info
sample_info <- read.table(sample.file, sep = "\t", header = TRUE, quote = "", fill = FALSE) %>% filter(Sample == sample)
tumor_doid <- sample_info[1, "TumorDOID"]
tumor_topnode_doid <- sample_info[1, "TopNodeDOID"]
tumor_whitelist <- sample_info[1, "Whitelist"]

# Load RUN annotation
run_annot <- read.table(annotation.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)

# Load consensus variants
consensus_variants <- read.table(consensus.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
	mutate(VAR = as.character(VAR))
rownames(consensus_variants) <- consensus_variants$VAR
vars <- consensus_variants %>% select(RUN:VAR_TVC, PASS:VAR_HG19, -CALLER)

# Genes
target_genes <- (read.table(genes.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>%
	mutate(NEW_GENE = ifelse(!is.na(NEW_GENE), NEW_GENE, OLD_GENE)))$NEW_GENE
genes_whitelist <- read.table(whitelist.file, header = TRUE, sep = ";", fill = TRUE, na.strings = c("NA", ".", ""))
genes_whitelist <- as.character(na.omit(unique(genes_whitelist[[tumor_whitelist]])))
genes_whitelist <- genes_whitelist[genes_whitelist %in% target_genes]

# MANE data (only target genes)
mane_TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(mane.gtf, format = "gtf")))
mane_gene <- fread(mane_gene.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% target_genes)
mane_gene_names <- mane_gene$GENE_ENSEMBL; names(mane_gene_names) <- mane_gene$GENE_SYMBOL
mane_exon <- makeGRangesFromDataFrame(
	fread(mane_exon.file, header = TRUE, sep = "\t", quote = "", nThread = 1) %>% filter(GENE_SYMBOL %in% target_genes), 
	keep.extra.columns = TRUE)

# Load GENIE oncogenic mutations
if(file.exists(genie.file)){
	genie_mut <- read.table(genie.file, header = TRUE, sep = "\t", quote = "", fill = FALSE) %>% 
		arrange(ONCO_POINTS, GENIE_CNT) %>% select(MUT_ID, GENIE_CNT, GENE, CHROM, POS, ONCOGENICITY)
} else {
	genie_mut <- data.frame(MUT_ID = NA, GENIE_CNT = NA, GENE = NA, CHROM = NA, POS = NA, ONCOGENICITY = NA)
}

# CIViC clinical data
civic_data <- read.table(civic.file, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = FALSE)
colnames(civic_data) <- paste0("CIViC_", colnames(civic_data))



# Data processing ---------------------------------------------------------

# Join the annotation info from the annotated run variants
# Add variant FLAGS:
#	- OK: No flags
#	- NoCall (LowCallers, LowAD, LowDP, LowVAF, Blacklist)
#	- Germline: (gnomAD_MAX_AF > germlinepVAF | AF > germlineVAF) & (!(CANCER_HOTSPOT | PANEL_HOTSPOT))
#	- ProblematicRegion: found in a problematic region of the genome and NOT cancer HOTSPOT
#	- Recurrent: VAR found in the list of panel recurrent mutations and NOT cancer HOTSPOT
annot_vars <- left_join(vars, run_annot, by = "VAR") %>% 
	mutate(
		GENE_WHITELIST = GENE_SYMBOL %in% genes_whitelist, 
		GERMLINE = (!(CANCER_HOTSPOT | PANEL_HOTSPOT)) & (gnomAD_MAX_AF > somatic_max_pvaf | AF > somatic_max_vaf),
		NO_FLAGS = !(
			!PASS | GERMLINE | (!CTR_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT)) | 
			(PROBLEMATIC_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT)) | RECURRENT)) %>%
	rowwise() %>%
	mutate(FLAGS = paste(var_flags[c(
		NO_FLAGS, !PASS, GERMLINE, (!CTR_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT)), 
		(PROBLEMATIC_REGION & !(CANCER_HOTSPOT | PANEL_HOTSPOT)), RECURRENT)], collapse = ";")) %>% 
	as.data.frame()

# Clinical classification (AMP/ASCO/CAP guidelines)
#	- Predictive, prognostic & diagnostic evidences in CIViC resource
# 	- Tier IA: CIViC evidence A (guidelines) in the SAME tumor
# 	- Tier IB: CIViC evidence B (trials) in the SAME tumor
# 	- Tier IIC: CIViC evidence A/B in OTHER tumor OR CIViC evidence C
# 	- Tier IID: CIViC evidence D
# 	- Tier III: NOT (Likely) Benign based on oncogenicity
# 	- Tier IV: (Likely) Benign based on oncogenicity
vars_civic <- annot_vars %>% 
	filter(CIViC_VARIANT_ID %in% civic_data$CIViC_VARIANT_ID | MUT_ID %in% civic_data$CIViC_ALTERATION | OTHER_ID %in% civic_data$CIViC_ALTERATION)
civic_annot <- civic_data %>% mutate(VAR = ".", CLIN_STATUS = ".") %>% filter(CIViC_ALTERATION == ".") %>% select(-CIViC_VARIANT_ID)
for(n in seq_len(nrow(vars_civic))){
	aux_var <- vars_civic[n, ]
	aux_civic <- civic_data %>% 
		filter(CIViC_VARIANT_ID %in% aux_var$CIViC_VARIANT_ID | CIViC_ALTERATION %in% c(aux_var$MUT_ID, aux_var$OTHER_ID)) %>%
		mutate(
			VAR = aux_var$VAR,
			isTumor = CIViC_TUMOR_DOID %in% tumor_doid | CIViC_TOPNODE_DOID %in% tumor_topnode_doid | CIViC_TUMOR_NAME %in% "cancer",
			CLIN_STATUS = ifelse(
				isTumor & CIViC_EVIDENCE_LEVEL %in% "A", "Tier IA", ifelse(
					isTumor & CIViC_EVIDENCE_LEVEL %in% "B", "Tier IB", ifelse(
						(CIViC_EVIDENCE_LEVEL %in% c("A", "B") & !isTumor) | CIViC_EVIDENCE_LEVEL %in% "C", "Tier IIC", ifelse(
							CIViC_EVIDENCE_LEVEL %in% "D", "Tier IID", NA))))) %>%
		select(-isTumor, -CIViC_VARIANT_ID) %>%
		arrange(CLIN_STATUS) %>%
		head(1)
	civic_annot <- rbind(civic_annot, aux_civic)
	rm(aux_var, aux_civic)
}
rm(n)
annot_vars <- annot_vars %>% left_join(civic_annot, by = "VAR") %>% select(-CIViC_GENE, -CIViC_MUT) %>% 
	mutate(
		CLIN_STATUS = ifelse(!is.na(CLIN_STATUS), CLIN_STATUS, ifelse(
			ONCOGENICITY %in% c("Oncogenic", "Likely Oncogenic", "VUS"), "Tier III", "Tier IV")),
		CLIN_RELEVANT = CLIN_STATUS %in% c("Tier IA", "Tier IB", "Tier IIC", "Tier IID"))

# Create final table of results
annot_vars <- annot_vars %>% select(
	RUN, SAMPLE, TYPE, CLASS, VAR, AD_ALT, AF, GENE_SYMBOL, MUTATION, FLAGS, FILTER, PANEL_HOTSPOT,
	CLIN_RELEVANT, CLIN_STATUS, ONCOGENICITY, ONCO_POINTS, ONCO_CODES, 
	CONSEQUENCE, IMPACT, EXON, INTRON, DP, AD_REF, OVERLAP_CALLERS, MATCH_CALLERS, CALLERS,
	HGVSg, HGVSc, HGVSp, HGVSp_SHORT, HGVSc_ENSEMBL, HGVSc_REFSEQ, HGVSp_ENSEMBL, 

	dbSNP_ID, ClinVar_ID, EXISTING_VARIATION, gnomAD_MAX_AF,
	NMD_ESCAPE, GENIE_CNT, SOMATIC_WHITELIST, HOTSPOT_MUT_CNT, HOTSPOT_POS_CNT, BENIGN, 

	GENE_WHITELIST, MMR_GENE, ONCOGENE, TSG, DRIVER, CANONICAL_DRIVER, GENE_ENSEMBL, GENE_HGNC, BIOTYPE, CANONICAL, 
	MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, TSL, TRANSCRIPT_ENSEMBL, TRANSCRIPT_REFSEQ, 
	cDNA_POS, CDS_POS, CCDS, PROTEIN_ENSEMBL, PROTEIN_POS, AA, CODONS, CODING, NONCODING, 

	contains("gnomAD"),

	POLYPHEN_TERM, POLYPHEN_SCORE, SIFT_TERM, SIFT_SCORE, REVEL_TERM, REVEL_SCORE, AlphaMissense_TERM, AlphaMissense_SCORE, 
	CADD_TERM, CADD_SCORE, PREDICTOR_PATHOGENIC, PREDICTOR_BENIGN,

	contains("ClinVar_"),

	PASS, LowCallers, LowAD, LowDP, LowVAF, PANEL_BLACKLIST,

	NO_FLAGS, CANCER_HOTSPOT, GERMLINE, CTR_REGION, PROBLEMATIC_REGION, RECURRENT,
	
	NULL_VARIANT, STOPLOSS_VARIANT, INFRAME_INDEL, MISSENSE_VARIANT, SILENT_VARIANT,
	ONCOGENIC_SOP_MUT, ONCOGENIC_SOP_POS, ONCOGENIC_VALID_VAR, ONCOGENIC_VALID_MUT, SBVS1:OVS1,

	MUT_ID, CIViC_ALTERATION, contains("CIViC_"),

	CHROM, START, END, REF, ALT, STRAND, VAR_HG19, VAR_MUTECT2, VAR_OCTOPUS, VAR_PISCES, VAR_VARDICT, VAR_TVC
)
annot_vars[is.na(annot_vars) | annot_vars == "NA" | annot_vars == ""] <- "."

# Generate a VCF file
aux_header <- readLines(headerVcf.file)
aux_header[grep("##fileDate", aux_header)] <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
aux_header[grep("##reference", aux_header)] <- paste0("##reference=file:", fasta_name)
if(nrow(annot_vars)){
	aux_vars <- annot_vars %>% mutate(AD = paste0(AD_REF, ",", AD_ALT), GT = ifelse(AD_REF != 0, "0/1", "1/1"))
	aux_info <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c(
			"VAR", "RUN", "SAMPLE", "TYPE", "AF", "DP", "AD_REF", "AD_ALT", "END", "VAR_HG19", "MATCH_CALLERS", "OVERLAP_CALLERS", 
			"CLASS", "GENE_SYMBOL", "MUTATION", "CLIN_STATUS", "ONCOGENICITY", "ONCO_POINTS", "ONCO_CODES",
			"CONSEQUENCE", "IMPACT", "EXON", "INTRON", "HGVSg", "HGVSc", "HGVSp", "HGVSp_SHORT", "HGVSc_ENSEMBL", "HGVSc_REFSEQ", "HGVSp_ENSEMBL", 
			
			"dbSNP_ID", "ClinVar_ID", "EXISTING_VARIATION", "gnomAD_MAX_AF",
			"NMD_ESCAPE", "GENIE_CNT", "HOTSPOT_MUT_CNT", "HOTSPOT_POS_CNT",

			"GENE_ENSEMBL", "GENE_HGNC", "BIOTYPE", "APPRIS", "TSL", "TRANSCRIPT_ENSEMBL", "TRANSCRIPT_REFSEQ", 
			"cDNA_POS", "CDS_POS", "CCDS", "PROTEIN_ENSEMBL", "PROTEIN_POS", "AA", "CODONS", 

			"gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF", "gnomADe_NFE_AF", 
			"gnomADe_OTH_AF", "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
			"gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_OTH_AF", "gnomADg_SAS_AF",

			"POLYPHEN_TERM", "POLYPHEN_SCORE", "SIFT_TERM", "SIFT_SCORE", "REVEL_TERM", "REVEL_SCORE", 
			"AlphaMissense_TERM", "AlphaMissense_SCORE", "CADD_SCORE", "CADD_TERM",

			"ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN", "ClinVar_ORIGIN",

			"MUT_ID", "CIViC_ALTERATION", "CIViC_VARIANT_ID", "CIViC_EVIDENCE_LEVEL", "CIViC_EVIDENCE_RATING", "CIViC_EVIDENCE_TYPE", 
			"CIViC_EFFECT", "CIViC_DRUG", "CIViC_TUMOR_NAME", "CIViC_TUMOR_DOID", "CIViC_TUMOR_CODE", "CIViC_TOPNODE_NAME", 
			"CIViC_TOPNODE_DOID", "CIViC_MOLECULAR_ID", "CIViC_VARIANT_GROUP", "CIViC_EVIDENCE_SCORE", "CIViC_ORIGIN",
			
			"CHROM", "START", "REF", "ALT", "STRAND", "VAR_MUTECT2", "VAR_OCTOPUS", "VAR_PISCES", "VAR_VARDICT"),
		function(x) paste0(x, "=", aux_vars[n, x]) )), collapse = ";")  }))
	aux_flags <- paste0(
		aux_vars$CALLERS, ";", aux_vars$FLAGS, ";",
		unlist(lapply(seq_len(nrow(aux_vars)), function(n){
			paste(unlist(lapply(c(
				"PANEL_HOTSPOT", "PANEL_BLACKLIST", "SOMATIC_WHITELIST", "CANCER_HOTSPOT", 
				"GENE_WHITELIST", "MMR_GENE", "ONCOGENE", "TSG", "DRIVER", "CANONICAL_DRIVER",
				"CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CODING", "NONCODING", "CTR_REGION",
				"PREDICTOR_PATHOGENIC", "PREDICTOR_BENIGN", "ClinVar_PATHOGENIC", "PROBLEMATIC_REGION",
				"NULL_VARIANT", "STOPLOSS_VARIANT", "INFRAME_INDEL", "MISSENSE_VARIANT", "SILENT_VARIANT",
				"ONCOGENIC_SOP_MUT", "ONCOGENIC_SOP_POS", "ONCOGENIC_VALID_VAR", "ONCOGENIC_VALID_MUT",
				"SBVS1", "SBS1", "SBP2", "SBP1", "OP4", "OP3", "OP1", "OM3", "OM4", "OM2", "OS3", "OS2", "OS1", "OVS1", "CLIN_RELEVANT"), 
			function(x) ifelse(aux_vars[n, x], paste0(x, ";"), "")  )), collapse = "")  })))
	aux_info <- paste0(aux_info, ";", aux_flags)
	aux_format <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("GT", "AD", "AF", "DP"), function(x) paste0(aux_vars[n, x]) )), collapse = ":")  }))
	aux_df <- aux_vars %>% mutate(
		POS = START, ID = ".", QUAL = ".", FILTER = FILTER, INFO = aux_info,
		FORMAT = paste(c("GT", "AD", "AF", "DP"), collapse = ":"), SAMPLE = aux_format) %>%
		select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)
	vcf_text <- c(
		aux_header,
		paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"),
		unlist(lapply(seq_len(nrow(aux_df)), function(x) paste(aux_df[x, ], collapse = "\t"))))

	rm(aux_header, aux_vars, aux_info, aux_flags, aux_format, aux_df)
} else {
	vcf_text <- c(aux_header, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_id), collapse = "\t"))
}



# Save data ---------------------------------------------------------------

# Set factors
# Sort by CLIN_RELEVANT, FLAGS, FILTER, SOMATIC, CLIN_STATUS, PANEL_HOTSPOT, ONCOGENICITY, WHITELIST, GENE_SYMBOL, AF
annot_vars <- annot_vars %>% 
	mutate(ONCOGENICITY = factor(ONCOGENICITY, levels = oncogenicity_levels)) %>%
	arrange(
		desc(CLIN_RELEVANT), desc(NO_FLAGS), desc(PASS), GERMLINE, CLIN_STATUS, desc(PANEL_HOTSPOT), 
		ONCOGENICITY, desc(ONCO_POINTS), desc(GENE_WHITELIST), GENE_SYMBOL, desc(AF))

# Kepp OK (no flags), somatic and germline variants
annot_ok <- annot_vars %>% filter(NO_FLAGS)
annot_somatic <- annot_vars %>% filter(!GERMLINE)
annot_germline <- annot_vars %>% filter(GERMLINE)

# Save the annotated variants
write.table(annot_vars, file = paste0(sample_id, "_SmallVariant_annot.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save the VCF file
writeLines(con = paste0(sample_id, "_SmallVariant_annot.vcf"), text = vcf_text)

# An excel file with the results
addWorksheet(results_excel, "OK")
writeDataTable(results_excel, "OK", annot_ok)
addWorksheet(results_excel, "Somatic")
writeDataTable(results_excel, "Somatic", annot_somatic)
addWorksheet(results_excel, "Germline")
writeDataTable(results_excel, "Germline", annot_germline)
addWorksheet(results_excel, "Consensus")
consensus_variants <- consensus_variants[annot_vars$VAR, ]
writeDataTable(results_excel, "Consensus", consensus_variants)
saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)


# Data visualization ---------------------------------------------------------------

# Plot genes with top detected variants
#	- Tier I, Tier II and Oncogenic variants always considered
# 	- OK non-benign variants (If there is a panel hotspot focus on that)
top_vars <- annot_vars %>% filter(
	CLIN_RELEVANT | ONCOGENICITY %in% "Oncogenic" | (NO_FLAGS & !BENIGN & ((!is_panel_hotspots) | (is_panel_hotspots & PANEL_HOTSPOT))))
plot_genes <- unique(top_vars$GENE_SYMBOL)

# Set ranges of somatic variants for plotting in each selected gene
#	- Only genes with top variants are selected
#	- OK or selected top variants
plot_somatic_vars <- annot_vars %>% filter(GENE_SYMBOL %in% plot_genes & (NO_FLAGS | VAR %in% top_vars$VAR))
plot_somatic_vars_ranges <- makeGRangesFromDataFrame(plot_somatic_vars, keep.extra.columns = TRUE)

# Set ranges of GENIE pathogenic mutations
genie_mut_ranges <- makeGRangesFromDataFrame(
	start.field = "POS", end.field = "POS", keep.extra.columns = TRUE,
	genie_mut %>% 
		mutate(DETECTED = MUT_ID %in% c(plot_somatic_vars$MUT_ID, plot_somatic_vars$VAR)) %>% 
		left_join(plot_somatic_vars, "MUT_ID") %>% 
		mutate(CHROM = CHROM.x, GENIE_CNT = GENIE_CNT.x, ONCOGENICITY = ONCOGENICITY.x) %>%
		select(CHROM, POS, GENE, GENIE_CNT, DETECTED, CLIN_RELEVANT, CLIN_STATUS, ONCOGENICITY) %>% arrange(DETECTED, desc(CLIN_STATUS), ONCOGENICITY)
)

if(nrow(top_vars)){
	
	addWorksheet(results_excel, "Plot_ByGene")

	# Setup karyoplot parameters
	plot.params <- getDefaultPlotParams(plot.type = 3)
	plot.params$leftmargin <- 0.06
	plot.params$rightmargin <- 0.1
	plot.params$topmargin <- 35
	plot.params$bottommargin <- 15
	plot.params$ideogramheight <- 7
	plot.params$data1inmargin <- 0
	plot.params$data1height <- 350
	plot.params$data2height <- 0

	excel_startRow <- 1

	for(gene in plot_genes){

		aux_vars <- plot_somatic_vars_ranges[plot_somatic_vars_ranges$GENE_SYMBOL %in% gene]
		aux_vars <- aux_vars[rev(seq_len(length(aux_vars)))]
		aux_vars_colors <- unlist(lapply(aux_vars$ONCOGENICITY, function(x) oncogenicity_colors[x]))
		aux_vars_shape1 <- ifelse(aux_vars$CLIN_RELEVANT, 17, 16)
		aux_vars_shape2 <- ifelse(aux_vars$CLIN_RELEVANT, 2, 1)

		aux_top_vars <- aux_vars[aux_vars$VAR %in% top_vars$VAR]

		aux_genie <- genie_mut_ranges[genie_mut_ranges$GENE %in% gene]
		aux_genie_notdetected <- aux_genie[!aux_genie$DETECTED]
		aux_genie_notdetected_colors <- ifelse(aux_genie_notdetected$ONCOGENICITY %in% "Oncogenic", "grey30", "grey70")
		aux_genie_detected <- aux_genie[aux_genie$DETECTED]
		aux_genie_detected_colors <- unlist(lapply(aux_genie_detected$ONCOGENICITY, function(x) oncogenicity_colors[x]))
		aux_genie_detected_shape1 <- ifelse(aux_genie_detected$CLIN_RELEVANT, 17, 16)
		aux_genie_detected_shape2 <- ifelse(aux_genie_detected$CLIN_RELEVANT, 2, 1)
		
		aux_exons <- c(subsetByOverlaps(mane_exon, aux_vars)$EXON, subsetByOverlaps(mane_exon, aux_genie)$EXON)
		if(!length(aux_exons)){next}
		aux_exons <- min(aux_exons):max(aux_exons)
		aux_exon_ranges <- mane_exon[mane_exon$GENE_SYMBOL %in% gene & mane_exon$EXON %in% aux_exons]
		if(!length(aux_exon_ranges)){next}
		aux_gene <- mane_gene_names[unique(aux_exon_ranges$GENE_SYMBOL)]
		aux_transcript <- unique(aux_exon_ranges$TRANSCRIPT_ENSEMBL)
		
		aux_region <- toGRanges(unique(seqnames(aux_exon_ranges)), min(start(aux_exon_ranges)), max(end(aux_exon_ranges)))+200
		width <- width(aux_region)
		zoom <- aux_region + width*0.02

		png(paste0(sample_id, "_SmallVariant_", gene, ".png"), width = 3000, height = 1200, res = 200)

		kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, cex = 1.8, plot.params = plot.params, zoom = zoom)
		kpAddBaseNumbers(kp, cex = 1, minor.ticks = FALSE, tick.len = 5, add.units = TRUE, digits = 3, tick.dist=round(width(zoom)/5, -2))
		kpAddChromosomeNames(kp, cex = 1.2, yoffset = -360, xoffset = -0.44)
		kpAddMainTitle(kp, main = paste0(gene, " (", aux_transcript, ")"), cex = 1.7)

		genes.data <- makeGenesDataFromTxDb(mane_TxDb, kp, plot.transcripts = TRUE, plot.transcripts.structure = TRUE)
		is_gene <- names(genes.data$genes) %in% aux_gene
		genes.data$genes <- genes.data$genes[is_gene]
		genes.data$transcripts <- genes.data$transcripts[is_gene]
		is_transcript <- genes.data[["transcripts"]][[aux_gene]]$tx_name %in% aux_transcript
		genes.data[["transcripts"]][[aux_gene]] <- genes.data[["transcripts"]][[aux_gene]][is_transcript]
		genes.data[["coding.exons"]] <- genes.data[["coding.exons"]][is_transcript]
		genes.data[["non.coding.exons"]] <- genes.data[["non.coding.exons"]][is_transcript]

		kpPlotMarkers(
			kp, aux_exon_ranges, labels = aux_exon_ranges$EXON, adjust.label.position = TRUE, label.dist = 0.0003, 
			lwd = 0.4, cex = 0.8, text.orientation = "horizontal", marker.parts = c(0.95, 0.045, 0.005), data.panel = 1, r0 = 0, r1 = 1.3, 
			line.color = "grey75", label.color = "black", ignore.chromosome.ends = TRUE)
		suppressWarnings(kpPlotGenes(
			kp, genes.data, add.gene.names=FALSE, col="black", data.panel = 1, r0 = 0.5, r1 = 0.55, add.strand.marks=TRUE, mark.height = 0.25,
			coding.exons.col="black", coding.exons.border.col="black", non.coding.exons.col="grey50", non.coding.exons.border.col="grey50"))			
		
		kpAddLabels(kp, labels = "GENIE oncogenic variants (counts)", cex = 0.8, label.margin = 0.04, srt = 90, pos = 1, data.panel = 1, r0 = 0.5, r1 = 0.05)
		if(length(aux_genie)){
			ymax <- plyr::round_any(max(aux_genie$GENIE_CNT), 100, f = ceiling)
			ymax_sep <- ifelse(ymax == 100, plyr::round_any(ymax/5, 10, f = ceiling), plyr::round_any(ymax/5, 100, f = ceiling))
			kpAxis(
				kp, ymin = 0, ymax = ymax, cex = 0.7, data.panel = 1, r0 = 0.5, r1 = 0.05, tick.pos = seq(0, ymax, ymax_sep), tick.len = width(zoom)/200)
			if(length(aux_genie_notdetected)){
				kpPoints(
					kp, aux_genie_notdetected, y = aux_genie_notdetected$GENIE_CNT, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.5, r1 = 0.05, 
					cex = 0.6, pch = 1, col = aux_genie_notdetected_colors)
				kpPoints(
					kp, aux_genie_notdetected, y = aux_genie_notdetected$GENIE_CNT, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.5, r1 = 0.05, 
					cex = 0.6, pch = 16, col = aux_genie_notdetected_colors)
			}
			if(length(aux_genie_detected)){
				kpPoints(
					kp, aux_genie_detected, y = aux_genie_detected$GENIE_CNT, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.5, r1 = 0.05, 
					cex = 0.8, pch = aux_genie_detected_shape1, col = aux_genie_detected_colors)
				kpPoints(
					kp, aux_genie_detected, y = aux_genie_detected$GENIE_CNT, ymin = 0, ymax = ymax, data.panel = 1, r0 = 0.5, r1 = 0.05, 
					cex = 0.8, pch = aux_genie_detected_shape2, col = "black")
			}
		}

		kpAxis(kp, ymin = 0, ymax = 1, cex = 0.7, data.panel = 1, r0 = 0.535, r1 = 0.95, numticks = 5, tick.len = width(zoom)/200)
		kpAddLabels(kp, labels = "Allele frequency", cex = 0.8, label.margin = 0.04, srt = 90, pos = 1, data.panel = 1, r0 = 0.535, r1 = 0.95)
		kpAbline(kp, h = min_vaf, data.panel = 1, r0 = 0.535, r1 = 0.95, col = "red", ymin = 0, ymax = 1, lwd = 0.3)
		suppressWarnings(kpPlotMarkers(
			kp, aux_top_vars, labels = aux_top_vars$MUTATION, y = aux_top_vars$AF+0.05, cex = 0.6, text.orientation = "vertical",
			adjust.label.position = TRUE, data.panel = 1, r0 = 0.535, r1 = 0.95, label.margin = 0.01, label.dist = 0.001, 
			marker.parts = c(0.94, 0.055, 0.005), line.color = "black", label.color = "black", ignore.chromosome.ends = TRUE))
		kpBars(
			kp, aux_vars, y1 = aux_vars$AF, ymin = 0, ymax = 1, data.panel = 1, r0 = 0.535, r1 = 0.95, cex = 0.8, 
			col = aux_vars_colors, border = aux_vars_colors)
		kpPoints(
			kp, aux_vars, y = aux_vars$AF, ymin = 0, ymax = 1, data.panel = 1, r0 = 0.535, r1 = 0.95, cex = 0.8, 
			pch = aux_vars_shape1, col = aux_vars_colors)
		kpPoints(
			kp, aux_vars, y = aux_vars$AF, ymin = 0, ymax = 1, data.panel = 1, r0 = 0.535, r1 = 0.95, cex = 0.8, 
			pch = aux_vars_shape2, col = "black")

		legend(
			0.92, 0.97, legend = c("Tier I/II", "Tier III/IV"), col = "black", title = "AMP/ASCO/CAP", 
			pch = c(17, 16), cex = 0.6, pt.cex = 0.7)
		legend(
			0.92, 0.8, legend = names(oncogenicity_colors), col = oncogenicity_colors, title = "ClinGen/CGC/VICC", 
			pch = 16, cex = 0.6, pt.cex = 0.7)
		legend(
			0.92, 0.40, legend = c("Oncogenic", "Likely Oncogenic"), col = c("grey30", "grey70"), title = "ClinGen/CGC/VICC", 
			pch = 16, cex = 0.6, pt.cex = 0.7)

		dev.off()

		insertImage(
			results_excel, "Plot_ByGene", paste0(sample_id, "_SmallVariant_", gene, ".png"), 
			units = "px", width = 3000, height = 1200, dpi = 200, startRow = excel_startRow)
		saveWorkbook(results_excel, results_excel.file, overwrite = TRUE)
		
		excel_startRow <- excel_startRow + 31
		
	}
	rm(gene, excel_startRow)
}

