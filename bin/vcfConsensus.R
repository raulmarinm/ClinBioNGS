#!/usr/bin/env Rscript

# ============================================================================ #
# What: VCF consensus from multiple callers
# ============================================================================ #

# General settings --------------------------------------------------------

# Set default R options
.libPaths("opt/R/lib/R/library")
options(dplyr.summarise.inform = FALSE)

# Load libraries
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
suppressPackageStartupMessages(suppressWarnings(library(VariantAnnotation)))
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
chain <- args_list[["chain"]]
panel_hotspots.bed <- args_list[["panelHotspotsBed"]]
panel_blacklist.bed <- args_list[["panelBlacklistBed"]]
genome <- args_list[["genome"]]
vcfHeader <- args_list[["vcfHeader"]]
caller_names <- unlist(strsplit(args_list[["callers"]], ","))
min_callers <- as.numeric(args_list[["minCallers"]])
min_ad <- as.numeric(args_list[["minAD"]])
min_dp <- as.numeric(args_list[["minDP"]])
min_vaf <- as.numeric(args_list[["minVAF"]])
min_ad_hotspot <- as.numeric(args_list[["minADHotspot"]])
min_dp_hotspot <- as.numeric(args_list[["minDPHotspot"]])
min_vaf_hotspot <- as.numeric(args_list[["minVAFHotspot"]])

# Set other variables
sample_id <- paste0(sample, "_", type)
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
callers_AD_REF <- c("pisces" = "AD_1", "mutect2" = "AD_1", "vardict" = "AD_1", "octopus" = "AD_1", "tvc" = "RO_1")
callers_AD_ALT <- c("pisces" = "AD_2", "mutect2" = "AD_2", "vardict" = "AD_2", "octopus" = "AD_2", "tvc" = "AO_1")
callers_DP <- c("pisces" = "DP_1", "mutect2" = "DP_1", "vardict" = "DP_1", "octopus" = "ADP_1", "tvc" = "DP_1")
callers_AF <- c("pisces" = "VF_1", "mutect2" = "AF_1", "vardict" = "AF_1", "octopus" = "AF_2")
caller_flag <- c("Mutect2", "Octopus", "Pisces", "TVC", "VarDict")
filter_flags <- c("PASS", "LowCallers", "LowAD", "LowDP", "LowVAF", "Blacklist")

# Set VCF files
callers.vcf <- list.files(pattern = "\\.vcf$")


# Load data ---------------------------------------------------------------

# Panel hotspots regions
if(file.exists(panel_hotspots.bed)){
	panel_hotspots <- read.table(panel_hotspots.bed, header = FALSE, sep = "\t", quote = "", fill = FALSE)[, 1:3]
	colnames(panel_hotspots) <- c("seqnames", "start", "end")
	panel_hotspots_ranges <- makeGRangesFromDataFrame(
		makeGRangesFromDataFrame(panel_hotspots, starts.in.df.are.0based = TRUE) %>%
			as.data.frame() %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
			mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.frame())
} else {panel_hotspots_ranges <- GRanges()}

# Panel blacklist regions
if(file.exists(panel_blacklist.bed)){
	panel_blacklist <- read.table(panel_blacklist.bed, header = FALSE, sep = "\t", quote = "", fill = FALSE)[, 1:3]
	colnames(panel_blacklist) <- c("seqnames", "start", "end")
	panel_blacklist_ranges <- makeGRangesFromDataFrame(
		makeGRangesFromDataFrame(panel_blacklist, starts.in.df.are.0based = TRUE) %>%
			as.data.frame() %>% mutate(ID = 1:n()) %>% uncount(width) %>% group_by(ID) %>% 
			mutate(pos = 1:n(), end = start+pos-1, start=end) %>% as.data.frame())
} else {panel_blacklist_ranges <- GRanges()}

# Collect variants from each caller
# 	- Make an intra-caller consensus of multiallelic variants
for (caller in caller_names){

	## Read VCF file
	vars.vcf <- suppressWarnings(readVcf(grep(caller, callers.vcf, value = TRUE)))

	# Keep unique variants by position based on the AF
	if(length(vars.vcf)){
		if(caller == "tvc"){
			aux_ad <- unlist(strsplit(callers_AD_ALT[caller], "_"))
			aux_dp <- unlist(strsplit(callers_DP[caller], "_"))
			aux_af <- sapply(geno(vars.vcf)[[aux_ad[1]]], `[`, as.integer(aux_ad[2])) / sapply(geno(vars.vcf)[[aux_dp[1]]], `[`, as.integer(aux_dp[2]))
			rm(aux_ad, aux_dp)
		} else {
			aux_af <- sapply(geno(vars.vcf)[[unlist(strsplit(callers_AF[caller], "_"))[1]]], `[`, as.integer(unlist(strsplit(callers_AF[caller], "_"))[2]))
		}

		consensus_df <- data.frame(VAR = rownames(vars.vcf), POS = paste0(as.character(seqnames(vars.vcf)), ":", start(vars.vcf)) , AF = aux_af) %>% 
			group_by(POS) %>% filter(AF == max(AF)) %>% ungroup() %>% distinct(POS, .keep_all = TRUE)
		
		rm(aux_af)
	} else {
		consensus_df <- data.frame(VAR = NULL)
	}
	unique.vcf <- vars.vcf[consensus_df$VAR]

	## Assign variants to each caller
	assign(paste0(caller, ".vcf"), unique.vcf)
	assign(paste0(caller, "_names"), rownames(unique.vcf))
	assign(paste0(caller, ".ranges"), rowRanges(unique.vcf))

	rm(vars.vcf, unique.vcf)
	
}
rm(caller)



# Data processing ---------------------------------------------------------

# Get MATCH_CALLERS and OVERLAP_CALLERS info
for (caller in caller_names){
	
	other_callers <- grep(caller, caller_names, value = TRUE, invert = TRUE)

	vars.vcf <- get(paste0(caller, ".vcf"))
	vars_names <- get(paste0(caller, "_names"))
	vars.ranges <- get(paste0(caller, ".ranges"))

	vars_match <- lapply(other_callers, function(x) vars_names %in% get(paste0(x, "_names")))
	names(vars_match) <- unlist(lapply(other_callers, function(x) paste0("MATCH_", x)))
	vars_match <- as.data.frame(vars_match)
	rownames(vars_match) <- vars_names
	vars_match$MATCH_CALLERS <- rowSums(vars_match)+1

	vars_overlap <- lapply(other_callers, function(x) overlapsAny(vars.ranges, get(paste0(x, ".ranges"))))
	names(vars_overlap) <- unlist(lapply(other_callers, function(x) paste0("OVERLAP_", x)))
	vars_overlap <- as.data.frame(vars_overlap)
	rownames(vars_overlap) <- vars_names
	vars_overlap$OVERLAP_CALLERS <- rowSums(vars_overlap)+1

	vars_info <- cbind(vars_match, vars_overlap)
	vars_info$VAR <- rownames(vars_info)

	assign(paste0(caller, "_info"), vars_info)

	rm(other_callers, vars.vcf, vars_names, vars.ranges, vars_match, vars_overlap, vars_info)
}
rm(caller)


# Convert Hg38 to Hg19 coordinates
for (caller in caller_names){

	vars.vcf <- get(paste0(caller, ".vcf"))

	if(length(vars.vcf)){
		if (!seqlevelsStyle(vars.vcf) == "UCSC") { stop(paste0("Genome style must be UCSC")) }

		vars_hg19_ranges <- unlist(liftOver(vars.vcf, import.chain(chain)))
		vars_hg19_ranges <- vars_hg19_ranges[!duplicated(names(vars_hg19_ranges))]
		vars_hg19 <- data.frame(
			VAR = names(vars_hg19_ranges), CHROM = seqnames(vars_hg19_ranges),
			START_HG19 = start(vars_hg19_ranges), END_HG19 = end(vars_hg19_ranges),
			REF = vars_hg19_ranges$REF, ALT = sapply(vars_hg19_ranges$ALT, function(x) as.character(x[[1]]))) %>%
			mutate(VAR_HG19 = paste0(CHROM, ":", START_HG19, "_", REF, "/", ALT)) %>%
			select(VAR, VAR_HG19, START_HG19, END_HG19)

		assign(paste0(caller, "_hg19"), vars_hg19)

		rm(vars.vcf, vars_hg19_ranges, vars_hg19)
	}
}
rm(caller)


# Final table of variants
for (caller in caller_names){

	vars.vcf <- get(paste0(caller, ".vcf"))

	if(length(vars.vcf)){

		vars_ranges <- rowRanges(vars.vcf)
		vars_hg19 <- get(paste0(caller, "_hg19"))
		vars_info <- get(paste0(caller, "_info"))

		aux_ad_ref <- sapply(geno(vars.vcf)[[unlist(strsplit(callers_AD_REF[caller], "_"))[1]]], `[`, as.integer(unlist(strsplit(callers_AD_REF[caller], "_"))[2]))
		aux_ad_alt <- sapply(geno(vars.vcf)[[unlist(strsplit(callers_AD_ALT[caller], "_"))[1]]], `[`, as.integer(unlist(strsplit(callers_AD_ALT[caller], "_"))[2]))
		aux_dp <- sapply(geno(vars.vcf)[[unlist(strsplit(callers_DP[caller], "_"))[1]]], `[`, as.integer(unlist(strsplit(callers_DP[caller], "_"))[2]))
		if(caller == "tvc"){
			aux_af <- aux_ad_alt / aux_dp
		} else {
			aux_af <- sapply(geno(vars.vcf)[[unlist(strsplit(callers_AF[caller], "_"))[1]]], `[`, as.integer(unlist(strsplit(callers_AF[caller], "_"))[2]))
		}

		variants <- data.frame(
			RUN = run, SAMPLE = sample, CALLER = caller, TYPE = ifelse(isSNV(vars.vcf), "SNV", "INDEL"), VAR = rownames(vars.vcf),
			AF = round(aux_af, 4), DP = aux_dp, AD_REF = aux_ad_ref, AD_ALT = aux_ad_alt,
			OVERLAP_CALLERS = vars_info[rownames(vars.vcf), "OVERLAP_CALLERS"], MATCH_CALLERS = vars_info[rownames(vars.vcf), "MATCH_CALLERS"],
			CHROM = seqnames(vars_ranges), START = start(vars_ranges), END = end(vars_ranges),
			REF = vars_ranges$REF, ALT = sapply(vars_ranges$ALT, function(x) as.character(x[[1]]))) %>%
			left_join(vars_hg19, by = "VAR")
		
		variants[is.na(variants)] <- "."


		rm(vars.vcf, vars_ranges, vars_hg19, vars_info, aux_ad_ref, aux_ad_alt, aux_dp, aux_af)

	} else {
		variants <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 19))
		colnames(variants) <- c(
			"RUN", "SAMPLE", "CALLER", "TYPE", "VAR", "AF", "DP", "AD_REF", "AD_ALT", "OVERLAP_CALLERS", 
			"MATCH_CALLERS", "CHROM", "START", "END", "REF", "ALT", "VAR_HG19", "START_HG19", "END_HG19")
	}
	
	assign(paste0(caller, "_variants"), variants)
	rm(variants)
	
}
rm(caller)


# Join all variants from each caller
# Sort correctly by chr (chr1, chr2, ...)
all_variants <- do.call(rbind, mget(unlist(lapply(caller_names, function(x) paste0(x, "_variants")))))
all_variants$CHROM <- factor(all_variants$CHROM, levels = chr_levels)
all_variants <- all_variants %>% arrange(CHROM, START, END, REF, ALT)
if (nrow(all_variants)){
	all_ranges <- makeGRangesFromDataFrame(all_variants, start.field = "START", end.field = "END")
} else { all_ranges <- GRanges() }


# Inter-caller consensus of overlapping variants
# 	- Keep the most recurrently detected variant (MATCH_CALLERS)
# 	- Select the caller with higher AF to use its variant metrics and its corresponding variants
aux_consensus_vars <- all_variants %>% mutate(
	VAR_PISCES = "", VAR_MUTECT2 = "", VAR_VARDICT = "", VAR_OCTOPUS = "", VAR_TVC = "", 
	AF_PISCES = "", AF_MUTECT2 = "", AF_VARDICT = "", AF_OCTOPUS = "", AF_TVC = "", 
	DP_PISCES = "", DP_MUTECT2 = "", DP_VARDICT = "", DP_OCTOPUS = "", DP_TVC = "", 
	AD_REF_PISCES = "", AD_REF_MUTECT2 = "", AD_REF_VARDICT = "", AD_REF_OCTOPUS = "", AD_REF_TVC = "", 
	AD_ALT_PISCES = "", AD_ALT_MUTECT2 = "", AD_ALT_VARDICT = "", AD_ALT_OCTOPUS = "", AD_ALT_TVC = ""
)
if (nrow(aux_consensus_vars)){

	## Skip variants that do NOT overlap with any caller
	vars_no_overlap <- aux_consensus_vars %>% filter(OVERLAP_CALLERS == 1)
	for (caller in unique(vars_no_overlap$CALLER)){
		aux_vars <- vars_no_overlap %>% filter(CALLER == caller)
		vars_no_overlap[rownames(aux_vars), paste0("VAR_", toupper(caller))] <- aux_vars$VAR
		vars_no_overlap[rownames(aux_vars), paste0("AF_", toupper(caller))] <- aux_vars$AF
		vars_no_overlap[rownames(aux_vars), paste0("DP_", toupper(caller))] <- aux_vars$DP
		vars_no_overlap[rownames(aux_vars), paste0("AD_REF_", toupper(caller))] <- aux_vars$AD_REF
		vars_no_overlap[rownames(aux_vars), paste0("AD_ALT_", toupper(caller))] <- aux_vars$AD_ALT
		rm(aux_vars)
	}
	rm(caller)

	## Inter-caller consensus
	vars_overlap <- aux_consensus_vars %>% filter(OVERLAP_CALLERS > 1)
	vars_overlap_ranges <- all_ranges[rownames(vars_overlap)]
	aux_overlap <- vars_overlap %>% filter(FALSE)
	vars_processed <- c()
	for (nvar in seq_len(nrow(vars_overlap))){

		if (vars_overlap[nvar, "VAR"] %in% vars_processed){next}

		last_var <- nvar
		while(tryCatch(
			any(overlapsAny(vars_overlap_ranges[nvar:last_var], vars_overlap_ranges[last_var+1])),
			error=function(e){FALSE})){ last_var <- last_var+1 }  
		
		aux_vars <- vars_overlap[nvar:last_var, ]

		aux_vars_unique <- aux_vars %>% group_by(VAR) %>% 
			summarise(
				CALLER = CALLER[AF == max(AF)][1], AD_ALT = AD_ALT[AF == max(AF)][1], 
				AF = max(AF), MATCH_CALLERS = unique(MATCH_CALLERS)) %>% 
			as.data.frame()
		
		aux_consensus <- aux_vars_unique %>% filter(MATCH_CALLERS == max(aux_vars_unique$MATCH_CALLERS))
		if(nrow(aux_consensus) > 1){aux_consensus <- aux_consensus %>% filter(AF == max(aux_consensus$AF))}
		if(nrow(aux_consensus) > 1){aux_consensus <- aux_consensus %>% filter(AD_ALT == max(aux_consensus$AD_ALT))}
		if(nrow(aux_consensus) > 1){aux_consensus <- aux_consensus[1,]}

		aux_df <- aux_vars %>% filter(CALLER == aux_consensus$CALLER)
		for (caller in unique(aux_vars$CALLER)){
			aux_caller <- aux_vars %>% filter(CALLER == caller)
			aux_df[, paste0("VAR_", toupper(caller))] <- paste(aux_caller$VAR, collapse = ",")
			aux_df[, paste0("AF_", toupper(caller))] <- paste(aux_caller$AF, collapse = ",")
			aux_df[, paste0("DP_", toupper(caller))] <- paste(aux_caller$DP, collapse = ",")
			aux_df[, paste0("AD_REF_", toupper(caller))] <- paste(aux_caller$AD_REF, collapse = ",")
			aux_df[, paste0("AD_ALT_", toupper(caller))] <- paste(aux_caller$AD_ALT, collapse = ",")
			rm(aux_caller)
		}
		rm(caller)

		vars_processed <- c(vars_processed, aux_vars_unique$VAR)
		aux_overlap <- rbind(aux_overlap, aux_df)

		rm(last_var, aux_vars, aux_vars_unique, aux_consensus, aux_df)

	}
	rm(nvar, vars_processed)

	consensus_variants <- rbind(aux_overlap, vars_no_overlap)
	consensus_variants[consensus_variants == ""] <- "."
	consensus_ranges <- makeGRangesFromDataFrame(consensus_variants, start.field = "START", end.field = "END")

} else { 
	consensus_variants <- aux_consensus_vars
	consensus_ranges <- GRanges()
}


# CALLERS column
# LowCallers column: OVERLAP_CALLERS < minCallers
# Add CALLING flags:
#	- PASS: No flags
#	- LowCallers, LowAD (AD_ALT < minAD), LowDP (DP < minDP), LowVAF (AF < minVAF), Blacklist
#	- Variants in the panel hotspots regions have specific thresholds and are NOT blacklisted
consensus_variants <- consensus_variants %>% 
	mutate(
		PANEL_HOTSPOT = overlapsAny(consensus_ranges, panel_hotspots_ranges),
		PANEL_BLACKLIST = overlapsAny(consensus_ranges, panel_blacklist_ranges),
		LowCallers = OVERLAP_CALLERS < min_callers,
		LowAD = (AD_ALT < min_ad & !PANEL_HOTSPOT) | (AD_ALT < min_ad_hotspot & PANEL_HOTSPOT), 
		LowDP = (DP < min_dp & !PANEL_HOTSPOT) | (DP < min_dp_hotspot & PANEL_HOTSPOT), 
		LowVAF = (AF < min_vaf & !PANEL_HOTSPOT) | (AF < min_vaf_hotspot & PANEL_HOTSPOT),
		PASS = !(LowCallers | LowAD | LowDP | LowVAF | (PANEL_BLACKLIST & !PANEL_HOTSPOT))) %>%
	rowwise() %>%
	mutate(
		CALLERS = paste(caller_flag[c(VAR_MUTECT2 != ".", VAR_OCTOPUS != ".", VAR_PISCES != ".", VAR_TVC != ".", VAR_VARDICT != ".")], collapse = ";"),
		FILTER = paste(filter_flags[c(PASS, LowCallers, LowAD, LowDP, LowVAF, (PANEL_BLACKLIST & !PANEL_HOTSPOT))], collapse = ";")) %>%
	as.data.frame() %>%
	arrange(CHROM, START, END, REF, ALT)


# Generate a VCF file
aux_header <- readLines(vcfHeader)
aux_header[grep("##fileDate", aux_header)] <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
aux_header[grep("##reference", aux_header)] <- paste0("##reference=file:", genome)
if(nrow(consensus_variants)){
	aux_vars <- consensus_variants %>% mutate(
		AD = paste0(AD_REF, ",", AD_ALT),
		CallMUT = VAR_MUTECT2 != ".", CallOCT = VAR_OCTOPUS != ".", CallPIS = VAR_PISCES != ".", 
		CallTVC = VAR_TVC != ".", CallVAR = VAR_VARDICT != ".", GT = ifelse(AD_REF != 0, "0/1", "1/1"),
		VC = ifelse(TYPE == "SNV", "SNV", ifelse(nchar(REF) > nchar(ALT), "DEL", ifelse(nchar(REF) < nchar(ALT), "INS", "COMPLEX")))
	)
	aux_info <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c(
			"VAR", "RUN", "SAMPLE", "CALLER", "TYPE", "AF", "DP", "AD_REF", "AD_ALT", "END", "MATCH_CALLERS", "OVERLAP_CALLERS", "VC", 
			"VAR_MUTECT2", "VAR_OCTOPUS", "VAR_PISCES", "VAR_VARDICT", "AF_MUTECT2", "AF_OCTOPUS", "AF_PISCES", "AF_VARDICT",
			"DP_MUTECT2", "DP_OCTOPUS", "DP_PISCES", "DP_VARDICT", "AD_REF_MUTECT2", "AD_REF_OCTOPUS", "AD_REF_PISCES", "AD_REF_VARDICT",
			"AD_ALT_MUTECT2", "AD_ALT_OCTOPUS", "AD_ALT_PISCES", "AD_ALT_VARDICT", "VAR_HG19"),
			function(x) paste0(x, "=", aux_vars[n, x]) )), collapse = ";")
	}))
	aux_flags <- paste0(aux_vars$CALLERS, ";", unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("PANEL_BLACKLIST", "PANEL_HOTSPOT"), function(x) ifelse(aux_vars[n, x], paste0(x, ";"), "")  )), collapse = "")})))
	aux_info <- paste0(aux_info, ";", aux_flags)
	aux_format <- unlist(lapply(seq_len(nrow(aux_vars)), function(n){
		paste(unlist(lapply(c("GT", "AD", "AF", "DP"), function(x) paste0(aux_vars[n, x]) )), collapse = ":")
	}))
	aux_df <- aux_vars %>% 
		mutate(
			CHROM = as.character(CHROM), POS = START, ID = ".", QUAL = ".", INFO = aux_info,
			FORMAT = paste(c("GT", "AD", "AF", "DP"), collapse = ":"),
			SAMPLE = aux_format) %>%
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

# Save the VCF file
writeLines(con = paste0(sample_id, "_SmallVariant_consensus.vcf"), text = vcf_text)

# Save the consensus table of variants
consensus_variants <- consensus_variants %>% select(
	RUN, SAMPLE, CALLER, TYPE, VAR, AF, DP, AD_REF, AD_ALT, FILTER, 
	PANEL_HOTSPOT, PANEL_BLACKLIST, OVERLAP_CALLERS, MATCH_CALLERS, CALLERS, 
	VAR_PISCES, VAR_MUTECT2, VAR_VARDICT, VAR_OCTOPUS, VAR_TVC,
	AF_PISCES, AF_MUTECT2, AF_VARDICT, AF_OCTOPUS, AF_TVC,
	DP_PISCES, DP_MUTECT2, DP_VARDICT, DP_OCTOPUS, DP_TVC,
	AD_REF_PISCES, AD_REF_MUTECT2, AD_REF_VARDICT, AD_REF_OCTOPUS, AD_REF_TVC,
	AD_ALT_PISCES, AD_ALT_MUTECT2, AD_ALT_VARDICT, AD_ALT_OCTOPUS, AD_ALT_TVC,
	PASS, LowCallers, LowAD, LowDP, LowVAF, 
	CHROM, START, END, REF, ALT, VAR_HG19, START_HG19, END_HG19
)
write.table(
	consensus_variants, file = paste0(sample_id, "_SmallVariant_consensus.txt"), 
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(
	list("Consensus" = consensus_variants), asTable = TRUE, overwrite = TRUE,
	file = paste0(run, "_", sample_id, "_SmallVariant_", format(Sys.time(), "%Y%m%d"), ".xlsx"))

# Save all detected variants
all_variants <- all_variants %>% select(
	RUN, SAMPLE, CALLER, TYPE, VAR, AF, DP, AD_REF, AD_ALT, OVERLAP_CALLERS, MATCH_CALLERS, 
	CHROM, START, END, REF, ALT, VAR_HG19, START_HG19, END_HG19
)
write.table(all_variants, file = paste0(sample_id, "_SmallVariant_all.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Save overlapping and match info of each caller
for (caller in caller_names){
	write.table(
		get(paste0(caller, "_info")), file = paste0(sample_id, "_SmallVariant_", caller, ".txt"),
		col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

