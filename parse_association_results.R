#!/usr/bin/env Rscript --vanilla

# ==============================================================================
# Copyright 2019, Asthma Collaboratory
# authors:
# -- Kevin L. Keys
# -- Jennifer Elhawary
# -- Andrew M. Zeiger
# -- Annie Li
# -- Oona Risse-Adams 
# ==============================================================================


#==============================================================================
# parse options
#==============================================================================

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

option_list = list(
    make_option(
        c("-a", "--summary-stats-file"),
        type    = "character",
        default = NULL,
        help    = "File of summary statistics from ROH association analysis.",
        metavar = "character"
    ),  
    make_option(
        c("-b", "--results-file"),
        type    = "character",
        default = NULL,
        help    = "File of ROH association results. This should contain results for all chromosomes.",
        metavar = "character"
    ),  
    make_option(
        c("-c", "--significant-SNPs-file"),
        type    = "character",
        default = NULL, 
        help    = "Output file for significant SNPs.",
        metavar = "character"
    ),  
    make_option(
        c("-d", "--suggestive-SNPs-file"),
        type    = "character",
        default = NULL, 
        help    = "Output file for _suggestive_ SNPs.",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n")
print(opt)

summary.stats  = opt$summary_stats_file
results.file   = opt$results_file
sig.SNPs.file  = opt$significant_SNPs_file
sugg.SNPs.file = opt$suggestive_SNPs_file

# ==============================================================================
# read results 
# ==============================================================================
sumstats    = fread(summary.stats, sep = "\t")
roh.results = fread(results.file)


# ==============================================================================
# append significant SNPs in sig.SNPs.file 
# ==============================================================================
sig.SNPs = roh.results %>%
    dplyr::filter(p < sumstats$Significance_Threshold) %>%
    as.data.table
fwrite(sig.SNPs, file = sig.SNPs.file, append = TRUE)


# ==============================================================================
# append suggestive SNPs in sugg.SNPs.file 
# ==============================================================================
sugg.SNPs = roh.results %>%
    dplyr::filter(
        p >= sumstats$Significance_Threshold &
        p <  sumstats$Suggestive_Threshold
    ) %>%
    as.data.table
fwrite(sugg.SNPs, file = sugg.SNPs.file, append = TRUE)
