#!/usr/bin/env Rscript --vanilla

# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by Kevin L. Keys, Jennifer Elhawary, Andrew Zeiger, Annie Li, Oona Risse-Adams 
# ==============================================================================


#==============================================================================
# parse options
#==============================================================================

suppressMessages(library(optparse))
suppressMessages(library(stringr))

option_list = list(
    make_option(
        c("-a", "--phenotype-name"),
        type    = "character",
        default = NULL,
        help    = "Name of phenotype to analyze, as character string.",
        metavar = "character"
    ),  
    make_option(
        c("-b", "--covariates"),
        type    = "character",
        default = NULL,
        help    = "Comma-separated list of covariates to include in association analysis.",
        metavar = "character"
    ),  
    make_option(
        c("-c", "--source-code-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where project analysis function code is stored.",
        metavar = "character"
    ),  
    make_option(
        c("-d", "--population"),
        type    = "character",
        default = NULL, 
        help    = "Codename (e.g. 'PR' or 'MX') of current population to analyze.",
        metavar = "character"
    ),  
    make_option(
        c("-e", "--GALA-ROH-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where GALA ROH calls are stored.",
        metavar = "character"
    ),  
    make_option(
        c("-f", "--SAGE-ROH-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where SAGE ROH calls are stored.", 
        metavar = "character"
    ),  
    make_option(
        c("-g", "--scratch-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory with available space for temporary storage.",
        metavar = "character"
    ),  
    make_option(
        c("-i", "--output-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where output will be saved.",
        metavar = "character"
    ),   
    make_option(
        c("-j", "--GLM-type"),
        type    = "character",
        default = "gaussian", 
        help    = "Name of link function (e.g. 'gaussian') for GLM regression.", 
        metavar = "character"
    ),
    make_option(
        c("-k", "--num-cores"),
        type    = "integer",
        default = 10, 
        help    = "Number of CPU cores to use for parallel calculations.", 
        metavar = "integer"
    ),
    make_option(
        c("-l", "--minimum-samples-at-probe"),
        type    = "integer",
        default = 8, 
        help    = "Minimum number of ROH segments for analysis of ROH probes.", 
        metavar = "integer"
    ),
    make_option(
        c("-m", "--R-library-path"),
        type    = "character",
        default = NULL, 
        help    = "Path to project R library.", 
        metavar = "character"
    ),
    make_option(
        c("-n", "--output-suffix"),
        type    = "character",
        default = "ROH.R.out", 
        help    = "Suffix for association output files.", 
        metavar = "character"
    ),
    make_option(
        c("-o", "--phenotype-file"),
        type    = "character",
        default = "NULL", 
        help    = "Path to unified GALA/SAGE phenotype file for ROH analyses.", 
        metavar = "character"
    ),
    make_option(
        c("-p", "--input-prefix"),
        type    = "character",
        default = "NULL", 
        help    = "Prefix to file path where ROH segments are stored (split by chromosome).", 
        metavar = "character"
    ),
    make_option(
        c("-q", "--plot-type"),
        type    = "character",
        default = "png", 
        help    = "File extension for plots, e.g. 'png' or 'pdf'. (Default:default)", 
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n")
print(opt)

pheno.name           = opt$phenotype_name
covariate.list       = opt$covariates
source.dir           = opt$source_code_directory 
pop.code             = opt$population
gala.roh.dir         = opt$GALA_ROH_directory 
sage.roh.dir         = opt$SAGE_ROH_directory 
scratch.dir          = opt$scratch_directory
out.dir              = opt$output_directory
glm.type             = opt$GLM_type 
ncores               = as.integer(opt$num_cores)
min.samples.at.probe = as.integer(opt$minimum_samples_at_probe)
library.path         = opt$R_library_path 
out.suffix           = opt$output_suffix 
phenotype.filepath   = opt$phenotype_file
input.prefix         = opt$input_prefix
plot.type            = opt$plot_type

# use for tagging saved results
# this prints in YYYY-MM-DD format
the.date = Sys.Date()


# ==============================================================================
# source analysis functions from external codebase 
# ==============================================================================

# add project library path to .libPaths
.libPaths(c(library.path, .libPaths()))

# paths to R scripts required for analyses
postprocessing.routines.path = file.path(source.dir, "postprocessing_routines.R")
analysis.routines.path       = file.path(source.dir, "analysis_routines.R")
plotting.routines.path       = file.path(source.dir, "plotting_routines.R")
environment.path             = file.path(source.dir, "set_R_environment.R")

# set_R_environment.R loads the packages required for analyses
# the key is that the package libraries are sourced from a common project folder
# if the package is not yet installed, then this script will attempt to install it
source(environment.path)

# analysis_routine.R contains the analysis functions for association testing
# it merges the ROH probes with the phenotype data and then fits a generalized linear model to the data
source(analysis.routines.path)

# analysis_routines.R creates all diagnostic and association plots 
source(plotting.routines.path)

# postprocessing_routines.R merges the chromosome results into one data frame with the probes ordered by p-value
source(postprocessing.routines.path)
	

# ==============================================================================
# script variables, file paths, directories
# ==============================================================================

# Output file prefix has phenotype, population, and date
# ex: ${outdir}/results/${PHENOTYPE}.${POP}.${YEAR}-${MONTH}-${DAY}
out.pfx = file.path(out.dir, "results", str_c(pop.code, pheno.name, the.date, sep = "."))

# output file for *concatenated* results will have similar file path`
out.all.chr = str_c(out.pfx, "ALLCHR.txt", sep = ".")

cat("Setting variables complete.\n")


# ==============================================================================
# Load data
# ==============================================================================

cat("Loading data for population ", pop.code, "\n")

# load unified phenotype data table from file
phenotype.df = fread(phenotype.filepath, header = TRUE)

# subset phenotype.df to current pop
phenotype.df = phenotype.df>%
    dplyr::filter(Pop_Code == pop.code)

fwrite(phenotype.df, file = paste0(pop.code, ".txt"))

print("Loading data complete")


# ==============================================================================
# Specify model
# ==============================================================================

# This formula is VERY IMPORTANT
# In conjunction with "type", it specifies the linear model to be fitted
# For complicated reasons, the form should always be
#     outcome ~ one.roh.snp + ...(covariates)...
# Ensure that model.formula points to covariates in your phenotype data frame!
# Errors arising from typos or missing covariates can create cryptic errors
# The following code creates a formula of the following form:
# model.formula = "outcome ~ one.roh.snp + var1 + var2 + var3 + var4 + var5 + var6 + var7"
covariate.formula = str_c(unlist(str_split(covariate.list, ",")), collapse = " + ")
formula.rhs       = str_c("one.roh.snp", covariate.formula, sep = " + ")
model.formula     = str_c(pheno.name, formula.rhs, sep = " ~ ") 

cat("Model formula specified:\n")
print(model.formula)
cat("\n")


# ==============================================================================
# run association analysis
# ==============================================================================

cat("Running analysis for population ", pop.code, "\n")

# this runs association for each chromosome over one population
PerformAssociationAnalysis(
    input.prefix = input.prefix,
    output.prefix = out.pfx,
    phenotype.df = phenotype.df,
    model.formula = model.formula,
    suffix = out.suffix,
    type = glm.type,
    ncores = ncores,
    min.samples.at.probe = min.samples.at.probe,
    library.path = library.path
)

cat("Perform association analysis complete.\n")


# ==============================================================================
# post-processing and visualization
# ==============================================================================

# combine results from PerformAssociationAnalysis() into one dataframe
# output is data frame for later use for plotting 
# note: output of PerformAssociationAnalysis (which ROH probes merged with the phenotype data) becomes INPUT here
# note also: - (hyphen) is not allowed when naming variables
results.df = ConcatenateResults(out.pfx, out.all.chr)

cat("Concatenate Results complete.\n")


# ==============================================================================
# compute the significance threshold
# ==============================================================================

# use 'coda' to compute the effective number of tests
# use that number to determine a Bonferroni threshold
# this function returns a list with the following information:
# -- the number of independent tests
# -- the Bonferroni correction (which will be inserted as the threshold for the Manhattan plot)
# -- the adjusted significance threshold 
sig.info = ComputeSignificanceThreshold(out.all.chr)

# pull the significance threshold 
significance.threshold = sig.info$bon.corr

# the suggestive threshold is 1/(2* the number of independent tests), or 10 * the significance threshold
suggestive.threshold = sig.info$sugg.thresh

# setting `threshold` to the significance threshold means that only significant SNPs are labled
threshold = significance.threshold

cat("Significance threshold: ", significance.threshold, "\n")
cat("Suggestive threshold: ", suggestive.threshold, "\n")
cat("Computing significance threshold complete.\n")


# ==============================================================================
# create diagnostic plots
# ==============================================================================

# Define variables for the Manhattan and QQ plot input
# Define the output file path and the titles for each plot
manhattan.plot.filepath = file.path(out.dir, "figures", paste(pop.code, pheno.name, the.date, "manhattan", plot.type, sep = "."))
qq.plot.filepath = file.path(out.dir, "figures", paste(pop.code, pheno.name, the.date, "qq", plot.type, sep = "."))

manhattan.plot.title = paste("Manhattan Plot for", pheno.name, "in", pop.code, sep = " ")
qq.plot.title = paste("QQ Plot for", pheno.name, "in", pop.code, sep = " ")

# plot both Manhattan plots and QQ plots
# the function also calculates the genomic inflation lambda as well
## NOTE - you may need to change the axis limits of the QQ plots depending on the significance of your GWAS results
CreateDiagnosticPlots(out.all.chr,
    manhattan.plot.filepath,
    qq.plot.filepath,
    manhattan.plot.title = manhattan.plot.title,
    qq.plot.title = qq.plot.title,
    threshold = threshold,
    significance.threshold = significance.threshold,
    suggestive.threshold = suggestive.threshold
)

cat("Create diagnostic plots complete.\n")

# save a file for summary statistics
data.for.pheno = na.omit(subset(phenotype.df, select=c(pheno.name, unlist(str_split(covariate.list, ",")))))
reg.type = ifelse(glm.type == "gaussian", "linear", "logistic")
N = dim(data.for.pheno)[1]
if ( reg.type == "logistic") {
    category0_N = sum(data.for.pheno[, pheno.name] == 0, na.rm = TRUE)
    category1_N = sum(data.for.pheno[, pheno.name] == 1, na.rm = TRUE)
    my.mean    = NA
    my.median  = NA
    my.25thpc  = NA
    my.75thpc  = NA
    my.range   = NA 
} else {
    category0_N = NA
    category1_N = NA
    my.summary = summary(data.for.pheno[[pheno.name]], na.rm = TRUE)
    my.mean    = as.numeric(my.summary["Mean"])
    my.median  = as.numeric(my.summary["Median"])
    my.25thpc  = as.numeric(my.summary["1st Qu."])
    my.75thpc  = as.numeric(my.summary["3rd Qu."])
    my.range   = paste(my.summary["Min."], my.summary["Max."], sep = " - ")
}


summary.stats = data.table(
    "Phenotype" = pheno.name,
    "Regression_Type" = reg.type,
    "Population" = pop.code,
    "Sample_Size" = N,
    "Category0_N" = category0_N,
    "Category1_N" = category1_N,
    "Mean" = my.mean,
    "Median" = my.median,
    "25th_Percentile" = my.25thpc,
    "75th_Percentile" = my.75thpc,
    "Range" = my.range,
    "Covariates" = covariate.list,
    "Significance_Threshold" = significance.threshold,
    "Suggestive_Threshold" = suggestive.threshold
)

# save table to file
summary.filepath = file.path(out.dir, "results", paste(pop.code, "summarystats", "txt", sep = ".")) 
fwrite(x = summary.stats, file = summary.filepath, sep = "\t", quote = FALSE, na = "NA")

cat("Summary stats written to file\n")
