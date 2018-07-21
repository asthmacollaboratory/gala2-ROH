# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# authors:
# -- Kevin L. Keys
# -- Page C. Goddard
# -- Andy M. Zeiger
# -- Annie Li

# This script contains subroutines for plotting and analyzing GWAS results. 
# Source this script after setting the R environment, e.g.
#
#     source("set_R_environment.R")
#     source("postprocessing_routines.R")
#
# This ensures that all requisite libraries are loaded before defining the
# plotting functions contained herein.
#
# This script postprocessing_routines.R contains no executable code.
#
# It has no input or output.
# ==============================================================================

ConcatenateResults = function(input.prefix, output.path, input.suffix = "ROH.R.out.results", sort.result = FALSE) {
    # ConcatenateResults 
    #
    # This function merges the chromosome results into one data frame with probes ordered by P value.
    # The function will write the merged data frame to file and return it.
    #
    # Args:
    #	  input.prefix: the title of the output files. Should coincide with argument output.prefix from
    #         function PerformAssociationAnalysis(), e.g. ${ROH_DIRECTORY}/${FILENAME_BEFORE_CHROMOSOME}
    #     output.path: a complete file path where results will be stored.
    #     input.suffix: the suffix of the file after chromosome. Should coincide with argument
    #         suffix from function PerformAssociationAnalysis().
    #         Default: "ROH.R.out.results"
    #     sort.result: should final data.frame be sorted by p-value?
    #         Default: FALSE
    #
    # Output: One data frame with all results

    # read the results for chromosome 1 into memory
	chr = 1
	input.file.path = paste(input.prefix, chr, "ROH.R.out.results", sep = ".")
	results.df = fread(input.file.path, header = TRUE)
  
    # make a data frame for concatenating results files
    # will tack chromosome number as new leftmost column
    nROHs   = dim(results.df)[1]
    gwas.df = data.table(chr = rep(chr,nROHs), results.df)

    # repeat this process for remaining chromosomes
	for (chr in 2:22) {
        input.file.path = paste(input.prefix, chr, "ROH.R.out.results", sep = ".")
        results.df = fread(input.file.path, header = TRUE)
        nProbes    = dim(results.df)[1]
		gwas.temp  = data.table(chr = rep(chr, nProbes), results.df)
		gwas.df    = rbind(gwas.df, gwas.temp)

        # keep workspace tidy during loop
		rm(gwas.temp)
		rm(results.df)
    }

    # sort frame by p-value
    if (sort.result) {
        gwas.df = gwas.df[order(gwas.df$p),]
    }
  
    # write data frame to file
    fwrite(gwas.df, file = output.path)
    return(gwas.df)
}

CreateDiagnosticPlots = function(results.filepath, manhattan.plot.filepath, qq.plot.filepath, manhattan.plot.title = "Manhattan plot", threshold = 5e-8, highlight.SNPs = NULL, manhattan.ylims = c(0,8), color = c("black", "blue"), significance.threshold = 5e-8, suggestive.threshold = 1e-7, qq.plot.title = "QQ Plot", qq.plot.subtitle = NULL, qq.xlim = NULL, qq.ylim = NULL) {
    # CreateDiagnosticPlots
	#
	# This function merges the chromosome results into one file with probes ordered by P value
	#
	# Args:
	#	  results.filepath: path to concatenated GWAS results for one population
	#	  manhattan.plot.filepath: path where Manhattan plot will be saved
	#	  qq.plot.filepath: file name of QQ plot that will be written to working directory
	#	  qq.plot.title: title of QQ plot
    # 	  threshold: pvalue limit for labeling SNPs on the plot. SNPs with p-values greater than "threshold"
    #         are not plotted. Low values of "threshold" make plotting slow and may overlabel the plot. 
    #         Default: 5e-8 (genome-wide significance)
    #     highlight.SNPs: vector of SNP ids to highlight,
    #         e.g.  highlight.SNPs = c("rs12345", "rs90181294", "rs556782")
    #         Default: NULL (no SNPs to highlight)
    #     ylim: the Y-axis limits for the plot
    #         e.g. [c(min,max)]
    #         Default: c(0,8), which plots p-values up to 1e-8
    #     color = a vector of colors (min 2 colors)
    #         e.g. [color = c("color1", "color2", "color3")]
    #         Default: c("black", "blue")
    #     title: an informative plot title
    #         Default: "Manhattan plot"
    #     significance.threshold: the Bonferroni significance threshold.
    #         Default: 5e-8
    #     suggestive.threshold: the suggestive threshold of significance.
    #         Default: 1e-7
    #     qq.plot.subtitle: an informative subtitle for the QQ plot. Leave as NULL to record genomic lambda.
    #         Default: NULL (will actually output genomic lambda, e.g. "λ = 1.0")
    #     qq.xlim, qq.ylim: the X-axis Y-axis limits for the QQ plot
    #         e.g. [qq.xlim = c(xmin, xmax), qq.ylim = c(ymin, ymax)]
    #         Default: NULL (use default axis limits)
    #
    # Output: nothing

    # Load in GWAS results
	gwas.df = fread(results.filepath, header = TRUE)

	# Reformat dataframe to correct way for manhattan plot
	gwas.df.sub = subset(gwas.df, select = c("chr", "position", "Probe", "p"))
	colnames(gwas.df.sub) = c("CHR", "BP", "SNP", "P")

	# create Manhttan plot of ROH GWAS results
    manhattan.plot = CreateManhattanPlot(gwas.df.sub,
                        threshold = threshold,
                        highlight.SNPs = highlight.SNPs,
                        ylims = manhattan.ylims,
                        color = color,
                        title = manhattan.plot.title,
                        significance.threshold = significance.threshold,
                        suggestive.threshold = suggestive.threshold,
                        save.as = manhattan.plot.filepath)
    qq.plot = CreateQQPlot(gwas.df.sub,
                title = qq.plot.title,
                xlim = qq.xlim,
                ylim = qq.ylim,
                save.as = qq.plot.filepath)

    # DISCUSS: perhaps return the diagnostic plots?
    return()
}


ComputeSignificanceThreshold = function(input.file) {
    # ComputeSignificantThreshold
    #
	# The function uses coda to analyze the p-values from a GWAS.
    # It then computes the effective number of independent tests
    # and the corresponding significance threshold for your data
	# 
	# Args: 
	#   input.file: Your input file name. The expected format is
    #     the output from PerformAssociationAnalysis().
	#
	# Returns:
	#   A list with four entries:
    #     -- number of independent tests
    #     -- Bonferroni correction for multiple testing 
    #     -- the significance threshold (1 / num_ind_tests)
    #     -- the suggestive threshold
	
	# =================================================
	# Load input file
	# =================================================
    input = fread(input.file, fill = TRUE)

    # =================================================
    # Clean data frame
    # =================================================
    # Make copy of data frame because later the data frame will be rearranged by chromosome number and base position
    copy.input = copy(input)

    # Sort the data frame based on chromosome number and base position, rather than sorting the data by p-values
    # Because effectiveSize() adjusts for autocorrelation between p-values, so if the function is computed on a sorted p-value, you will get a really small number
    setorder(copy.input, chr, Probe)

	# =================================================
	# Calculate the significance threshold
	# =================================================
	# coda: calculates adjusted significance threshold with Bonferroni correction
	# effectiveSize() adjusts for autocorrelation (correlation between values)
	# effectiveSize() basically calculates the number of independent tests
	#   Need to add -log10 when calculating effectiveSize since the function requires large numbers not small numbers(ex:10^-6 will be converted to 6) 
    total_eff = effectiveSize(-log10(copy.input$p))
	# Output is :
	    # var 1
	    # Number
	# Since we only need the number for the threshold, as.numeric() will fix this and return just a numeric value, and this numeric value is carried throughout the rest of the code
    total_eff = as.numeric(total_eff)

    # Bonferroni correction: divide alpha level(0.05) by the number of independent test
    # <<- makes p.adj variable global, so that it can be called in other scripts
    p.adj <<- 0.05/total_eff

    # After -log10 transformation
    transformation = -log10(p.adj)

    # Compute the suggestive threshold
    suggestive <<- 1/(2*total_eff)

    # List out number of independent tests, the adjusted threshold, and the threshold after -log10 transformation
    return(list("n.indep.test" = total_eff, "bon.corr" = p.adj, "transformation" = transformation, "sugg.thresh" = suggestive))
}
