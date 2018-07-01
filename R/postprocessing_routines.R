ConcatenateResults = function(input.prefix, output.path, input.suffix = "ROH.R.out.results") {
    # results_processing
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
	  gwas.df = gwas.df[order(gwas.df$p),]
  
    # write data frame to file
	  fwrite(gwas.df, file = output.path)
    return(gwas.df)
}

CreateDiagnosticPlots = function(results.filepath, manhattan.plot.filepath, qq.plot.filepath, manhattan.plot.title = "Manhattan plot", threshold = 1e-2, highlight.SNPs = NULL, manhattan.ylims = c(0,8), color = c("black", "blue"), significance.threshold = 5e-8, suggestive.threshold = 1e-7, qq.plot.title = "QQ Plot", qq.plot.subtitle = NULL, qq.xlim = NULL, qq.ylim = NULL) {
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
    #         are not plotted. High values of "threshold" make plotting speedy but make plots look strange.
    #         Default: 1e-2
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
	  gwas.df$BP = gwas.df$Probe
	  gwas.df.sub = subset(gwas.df, select = c("chr", "Probe", "p", "BP"))
	  colnames(gwas.df.sub) = c("CHR", "SNP", "P", "BP")

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