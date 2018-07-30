# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by Kevin L. Keys, Pagé C. Goddard, Andrew Zeiger, Annie Li, Oona Risse-Adams
#
# This script contains several subroutines to facilitate plotting.
# It relies heavily on ggplot2 and other tidyverse packages.
# Source this script after setting the R environment, e.g.
#
#     source("set_R_environment.R")
#     source("plotting_routines.R")
#
# This ensures that all requisite libraries are loaded before defining the
# plotting functions contained herein.
#
# This script plotting_routines.R contains no executable code.
#
# It has no input or output.
# ==============================================================================

# ==============================================================================
# define functions here
# ==============================================================================

CreateManhattanPlot = function(df, threshold = 5e-8, highlight.SNPs = NULL, ylims = c(0,8), color = c("black", "blue"),
    title = "Manhattan plot", significance.threshold = 5e-8, suggestive.threshold = 1e-7, save.as = NULL){
    # Create a Manhattan Plot
    #
    # This function creates a Manhattan plot. It expects a data frame with the following four labeled columns:
    #     CHR: the chromosome to plot
    #     SNP: the single nucleotide polymorphisms (one per row)
    #     BP:  the base pair (position) of each SNP to plot
    #     P:   the p-value from the association test, one per SNP
    #
    # Args:
    #     df: data frame with colnames [SNP, CHR, BP, P]
    # 	  threshold: pvalue limit for labeling SNPs on the plot. SNPs with p-values greater than "threshold"
    #         are not plotted. BEWARE: high values of "threshold" can potentially lead to many SNPs
    #         being labeled and thereby make make plots unreadable!
    #         Default: 5e-8 (label any SNP with "standard" Bonferroni-corrected genome-wide significance)
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
    #         Default: 5e-8 ("standard" genome-wide significance)
    #     suggestive.threshold: the suggestive threshold of significance.
    #         Default: 1e-7
    #     save.as: a filepath for saving the plot to file
    #         Default: NULL (do not save to file)
    # Outputs:
    #    g is a ggplot object containing the Manhattan plot

	# format df with (complicated) dplyr filtering
	df.tmp = df %>%

		# Compute chromosome size
		group_by(CHR) %>%
		summarise(chr_len = max(BP)) %>%

		# Calculate cumulative position of each chromosome
		mutate(tot = cumsum(chr_len) - chr_len) %>%
		select(-chr_len) %>%

		# Add this info to the initial dataset
		left_join(df, ., by = c("CHR" = "CHR")) %>%

		# Add a cumulative position of each SNP
		arrange(CHR, BP) %>%
		mutate(BPcum = BP + tot) %>%

		# Add highlight and annotation information
		mutate(is_highlight = ifelse(SNP %in% highlight.SNPs, "yes", "no")) %>%
		mutate(is_annotate = ifelse(P < threshold, "yes", "no"))  ### done filtering!

	# get chromosome center positions for x-axis
	axisdf = df.tmp %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2 )

	# plot with filtered data frame
	# we will construct this ggplot stepwise
	g = ggplot(df.tmp, aes(x = BPcum, y = -log10(P)))

	# Show all points
	g = g + geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
	scale_color_manual(values = rep(color, 22))

	# custom X axis
	# note: expand = c(0, 0) removes space between plot area and x axis
	g = g + scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
	scale_y_continuous(expand = c(0, 0), limits = ylims)

	# add plot and axis titles
	g = g + ggtitle(paste0(title)) + labs(x = "Chromosome")

	# add genome-wide significant.threshold and suggestive.threshold lines
	g = g + geom_hline(yintercept = -log10(significance.threshold), color = "red") +
	geom_hline(yintercept = -log10(suggestive.threshold), linetype = "dashed", color = "blue")

	# add highlighted points
	g = g + geom_point(data = subset(df.tmp, is_highlight == "yes"), color = "orange", size = 2)

	# add label using ggrepel to avoid overlapping
    df.label = df.tmp[df.tmp$is_annotate == "yes",]
    if (dim(df.label)[1] > 0) {
        g = g + geom_label_repel(
            data = df.label,
            aes(label = as.factor(SNP), alpha = 0.7),
            size = 5,
            force = 1.3
        )
    }

	# custom the theme
	g = g + theme_bw(base_size = 22) +
		theme(
			plot.title = element_text(hjust = 0.5),
			legend.position = "none",
			panel.border = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank()
	)

    # save to file?
    if (!is.null(save.as)) {
        ggsave(save.as, plot = g)
    }

	return(g)
}

CreateQQPlot = function(df, title = "QQ Plot", subtitle = NULL, xlim = NULL, ylim = NULL, save.as = NULL) {
    # Create a Quantile-Quantile Plot
    #
    # This function creates a QQ plot of GWAS p-values. It expects a data frame with the following four labeled columns:
    #     CHR: the chromosome to plot
    #     SNP: the single nucleotide polymorphisms (one per row)
    #     BP:  the base pair (position) of each SNP to plot
    #     P:   the p-value from the association test, one per SNP
    #
    #
    # Args:
    #     df: data frame with colnames [SNP, CHR, BP, P]
    #     title: an informative plot title
    #         Default: "QQ Plot"
    #     subtitle: an informative plot subtitle
    #         Default: NULL (will actually output genomic lambda, e.g. "λ = 1.0")
    #     xlim, ylim: the X-axis Y-axis limits for the plot
    #         e.g. [xlim = c(xmin, xmax), ylim = c(ymin, ymax)]
    #         Default: NULL (use default axis limits)
    #     save.as: a filepath for saving the plot to file
    #         Default: NULL (do not save to file)
    # Outputs:
    #    g is a ggplot object containing the Manhattan plot

    # could theoretically use new ggplot2 3.0.0 functions for this
    # unfortunately, they do not handle -log10 transformed p-values nicely
    # code below is appropriate for untransformed p-values
    # note: the stat_qq() layer expects an aesthetic "sample"
    #g = ggplot(df, aes(sample = P)) + geom_qq(distribution = stats::qunif)
    #g = g + geom_qq_line(distribution = stats::qunif)

    # make a transformed copy of the p-values for plotting
    df.copy = data.frame("observed" = -log10(sort(df$P)), "expected" = -log10(ppoints(length(df$P))))
    g = ggplot(df.copy, aes(x = expected, y = observed)) + geom_point()

    # draw the qq-line
    g = g + geom_abline(intercept = 0, slope = 1, color = "red")

    # prettify the axis labels
    g = g + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")

    # adjust the axis limits
    g = g + coord_cartesian(xlim = xlim, ylim = ylim)

    # compute genomic control factor
    # note that we operate on the original and untransformed p-values!
    # when adding title, also add subtitle with genomic lambda included
    qchi = qchisq(df$P, 1, lower.tail = FALSE)
    genomic.control.factor = median(qchi[!is.nan(qchi)], na.rm = TRUE) / qchisq(0.5, 1)
    g = g + labs(title = title, subtitle = bquote(lambda == .(genomic.control.factor)))

    # save to file?
    if (!is.null(save.as)) {
        ggsave(save.as, plot = g)
    }
    return(g)
}


DensityAndHistogramPlot = function(input.pheno, output.file, phenotype, title, subtitle = NULL, hist.binwidth = 1, hist.border.color = "blue1", hist.fill.color = "blue", alpha = 0.6, line.color = "blue", size = 3) {
    # DensityAndHistogramPlot()
    #
    # This function plots a Kernel density plot that overlays a histogram to show the distribution of a continuous variable
    #   The plot can display the distribution of the phenotype and can help you decide if you may need to transform your data set
    # Written by Annie Li
    #
    # Args:
    #     input.pheno: The input phenotype file containing the values for the continous phenotype
    #     output.pheno: The output file path containing the full path to the output destination and the name of the saved plot (Make sure it ends in a file extension such as .png or .pdf)
    #     phenotype: Name of column where your phenotype values are located
    #     title: Title of your plot
    #     subtitle: Subtitle of your plot (default is NULL)
    #     hist.binwidth: The binwidth size of the histogram (default is 1)
    #     hist.border.color: The color for the borders of the histogram (default color is "blue1")
    #     hist.fill.color: The color for the histogram body (fill) (default color is "blue")
    #     alpha: The transparency level of the histogram (default is 0.6)
    #     line.color: The color for the Kernel density line (default color is "blue")
    #     size: The thickness of the Kernel density line (default size is 3)
    #
    # The color names can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    #   Please use one color scheme for each phenotype
    #
    # PUT QUOTES "" AROUND ALL YOUR INPUT NAMES
    #
    # Output: a kernel density plot overlaying a histogram for your continuous phenotype


    # Load the phenotype file
    pheno = fread(input.pheno)

    # Plot the phenotype data
    # Note:
    #   ggplot(<data>, aes(x = <data on x-axis>)) +         *aes_string() if the input is a character string or contains quotes
    #       geom_histogram(aes(y = ..density..))            Specify y = ..density.. so that you can plot the Kernel density line
    pheno.plot = ggplot(pheno, aes_string(x = phenotype)) +
        geom_histogram(aes(y = ..density..), binwidth = hist.binwidth, color = hist.border.color, fill = hist.fill.color, alpha = alpha) +
        geom_density(color = line.color, size = size) +
        labs(title = title, subtitle = subtitle)

    # Save the plot
    ggsave(filename = output.file, plot = pheno.plot)

    # return the plot
    return(pheno.plot)
}


RohPlot = function(input.pheno, phenotype, input.dosage, probe.min, probe.max, size = 20, roh.title = "ROH vs. phenotype", roh.subtitle = NULL, output.file) {
    # RohPlot
    #
    # This function produces a plot displaying the ROH data for a specific phenotype and population; the plot has Probes on the x-axis and the phenotype on the y-axis, each horizontal line represents the ROH for an individual in the population
    #
    # Arguments:
    #   input.pheno: Input phenotype file produced by the phenotype script; the input file should contain the SubjectIDs and the phenotype under investigation
    #   phenotype: The column name for your phenotype in your phenotype file
    #   input.dosage: Input dosage file containing the ROH data
    #       Note that you must specify the population and chromosome of interest
    #       The file has SubjectID as rows, and probes as columns (saved as chr and pos; however, chromosome is same throughout, so probes depend on pos (position)), and for each Subject and probe there is either the value 0 or 1, meaning is the individual's probe located inside an ROH segment (0 = probe not in ROH, 1 = probe in ROH)
    #   probe.min: The minimum probe position; this should match the range you are using in the Locus Zoom plots (For example: I looked at the probes between the range of 116000000 to 117000000, so the minimum probe position is 116000000)
    #   probe.max: The maximum probe position; this should match the range you are using in the Locus Zoom plots (For example: I looked at the probes between the range of 116000000 to 117000000, so the maximum probe position is 117000000)
    #   size: The size of each point on the graph, the default in 20
    #   roh.title: The title of the ROH plot (default is "ROH vs. phenotype")
    #   roh.subtitle: The subtitle of the ROH plot (default is NULL)
    #       Note that adding a subtitle compresses the graph
    #   output.file: Output file name, put the full file path to where you want to save your plot to
    #       *** MAKE SURE YOU END THE FILE NAME WITH .png OR .pdf
    #
    # Output: A plot displaying the ROH data for a specific phenotype and population; the plot has Probes on the x-axis and the phenotype on the y-axis, each horizontal line represents the ROH for an individual in the population


    # =================================================================
    # Load data
    # =================================================================

    # Load phenotype file
    pheno = read.table(input.pheno, header = T, sep = "\t")

    # Load dosage data on ROH
    dosage = read.table(input.dosage, header = T)

    # =================================================================
    # Subset, transpose, and melt data frame
    # =================================================================

    # Subset dosage data for only the chromosome segment of interset (this segment/section should match the segment/section used in your locus zoom plot
    dosage.segment = subset(dosage, pos > probe.min & pos < probe.max)

    # Transpose the dosage data so that the columns are now probes and the SubjectIDs are now rows
    #   This allows for us to melt our data
    dosage.transpose = t(dosage.segment)

    # Since the column names are just the column number, rename the column names after the probe position, which is in the second row, this will be important when we melt our data, which looks at the column names
    colnames(dosage.transpose) = dosage.transpose[2, ]

    # Melt the data so that you will have a column for SubjectID and ROH data
    #   Essentially melt() will collapse all the columns into one column
    # An example of how your data may look after melt()
    #   Var1 is the SubjectIDs, Var2 is the probe, and value is whether the individual has the specified probe in an ROH segment
    #        Var1     Var2 value
	#    1 HR1809 45001031     0
	#    2 HR1755 45001031     0
	#    3 HR1734 45001031     0
	#    4 HR1661 45001031     0
	#    5 HR1727 45001031     0
	#    6 HR1832 45001031     0
    # Remove the chromosome and probe rows (first and second row), since we already have probes set as column names and all the probes are in the same chromosome
    dosage.melt = melt(dosage.transpose[-c(1:2), ])

    # Change the column names of the melted data, so that Var1 will be "SubjectID", Var2 will be "Probe", and value will be "ROH"
    colnames(dosage.melt) = c("SubjectID", "Probe", "ROH")

    # Subset the phenotype file for only the SubjectID and phenotype data (the covariates and principal components are not needed for plotting ROH)
    pheno.subset = subset(pheno, select = c("SubjectID", phenotype))

    # Merge the dosage data with the phenotype data by the column SubjectID
    merge.df = merge(dosage.melt, pheno.subset, by = "SubjectID", all.y = T)

    # Order the data by probe position and phenotype
    merge.df = setorderv(merge.df, cols = c(phenotype), order = 1)

    # =====================================================================
    # Plotting ROH data
    # =====================================================================

    # Set the colors for the plot to be white for 0 and black for 1
    set.colors = c("0" = "white", "1" = "black")

    # Plot the data using ggplot()
    # The probes will be on the x-axis, while the sorted phenotype will be on the y-axis
    # Note that get() was used since the phenotype input is a character string, so get() will return the value of the character string, which is the column containing the phenotype values
    # The plot will be colored based on the ROH values 0 or 1, we will later use scale_color_manual() to reset the colors to white and black
    #   Note in aes() color must be a character value, so convert the ROH values to character values
    # geom_point() plots the points
    #   shape 15 is a square, so that ROH segments (with consecutive 1s) will form a line
    #   size adjusts for size of the point (in our case it is a square)
    # labs() sets the labels of the title, subtitle, x-axis, y-axis, and legend
    roh.plot = ggplot(data = merge.df, aes(x = Probe, y = get(phenotype))) +
        geom_point(shape = 15, size = size, aes(color = as.character(ROH))) +
        scale_color_manual(values = set.colors) +
        labs(title = roh.title, subtitle = roh.subtitle, x = "Probes", y = phenotype, color = "Probe in ROH")

    # Use ggsave() to save ggplots
    ggsave(roh.plot, file = output.file)

    # return the plot
    return(roh.plot)
}
