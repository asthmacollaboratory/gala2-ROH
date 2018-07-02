# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by Kevin L. Keys, Pagé C. Goddard, Andrew Zeiger
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

    # the stat_qq() layer expects an aesthetic "sample"
    g = ggplot(df, aes(sample = P)) + geom_qq()

    # note: we draw the qq-line manually
    # in future version of ggplot2, this can be done within ggplot itself; see here
    # https://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
    g = g + geom_abline(intercept = mean(df$P), slope = sd(df$P))

    # prettify the axis labels
    g = g + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")

    # adjust the axis limits
    g = g + coord_cartesian(xlim = xlim, ylim = ylim)

    # compute genomic control factor
    # when adding title, also add subtitle with genomic lambda included
    genomic.control.factor = median(qchisq(1 - df$P, 1)) / qchisq(0.5, 1)
    g = g + labs(title = title, subtitle = bquote(lambda == .(genomic.control.factor)))

    # save to file?
    if (!is.null(save.as)) {
        ggsave(save.as, plot = g)
    }
    return(g)
}
