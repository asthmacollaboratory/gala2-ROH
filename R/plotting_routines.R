# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by Kevin L. Keys, Pagé C. Goddard 
#   
# This script contains several subroutines to facilitate plotting.
# It relies heavily on ggplot2.
# Source this script after setting the R environment, e.g.
#
#     source("set_R_environment.R")
#     source("plotting_routines.R")
#
# This ensures that all requisite libraries are loaded before defining the
# plotting functions contained herein.
# plotting_routines.R contains no executable code. 
# It has no input or output.
# ==============================================================================

# ==============================================================================
# define functions here 
# ==============================================================================

ManhattanPlot = function(df, significance.threshold = 5e-8, suggestive.threshold = 1e-6,
	threshold, hlight, ylims, col, title = "Manhattan plot"){
    # Plot a Manhattan Plot 
    # 
    # This function creates a Manhattan plot.
    #
    # Args:
    #     df: data frame with colnames [SNP, CHR, BP, P]
    #         The colnames stand for (SNP, CHRomosome, Base Pair, P-value).
    # 	  threshold: pvalue limit for labeling SNPs on the plot
    #     highlight.SNPs: vector of SNP ids to highlight,
    #         e.g.  highlight.SNPs = c("rs12345", "rs90181294", "rs556782")
    #     ylim: the Y-axis limits for the plot
    #         e.g. [c(min,max)]
    #     col = name for vector of colors (min 2 colors)
    #         e.g. [blues = c("color1", "color2", "color3")]
    #     significance.threshold: the Bonferroni significance threshold.
    #         Default: 5e-8
    #     suggestive.threshold: the suggestive threshold of significance.
    #         Default: 1e-6
    #     title: an informative plot title 
    #
    # Outputs:
    #    g is a ggplot object containing the Manhattan plot
	# TODO(Andy): test this with initial results 

	### format df with (complicated) dplyr filtering
	df.tmp = df %>% 

		# Compute chromosome size
		group_by(CHR) %>% 
		summarise(chr_len = max(BP)) %>% 

		# Calculate cumulative position of each chromosome
		mutate(tot = cumsum(chr_len) - chr_len) %>%
		select(-chr_len) %>%

		# Add this info to the initial dataset
		left_join(df, ., by = c("CHR " ="CHR")) %>%

		# Add a cumulative position of each SNP
		arrange(CHR, BP) %>%
		mutate(BPcum = BP + tot) %>%

		# Add highlight and annotation information
		mutate(is_highlight = ifelse(SNP %in% highlight.SNPs, "yes", "no")) %>%
		mutate(is_annotate = ifelse(P < threshold, "yes", "no"))  # done filtering
  
	# get chromosome center positions for x-axis
	axisdf = df.tmp %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2 )

	### plot with filtered data frame 
	# we will construct this ggplot stepwise
	g = ggplot(df.tmp, aes(x = BPcum, y = -log10(P)))

	# Show all points
	g = g + geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
	scale_color_manual(values = rep(col, 22))

	# custom X axis
	# note: expand = c(0, 0) removes space between plot area and x axis 
	g = g + scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
	scale_y_continuous(expand = c(0, 0), limits = ylims)

	# add plot and axis titles
	g = g + ggtitle(paste0(title)) +
	labs(x = "Chromosome")

	# add genome-wide significant.threshold and suggestive.threshold lines
	g = g + geom_hline(yintercept = -log10(significance.threshold)) +
	geom_hline(yintercept = -log10(suggestive.threshold), linetype="dashed")

	# Add highlighted points
	g = g + geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2)

	# Add label using ggrepel to avoid overlapping
	g = g + geom_label_repel(
		data = df.tmp[df.tmp$is_annotate == "yes",],
		aes(label = as.factor(SNP), alpha = 0.7),
		size = 5,
		force = 1.3
	)

	# Custom the theme:
	g = g + theme_bw(base_size = 22) +
		theme( 
			plot.title = element_text(hjust = 0.5),
			legend.position = "none",
			panel.border = element_blank(),
			panel.grid.major.x = element_blank(),
			panel.grid.minor.x = element_blank()
	)

	return(g)
}
