# ==============================================
# This script plots a Kernel density plot over a histogram for a continuous phenotype
#   The plot can display the distribution of the phenotype and can help you decide if you may need to transform your data set
# Written by Annie Li
# July 27, 2018
# ===============================================

DensityAndHistogramPlot = function(input.pheno, output.file, phenotype, title, subtitle = NULL, hist.binwidth = 1, hist.border.color = "blue1", hist.fill.color = "blue", alpha = 0.6, line.color = "blue", size = 3) {
    # DensityAndHistogramPlot()
    # 
    # This function plots a Kernel density plot that overlays a histogram to show the distribution of a continuous variable
    # 
    # Args: 
    #     input.pheno: The input phenotype file containing the values for the continous phenotype
    #     output.pheno: The output file path containing the full path to the output destination and the name of the saved plot (Make sure it end in .png)
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
    # PLEASE PUT HAVE QUOTES "" AROUND ALL YOU INPUT NAMES
    # 
    # Output: Kernel density plot overlaying a histogram for your continuous phenotype


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
}
