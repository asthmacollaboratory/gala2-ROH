# =====================================================================
# This scripts contains the RohPlot() function, which plots the ROH data for the same segment you looked at in Locus Zoom
#   Note that Zachary Szpiech's PERL script was used to extract the segment of data used for the local ancestry plot
# Assistance provided by Zachary Szpiech, Kevin L. Keys, Jennifer Elhawary, Andy Zeiger, and Oona Risse-Adams
# Written by Annie Li
# =====================================================================

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
    #       *** MAKE SURE YOU END THE FILE NAME WITH .png
    #
    # Output: A plot displaying the ROH data for a specific phenotype and population; the plot has Probes on the x-axis and the phenotype on the y-axis, each horizontal line represents the ROH for an individual in the population

    # =================================================================
    # Load packages needed
    # =================================================================
    
    library(ggplot2)
    library(reshape2)
    
    # =================================================================
    # Load data
    # =================================================================
    
    # Load phenotype file
    pheno = read.table(input.pheno, header = T, sep = "\t")

    print(paste("Finished loading", input.pheno, sep = " "))

    # Load dosage data on ROH
    dosage = read.table(input.dosage, header = T)   

    print(paste("Finished loading", input.dosage, sep = " ")) 

    # =================================================================
    # Subset, transpose, and melt data frame
    # =================================================================
    
    # Subset dosage data for only the chromosome segment of interset (this segment/section should match the segment/section used in your locus zoom plot
    dosage.segment = subset(dosage, pos > probe.min & pos < probe.max)

    print("Finished subsetting dosage data")
    
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
    
    print("Melted dosage data")

    # Change the column names of the melted data, so that Var1 will be "SubjectID", Var2 will be "Probe", and value will be "ROH"
    colnames(dosage.melt) = c("SubjectID", "Probe", "ROH")

    # Subset the phenotype file for only the SubjectID and phenotype data (the covariates and principal components are not needed for plotting ROH)
    pheno.subset = subset(pheno, select = c("SubjectID", phenotype))

    print("Finished subsetting phenotype data")

    # Merge the dosage data with the phenotype data by the column SubjectID
    merge.df = merge(dosage.melt, pheno.subset, by = "SubjectID", all.y = T)

    print("Finished merging dosage data with phenotype data")

    # Order the data by probe position and phenotype 
    merge.df = setorderv(merge.df, cols = c(phenotype), order = 1)
    
    print("Ready to plot")

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
}
