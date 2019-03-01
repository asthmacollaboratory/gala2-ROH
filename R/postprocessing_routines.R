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
	input.file.path = paste(input.prefix, chr, input.suffix, sep = ".")
	results.df = fread(input.file.path, header = TRUE)

    # make a data frame for concatenating results files
    # will tack chromosome number as new leftmost column
    nROHs   = dim(results.df)[1]
    gwas.df = data.table(chr = rep(chr,nROHs), results.df)

    # repeat this process for remaining chromosomes
	for (chr in 2:22) {
        input.file.path = paste(input.prefix, chr, input.suffix, sep = ".")
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
    p.adj = 0.05/total_eff

    # After -log10 transformation
    transformation = -log10(p.adj)

    # Compute the suggestive threshold
    suggestive = 1/(2*total_eff)

    # List out number of independent tests, the adjusted threshold, and the threshold after -log10 transformation
    return(list("n.indep.test" = total_eff, "bon.corr" = p.adj, "transformation" = transformation, "sugg.thresh" = suggestive))
}

RohLancPlot = function(input.pheno, phenotype, input.dosage, input.ancestry, probe.min, probe.max, colorblind.friendly = FALSE, size = 0.05, sig.probe, roh.title = "ROH Plot", roh.subtitle = NULL, ancestry.title = "Local Ancestry Plot", ancestry.subtitle = NULL, sum.title = "Local Ancestry Summary Plot", sum.subtitle = NULL, roh.lanc.title = "ROH Local Ancestry", roh.lanc.subtitle = NULL, output.dir, output.pfx) {
    # RohLancPlot
    #
    # This function cleans the data then uses the functions RohPlot() and LancPlot() to plot the ROH data and the local ancestry data
    # 
    # Args:
    #     input.pheno: Input file containing phenotype data
    #     phenotype:  The column name that contains the phenotype values
    #     input.dosage:  Input file containing the ROH and probe data
    #     input.ancestry: Input file containing the local ancestry data of a specific region
    #                   The input file is created using the lanc_overlap.pl script created by Zachary Szpiech
    #     probe.min: The minimum probe position; this should match the range you are using in the Locus Zoom plots (For example: I looked at the probes between the range of 116000000 to 117000000, so the minimum probe position is 116000000)
    #     probe.max: The maximum probe position; this should match the range you are using in the Locus Zoom plots (For example: I looked at the probes between the range of 116000000 to 117000000, so the maximum probe position is 117000000)
    #     colorblind.friendly: If TRUE, a color blind friendly palette is used instead of the normal red, blue, yellow scale for the ROH local ancestry plot (default is FALSE)
    #     size: The size of each point on the graph (default is 0.05)
    #     sig.probe: The marker for the probe of interest. A vertical line will be plotted at the probe of interest
    #     roh.title: The title of the ROH plot (default is "ROH Plot")
    #     roh.subtitle: The subtitle of the ROH plot (default is NULL)
    #     ancestry.title: The title of the local ancestry plot (default is "Local Ancestry Plot")
    #                   Note that adding a subtitle compresses the graph
    #     ancestry.subtitle: The subtitle of the local ancestry plot (default is NULL)
    #                   Note that adding a subtitle compresses the graph
    #     sum.title: The title of the summary local ancestry plot (default is "Local Ancestry Summary Plot")
    #     sum.subtitle: The subtitle of the summary local ancestry plot (default is NULL)
    #                   Note that adding a subtitle compresses the graph
    #     roh.lanc.title: The title of the ROH plot colored by ancestry (default is "ROH Local Ancestry Plot")
    #     roh.lanc.subtitle: The subtitle of the ROH plot colored by ancestry (default is NULL)
    #                   Note that adding a subtitle compresses the graph
    #     output.dir: The output directory
    #     output.pfx: The prefix of the saved plot names
    #                 The following suffix will be added:
    #                   ".lanc.png" for the local ancestry plot
    #                   ".sum.lanc.png" for the local ancestry summary plot
    #                   ".roh.lanc.png" for the ROH local ancestry plot
    # Output: 3 plot displaying the local ancestry data for a specific phenotype, population, and region on a chromosome.
    #         All plots has probes on the x-axis and the phenotype on the y-axis


    # ====================================================================================
    # Load the data
    # ====================================================================================
    
    # Load the phenotype data
    pheno =  fread(input.pheno)

    # Load the ROH dosage data
    dosage = read.table(input.dosage, header = T)

    # Load the local ancestry data
    local.ancestry = read.table(input.ancestry, header = T)
    
    # ====================================================================================
    # Clean phenotype and ROH dosage data
    # ====================================================================================
    
    # For the phenotype data, subset for the SubjectIDs and the phenotype columns
    pheno.subset = subset(pheno, select = c("SubjectID", phenotype))

    # For the dosage data, subset for the region of interest (this segment should match the segment used in your locus zoom plot
    dosage.segment = subset(dosage, pos >= probe.min & pos <= probe.max)

    # Transpose the dosage data so that the columns are now probes and the SubjectIDs are now rows
    #   This allows for us to melt our data 
    dosage.transpose = t(dosage.segment)

    # Since the column names are just the column number, rename the column names after the probe position, which is in the second row
    # This will be important when we melt our data, which looks at the column names
    colnames(dosage.transpose) = dosage.transpose[2, ]

    # Melt the data so that you will have a column for SubjectID and ROH data
    #   Essentially melt() will collapse all the columns into one column
    # An example of how your data may look after melt()
    #   Var1 is the SubjectIDs, Var2 is the probe, and value is whether the individual has the specified probe in an ROH segment
    #        Var1     Var2 value
	#    1 A      45001031     0
	#    2 B      45001031     0
	#    3 C      45001031     0
	#    4 D      45001031     0
	#    5 E      45001031     0
	#    6 F      45001031     0
    # Remove the chromosome and probe rows (first and second row), since we already have probes set as column names and all the probes are in the same chromosome
    dosage.melt = melt(dosage.transpose[-c(1:2), ])
    
    print("Melted dosage data")

    # Change the column names of the melted data, so that Var1 will be "SubjectID", Var2 will be "Probe", and value will be "ROH"
    colnames(dosage.melt) = c("SubjectID", "Probe", "ROH")

    # Subset the phenotype file for only the SubjectID and phenotype data (the covariates and principal components are not needed for plotting ROH)
    pheno.subset = subset(pheno, select = c("SubjectID", phenotype))

    print("Finished subsetting phenotype data")

    # Merge the dosage data with the phenotype data by the column SubjectID
    #   Set all.x = T, so that all the phenotype data/values will be there after the merge
    merge.roh.pheno.df = merge(pheno.subset, dosage.melt, by = "SubjectID", all.x = T)

    print("Finished merging ROH data and phenotype data")

    # First we will check how many rows have missing ROH data
    print("Number of rows with missing ROH data")
    print(sum(is.na(merge.roh.pheno.df$ROH)))

    # Check how many individuals have missing ROH data
    print("Number of individuals with missing ROH data")
    print(length(unique(merge.roh.pheno.df$SubjectID[is.na(merge.roh.pheno.df$ROH)])))

    # Check which individuals have missing ROH data
    print("Individuals with missing ROH data")
    print(unique(merge.roh.pheno.df$SubjectID[is.na(merge.roh.pheno.df$ROH)]))

    # Print all the rows of individuals with missing ROH data
    print("Rows with missing ROH data")
    print(merge.roh.pheno.df[is.na(merge.roh.pheno.df$ROH), ])
    print(sum(is.na(merge.roh.pheno.df$ROH)))

    # Remove all the rows where ROH has NAs 
    #   The phenotype data has individuals not contained in the local ancestry data
    merge.roh.pheno.df = merge.roh.pheno.df[!is.na(merge.roh.pheno.df$ROH), ]

    # Check all the missing ROH data are removed
    print("Number of rows with missing ROH data")
    print(sum(is.na(merge.roh.pheno.df$ROH)))

    # Reorder the data frame
    merge.roh.pheno.df = setorderv(merge.roh.pheno.df, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))

    # Since there are missing probe data, the plots contain vertical gaps
    # To solve this problem, we will assign each unique probe a number
    # First we will create a new data frame with one column containing all the unique probes, and another column containing the unique number corresponding with the probe
    roh.probe.df = data.table(x = unique(merge.roh.pheno.df$Probe), y = order(unique(sort(merge.roh.pheno.df$Probe))))
    
    # Rename the columns "Probe" and "probe.order"
    colnames(roh.probe.df) = c("Probe", "probe.order")

    # Merge the ROH data with the probe order
    roh.probe.merge = merge(merge.roh.pheno.df, roh.probe.df, by = "Probe")

    # Create a new data frame that consists of the SubjectID and phenotype columns
    roh.pheno.df = subset(merge.roh.pheno.df, select = c("SubjectID", phenotype))    
    
    # Reorder the data frame to ensure that it is ordered by phenotype and SubjectID
    roh.pheno.df = setorderv(roh.pheno.df, cols = c(phenotype, "SubjectID"), order = c(1, 1))
    
    # Create a new data frame which contains the unique SubjectID and phenotype rows
    roh.unique.rows.df = unique(roh.pheno.df)

    # Create a new column in the data frame pheno.order, which contains the unique order number for each individual and their phenotype value
    roh.unique.rows.df$pheno.order = order(roh.unique.rows.df[, get(phenotype)])
    
    # Merge the ROH data ,with the ordered probe, with the newly created phenotype order numbers
    roh.final.merge = merge(roh.probe.merge, roh.unique.rows.df, by = c("SubjectID", phenotype))
    
    # Reorder the data frame
    roh.final.merge = setorderv(roh.final.merge, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))
    
    print("Finished creating creating probe and phenotype order for ROH data, now ready to plot ROH")

    # ===============================================================================
    # Plot ROH data
    # ===============================================================================

    # We will place a vertical line at the probe of interest, so that we will be able to identify the ROH data at that probe
    # If the probe of interest is not in the data set, the nearest probe will be used instead
    significant.probe.index = which.min(abs(roh.final.merge$Probe - sig.probe)) 
    significant.probe = roh.final.merge[significant.probe.index, ]$Probe
    significant.probe.order = roh.final.merge[significant.probe.index, ]$probe.order
    
    print(paste0("Vertical line plotted at:", significant.probe))

    # Plot the ROH data using the RohPlot() function, which is defined in the plotting_routines.R script
    RohPlot(input.df = roh.final.merge, 
            phenotype = phenotype, 
            size = size, 
            roh.title = roh.title, 
            roh.subtitle = roh.subtitle, 
            sig.probe = significant.probe.order, 
            output.dir = output.dir, 
            output.pfx = output.pfx)
    
    # ===============================================================================
    # Clean local ancestry data
    # ===============================================================================
    # For the local ancestry data, subset for the ancOverlap and the ancTypes columns
    local.ancestry.subset = subset(local.ancestry, select = c("ancOverlap", "ancTypes"))
    
    # For the local ancestry data, the SubjectIDs are set as rownames
    # In order to merge the local ancestry data with the already merged data frame (ROH dosage data and phenotype data), we will create a column called "SubjectID" containing all the SubjectIDs
    local.ancestry.subset$SubjectID = rownames(local.ancestry.subset)

    # Create a for loop that loops over each row of the local ancestry data and splits the values in ancOverlap that has two overlap areas
    for (i in 1:nrow(local.ancestry.subset)) {
        # strsplit() seperates values where ancOverlap has 2 values into a list containing 2 values, for example "123, 456", becomes "123" "456"
        # unlist() allows you to access the strsplit values in the list
        # Create two new columns called ancOverlap1 and ancOverlap2
        #   ancOverlap1 contains the first overlap region
        #   ancOverlap2 contains the second overlap region
        # Example: ancOverlap -> ancOverlap1    ancOverlap2
        #            123,456        123             456
        two.overlap  = unlist(strsplit(local.ancestry.subset$ancOverlap[i], ","))
        local.ancestry.subset$ancOverlap1[i] = two.overlap[1]
        local.ancestry.subset$ancOverlap2[i] = two.overlap[2]

        # Repeat for ancTypes because there are some rows that contain 2 ancestry type values
        #   Example: "EA, EE" will be split to "EA" "EE"
        # Create two new columns called ancTypes1 and ancTypes2
        #   ancTypes1 contains the first ancestry type
        #   ancTypes2 contains the second ancestry type
        # Example: ancTypes -> ancTypes1    ancTypes2
        #           EE,EA          EE           EA
        two.ancestry = unlist(strsplit(local.ancestry.subset$ancTypes[i], ","))
        local.ancestry.subset$ancTypes1[i] = two.ancestry[1]
        local.ancestry.subset$ancTypes2[i] = two.ancestry[2]

    }
    print("Finished creating ancOverlap1, ancOverlap2, ancTypes1, and ancTypes2 columns in local ancestry data")

    # Merge the local ancestry data with the already merged ROH and phenotype data
    #   Add all.x = T to ensure that all the phenotype data are there
    merge.all.df = merge(merge.roh.pheno.df, local.ancestry.subset, by = "SubjectID", all.x = T)
    
    # First we will check how many rows have missing ancestry data
    print("Number of rows with missing ancestry data")
    print(sum(is.na(merge.all.df$ancOverlap)))

    # Check how many individuals have missing ancestry data
    print("Number of individuals with missing ancestry data")
    print(length(unique(merge.all.df$SubjectID[is.na(merge.all.df$ancOverlap)])))

    # Check which individuals have missing ancestry data
    print("Individuals with missing ancestry data")
    print(unique(merge.all.df$SubjectID[is.na(merge.all.df$ancOverlap)]))

    # Print all the rows of individuals with missing ancestry data
    print("Rows with missing ancestry data")
    print(merge.all.df[is.na(merge.all.df$ancOverlap), ])

    # Remove all the rows where ancOverlap has NAs 
    #   The phenotype data has individuals not contained in the local ancestry data
    merge.all.df = merge.all.df[!is.na(merge.all.df$ancOverlap), ]

    # Check all the missing ancestry data are removed
    print("Number of rows with missing ancestry data")
    print(sum(is.na(merge.all.df$ancOverlap)))

    print("Finished merging local ancestry data with ROH data and phenotype data")

    # Create two new columns, called overlap1 and overlap2, where they contain the values TRUE or FALSE
    #   overlap1 is T if a probe is between the minimum probe value and the minimum probe + the first overlap value (probe.min <= probe <= probe.min + first overlap)
    #   overlap2 is T if a probe is between the (minimum probe + first overlap) and the maximum probe (probe.min + first overlap <= probe <= probe.max))
    # as.numeric() is used to convert the values from character values to numeric values
    merge.all.df$overlap1 = ifelse(as.numeric(merge.all.df$Probe) >= as.numeric(probe.min) & as.numeric(merge.all.df$Probe) <= as.numeric(probe.min) + as.numeric(merge.all.df$ancOverlap1), T, F)

    merge.all.df$overlap2 = ifelse(as.numeric(merge.all.df$Probe) <= as.numeric(probe.max) & as.numeric(merge.all.df$Probe) >= as.numeric(probe.min) + as.numeric(merge.all.df$ancOverlap1), T, F)

    # Create a new data frame for the probes that are between the minimum probe value and the minimum probe + the first overlap value 
    overlap1.df = subset(merge.all.df, overlap1 == "TRUE" & overlap2 == "FALSE", select = c("SubjectID", "Probe", "ROH", phenotype, "ancOverlap1", "ancTypes1"))
    
    # Create a new data frame for the probes that are bewteen the (minimum probe + first overlap) and the maximum probe 
    overlap2.df = subset(merge.all.df, overlap2 == "TRUE" & overlap1 == "FALSE", select = c("SubjectID", "Probe", "ROH", phenotype, "ancOverlap2", "ancTypes2"))
    
    print("Finished creating new data frames overlap1.df and overlap2.df")

    # Rename the columns for the newly created data frames so that we can later bind the two data frames together
    colnames(overlap1.df) = c("SubjectID", "Probe", "ROH", phenotype, "ancOverlap", "ancTypes")
    colnames(overlap2.df) = c("SubjectID", "Probe", "ROH", phenotype, "ancOverlap", "ancTypes")

    # Use rbind() to bind the two data frames together
    # This newly created data frame will all have rows with only one ancOverlap and one ancTypes value 
    clean.df = rbind(overlap1.df, overlap2.df)
    
    print("Finished binding overlap1.df and overlap2.df")

    # Currently we have the following in our data frame 
    #         pop  chr  queryAnc  queryEurAlleles  queryAfrAlleles  queryNamAlleles    start       end      ancOverlap  ancTypes
    # ID      PR  chr11   EN            1                0                 1         116000000   117000000   1000001       EN
    # However, we want ancTypes to have only one ancestry per row, so we will create 2 new columns called ancestry1 and ancestry2
    #   The column ancestry1 will contain the first ancestry
    #   The column ancestry2 will contain the second ancestry
    #     Example: ancTypes -> ancestry1  ancestry2
    #                 EN           E          N
    clean.df$ancestry1 = unlist(lapply(clean.df$ancTypes, substring, 1, 1))
    clean.df$ancestry2 = unlist(lapply(clean.df$ancTypes, substring, 2, 2))
    
    print("Finished creating additional ancestry columns")

    # Reorder data frame
    clean.df = setorderv(clean.df, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))

    # ====================================================================
    # Make local ancestry data frame
    # ====================================================================

    # Since there are missing probe data, the plots contain vertical gaps
    # To solve this problem, we will assign each unique probe a number
    # First we will create a new data frame called probe.df, with one column containing all the unique probes, and another column containing the unique number corresponding with the probe
    probe.df = data.table(x = unique(clean.df$Probe), y = order(unique(sort(clean.df$Probe))))
    
    # Rename the columns "Probe" and "probe.order"
    colnames(probe.df) = c("Probe", "probe.order")
    
    print("Finished creating probe.order columns")

    # Merge the clean.df data frame with the newly created probe.df data frame
    merge.probe.df = merge(clean.df, probe.df, by = "Probe")
    
    # Reorder the merged data frame by phenotype, SubjectID, and Probe
    merge.probe.df = setorderv(merge.probe.df, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))
    
    print("Finished merging probe.order column to data frame")

    # Create a new data frame called pheno.df that consists of the SubjectID and phenotype columns
    pheno.df = subset(clean.df, select = c("SubjectID", phenotype))    
    
    # Reorder the data frame to ensure that it is ordered by phenotype and SubjectID
    pheno.df = setorderv(pheno.df, cols = c(phenotype, "SubjectID"), order = c(1, 1))
    
    # Create a new data frame called unique.rows.df which contains the unique SubjectID and phenotype rows
    unique.rows.df = unique(pheno.df)

    print("Finished creating pheno.order columns")

    # Create a new data frame called unique.rows, which contains all the unique rows (unique SubjectID and phenotype)
    # Create 2 new columns in the data frame unique.rows called pheno.order1 and 2, which contains the unique order number for each individual and their phenotype value
    # pheno.order1 contains the order number for the first ancestry call, and will be all odd numbers
    # pheno.order2 contains the order number for the second ancestry call, and will be all even numbers
    #   This way, each individual will have 2 lines, one for the first ancestry call and the other for the second ancestry call
    #     In addition, having odd and even rows ensure that the ancestry calls for an individual will stay together and will not be overlapping or overrided
    unique.rows.df$pheno.order1 = order(unique.rows.df[, get(phenotype)]) *2 - 1
    unique.rows.df$pheno.order2 = order(unique.rows.df[, get(phenotype)]) *2
    
    # Merge the ordered phenotype numbers with the ordered probe numbers (which is in merge.probe.df)
    lanc.merge.df = merge(merge.probe.df, unique.rows.df, by = c("SubjectID", phenotype))

    # Resort the merged data frame
    lanc.merge.df = setorderv(lanc.merge.df, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))
    
    print("Finished creating phenotype order for local ancestry plot, now ready to plot local ancestry")

    # ===================================================================
    # Make summary local ancestry data frame
    # ===================================================================
    
    # Melt the columns ancestry1 and ancestry2 together, so you have a combined row of haploid ancestry
    melt.df = melt(clean.df, measure.vars = c("ancestry1", "ancestry2"), variable.name = "ancestry1,2", value.name = "ancestry")    

    # Create a new data frame called summarise.df, where the mean of the phenotype values for each probe and ancestry is calculated
    # So for each probe, there should be a total of 3 means calculated, one for each ancestry
    summarise.df = melt.df %>% group_by(ancestry, Probe) %>% summarise(mean = mean(get(phenotype)))
    
    print("Finished computing mean phenotype value for each probe and ancestry, now ready to plot summary of local ancestry data")

    # ===================================================================
    # Make ROH and local ancestry data frame
    # ===================================================================

    # Since there are missing probe data, the plots contain vertical gaps
    # To solve this problem, we will assign each unique probe a number for the ROH Local Ancestry plot
    # First we will create a new data frame called probe.df1, with one column containing all the unique probes, and another column containing the unique number corresponding with the probe
    probe.df1 = data.table(x = unique(clean.df$Probe), y = order(unique(sort(clean.df$Probe))))

    # Rename the columns "Probe" and "probe.order"
    colnames(probe.df1) = c("Probe", "probe.order")

    # Merge the subsetted ROH data set with the probe order
    roh.lanc.probe.merge = merge(clean.df, probe.df1, by = "Probe")

    # Create a new data frame called pheno.df1 that consists of the SubjectID and phenotype columns
    pheno.df1 = subset(clean.df, select = c("SubjectID", phenotype))    
    
    # Reorder the data frame to ensure that it is ordered by phenotype and SubjectID
    pheno.df1 = setorderv(pheno.df1, cols = c(phenotype, "SubjectID"), order = c(1, 1))

    # Create a new data frame called unique.rows.df1 which contains the unique SubjectID and phenotype rows
    unique.rows.df1 = unique(pheno.df1)

    # Create a new column in the data frame unique.rows.df1 called pheno.order, which contains the unique order number for each individual and their phenotype value
    
    unique.rows.df1$pheno.order = order(unique.rows.df1[, get(phenotype)])

    # Merge the roh.lanc.df data frame (Contains subsetted ROH data and probe order) with the newly created phenotype order numbers
    roh.lanc.merge = merge(roh.lanc.probe.merge, unique.rows.df1, by = c("SubjectID", phenotype))
    
    # Reorder the data frame
    roh.lanc.merge = setorderv(roh.lanc.merge, cols = c(phenotype, "SubjectID", "Probe"), order = c(1, 1, 1))

    # Select for rows where ROH has the value of 1, which means that the probe is in an ROH
    roh.lanc.df = subset(roh.lanc.merge, ROH == 1)

    # Reorder the data frame by phenotype, SubjectID, and Probe in ascending order

    print("Finished creating phenotype order for ROH local ancestry plot and subsetting for probes in ROH; now ready to plot local ancestry for probes in ROH")

    # ===================================================================
    # Plot local ancestry plots
    # ===================================================================

    # We will place a vertical line at the probe of interest, so that we will be able to identify the ROH data at that probe
    # If the probe of interest is not in the data set, the nearest probe will be used instead
    significant.probe.index1 = which.min(abs(lanc.merge.df$Probe - sig.probe)) 
    significant.probe1 = lanc.merge.df[significant.probe.index1, ]$Probe
    significant.probe.order1 = lanc.merge.df[significant.probe.index1, ]$probe.order
    
    print(paste0("Local Ancestry: Vertical line plotted at:", significant.probe1))

    # Repeat for the summary local ancestry plot
    significant.probe.index2 = which.min(abs(summarise.df$Probe - sig.probe))
    significant.probe2 = summarise.df[significant.probe.index2, ]$Probe
    
    print(paste0("Summary Local Ancestry: Vertical line plotted at:", significant.probe2))
    
    # Repeat for the ROH local ancestry plot
    significant.probe.index3 = which.min(abs(roh.lanc.df$Probe - sig.probe))
    significant.probe3 = roh.lanc.df[significant.probe.index3, ]$Probe
    significant.probe.order3 = roh.lanc.df[significant.probe.index3, ]$probe.order
    
    print(paste0("ROH Local Ancestry: Vertical line plotted at:", significant.probe3))

    # Plot all three local ancestry plots using the function LancPlot() which is defined in the plotting_routines.R script
    LancPlot(input.lanc = lanc.merge.df, 
             input.sum.lanc = summarise.df, 
             input.roh.lanc = roh.lanc.df, 
             phenotype = phenotype, 
             colorblind.friendly = colorblind.friendly,
             size = size, 
             lanc.sig.probe = significant.probe.order1, 
             sum.lanc.sig.probe = significant.probe2, 
             roh.lanc.sig.probe = significant.probe.order3, 
             ancestry.title = "Local Ancestry Plot", 
             ancestry.subtitle = ancestry.subtitle, 
             sum.title = sum.title, 
             sum.subtitle = sum.subtitle, 
             roh.lanc.title = roh.lanc.title, 
             roh.lanc.subtitle = roh.lanc.subtitle, 
             output.dir = output.dir, 
             output.pfx = output.pfx)
}
