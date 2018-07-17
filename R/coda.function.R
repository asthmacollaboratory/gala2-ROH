# This script contains the function ComputeSignificanceThreshold()
# It will compute the adjusted significance threshold 
# Written by Annie Li

ComputeSignificanceThreshold = function(input.file) {
	# The function computes the significance threshold for your data
	# 
	# Args: 
	#   input.file: Your input file name
	#
	# Returns:
	#   The number of independent tests, Bonferroni correction, the significance threshold, and the suggestive threshold
	
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
