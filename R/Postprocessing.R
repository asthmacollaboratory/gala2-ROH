# Postprocessing

results_processing = function(wd, out_name) {
	# results_processing
	#
	# This function merges the chromosome results into one file with probes ordered by P value
	#
	# Args:
	#	wd: set working directory
	#	out_name: title of your output files from association testing
	#
	# Output: Results file
	#
	setwd(wd)
	chr = 1
	file = paste(out_name, chr, '.ROH.R.out.results', sep = "")
	data = read.table(file, header = T)
	gwas = data.frame(chr = rep(chr,dim(data)[1]),data)

	for(chr in 2:22){
		file = paste(out_name, chr, '.ROH.R.out.results', sep = "")
		data = read.table(file,header = T)
		df = data.frame(chr = rep(chr,dim(data)[1]),data)
		gwas = rbind(gwas,df)
		rm(df)
		rm(data)
		}
	x = gwas[order(gwas$p),]
	write.table(x, file = "ordered_results.txt")	
	}

diagnostic_plots = function(wd, results, manhattan_file_name, manhattan_title, qq_file_name, qq_title) {
	# results_processing
	#
	# This function merges the chromosome results into one file with probes ordered by P value
	#
	# Args:
	#	wd: set working directory
	#	results: name of ordered results from results processing function
	#	out_name: title of your output files from association testing
	#	manhattan_file_name: file name of Manhattan plot that will be written to working directory
	#	manhattan_title: title of Manhattan plot
	#	qq_file_name: file name of QQ plot that will be written to working directory
	#	qq_title: title of QQ plot
	#
	# Output: Results file
	#	
	# Load in GWAS results
	setwd(wd)
	gwas = read.table(results, header=T)

	# Reformat dataframe to correct way for manhattan plot
	gwas$BP = gwas$Probe
	gwas = subset(gwas, select = -c(3,4,5) )
	colnames(gwas) = c("CHR", "SNP", "P", "BP")

	# Visualize regression results
	options(bitmapType = 'cairo')
	png(file = manhattan_file_name)
	manhattan(gwas, main = manhattan_title)
	dev.off()

	# Q-Q plot
	#Calculate lambda
	    pvalue1 = gwas$P
	    chisq1 = qchisq(1-pvalue1,1)
	    lambda1 = median(chisq1)/qchisq(0.5,1)
	options(bitmapType = 'cairo')
	png(filename = qq_file_name)
	qq(gwas$P, main = qq_title, xlim = c(0, 8), ylim = c(0, 8))
	dev.off()
	print(lambda1)
	}