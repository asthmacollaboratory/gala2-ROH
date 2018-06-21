# Post-processing function

results_processing = function(results, results_out) {
	chr = 1
	file = paste(results, chr, results_out, sep = "")
	data = read.table(file, header = T)
	gwas = data.frame(chr = rep(chr,dim(data)[1]),data)

	for(chr in 2:22){
		file = paste(results, chr, results_out, sep = "")
		data = read.table(file,header = T)
		df = data.frame(chr = rep(chr,dim(data)[1]),data)
		gwas = rbind(gwas,df)
		rm(df)
		rm(data)
		}
	x = gwas[order(gwas$p),]
	write.table(x, file = "ordered_probes.txt")	
	}

results_processing(results, results_out)


diagnostic_plots = function(results, manhattan_file_name, manhattan_title, qq_file_name) {
	# Load in GWAS results
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
	qq(gwas$P, xlim = c(0, 8), ylim = c(0, 8))
	dev.off()
	print(lambda1)
	}

diagnostic_plots(results, manhattan_file_name, manhattan_title, qq_file_name)
