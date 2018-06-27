# association.function

Association.function = function(wd, out_name, pop, ROH.wd, in.pheno, suffix, type) {
	# Regression_Model
	#
	# This function merges ROH probes with phenotype data and then runs a linear or logistic regression
	#
	# Args:
	#	wd: specify the working directory to write results to
	#	out_name: this will be the title of your output files
	#	pop: choose which population you are working with.
	#	ROH.wd: ROH calls are in different working directories, specify which directory to call
	#	in.pheno: specify your phenotype dataframe
	#	suffix: when ROH dosage was produced, each cohort had a different name. SAGE = ".ROH.R.out.gz", GALA = ".ROH.R.out"
	#	type: type refers to the GLM family. If binary outcome (i.e. asthma status), include argument: type = "binomial". If continuous outcome (i.e. BDR), include argument: type = "gaussian"
	#
	# Output: Coefficients of linear model
	#	
	  foreach (i = 1:22) %dopar%{

	  inNameLSA = paste(pop,i,suffix, sep = "")
	  outName = paste(out_name,i,'.ROH.R.out.results',sep = "")

	  setwd(ROH.wd)
	  lsa = read.table(inNameLSA, header = T, check.names = F, stringsAsFactors = F)

	  new.lsa = lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,3:dim(lsa)[2]]
	  new.lsa = new.lsa[,names(new.lsa) %in% in.pheno$SubjectID]
	  name.order = order(colnames(new.lsa))
	  new.lsa = new.lsa[,name.order]

	  # Removes duplicates
		exclude = in.pheno$SubjectID[!(in.pheno$SubjectID %in% names(new.lsa))]
		in.pheno = in.pheno[ ! in.pheno$SubjectID %in% exclude, ]

	  pos = lsa[rowSums(lsa[,3:dim(lsa)[2]]) > 0,2]
	  result = t(apply(new.lsa, 1, Regression_Model, in.pheno, type))

	  result.final = cbind(pos,result)
	  colnames(result.final) = c('Probe','beta','se','t','p')
	  setwd(wd)
	  write.table(result.final,outName,quote = F,sep = '\t',row.names = F)
	   } 
	}
