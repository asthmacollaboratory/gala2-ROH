# Regression_Model.R

Regression_Model = function(new.lsa, in.pheno, type) {
	# Regression_Model
	#
	# This function regresses a phenotype of interest against covariates of interest
	#
	# Args:
	#	new.lsa: ROH probes
	#	pheno: dataset including relevent phenotype variables
	#	type: specify "binomial" for binary outcome, specify "gaussian" for continous outcome
	#
	# Output: Coefficients of linear model
	#
	     form <- formula(paste(outcome,predictor,covars))
	     tempresult = glm(form, data = in.pheno, family = type)
	     coef = summary(tempresult)$coef[2,]
	     return(coef)
}
