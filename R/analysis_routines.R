# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by Kevin L. Keys, Andrew Zeiger, Zachary Szpiech
#  
# This script contains analysis functions for association testing. 
# Source this script after setting the R environment, e.g.
#
#     source("set_R_environment.R")
#     source("analysis_routines.R")
#
# This ensures that all requisite libraries are loaded before defining the
# plotting functions contained herein.
#
# This script analysis_routines.R contains no executable code.
#
# It has no input or output.
# ==============================================================================

# ==============================================================================
# define functions here
# ==============================================================================


PerformAssociationAnalysis = function(input.prefix, output.prefix, phenotype.df, model.formula, suffix = "ROH.R.out", type = "gaussian", ncores = 22) {
    # PerformAssociationAnalysis
    #
    # This function merges ROH probes with phenotype data.
    # It then fits a generalized linear model to the data.
    # The results are written to file. This function returns nothing.
    # NOTA BENE: PerformAssociationAnalysis makes explicit use of %dopar% from the doParallel package.
    # PerformAssociationAnalysis chokes if no parallel backend has been configured before being called.
    #
    # Args:
    #     input.prefix: an incomplete file path to the ROH data, e.g.
    #         "${ROH_DIRECTORY}/${FILENAME_BEFORE_CHROMOSOME}"
    #	  output.prefix: the title of the output files.
    #         The full output is "${output.prefix}.ROH.R.out.results"
    #	  phenotype.df: the dataframe with phenotype and covariates
    #     model.formula: a formula object used for the linear model fit. Can produce with, e.g.
    #         formula("y ~ x + z") where y is outcome, x is ROH segment, and z is covariate.
    #	  suffix: when ROH dosage was produced, each cohort had a different file suffix.
    #         SAGE = "ROH.R.out.gz"
    #         GALA = "ROH.R.out"
    #         Default: "ROH.R.out"
    #	  type: specify the GLM family. Must be a valid argument to "family" in glm().
    #         binary outcome --> logistic regression --> type = "binomial"
    #         continuous outcome --> linear regression --> type = "gaussian"
    #         Default: "gaussian"
    #     ncores: how many cores should be used for parallel calculations?
    #         Defatul: 22 (one for each chromosome)
    #
    # Output: Nothing

    # check arguments
    ncores >= 1 || stop("Argument ncores must be at least 1")

    # register a cluster to use for parallel execution
    cl = makeCluster(ncores)
    registerDoParallel(cl)

    # run association analysis for each chromosome separately
    foreach (chr = 1:22) %dopar% {

        input.file.path = paste(input.prefix, chr, suffix, sep = ".")
        output.file.path = paste(output.prefix, chr, "ROH.R.out.results", sep = ".")

        # read the ROH calls from file
        roh.data = read.table(input.file.path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

        # must subset the ROH data to included all samples with 1+ ROH segments
        # need two pieces:
        # (1) the number of SNPs in a ROH for each sample, and
        # (2) the column dimension of roh.data
        # the SNPs start in the 3rd column of roh.data
        nSNPs.in.ROH = rowSums(roh.data[,3:dim(roh.data)[2]])
        nCols        = dim(roh.data)[2]

        # subset roh.data for samples with at least one ROH in it
        roh.data.sub = roh.data[nSNPs.in.ROH > 0, 3:nCols]

        # find the subjects in
        roh.data.sub = roh.data.sub[, names(roh.data.sub) %in% phenotype.df$SubjectID]
        name.order   = order(colnames(roh.data.sub))
        roh.data.sub = roh.data.sub[,name.order]

        # remove duplicates
        included.ids = phenotype.df$SubjectID %in% names(roh.data.sub)
        #ids.to.exclude = phenotype.df$SubjectID[!included.ids]
        #phenotype.df.nodup = phenotype.df[ !phenotype.df$SubjectID %in% ids.to.exclude, ]
        phenotype.df.nodup = phenotype.df[included.ids, ]

        # define a regression kernel here
        # ideally we would define it outside of the function
        # however, in that case, %dopar% fails since function is defined on master process,
        # and it is NOT propagated to slave processes, so the regression fails
        # in contrast, defining it within the %dopar loop ensures that all processes can see it
        FitLinearModel = function(one.roh.snp, phenotype.df, model.formula, type = "gaussian") {
            # FitLinearModel
            #
            # This function regresses a phenotype of interest against one ROH probe
            # and any covariates included in the phenotype
            #
            # Args:
            #	  one.roh.snp: a single SNP ROH probes
            #	  phenotype.df: the dataframe with phenotype and covariates
            #     model.formula: the formula used in the linear model.
            #	  type: specify the GLM family. Must be a valid argument to "family" in glm().
            #         Default: "gaussian"
            #
            # Output: Coefficients of linear model

            temp.df = glm(model.formula, data = phenotype.df, family = type)
            model.coefficients = summary(temp.df)$coef[2,]
            return(model.coefficients)
        }

        # apply FitLinearModel as a kernel to the rows of roh.data.sub
        # this fits a linear model for each of the SNPs
        # remember, each SNP is coded 0 (not in ROH) or 1 (in ROH)
        # the anonymous function wraps FitLinearModel to clarify what argument roh.data.sub occupies in FitLinearModel
        result = t(apply(roh.data.sub, 1, function(z) FitLinearModel(z, phenotype.df.nodup, model.formula, type)))

        # add position as column to left of results
        snp.positions = roh.data[nSNPs.in.ROH > 0, 2]
        result.final = cbind(snp.positions, result)
        colnames(result.final) = c("Probe", "beta", "stderr", "t", "p")
        write.table(result.final, file = output.file.path, quote = FALSE, sep = "\t", row.names = FALSE)
    }

    # close the cluster
    stopCluster(cl)
    return()
}
