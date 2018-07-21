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



PerformAssociationAnalysis = function(input.prefix, output.prefix, phenotype.df, model.formula, suffix = "ROH.R.out", type = "gaussian", ncores = 22, min.samples.at.probe = 2, library.path = "/media/BurchardRaid01/LabShare/Data/share_data_projectInProgress/ROH_project/R_libraries") {
    # PerformAssociationAnalysis
    #
    # This function merges ROH probes with phenotype data. It then fits a generalized linear model to the data.
    # The results are written to file. This function returns nothing.
    # NOTA BENE: PerformAssociationAnalysis makes explicit use of %dopar% from the doParallel package.
    # PerformAssociationAnalysis will set up a parallel backend, farm out computations to slave cores,
    # and then shut down the backend when finished. 
    #
    # Args:
    #     input.prefix: an incomplete file path to the ROH data, e.g.
    #         "${ROH_DIRECTORY}/${FILENAME_BEFORE_CHROMOSOME}"
    #	  output.prefix: the title of the output files.
    #         The full output is "${output.prefix}.ROH.R.out.results"
    #	  phenotype.df: the dataframe with phenotype and covariates

    #     model.formula: a string to create a formula object used for the linear model fit. Can produce as, e.g.,
    #         "y ~ x + z" where y is outcome, x is ROH segment, and z is covariate.
    #	  suffix: when ROH dosage was produced, each cohort had a different file suffix.
    #         SAGE = "ROH.R.out.gz"
    #         GALA = "ROH.R.out"
    #         Default: "ROH.R.out"
    #	  type: specify the GLM family. Must be a valid argument to "family" in glm().
    #         binary outcome --> logistic regression --> type = "binomial"
    #         continuous outcome --> linear regression --> type = "gaussian"
    #         Default: "gaussian"
    #     ncores: how many cores should be used for parallel calculations?
    #         Default: 22 (one for each chromosome)
    #     min.samples.at.probe: at minimum, how many samples must have ROH at a probe to merit analysis?
    #         Default: 2 (probes with 0 or 1 samples with ROH discarded) 
    #
    # Output: Nothing

    # check arguments
    ncores >= 1 || stop("Argument ncores must be at least 1")
    min.samples.at.probe >= 0 || stop("Argument min.samples.at.probe must be at least 0")

    # register a cluster to use for parallel execution
    cl = makeCluster(ncores)
    registerDoParallel(cl)

    # run association analysis for each chromosome separately
    foreach (chr = 1:22) %dopar% {

        # this line ensures that all of the parallel workers look to the same library folder
        # library.path is defined in set_R_environment.R
        .libPaths(c(library.path, .libPaths()))

        input.file.path = paste(input.prefix, chr, suffix, sep = ".")
        output.file.path = paste(output.prefix, chr, "ROH.R.out.results", sep = ".")

        # read the ROH calls from file
        # use read.table in lieu of fread since the latter cannot easily handle gzipped files
        roh.data = read.table(input.file.path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

        # must subset the ROH data to included all samples with a minimum number of probes in ROH segments
        # need two pieces:
        # (1) the number of SNPs in a ROH for each sample, and
        # (2) the column dimension of roh.data
        # the SNPs start in the 3rd column of roh.data
        nCols        = dim(roh.data)[2]
        nSNPs.in.ROH = rowSums(roh.data[,3:nCols])

        # subset roh.data for samples with number of ROHs at least n = min.samples.at.probe
        probes.to.analyze = nSNPs.in.ROH >= min.samples.at.probe
        roh.data.sub = roh.data[probes.to.analyze, 3:nCols]

        # find the subjects in common between roh.data.sub and phenotype.df
        # discard subjects not in both data frames
        # then reorder the columns in increasing alphanumeric order 
        subjects.in.common = names(roh.data.sub) %in% phenotype.df$SubjectID
        roh.data.sub = roh.data.sub[, subjects.in.common] 
        name.order   = order(colnames(roh.data.sub))
        roh.data.sub = roh.data.sub[,name.order]

        # now remove duplicate subjects
        included.ids = phenotype.df$SubjectID %in% names(roh.data.sub)
        phenotype.df.nodup = phenotype.df[included.ids, ]

        # define a regression kernel here
        # ideally we would define it outside of the function. however, in that case, %dopar% fails:
        # the function is defined on master process, and it does NOT propagate to slave processes. 
        # in contrast, defining it within the %dopar% loop ensures that all processes can see it
        FitLinearModel = function(one.roh.snp, phenotype.df, model.formula, type = "gaussian") {
            # FitLinearModel
            #
            # This function regresses a phenotype of interest against one ROH probe
            # and any covariates included in the phenotype
            #
            # Args:
            #	  one.roh.snp: a single SNP ROH probe
            #	  phenotype.df: the dataframe with phenotype and covariates
            #     model.formula: the formula used in the linear model.
            #	  type: the GLM family for regression, passed to argument "family" in glm().
            #         Default: "gaussian"
            #
            # Output: Coefficients of linear model

            temp.df = glm(formula(model.formula), data = phenotype.df, family = type)
            model.coefficients = summary(temp.df)$coef[2,]
            return(model.coefficients)
        }

        # apply FitLinearModel as a kernel to the rows of roh.data.sub
        # this fits a linear model for each of the SNPs
        # remember, each SNP is coded 0 (not in ROH) or 1 (in ROH)
        # the anonymous function wraps FitLinearModel to clarify what argument roh.data.sub occupies in FitLinearModel
        # TODO (Kevin): check if this can be rotated so that apply() works columnwise on roh.data.sub (potential performance boost!)
        result = t(apply(roh.data.sub, 1, function(z) FitLinearModel(z, phenotype.df.nodup, model.formula, type)))

        # add position as column to left of results
        snp.positions = roh.data[probes.to.analyze, 2]
        snp.labels = paste(chr, roh.data[probes.to.analyze, 2], sep = ":")
        result.final  = data.frame(cbind(snp.positions, snp.labels, result))

        # add number of samples analyzed to right of results
        num.samples  = apply(roh.data[,-c(1:3)], 1, function(z) sum(!is.na(z)))[probes.to.analyze]
        #result.final = data.frame(cbind(result, num.samples))
        result.final$nsamples = num.samples

        # also add dummy-coded alleles (G --> ROH = 1, C --> ROH = 0) to left
        result.final$eff_allele = "G"
        result.final$alt_allele = "C"

        # add number of ROH segments observed at the probe
        result.final$nROH = nSNPs.in.ROH[probes.to.analyze]

        # name columns and save to file
        colnames(result.final) = c("position", "Probe", "beta", "stderr", "t", "p", "nsamples", "eff_allele", "alt_allele", "nROH")
        write.table(result.final, file = output.file.path, quote = FALSE, sep = "\t", row.names = FALSE)
    }

    # shut down the cluster
    stopCluster(cl)
    return()
}
