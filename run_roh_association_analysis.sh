#!/usr/bin/env bash

# load environment variables from .env
source .env.sh

# binaries
RSCRIPT=${RSCRIPT}

# directories
ROHdir=${ROHdir}
scratchdir=${scratchdir}
repodir=${repodir}
codedir=${codedir}
gala_roh_dir=${gala_roh_dir}
sage_roh_dir=${sage_roh_dir}
phenofile=${phenofile}

# scripts
R_environment_script="${repodir}/R/set_R_environment.R"
R_parse_phenotypes="${repodir}/parse_phenotype_data.R"
R_run_analysis="${repodir}/run_roh_association_analysis.R"

# variables
pops=${pops}
phenos=${phenos}
glm_types=${glm_types}
covariate_lists=${covariate_lists}
outdirpfx=${outdirpfx}
outdirs=${outdirs}


# ---
# These are variables for association analysis
# these values correspond to defaults for function PerformAssociationAnalyses
# modify them here as needed 

# suffix for output files
outsfx=${outsfx}

# point R to the right directory when loading packages
library_path=${library_path}

# number of parallel cores to use in association analysis; sufficient to use 1 per chromosome
# BE VERY CAREFUL WITH THIS. with 3 pops, this can ask for ncores*3 processing cores 
# to avoid overburdening server when running multiple analyses, try ncores=10 
ncores=${ncores}

# how many samples are required before analyzing a probe?
# Rule of thumb, based on experience:
# To admit a probe for analysis, the minimum samples with an ROH at a probe should be set at >= 8
# probes where only 2-8 people with a ROH segment at that probe will yield unrealistic p-values
# probes with < 2 people with a ROH segument have no variance and can yield numeric errors 
min_samples_at_probe=${min_samples_at_probe}

# file prefixes for GALA, SAGE
galapfx=${galapfx}
sagepfx=${sagepfx}
# ---

# will loop over phenos + covariates
# requires checking lengths of these arrays
# throw error if they aren't the same length
# otherwise, proceed with analysis
npheno="${#phenos[@]}"
ncovariate_lists="${#covariate_lists[@]}"
nmodels="${#glm_types[@]}"
noutdir="${#outdirs[@]}"
if [[ ${npheno} -ne ${ncovariate_lists} ]] || [[ ${npheno} -ne ${nmodels} ]] || [[ ${npheno} -ne ${noutdir} ]]; then
    echo -e "Number of phenotypes = ${npheno}\nNumber of covariate lists = ${ncovariate_lists}\nNumber of glm types = ${nmodels}\nNumber of output directories = ${noutdir}\nAll of these numbers must match." 1>&2
    exit 1
fi

# loop over phenotypes and covariates
for i in $(seq 0 $((${npheno} - 1)) ); do 

    # phenotype name must match what is available in data files 
    pheno="${phenos[$i]}"

    # grab list of covariates for current phenotype
    covariates="${covariate_lists[$i]}"

    # this is the GLM class (e.g. "gaussian" for current $pheno)
    glm_type="${glm_types[$i]}"

    # make output directory if it doesn't already exist
    outdir="${outdirs[$i]}"
    mkdir -p ${outdir} 

    for pop in ${pops[@]}; do

        # $inputpfx changes depending on pop
        # it has one general form for GALA and another for SAGE
        inputpfx="${gala_roh_dir}/${galapfx}_${pop}"
        if [[ "${pop}" = "AA" ]]; then
            inputpfx="${sage_roh_dir}/${sagepfx}"
        fi

        # run analysis
        $RSCRIPT $R_run_analysis \
            --phenotype-name ${pheno} \
            --covariates ${covariates} \
            --source-code-directory ${codedir} \
            --population ${pop} \
            --GALA-ROH-directory ${gala_roh_dir} \
            --SAGE-ROH-directory ${sage_roh_dir} \
            --scratch-directory ${scratchdir} \
            --output-directory ${outdir} \
            --GLM-type ${glm_type} \
            --num-cores ${ncores} \
            --minimum-samples-at-probe ${min_samples_at_probe} \
            --R-library-path ${library_path} \
            --output-suffix ${outsfx} \
            --phenotype-file ${phenofile} \
            --input-prefix ${inputpfx}
    done
done
