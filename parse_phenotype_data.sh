#!/usr/bin/env bash

# source environment variables
source .env.sh

# binaries
RSCRIPT=${RSCRIPT}

# directories
MYHOME=${MYHOME}
ROHdir=${ROHdir}
outdir="${ROHdir}/analysis/data"

mkdir -p ${outdir}

# scripts
R_environment_script="${repodir}/R/set_R_environment.R"
R_parse_phenotypes="${repodir}/parse_phenotype_data.R"
R_run_analysis="${repodir}/run_roh_association_analysis.R"

# file paths
gala_data_file=$gala_data_file
sage_data_file=$sage_data_file
gala_bmi_file=$gala_bmi_file
gala_pca_file=$gala_pca_file
sage_pca_file=$sage_pca_file
gala_ancestry_file=$gala_ancestry_file
sage_ancestry_file=$sage_ancestry_file

$RSCRIPT $R_parse_phenotypes \
    --GALA-data-file ${gala_data_file} \
    --SAGE-data-file ${sage_data_file} \
    --R-environment-code ${R_environment_script} \
    --GALA-BMI-file ${gala_bmi_file} \
    --GALA-PCA-file ${gala_pca_file} \
    --SAGE-PCA-file ${sage_pca_file} \
    --GALA-ancestry-file ${gala_ancestry_file} \
    --SAGE-ancestry-file ${sage_ancestry_file} \
    --output-directory ${outdir}
