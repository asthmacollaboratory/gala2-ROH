#!/usr/bin/env bash
# =====================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys 
#
# This BASH script creates a METAL script for use with GWAS meta-analysis.
# The script needs the following command line parameters for input:
# -- input file 1
# -- input file 2
# -- input file 3
# -- output file
# -- path to METAL executable
#
# The order of the input files is irrelevant, but all three must appear 
# before the arguments to the output and the path to the METAL executable.
#
# Note: METAL will choke when producing output if the user lacks write permissions.
# Make sure to 
#
# Call:
#
#     bash metal_metaanalaysis.sh $INPUT1 $INPUT2 $INPUT3 $OUTPUTFILE $PATH_TO_METAL_EXECUTABLE 
# =====================================================================

set -e          # exit if any command has nonzero exit status
set -u          # exit if referencing any undefined variable
set -o pipefail # if any component of a pipe fails, then pass along nonzero exit status


# subroutine to output example call to this script
display_usage() {
    echo -e "This script performs meta-analysis of exactly three GWASes with METAL."
    echo -e "Usage: bash metal_metaanalaysis.sh \$INPUT1 \$INPUT2 \$INPUT3 \$OUTPUTFILE_PREFIX \$PATH_TO_METAL_EXECUTABLE\n" 
    echo -e "Example: bash metal_metaanalysis.sh input1.txt input2.txt input3.txt output /usr/bin/metal\n"
}


CMDLINEPARAM=5
# parse command line arguments
# if the user supplies less than 5 arguments, then output a usage suggestion and nonzero exit
# if the user asks for help, then also output the usage suggestion but with zero exit
# if the user supplies at least 5 arguments, then parse the first 5 and execute
if [[ $# -lt "${CMDLINEPARAM}" ]]; then
    display_usage
    exit 1
fi
if [[ ( $# == "--help") ||  $# == "-h" ]]
    display_usage
    exit 0
fi

if [[ $# -ge "${CMDLINEPARAM}" ]]; then
    input1=$1
    input2=$2
    input3=$3
    output=$4
    metal=$5
fi

# repeat arguments to user; helps with debugging
echo -e "Arguments given:"
echo -e "\tfirst input file: $input1"
echo -e "\tsecond input file: $input2"
echo -e "\tthird input file: $input3"
echo -e "\toutput file: $output"
echo -e "\tpath to METAL executable: $metal"
echo -e "All further arguments are ignored."

# make temporary working directory and file
$tmpdir=$(mktemp -d)
metalscript="${tmpdir}/script.metal"

# this here-document writes the METAL script using our BASH variables
# comments in METAL are also prefaced with a hash "#"
# commands are capitalized (always?)
(cat << SOMETAL
# enable genomic control correction 
GENOMICCONTROL ON

# use effect size and stderr analysis scheme
SCHEME STDERR

# declare the output file with no filtering conditions
OUTFILE ${output}_ .txt 

# describe and process the first input file 
MARKER Probe 
ALLELE eff_allele alt_allele 
EFFECT beta 
PVALUE p 
STDERR stderr
WEIGHT nsamples 
SEPARATOR COMMA

# process the first input file
PROCESS $input1 

# second and third input files have same format 
# can simply process them outright
PROCESS $input2 
PROCESS $input3 

# perform "metal"-analysis!
ANALYZE HETEROGENEITY 

# end of METAL script, followed by end of here-doc
QUIT
SOMETAL
) > $metalscript 

# run metal
$metal $metalscript

# cleanup script
rm -rf $tmpdir

# all done!
echo -e "\ndone."
