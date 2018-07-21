#!/usr/bin/env bash
#
# coded by Kevin L. Keys (2018)
# This BASH script creates a METAL script for use with GWAS meta-analysis.

# needs the following command line parameters for input:
# -- input file 1
# -- input file 2
# -- input file 3
# -- output file
# -- path to METAL executable
# the order of the input files is irrelevant

# parse command line arguments
CMDLINEPARAM=5
if [ $# -ge $CMDLINEPARAM ]; then
    input1=$1
    input2=$2
    input3=$3
    output=$4
    metal=$5
else
    input1="metal.input1.txt" 
    input2="metal.input2.txt" 
    input3="metal.input3.txt" 
    output="metal.output.txt" 
    metal="/usr/bin/metal"
fi

# make temporary working directory and file
tmpdir="/tmp/metal"
mkdir -p $tmpdir
metalscript="${tmpdir}/script.metal"

(cat << SOMETAL
# enable genomic control correction 
GENOMICCONTROL ON

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

# declare the output file with no filtering conditions
OUTFILE $output 

# use effect size and stderr analysis scheme
SCHEME STDERR

# perform "metal"-analysis!
ANALYZE HETEROGENEITY 

QUIT
SOMETAL
) > $metalscript 

# run metal
$metal $metalscript

# cleanup script
rm -rf $tmpdir
