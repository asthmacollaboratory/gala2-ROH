#!/usr/bin/env --bash
# =====================================================================
# copyright Asthma Collaboratory (2018)
# coded by A. Californian
#
# This script provides a template for shell scripts.
# It is meant to be run with Bourne Against Shell (BASH) only.
#
# This script has no inputs or outputs, apart from echoing "goodbye" upon exit.
# This happens from a bash exit trap; see script code below.
#
# Call:
#
#     bash example_bash_script.sh
# =====================================================================

# set 
set -e          # exit if any command has nonzero exit status
set -u          # exit if referencing any undefined variable
set -o pipefail # if any component of a pipe fails, then pass along nonzero exit status

# source any other shell scripts here like so:
# > source "${FILENAME}.sh"

# define any functions as needed
function finish {
    echo "goodbye"
}

# executable code goes here


# prefer $(COMMAND) versus backticking `COMMAND` 
files=$(ls -la $(pwd)) 

# for-loop syntax
for i in $(seq 1 3); do
    mystring="${i}";
done

# if-else statement syntax
condition1=0
condition2=1
condition3=0
if [[ "${condition1}" = "${condition2}" ]]; then
    condition3=2
elif [[ "${condition1}" = "${condition3}" ]]; then
    condition3=1
else
    condition3=-1
fi

# this is a bash exit trap
# the function finish would contain any cleanup code
# by trapping it, the code will execute when the script exits
trap finish EXIT
