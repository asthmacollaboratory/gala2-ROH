#!/usr/bin/env Rscript --vanilla
# ==============================================================================
# Copyright 2018, Asthma Collaboratory
# coded by A. Californian
#
# This script demonstrates example code in proper R style.
# It has no input or output.
# ==============================================================================

# ==============================================================================
# if reading code from another .R file, then do so here like so:
# > source("my_R_code.R")
# ==============================================================================

# ==============================================================================
# load libraries
# ==============================================================================
library(ggplot2)

# ==============================================================================
# define functions here
# ==============================================================================
SumTwoNumbers = function(x, y){
	# Computes the sum of two numbers.
    # Note that a nonfinite argument poisons the outcome.
	#
	# Args:
	#     x: A number.
	#     y: Another number.
	#
	# Returns:
	#     The sum x + y

	return(x + y)
}

# ==============================================================================
# executable code goes here
# ==============================================================================
x = 1
y = 2
z = SumTwoNumbers(x,y)
