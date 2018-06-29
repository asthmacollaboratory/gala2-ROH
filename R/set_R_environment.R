#!/usr/bin/env Rscript --vanilla
# ==============================================================================
# Copyright Asthma Collaboratory (2018) 
# Coded by Kevin L. Keys and Jennifer R. Liberto
#
# This script loads all of the libraries required for analyses.
# If a package is not installed, then the script will attempt to install it.
# Note that the installation code is somewhat brittle and unrobust.
# The content of this file is similar in spirit to an .Rprofile.
# However, it is not meant to replace an .Rprofile.
# Instead, load it at the beginning of GALA II ROH analysis scripts like so:
#
#     source("set_R_environment.R")
#
# Nota bene: all variables and subroutines are defined in the current workspace. 
# Consequently, all variables and functions are deletable with
#
#    rm(list = ls() )
#
# so be very careful when cleaning workspaces!
# ==============================================================================

# ==============================================================================
# function definitions
# ==============================================================================

LoadPackage = function(package.name, library.path){
	# LoadPackage 
	#
	# This function loads a package into the current workspace.
	# If the package is not installed, then LoadPackage will attempt to install it. 
	# The installed package is then silently loaded

	if(!require(package.name, lib.loc = library.path, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
		install.packages(package.name, lib = library.path)
		library(package.name, lib.loc = library.path, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
	}
	return()
}

AutoloadPackages = function(vector.of.package.names, library.path){
	# AutoloadPackages
	#
	# This function calls LoadPackage on a vector of package names.
	#
	# Args:
	#	list.of.package.names: a vector of package names, e.g. c("ggplot2", "MASS") 
	#
	# Output: NULL

	invisible(sapply(vector.of.package.names, LoadPackage, library.path))
	return()
}

# ==============================================================================
# environment variables
# ==============================================================================
cran.mirror     = "https://cran.cnr.berkeley.edu/"  # use UC Berkeley CRAN mirror
max.print.lines = 200  # default number of lines to print
editor          = "vim"  # default text editor
library.path    = "/media/BurchardRaid01/LabShare/Data/share_data_projectInProgress/ROH_project/R_libraries"
bitmap.type     = "cairo"  # needed for producing raster plots

# this vector should contain all packages required for analysis
auto.loads = c("ggplot2", "dplyr", "data.table", "ggpubr", "doParallel", "readr",
    "qqman", "ggrepel", "RColorBrewer", "grid", "gridExtra") 

# ==============================================================================
# executable code
# ==============================================================================

# set CRAN mirror to UC Berkeley
local({
    r = getOption("repos")
    r["CRAN"] = cran.mirror 
    options(repos=r)
})

# set group R library path
.libPaths(library.path)

# R allows strings as factors by default
# this is one of the most frustrating engineering decisions in the language
# make the default to *not* convert strings to factors
options(stringsAsFactors=FALSE)

# max number of printed lines hopefully set to something reasonable
options(max.print=max.print.lines)

# default text editor
options(editor=editor)

# turn off tk since we do not need it here
options(menu.graphics=FALSE)

# enable raster graphics
options(bitmapType = "cairo")

# autoload packages
AutoloadPackages(auto.loads, library.path)
