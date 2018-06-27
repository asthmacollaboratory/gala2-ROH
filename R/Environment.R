# Environment

cran.mirror     = "https://cran.cnr.berkeley.edu/"  # use UC Berkeley CRAN mirror
max.print.lines = 200  # default number of lines to print
editor          = "vim"  # default text editor

  # this vector should contain all packages required for analysis
  auto.loads = c("dplyr", "ggplot2", "data.table", "ggpubr", "doParallel", "readr",
    "qqman", "ggrepel", "RColorBrewer", "grid", "gridExtra", "coda") 

  AutoloadPackages(auto.loads)

  # Create 22 parellel cores to run association tests for each chromosome
  registerDoParallel(cores = 22)
