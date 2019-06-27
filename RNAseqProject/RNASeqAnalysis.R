#### Gene expression analysis of myleofibrosis ####

# This script includes examples of code 

# Notes on using RStudio:
#   Any lines starting with # are not code; these are notes for humans and R ignores them
#   Hold down the Control key and press enter to execute (run) the line or chunk of code where the cursor is located

#### Install and load software ####

# install tools for data manipulation and plotting
install.packages("tidyverse")

# install tools for genomic analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# if you reveice the following prompts:
#   "Update all/some/none", type a (for all) and hit Enter
#   "Do you want to install from sources the package which needs compilation? (Yes/no/cancel)" hit Enter 

# Note: package installation only needs to happen once per computer, but you'll need to load the packages every time you reopen RStudio

# load packages
library(tidyverse)
library(DESeq2)