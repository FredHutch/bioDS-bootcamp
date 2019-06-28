#### Gene expression analysis of myleofibrosis ####

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

#### Code above here is included in student script ####

#### Examine metadata ####

# import phenotype data
metadata_df <- read.csv("data/metadata.csv")

# examine import
colnames(metadata_df)
head(metadata_df)

# summarize phenotype data
metadata_df %>% 
  group_by(diagnosis, age_range, sex) %>% 
  summarize(n())
## summarize genotype data
metadata_df %>%
  group_by(genotype_jak2, genotype_calr) %>%
  summarize(n())

#### Create summarized experiment object ####

# import counts as matrix
counts_df <- read.csv("Data/RNAseqGeneCounts.csv", row.names=1)
counts_mat <- as.matrix(counts_df)

# create metadata matrix
metadata_df <- read.csv("Data/metadata.csv", row.names=1)
meta_mat <- as.matrix(subset(metadata_df, select=-molecular_id))
nrow(meta_mat) 
# use first column as row names
rownames(meta_mat)<- metadata_df$molecular_id

# match column names for counts with row names for metadata
meta_mat <- meta_mat[match(colnames(counts_mat),row.names(meta_mat)),]
head(meta_mat, n = 1)

# create deseq object
counts_deseq <- DESeqDataSetFromMatrix(counts_mat, meta_mat, design = ~ diagnosis)
# set base level to normal (control)
counts_deseq$diagnosis <- relevel(counts_deseq$diagnosis, "normal")

# remove low gene counts (must be greater than 1)
counts_deseq <- counts_deseq[rowSums(counts(counts_deseq)) > 1, ] 

#### Analyze data ####

# perform differential gene analysis
deseq <- DESeq(counts_deseq)

# view results 
deseq_results <- results(deseq)
head(deseq_results, n = 2)

# MA plot
plotMA(deseq)
plotMA(deseq, ylim=c(-2, 2))
# red points have p value less than 0.1
# triangles are points occurruing outside the window shown

# view count plots for single genes
DESeq2::plotCounts(counts_deseq, "CALR", intgroup = "diagnosis")
DESeq2::plotCounts(counts_deseq, "LINC00221", intgroup = "diagnosis")
DESeq2::plotCounts(counts_deseq, "ACOT6", intgroup = "diagnosis")
DESeq2::plotCounts(counts_deseq, "ADAMTS7", intgroup = "diagnosis")
# other possibilities: CALR, LINC00221, ACOT6, ADAMTS7
# JAK2, MPL, TET2
DESeq2::plotCounts(counts_deseq, "JAK2", intgroup = "diagnosis")
DESeq2::plotCounts(counts_deseq, "MPL", intgroup = "diagnosis")
DESeq2::plotCounts(counts_deseq, "TET2", intgroup = "diagnosis")

