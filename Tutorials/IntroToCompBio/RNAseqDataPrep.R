#### Preparation of RNAseq data ####

# This script includes some manipulations of the data that need to occur to have everything ready for analysis.

# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("rprojroot")

# load packages
library(DESeq2)
library(rprojroot)

# set top level directory
proj <- find_root_file(criterion = has_file(".git/index"))

#### Prepare metadata ####

# copy original data files to project directory
file.copy(file.path(proj, "DataSets", "RNASeqData", "phenotypeTable.csv"), file.path(proj, "RNASeqProject", "Data", "phenotypeTable.csv"))
file.copy(file.path(proj, "DataSets", "RNASeqData", "molecularDataSets.csv"), file.path(proj, "RNASeqProject", "data", "molecularDataSets.csv"))

# import phenotype data
df_pheno <- read.csv(file.path(proj, "DataSets", "RNASeqData", "phenotypeTable.csv"), 
                     header =TRUE)
# import information about molecular data
df_molecular <- read.csv(file.path(proj, "DataSets", "RNASeqData", "molecularDataSets.csv"), 
                         header =TRUE, sep=",")
# combine phenotype and molecular data into metadata
df_metadata <- dplyr::full_join(df_pheno, df_molecular, 
                                by = "assay_material_id")
# inspect output
head(df_metadata, n = 3)
# write metadata to file
write.csv(df_metadata, file = file.path(proj, "RNASeqProject", "Data", "metadata.csv"), row.names = FALSE)

#### Combine data files ####

# import metadata (if neceassary)
df_metadata <- read.csv(file.path(proj, "RNASeqProject", "Data", "metadata.csv"))
str(df_metadata)

# make list of directories containing data
RNADirectoryList <- list.dirs(path = file.path(proj, "DataSets", "RNASeqData"), recursive = FALSE)
# find all files for htseq
FileList_htseq1 <- sapply(RNADirectoryList, 
                         function(x){list.files(path = x, 
                                                full.names = TRUE, 
                                                pattern = "htseq.txt") }) 
# create data frame of file names
htseq_FileID_df <- as.data.frame(cbind(FileList_htseq1, 
                                       gsub("^.*-","", RNADirectoryList)), 
                                 stringsAsFactors = F)
colnames(htseq_FileID_df) <- c("Path", "molecular_id")
head(htseq_FileID_df, n = 3)
# create list of all data frames
listOf_alldf <- lapply(seq(1:nrow(htseq_FileID_df)),
                       function(i){ 
                         X <- read.delim(file = htseq_FileID_df$Path[i],
                                         header = FALSE);
                         colnames(X) <- c("Gene", htseq_FileID_df$molecular_id[i]);
                         return(X)
                       } )
# create single data frame of counts
htseq_genecounts_df <- plyr::join_all(listOf_alldf, by = NULL, 
                                       type = "full", match = "all")
# inspect output
head(htseq_genecounts_df, n = 1)

# save to file
write.csv(htseq_genecounts_df, file.path(proj, "RNASeqProject", "Data", "RNAseqGeneCounts.csv"), 
          row.names = FALSE)

#### Create summarized experiment object ####

# create count matrix
htseqCountsMat <- as.matrix(htseq_genecounts_df[,-1]); ncol(htseqCountsMat)  
# use first column as row names
rownames(htseqCountsMat)<- htseq_genecounts_df[,1] 
head(htseqCountsMat, n=1)

# create phenotype matrix
phenoMat <- as.matrix(subset(df_metadata, select=-molecular_id)); nrow(phenoMat) 
# use first column as row names
rownames(phenoMat)<- df_combined$molecular_id

# join count and phenotype matrices
phenoMat <- phenoMat[match(colnames(htseqCountsMat),row.names(phenoMat)),]
head(phenoMat, n = 1)

# create deseq object
dseq_set_htseq <- DESeqDataSetFromMatrix(htseqCountsMat,phenoMat, design = ~ diagnosis)
# set base level to normal (control)
dseq_set_htseq$diagnosis <- relevel(dseq_set_htseq$diagnosis, "normal")

# remove low gene counts (must be greater than 1)
dseq_set_htseq <- dseq_set_htseq[rowSums(counts(dseq_set_htseq)) > 1, ] 

#Make copy of htseq data set so original is not modified, multiple formulas can be applied
dseq_set_htseq_copy <- dseq_set_htseq
design(dseq_set_htseq_copy) <- formula(~ diagnosis)
dseq_set_htseq_copy <- DESeq(dseq_set_htseq_copy)

#### Analyze data ####

# view results 
res_dseq_set_htseq_copy <- results(dseq_set_htseq_copy)
head(res_dseq_set_htseq_copy, n = 2)

# MA plot
plotMA(dseq_set_htseq_copy, ylim=c(-2, 2))
# red points have p value less than 0.1. 
# triangles are points occurruing outside the window shown

# view count plots for single genes
DESeq2::plotCounts(dseq_set_htseq, "CALR", intgroup = "diagnosis")
# other possibilities: LINC00221, ACOT6, ADAMTS7
# JAK2, MPL, TET2
