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
pheno_df <- read.csv(file.path(proj, "DataSets", "RNASeqData", "phenotypeTable.csv"), 
                     header =TRUE)
# import information about molecular data
mol_df <- read.csv(file.path(proj, "DataSets", "RNASeqData", "molecularDataSets.csv"), 
                         header =TRUE, sep=",")
# combine phenotype and molecular data into metadata
metadata_df <- dplyr::full_join(pheno_df, mol_df, 
                                by = "assay_material_id")
# inspect output
head(metadata_df, n = 3)
# write metadata to file
write.csv(metadata_df, file = file.path(proj, "RNASeqProject", "Data", "metadata.csv"), row.names = FALSE)

#### Combine data files ####

# import metadata (if neceassary)
metadata_df <- read.csv(file.path(proj, "RNASeqProject", "Data", "metadata.csv"))
str(metadata_df)

# make list of directories containing data
RNADirectoryList <- list.dirs(path = file.path(proj, "DataSets", "RNASeqData"), 
                              recursive = FALSE)
# find all files for htseq
filelist <- sapply(RNADirectoryList, 
                         function(x){list.files(path = x, 
                                                full.names = TRUE, 
                                                pattern = "htseq.txt") }) 
# create data frame of file names
file_id_df <- as.data.frame(cbind(filelist, 
                                       gsub("^.*-","", RNADirectoryList)), 
                                 stringsAsFactors = F)
colnames(file_id_df) <- c("Path", "molecular_id")
head(file_id_df, n = 3)
# create data frame of all data frames
all_df <- lapply(seq(1:nrow(file_id_df)),
                       function(i){ 
                         X <- read.delim(file = file_id_df$Path[i],
                                         header = FALSE);
                         colnames(X) <- c("Gene", file_id_df$molecular_id[i]);
                         return(X)
                       } )
# create single data frame of counts
genecounts_df <- plyr::join_all(all_df, by = NULL, 
                                       type = "full", match = "all")
# inspect output
head(genecounts_df, n = 1)

# save to file
write.csv(genecounts_df, file.path(proj, "RNASeqProject", "Data", "RNAseqGeneCounts.csv"))
