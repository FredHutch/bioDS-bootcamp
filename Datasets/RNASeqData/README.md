# RNA sequencing Data sets for Differential Gene expression
This folder contains data sets from RNA sequencing of 12 peripheral blood samples taken from myelofibrosis patients and 12 from normal blood donors.  The same RNA seq raw data were analyzed via two different approaches, HiSat2/htseqcount and salmon to create gene-count lists.  These lists can then be used to identify both whether there is a difference in gene expression between these two groups of samples, as well as more technical analyses such as whether the bioinformatic tools used to process the data have large impacts on the results of the analysis, and if there are any experimental covariates that may confound the analysis of biological differences.    

## Load libraries and import data
Start a new Rmd, and choose which libraries to load in the first chunk. Then read in the phenoTable and molecularDataSets csv's.  Join the two data frames by a common column to generate a complete data frame of all covariates for the RNA sequencing data in the folder.  

## Read in and merge `htseq.txt` files
List directories in the RNASeqData folder, read in only files that end in `htseq.txt` from each directory into a list of data frames.  Merge the data frames in the list into one large data frame with the rownames being the gene name, and one column for each sample containing the RNA seq counts of reads associated with each gene.  Name the columns with a unique identifier for each data set that will correspond to a column in the covariate data frame you made above.

## Summarize the data sets
Use dplyr to group by and summarize and ggplot to visualize how many groups of samples we have, how many RNA seq data sets are there in each biological group, identify possible confounding variables, assess library size and plot to identify outliers.  

## Create a normalized data frame
Normalize each sample's counts data based on over all library size for each sample.

## Use DESeq2 to do differential gene expression analysis
Identify the differentially expressed genes in the myelofibrosis samples as compared to the normal control samples.  

## Read in and merge `salmon` files, merge and normalize
Copy code above to now list directories in the RNASeqData folder, and read in only files that end in `salmon.txt` from each directory into a list of data frames.  Merge the data frames in the same way as was done for the previous data sets. Normalize the data as was done previously.  

## Bioinformatics tool assessment QC plots
With normalized data tables from both `htseq` and `salmon`, use ggplot to create visualizations of the various aspects of the data from the same samples analyzed using each approach in order to highlight possible differences and relative merits.

## Use DESeq2 to do differential gene expression analysis and compare with previous results
Using the `salmon` data, perform the same differential expression analysis as above and then compare the results of the two approaches in order to see if there are differences in the overall results when different bioinformatic tools are used.  
