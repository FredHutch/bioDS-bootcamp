# DNA sequencing Data sets for Variant Analysis in a Disease
This folder contains data sets from DNA sequencing of 12 peripheral blood samples taken from myelofibrosis patients using a targeted DNA sequencing panel.

## Data
1. molecularDataSets.csv -- this file contains annotation information about the sequencing details for each molecular data set made from each assay material.


2. phenotypeTable.csv -- this file contains annotation information about the clinical and specimen processing details for each assay material.

3. Sequencing data folders (e.g. JAK2-1-1-R0631-M00000812) -- folders are organized by patient, assay material id, and molecular id.  Each folder contains:
  - GATK.vcf file -- Variant Call Format file (https://en.wikipedia.org/wiki/Variant_Call_Format)
  - GATK.ann.txt file -- list of annotated variants
  - QCStats.csv file -- quality metrics summary file

## Annotation Definitions
  To see the current definitions of the annotations included alongside the datasets in this repo, please visit this link and [click on the "TGDCC" tab.](https://translationalgenomics.shinyapps.io/FHOntologyBrowser/)
