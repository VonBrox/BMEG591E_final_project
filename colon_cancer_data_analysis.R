# Goal: process the RNAseq read counts data from the late-onset vs. early-onset colon cancer study, and generate figures 1 and 2 from the paper.
# Optional: Use the data on the GTEX site to generate figures 3 and 4

# RNAseq raw reads file was downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213092

# Store the RNAseq raw reads file as a dataframe
library(readr)
library(dplyr)
library(tidyverse)

counts_data <- read_tsv('GSE213092_HRCRC_rowsort.txt')
counts_data <- as.data.frame(counts_data)
#rownames(counts_data) <- counts_data[,1]

# Normalize read counts for every gene with the trimmed mean of M-values (TMM) method using the edgeR package in R
# Convert the counts data into a format that DGEList will accept: convert to matrix

dge_counts <- as.matrix(counts_data)


mode(dge_counts[,2]) <-"integer"
#Remove column with NA values (now we don't have a gene name label, must add later)
dge_counts<-dge_counts[,-1]
#Remove rows with NAs
dge_counts <- dge_counts[!rowSums(is.na(dge_counts)),]


# Normalize read counts, creating DGEList object. [[Note: we have to assign these to the groups early and late colon cancer!]]
library(edgeR)

sample_labels <- read.csv('sample_labels.csv')
sample_labels <- sample_labels[!rowSums(is.na(sample_labels)),]

dge <- DGEList(counts=dge_counts, samples = as.data.frame(colnames(dge_counts)), genes = counts_data$Gene, group = sample_labels$disease_state)
dge <- calcNormFactors(dge, method = "TMM")


#Below: steps from assignment 6 for the GLM model DE comparison, to modify for group labelling

design <- model.matrix(~sample_labels$disease_state) 
# we are using timeName here to make sure that time is treated as a categorical variable. Had we more time points it might make sense to treat time as a value.
dge = estimateDisp(dge, design)

fit = glmQLFit(dge, design)
qlf = glmQLFTest(fit, coef=2) 
topTags(qlf)

# To expand on the work of the authors by introducing an appropriate QC step, we made an MA plot:

plotQLDisp(fit)

#Figure 1 in the study is a heatmap.

