# Goal: process the RNAseq read counts data from the late-onset vs. early-onset colon cancer study, and generate figures 1 and 2 from the paper.
# Optional: Use the data on the GTEX site to generate figures 3 and 4

# RNAseq raw reads file was downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213092

library(readr)
library(dplyr)
library(tidyverse)

# Store the RNAseq raw reads file as a dataframe
counts_data <- read_tsv('GSE213092_HRCRC_rowsort.txt')
counts_data <- as.data.frame(counts_data)
counts_data <- counts_data[!rowSums(is.na(counts_data)),] #Remove NA values from counts data

# Convert the counts data into a format that DGEList will accept: integer matrix
dge_counts <- as.matrix(counts_data)
mode(dge_counts) <-"integer"
dge_counts<-dge_counts[,-1] #Remove column with NA values

#Read in a table of sample labels, so each sample can be matched to disease state (early or late onset colon cancer)
sample_labels <- read.csv('sample_labels.csv')
sample_labels <- sample_labels[!rowSums(is.na(sample_labels)),] #Remove NA values

# Normalize read counts for every gene with the trimmed mean of M-values (TMM) method, using the edgeR package to create a DGEList object.
# 2 disease state groups: early-onset and late-onset colon cancer
library(edgeR)
dge <- DGEList(counts=dge_counts, samples = as.data.frame(colnames(dge_counts)), genes = counts_data$Gene, group = sample_labels$disease_state)
dge <- calcNormFactors(dge, method = "TMM")
head(dge$samples) # Data frame containing the normalization factors

#Create GLM model for differential gene expression comparison
design <- model.matrix(~sample_labels$disease_state) 
dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design)
qlf = glmQLFTest(fit, coef=2)
topTags(qlf, n=10, adjust.method="BH", sort.by="logFC", p.value=0.01)

# To expand on the work of the authors by introducing an appropriate QC step, we made an MA plot:
plotQLDisp(fit)

# Calculating FDR values for all the P-values using the same method as topTags
dge_fdr <- as.data.frame(p.adjust(qlf$table$PValue, method = "BH"))
colnames(dge_fdr) <- c("FDR")

# Creating a data frame of the DGE data that can be filtered as specified in the paper
dge_fc <- c(qlf$genes, qlf$table, dge_fdr)
dge_fc <- data.frame(dge_fc)
head(dge_fc)

# Filter DGE data based on following criteria: 
# p < 0.01, |log2FC| >1, logCPM > 1
dge_fc_filtered <- dge_fc[which(dge_fc$PValue < 0.01),]
dge_fc_filtered <- dge_fc_filtered[which(dge_fc_filtered$logCPM > 1),]
dge_fc_filtered <- dge_fc_filtered[which(abs(dge_fc_filtered$logFC) > 1),]
#This narrowed the list of candidate genes from 33119 to 141

# Sort based on fold change values
dge_fc_filtered <- dge_fc_filtered[order(-dge_fc_filtered$logFC),]

#To make a heatmap as done in the paper, we need to have a dataframe of normalized counts for the filtered gene list.
dge_counts_normalized <- cpm(dge) # Data frame of read counts normalized based on the TMM normalization factors

# Subset df2 to match the order and labels of df1
dge_counts_normalized_topFC <- dge_counts_normalized[dge_counts_normalized$ %in% df1$id, ]
df2_subset <- df2_subset[match(df1$id, df2_subset$id), ]

# Generate heatmap (as per Figure 1 in the paper)b using complexheatmap package
# (Gu, Z. Complex Heatmap Visualization. iMeta 2022.)
library(ComplexHeatmap)

Heatmap(dge_fc_filtered)



