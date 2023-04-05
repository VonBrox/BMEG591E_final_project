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

#NOTE: There were duplicate rows for genes labelled "01-Mar" and "02-Mar" in the original data.
#These seemed to be mislabeled, and were causing issues with downstream steps, so we decided to remove these rows.
#Ideally, we would delve into why these genes were mislabelled and correctly identify them, but we do not have sufficient information to do that.
counts_data <- counts_data[counts_data$Gene != "01-Mar", ]
counts_data <- counts_data[counts_data$Gene != "02-Mar", ]

#Set the row names to be the gene names
row.names(counts_data) <- counts_data[,1]
counts_data <- counts_data[, -1]

# Convert the counts data into a format that DGEList will accept: integer matrix
dge_counts <- as.matrix(counts_data)
mode(dge_counts) <-"integer"

#Read in a table of sample labels, so each sample can be matched to disease state (early or late onset colon cancer)
sample_labels <- read.csv('sample_labels.csv')
sample_labels <- sample_labels[!rowSums(is.na(sample_labels)),] #Remove NA values
sample_labels_subset <- sample_labels[, c(8,16)] # Create list of sample names with early vs late onset, to be used later for the heatmap
sample_labels_subset <- sample_labels_subset[order(sample_labels_subset$disease_state),] # Sort the list
sample_labels_subset$EorL <- c(as.numeric(sample_labels_subset$disease_state == "Early-onset colorectal cancer")) #Used for sample labelling in the later heatmap

# Normalize read counts for every gene with the trimmed mean of M-values (TMM) method, using the edgeR package to create a DGEList object.
# 2 disease state groups: early-onset and late-onset colon cancer
library(edgeR)
dge <- DGEList(counts=dge_counts, samples = as.data.frame(colnames(dge_counts)), genes = rownames(counts_data), group = sample_labels$disease_state)
dge <- calcNormFactors(dge, method = "TMM")
#head(dge$samples) # Data frame containing the normalization factors
dge_counts_normalized <- cpm(dge, log=TRUE) # Data frame of read counts normalized based on CPM and the TMM normalization factors, converted to log2
dge_counts_normalized_z_scores <- scale(dge_counts_normalized) # Calculating Z-score for the counts

#Create GLM model for differential gene expression comparison
design <- model.matrix(~sample_labels$disease_state)
dge = estimateDisp(dge, design)#Check this function
fit = glmQLFit(dge, design)#Check this function
qlf = glmQLFTest(fit, coef=2) #Check this function
topTags(qlf, n=10, adjust.method="BH", sort.by="logFC", p.value=0.01)

# To expand on the work of the authors by introducing an appropriate QC step, we made an MA plot:
plotQLDisp(fit)

# Calculating FDR values for all the P-values using the same method as topTags
dge_fdr <- as.data.frame(p.adjust(qlf$table$PValue, method = "BH"))
colnames(dge_fdr) <- c("FDR")

# Creating a data frame of the DGE data that can be filtered as specified in the paper
dge_fc <- c(qlf$genes, qlf$table, dge_fdr)
dge_fc <- data.frame(dge_fc)
dge_fc[,2] <- dge_fc[,2] * -1 # The paper's supplementary table had fold-change values with opposite sign because FC comparison had been made in the opposite direction. 
head(dge_fc)

# Filter DGE data based on following criteria: 
# p < 0.01, |log2FC| >1, logCPM > 1
dge_fc_filtered <- dge_fc[which(dge_fc$PValue < 0.01),]# try to instead filter by FDR, after I get same results
dge_fc_filtered <- dge_fc_filtered[which(dge_fc_filtered$logCPM > 1),]
dge_fc_filtered <- dge_fc_filtered[which(abs(dge_fc_filtered$logFC) > 1),]
#This narrowed the list of candidate genes from 33119 to 142

# Sort based on fold change values
dge_fc_filtered <- dge_fc_filtered[order(-dge_fc_filtered$logFC),]

# The authors have made their filtered DGE gene list available as a supplementary appendix 1 (see our github folder).
# It appears that many of the genes match between our list and theirs, and the values in the corresponding columns are usually close.
# A notable difference is that the top 4 logFC candidates in our table (SFTA3, SFTPB, SFRP1, SLC5A8) and the bottom 4 (LYPD2, TSIX, ITLN2, and LOC102723453) are not in the authors' table.
# Some of these have very large |logFC| values (e.g. SFTA3: 6.796528, LYPD3: -8.039594)
# None of these should have been filtered out based on the filtering criteria, so we think they should be included.
# For comparison, we will make heat maps based on both our list and the authors'.

#To make a heatmap as in the paper, we need to have a matrix of z-scores for only the genes remaining in the filtered data
heatmap_data <- as.data.frame(cbind(rownames(dge_counts_normalized_z_scores), dge_counts_normalized_z_scores))
colnames(heatmap_data)[1] <- "Gene"
heatmap_data_subset <- heatmap_data[heatmap_data$Gene %in% dge_fc_filtered$genes, ] # Get the normalized counts of the genes from the filtered and sorted DGE list 
heatmap_data_subset <- heatmap_data_subset[match(dge_fc_filtered$genes, heatmap_data_subset$Gene), ] # Sort as per the DGE list
heatmap_data <- heatmap_data_subset
rm(heatmap_data_subset)
heatmap_data <- as.matrix(heatmap_data) #Convert to a matrix
heatmap_data <- heatmap_data[,-1] #No longer need genes column for matching
mode(heatmap_data) <-"numeric"

# Generate heatmap (as per Figure 1 in the paper) using complexheatmap package (Gu, Z. Complex Heatmap Visualization. iMeta 2022.)
library(ComplexHeatmap)
library(circlize)
my_colors <- colorRamp2(c(-4, 0, 4), c("green", "black", "red"))

onset_labels = HeatmapAnnotation(df = as.data.frame(c(rep(0, 49), rep(1, 50))), 
                                gp = gpar(col = "black"), show_legend = FALSE, 
                                which = "column")

heatmap_1 <- Heatmap(heatmap_data, heatmap_legend_param = list(title = ""), col = my_colors, 
        show_row_dend = FALSE, show_column_dend = FALSE,
        row_order = rownames(heatmap_data), column_order = sample_labels_subset$X.1,
        top_annotation = onset_labels)
heatmap_1

# There are some general similarities between our heatmap and the heatmap shown in the paper. 
# It appears that there are differences in expression between the early- and late-onset cohorts for these genes.
# Both our heatmap and the heatmap in the paper have more high expression visible in the top left and bottom left quadrants.
# However, there seem to be some genes in our list with high expression across the board that drown out the signal.
# To determine whether the issue is with the list of genes we used, we will make a new heatmap using the genes listed in the appendix table.

#Read in the supplementary table
paper_dge_fc <- read.csv('cam45675-sup-0001-appendixs1.csv')

#As before, for the heatmap input, we make a matrix of z scores, this time for only the genes in the authors' gene list
heatmap_2_data <- as.data.frame(cbind(rownames(dge_counts_normalized_z_scores), dge_counts_normalized_z_scores))
colnames(heatmap_2_data)[1] <- "Gene"
heatmap_2_data_subset <- heatmap_2_data[heatmap_2_data$Gene %in% paper_dge_fc$Gene, ] # Get the normalized counts of the genes from the paper's DGE list 
heatmap_2_data_subset <- heatmap_2_data_subset[match(paper_dge_fc$Gene, heatmap_2_data_subset$Gene), ] # Sort as per the DGE list
heatmap_2_data <- heatmap_2_data_subset
rm(heatmap_2_data_subset)
heatmap_2_data <- as.matrix(heatmap_2_data) #Convert to a matrix
heatmap_2_data <- heatmap_2_data[,-1] #No longer need genes column for matching
mode(heatmap_2_data) <-"numeric"

# For the paper's data, generate heatmap (as per Figure 1 in the paper) using complexheatmap package (Gu, Z. Complex Heatmap Visualization. iMeta 2022.)
library(ComplexHeatmap)
library(circlize)
my_colors <- colorRamp2(c(-4, 0, 4), c("green", "black", "red"))

# onset_labels = HeatmapAnnotation(df = as.data.frame(c(rep(0, 49), rep(1, 50))), 
#                                  gp = gpar(col = "black"), show_legend = FALSE, 
#                                  which = "column")

onset_labels = HeatmapAnnotation(bar = c(rep(0, 49), rep(1, 50)),
                                 col = list(bar = c("1" = "yellow", "0" = "blue")),
                                 gp = gpar(col = "black"), show_legend = FALSE)

heatmap_2 <- Heatmap(heatmap_2_data, heatmap_legend_param = list(title = ""), col = my_colors, 
                     show_row_dend = FALSE, show_column_dend = FALSE,
                     row_order = rownames(heatmap_2_data), column_order = sample_labels_subset$X.1,
                     top_annotation = onset_labels)
heatmap_2

# This heat map does not look much better - if anything, it looks worse.
# This suggests that the problem probably lies with our normalization method,rather than our gene list.
# The paper's description of the normalization method is quite brief, so 
