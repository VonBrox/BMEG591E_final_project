
# Introduction:

Our analysis is based on the work of Ha et al. (2023) entitled:
"Reduced expression of alanyl aminopeptidase is a robust biomarker of non-familial adenomatous polyposis and non-hereditary nonpolyposis colorectal cancer syndrome early-onset colorectal cancer".
This study aimed to identify differential gene expression markers to distinguish early-onset
colorectal cancer from late-onset colorectal cancer of the subtypes non-familial	adenomatous	polyposis	(FAP)	and	non-hereditary nonpolyposis.
The group performed RNA-seq of 49 early-onset and 50 late-onset colorectal cancer samples.
They analyzed differentially-expressed genes between these groups and filtered them based on lof2FC, logCPM, and P-values.
They validated the identified differentially-expressed genes using data from TCGA (the cancer genome atlas) database, based on similar expression profiles.
They found that the gene analyl aminopeptidase (ANPEP) was significantly downregulated in early-onset colorectal cancer patients in both their cohort and the TCGA cohort.
This association was further supported by methylation data and information from the GTEx and GSE196006 datasets.
They concluded that the gene ANPEP was significantly down-regulated in early-onset colorectal cancer, and could serve as a biomarker.


#Methods:

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
dge_counts_normalized_z_scores <- apply(dge_counts_normalized, 2, function(x) (x - mean(x)) / sd(x)) # Calculating Z-score for the counts

#Create GLM model for differential gene expression comparison
design <- model.matrix(~sample_labels$disease_state)
dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design)
qlf = glmQLFTest(fit, coef=2) 
topTags(qlf, n=10, adjust.method="BH", sort.by="logFC", p.value=0.01) #Checking output

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

# The authors have made their filtered DGE gene list available in supplementary appendix 1.
paper_dge_fc <- read.csv('cam45675-sup-0001-appendixs1.csv')
paper_dge_fc #Authors' filtered gene list
dge_fc_filtered #Our filtered gene list for comparison.

# It appears that most genes match between our list and theirs, and the values in the corresponding columns are usually close.
# A notable difference is that the top 4 logFC candidates in our table (SFTA3, SFTPB, SFRP1, SLC5A8) and the bottom 4 (LYPD2, TSIX, ITLN2, and LOC102723453) are not in the authors' table.
# Some of these have very large |logFC| values (e.g. SFTA3: 6.796528, LYPD3: -8.039594)
# None of these should have been excluded based on the filtering criteria, so we think they should be part of the analysis.
# Furthermore, these genes have very low FDR values, which is further evidence in favor of their inclusion.
# If the authors had some other reason to exclude these genes, we think they should have mentioned this in the paper.
# For comparison, we will make heat maps based on both our gene list and the authors'.

#Heat map 1: based on our gene list
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

heatmap_ourList <- Heatmap(heatmap_data, heatmap_legend_param = list(title = ""), col = my_colors, 
        show_row_dend = FALSE, show_column_dend = FALSE,
        row_order = rownames(heatmap_data), column_order = sample_labels_subset$X.1,
        top_annotation = onset_labels)
heatmap_ourList

# There are some general similarities between our heatmap and the heatmap shown in the paper. 
# It appears that there are differences in expression between the early- and late-onset cohorts for these genes.
# Both our heatmap and the heatmap in the paper seem to have higher expression levels visible in the top left and bottom right quadrants.
# However, there seem to be some genes in our list with high expression across the board that overwhelm the signal.
# To determine whether the issue is with the list of genes we used, we will make a new heatmap using the genes listed in the appendix table.

#Heat map 2, based on authors' gene list:
#As before, for the heatmap input, we make a matrix of z scores, this time for only the genes in the authors' gene list
heatmap_authorList_data <- as.data.frame(cbind(rownames(dge_counts_normalized_z_scores), dge_counts_normalized_z_scores))
colnames(heatmap_authorList_data)[1] <- "Gene"
heatmap_authorList_data_subset <- heatmap_authorList_data[heatmap_authorList_data$Gene %in% paper_dge_fc$Gene, ] # Get the normalized counts of the genes from the paper's DGE list 
heatmap_authorList_data_subset <- heatmap_authorList_data_subset[match(paper_dge_fc$Gene, heatmap_authorList_data_subset$Gene), ] # Sort as per the DGE list
heatmap_authorList_data <- heatmap_authorList_data_subset
rm(heatmap_authorList_data_subset)
heatmap_authorList_data <- as.matrix(heatmap_authorList_data) #Convert to a matrix
heatmap_authorList_data <- heatmap_authorList_data[,-1] #No longer need genes column for matching
mode(heatmap_authorList_data) <-"numeric"

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

heatmap_authorList <- Heatmap(heatmap_authorList_data, heatmap_legend_param = list(title = ""), col = my_colors, 
                     show_row_dend = FALSE, show_column_dend = FALSE,
                     row_order = rownames(heatmap_authorList_data), column_order = sample_labels_subset$X.1,
                     top_annotation = onset_labels)
heatmap_authorList

# This heat map also has some visible differences between the early onset and late onset groups, but it has the same issue with visual noise as the first heatmap.
# This suggests that the problem with our method is primarily with our normalization steps, rather than our gene list, since this gene list was the same as the one used in the paper.
# The paper's description of the normalization method used to display the heatmap is brief and lacks detail (for example, z-score normalization was not mentioned), so the order in which steps were carried out is unclear.
# Possible further approaches to replicate the paper's results might include:
#   - Calculating the z scores before calculating cpm, or only calculating z scores without cpm.
#   - Removing high-signal genes and re-scaling (say, removing the top 6 genes with the highest sum of values across rows).
#   - Performing one or more of the normalization steps (TMM, cpm, and/or z-scores) after first filtering the genes.

# In the interest of time, we will move on to figure 2 in the paper: a volcano plot of DEGs using the EnhancedVolcano package.
# (Blighe K, Rana S, Lewis M (2022). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling. R package version 1.16.0, https://github.com/kevinblighe/EnhancedVolcano.)

#Volcano plot based on our unfiltered DGE gene list:
library(EnhancedVolcano)
volcano_ourList_unfiltered <- EnhancedVolcano(dge_fc,
                             lab = dge_fc$genes,
                             x = 'logFC',
                             y = 'PValue',
                             pCutoff = 0.01,
                             FCcutoff = 1,
                             pointSize = 3.0,
                             labSize = 2.0,
                             title = "",
                             caption = "")
volcano_ourList_unfiltered

# Volcano plot based on the authors' filtered DGE gene list:
library(EnhancedVolcano)
volcano_authorList_filtered <- EnhancedVolcano(paper_dge_fc,
                lab = paper_dge_fc$Gene,
                x = 'logFC',
                y = 'Pvalue',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 2.0,
                title = "",
                caption = "")
volcano_authorList_filtered

#Volcano plot based on our filtered DGE gene list:
library(EnhancedVolcano)
volcano_ourList_filtered <- EnhancedVolcano(dge_fc_filtered,
                             lab = dge_fc_filtered$genes,
                             x = 'logFC',
                             y = 'PValue',
                             pCutoff = 0.01,
                             FCcutoff = 1,
                             pointSize = 3.0,
                             labSize = 2.0,
                             title = "",
                             caption = "")
volcano_ourList_filtered

# These plots display the genes that meet the authors' defined threshold values for |log2FC| and -log10 P-value in the top right and left squares.
# The gene NTS appears to be at or near the same position in all 3 volcano plots as well as the Figure 2 plot shown in the paper,
# which is supporting evidence that our analysis replicated correctly.
# A zoomed-in image of the volcano_authorList_filtered plot shows that the positions of DKK4, CCL19, and PRSS33 match
# between this plot and the Figure 2 plot in the paper, which further supports that our replication of Figure 2 is successful.
# Other genes are not labelled in the paper but appear in the same positions between the 3 plots (e.g. DHRS1, CCDC88B, CDK12).

# Since the filtered volcano plots only include the genes that already met the threshold values,
# they only consist of genes that fall into the top right and left squares of the volcano plot. 
# The authors only included their filtered gene list in the appendix, so we cannot compare their unfiltered list to ours.

# When comparing our filtered list to the authors' filtered list, it is apparent that ours includes genes 
# with much higher |log2FC| and -log10P values than those identified by the authors.
# It remains unclear why the authors excluded these genes from the analysis given their low FDR values, and this bears further investigation.

# Conclusion:

# References:
# 1. Ha, Ye Jin, Yun Jae Shin, Ka Hee Tak, Jong Lyul Park, Jeong Hwan Kim, Jong Lyul Lee, Yong Sik Yoon, Chan Wook Kim, Seon Young Kim, and Jin Cheon Kim. “Reduced Expression of Alanyl Aminopeptidase Is a Robust Biomarker of Non-Familial Adenomatous Polyposis and Non-Hereditary Nonpolyposis Colorectal Cancer Syndrome Early-Onset Colorectal Cancer.” Cancer Medicine n/a, no. n/a. Accessed April 5, 2023. https://doi.org/10.1002/cam4.5675.
