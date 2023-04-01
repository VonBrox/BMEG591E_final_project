# Goal: process the RNAseq read counts data from the late-onset vs. early-onset colon cancer study, and generate figures 1 and 2 from the paper.
# Optional: Use the data on the GTEX site to generate figures 3 and 4

# RNAseq raw reads file was downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213092

# Store the RNAseq raw reads file as a dataframe
library(readr)
library(dplyr)

counts_data <- read_tsv('GSE213092_HRCRC.txt')
counts_data <- as.data.frame(counts_data)

counts_data <- counts_data %>%
  column_to_rownames("Gene")

counts_data %>%
  mutate_all(as.numeric)

# Normalize read counts for every gene with the trimmed mean of M-values (TMM) method using the edgeR package in R
# Convert the counts data into a format that DGEList will accept: convert to matrix

dge_counts <- as.matrix(counts_data)

#Will try to see if removing the row names and column names avoids the NA error
# row.names(dge_counts) <- NULL
# colnames(dge_counts) <- NULL

#did not work
which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(badNum & !is.na(x))
}
lapply(counts_data, which.nonnum)

#dge_counts <- subset(dge_counts, select = -c(1))
#colnames(dge_counts)<-NULL
#dge_counts[] <- lapply(dge_counts, as.numeric)


# Normalize read counts, creating DGEList object
library(edgeR)
dge <- DGEList(counts=dge_counts, samples = as.data.frame(colnames(counts_data)), genes = as.data.frame(row.names(counts_data)))
dge <- calcNormFactors(dge, method = "TMM")

#Figure 1 in the study is a heatmap.

