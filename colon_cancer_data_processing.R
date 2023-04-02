# Goal: Download fastq files from the late-onset vs. early-onset colon cancer study, and process them as performed in the paper.

# File SRA numbers: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP396591&o=acc_s%3Aa

# SRA tools tutorials: https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/#download-sequence-data-files-usi
# https://erilu.github.io/python-fastq-downloader/

# From the paper: "The sequencing reads had a length of 100 bp and were
# paired-end.
# These reads were aligned to the hg38 human
# reference genome19 using HISAT2 aligner version 2.1.0.20
# The aligned reads were counted using featureCounts in the
# Subread 2.0.3 package.21"

# Step 1: Perform fastqc quality check of the fastq files (not mentioned in original paper, but an important step to consider)
mkdir fastQC_output  
fastqc -o ./fastQC_output *.fastq.gz
