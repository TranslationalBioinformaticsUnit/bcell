#******************************************************#
# Aggregated archetype score - with chromVAR           #
# script: 17_archetype_chromVAR.R                      #
#******************************************************#


# 17. Aggregated archetype score - chromVAR (Archetypes from Viestra)



# Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVAR)
library(TFBSTools)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pals)
library(BuenColors)
library(circlize)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load internal function
source("chromVAR_internal_functions.R")



##>>>>>>>>>>>>>>>>>>>>> Inputs required

# 1- Set of OCR
peaks_file="consensus_all_cell_types.bed"

# Read consensus peaks
peaks <- read.table(peaks_file) #105,530
rownames(peaks) <- paste(peaks[,1],peaks[,2],peaks[,3],sep="_")
colnames(peaks) <- c("chr","start","end")

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks)


# 2- OCR counts and metadata
# counts
counts_file <- "count_data_atac.txt"
my_counts_matrix <- as.matrix(read.delim(counts_file, check.names = FALSE))
#remove outlier sample
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix)=="ImmatureB_ATAC_D215_S5")]

# metadata
atac_metadata <- read.table("metadata_atac.txt", sep="\t", header=TRUE)
rownames(atac_metadata) <- atac_metadata$SampleID
atac_metadata2 <- atac_metadata[colnames(my_counts_matrix),]
#order cell type labels
atac_metadata2$Tissue <- factor(as.character(atac_metadata2$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))


# 3 - OCR-archetype (motif) binary matrix
ocr_fp_arch_matrix <- read.table("ocr_arch_matrix.txt", sep="\t", header=TRUE, check.names = FALSE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Run chromVAR with internal function chromVAR_assay
chromVAR_assay(ocr_granges=peak_granges, ocr_counts=my_counts_matrix, ocr_metadata=atac_metadata2, ocr_motif_anno=ocr_fp_arch_matrix, assay_name="OCR_archetype", target_tfs=c("PAX","GATA","Ebox","RUNX","EBF1"))
#-------------------------------------
