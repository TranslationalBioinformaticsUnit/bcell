#************************************************#
# LINK OCRs TO GENES                             #
# Validation approach with public ChiP-seq data  #
# script: overlap_target_ocrs_vs_histone_marks.R #
#************************************************#

# Validation of candidate regulatory elements using ChiP data



# Required libraries                                       
library(GenomicRanges)
library(rtracklayer)



##>>>>>>>>>>>>>>>>>>>>> Input data
#Read candidate regulatory OCR
ocr_anno <- readRDS("atac_OCR_annotation.rds")
#Generate the Granges object
ocr_anno_gr <- makeGRangesFromDataFrame(ocr_anno, keep.extra.columns = TRUE) 
names(ocr_anno_gr) <- rownames(ocr_anno)
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> Overlap with histone marks
#Read public histone information from encode, create a Granges and compute the overlap with our OCRs

# H3K4me3
h3k4me3 <- read.table("H3K4me3_regions.bed", sep="\t")
names(h3k4me3)[1:4] <- c("chr","start","end","name")
h3k4me3_gr <- makeGRangesFromDataFrame(h3k4me3, keep.extra.columns = TRUE) #37,952

overlapping <- findOverlaps(ocr_anno_gr,h3k4me3_gr, minoverlap = 6)
ocr_gene$h3k4me3 <- rep("no",nrow(ocr_anno))
ocr_gene$h3k4me3[unique(overlapping@from)] <- "yes"

# H3K4me1
h3k4me1 <- read.table("H3K4me1_regions.bed", sep="\t")
names(h3k4me1)[1:4] <- c("chr","start","end","name")
h3k4me1_gr <- makeGRangesFromDataFrame(h3k4me1, keep.extra.columns = TRUE) #166,295

overlapping <- findOverlaps(ocr_anno_gr,h3k4me1_gr, minoverlap = 6)
ocr_gene$h3k4me1 <- rep("no",nrow(ocr_anno))
ocr_gene$h3k4me1[unique(overlapping@from)] <- "yes"

# H3K27ac
h3k27ac <- read.table("H3K27ac_regions.bed", sep="\t")
names(h3k27ac)[1:4] <- c("chr","start","end","name")
h3k27ac_gr <- makeGRangesFromDataFrame(h3k27ac, keep.extra.columns = TRUE) #159,079

overlapping <- findOverlaps(ocr_anno_gr,h3k27ac_gr, minoverlap = 6)
ocr_gene$h3k27ac <- rep("no",nrow(ocr_anno)) 
ocr_gene$h3k27ac[unique(overlapping@from)] <- "yes"
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> Overlap with histone marks
#Identify those OCRs link to active promoter (H3K4me3 + H3K27ac)
ocr_anno$active_promoter_histone_mark <- rep("no",nrow(ocr_anno))
ocr_anno$active_promoter_histone_mark[which(ocr_anno$h3k4me3=="yes" & ocr_anno$h3k27ac=="yes")] <- "yes"
#Identify those OCRs link to active promoter (H3K4me1 + H3K27ac)
ocr_anno$active_enhancer_histone_mark <- rep("no",nrow(ocr_anno))
ocr_anno$active_enhancer_histone_mark[which(ocr_anno$h3k4me1=="yes" & ocr_anno$h3k27ac=="yes")] <- "yes"
#-----------------------------------



saveRDS(ocr_anno, paste0(format(Sys.time(), "%Y%m%d"),"_atac_OCR_annotation_extended_histone.rds"))
#-----------------------------------
