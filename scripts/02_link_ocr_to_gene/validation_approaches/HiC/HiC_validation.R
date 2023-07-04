#*********************************************#
# LINK OCRs TO GENES                          #
# Validation approach with public HiC data    #
# script: HiC_validation.R                    #
#*********************************************#

# Validation of candidate regulatory elements using HiC from datasets associated 
# with the paper Vilarrasa-Blasi R, "Dynamics of genome architecture and chromatin function during human B cell differentiation and neoplastic transformation". Nat Commun 12:2021 651
# All data available at: https://ega-archive.org/datasets/EGAD00001006486


# Required libraries
library(GenomicRanges)

##>>>>>>>>>>>>>>>>>>>>> Input data
#Read candidate regulatory OCR
ocr_anno <- readRDS("atac_OCR_annotation.rds")
#Generate the Granges object
ocr_anno_gr <- makeGRangesFromDataFrame(ocr_anno, keep.extra.columns = TRUE) 
names(ocr_anno_gr) <- rownames(ocr_anno)
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> HiC compartment distribution; overlap with out OCRs

#Read HiC compartments
HiC_compartments <- read.table("HiC_compartments.txt", sep="\t", header=TRUE)
HiC_compartments_gr <- makeGRangesFromDataFrame(HiC_compartments, keep.extra.columns = TRUE) #27,881

#Compare with our ATAC--Seq regions (if overlap in two regions keep the compartment assignment from region with max overlap)
overlap_output <- findOverlaps(ocr_anno_gr,HiC_compartments_gr)

unique_regions <- unique(overlap_output@from)
ocr_anno$HiC_compartment <- rep(NA, nrow(ocr_anno))
for(i in 1:length(unique_regions)){
  hic_regions <- overlap_output@to[overlap_output@from==unique_regions[i]]
  if(length(hic_regions)>1){
    length_overlap <- vector()
    for(j in 1:length(hic_regions)){
      ii <- pintersect(ocr_anno_gr[unique_regions[i]],HiC_compartments_gr[hic_regions[j]])
      length_overlap <- c(length_overlap,end(ii)-start(ii))
    }
    hic_regions <- hic_regions[length_overlap==max(length_overlap)]
  }
  ocr_anno$HiC_compartment[unique_regions[i]] <- HiC_compartments_gr[hic_regions]$Compartment
}

table(ocr_anno$HiC_compartment)

saveRDS(ocr_anno, paste0(format(Sys.time(), "%Y%m%d"),"_atac_OCR_annotation_extended_HiC.rds"))
#-----------------------------------
