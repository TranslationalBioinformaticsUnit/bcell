#******************************************************#
# FOOTPRINT TO TF MOTIF                                #
# script: 16_unified_consensus_fp_ocr_motif_viestra.R  #
#******************************************************#


# 16. FP - Motifs (Archetypes from Viestra)


# Load required libraries
library(GenomicRanges)


##>>>>>>>>>>>>>>>>>>>>> Get list motifs from Vierstra
vierstra_annotation <- read.table("reference_vierstra_motifs.txt")
#30,258,213 x 8

#Generate the GRanges with target vierstra db annotation
vierstra_gr <- makeGRangesFromDataFrame(vierstra_annotation,
                                        keep.extra.columns=TRUE)
#30,258,213 ranges
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Generate a GRanges from our consensus footprint
#Get the GRanges of footprint
fp <- read.table("consensus_unified_footprint_all_cell_types.bed")
colnames(fp)[1:3] <- c("chr","start","end")
rownames(fp) <- paste(fp[,1],fp[,2],fp[,3], sep="_")
#keep just those FP with a minimum length of 6
fp_L6 <- fp[(fp[,3]-fp[,2])>=6,]
#380,922      6
FP_gr <- makeGRangesFromDataFrame(fp_L6)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Relate each FP to motif
archetypes <- unique(vierstra_gr$Motif)
archetype_info <- read.table("motif_annotations_archetypes.txt", sep="\t", header=TRUE)

fp_arch_matrix <- matrix(0,nrow=nrow(fp_L6),ncol=length(archetypes))
colnames(fp_arch_matrix) <- archetypes
rownames(fp_arch_matrix) <- names(FP_gr)

for(i in 1:length(archetypes)){
  cat(i," - ")
  #select granges for target motif
  vierstra_gr_motif <- vierstra_gr[vierstra_gr$Motif==archetypes[i]]
  #identify motif length
  archetype_width <- min(end(vierstra_gr_motif)-start(vierstra_gr_motif))
  #ranges overlaping (min overlap 90%)
  fp_arch_over90_association <- findOverlaps(vierstra_gr_motif,FP_gr,minoverlap=max(round(archetype_width*0.9),6))
  if(length(fp_arch_over90_association)>0){
    fp_arch_matrix[unique(fp_arch_over90_association@to),i] <- 1    
  }  
}

##Region without footprint-motif annotation:
cat("Regions without footprint-arch annotation:\n")
cat(length(which(rowSums(fp_arch_matrix)==0)),"\n") #146,134
##Peaks with footprint-motif annotation:
cat("Regions wit footprint-arch annotation:\n")
cat(length(which(rowSums(fp_arch_matrix)>0)),"\n") #234,788
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Save the output
save(list=c("fp_arch_matrix"),file=paste(format(Sys.Date(),"%Y%m%d"),"_unified_consensus_fp_arch_vierstra.RData",sep=""))
#-------------------------------------


  
##>>>>>>>>>>>>>>>>>>>>> FP to motif exploration
pdf(paste(format(Sys.Date(),"%Y%m%d"),"_histogram_unified_consensus_fp_arch.pdf",sep=""))
    hist(rowSums(fp_arch_matrix), breaks=30, main="Histogram: number of archetypes in a FP", xlab="Number of Archetypes")
dev.off()
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Generate a GRanges from our consensus OCRs
#Get the GRanges of OCRs
ocr <- read.table("consensus_all_cell_types.bed")
colnames(ocr)[1:3] <- c("chr","start","end")
rownames(ocr) <- paste(ocr[,1],ocr[,2],ocr[,3],sep="_")
#105,530 ocr x 6

OCR_gr <- makeGRangesFromDataFrame(ocr)
#105,530 ranges
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Compute the overlap between OCR and FP
##Overlap
findOverlaps(FP_gr,OCR_gr, type="within") #380,693
fp_by_ocr <- countOverlaps(OCR_gr,FP_gr)

cat("OCRs without footprint:",sum(fp_by_ocr==0),"\n") #32,759
cat("OCRs with footprint:",sum(fp_by_ocr>0),"\n") #72,771
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> FP to OCR exploration
pdf(paste(format(Sys.Date(),"%Y%m%d"),"_histogram_fp_to_ocr.pdf",sep=""))
hist(fp_by_ocr, breaks=100, main="Histogram: number of FP in an OCR", xlab="Number of footprints")
dev.off()
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Relate each FP-Arch to OCR
ocr_fp_association <- findOverlaps(FP_gr,OCR_gr, type="within")
#380,693 of links

ocr_fp_arch_matrix <- matrix(0,nrow=nrow(ocr),ncol=ncol(fp_arch_matrix))
colnames(ocr_fp_arch_matrix) <- colnames(fp_arch_matrix)
rownames(ocr_fp_arch_matrix) <- rownames(ocr)

for(i in 1:nrow(ocr_fp_arch_matrix)){
  if(i%in%ocr_fp_association@to){
    tmp_ocr_matrix <- fp_arch_matrix[ocr_fp_association@from[ocr_fp_association@to==i],]
    if(is.matrix(tmp_ocr_matrix)){
      ocr_fp_arch_matrix[i,which(colSums(tmp_ocr_matrix)>0)] <- 1
    } else {
      ocr_fp_arch_matrix[i,] <- tmp_ocr_matrix
    }
  }
} 

cat("OCRs without ocr-fp-arch annotation:",sum(rowSums(ocr_fp_arch_matrix)==0),"\n") #43,339
cat("OCRs with ocr-fp-arch annotation:",sum(rowSums(ocr_fp_arch_matrix)>0),"\n") #62,191
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Save the output
save(list=c("ocr_fp_arch_matrix"),file=paste(format(Sys.Date(),"%Y%m%d"),"_unified_consensus_ocr_fp_arch_vierstra.RData",sep=""))
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> OCR to motif exploration
pdf(paste(format(Sys.Date(),"%Y%m%d"),"_histogram_OCR_to_archetype.pdf",sep=""))
   hist(rowSums(ocr_fp_arch_matrix), breaks=40, main="Histogram: number of Archetypes in a OCR", xlab="Number of Archetypes")
dev.off()
#-------------------------------------