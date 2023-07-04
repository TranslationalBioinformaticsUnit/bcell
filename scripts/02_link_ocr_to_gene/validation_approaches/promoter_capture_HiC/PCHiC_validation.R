#*********************************************#
# LINK OCRs TO GENES                          #
# Validation approach with public PCHiC data  #
# script: PCHiC_validation.R                  #
#*********************************************#

# Validation of candidate regulatory elements using PCHiC from datasets associated 
# with the paper Javierre* / Burren* / Wilder* / Kreuzhuber* /Hill* et al. (2016). Genomic regulatory architecture links enhancers and disease variants to target gene promoters. Cell 167.
# All data available at: https://osf.io/u8tzp/ 


# Required libraries
require(rtracklayer)
require(liftOver)
library(GenomicRanges)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pals)
library(BuenColors)
library(circlize)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)



##>>>>>>>>>>>>>>>>>>>>> Input data
#Read candidate regulatory OCR
ocr_gene <- readRDS("candidate_CREs.rds")
candidate_peaks <- matrix(unlist(strsplit(as.character(ocr_gene$peak),"_")),ncol=3,byrow = TRUE)
candidate_peaks <- data.frame(chr=candidate_peaks[,1],start=candidate_peaks[,2],end=candidate_peaks[,3],peak_type=ocr_gene$peak_type)
#Generate the Granges object
peak_regulation_gr <- makeGRangesFromDataFrame(candidate_peaks, keep.extra.columns = TRUE) #62,559
names(peak_regulation_gr) <- ocr_gene$peak
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> Public interaction regions
# Load public data
data <- read.table("PCHiC_peak_matrix_cutoff5.tsv", header=TRUE) #728,838 x 30

# Keep only those rows where nB are >=5 (Naive B cell significant interactions)
Bdata <- data[which(data$nB>=5),-c(12:26)] #197,442 x 15

# Identify links promoter-promoter regions and links larger than 1Mb to do it comparable with our approach
Bdata$PPI <- rep("no",nrow(Bdata)) #promoter-promoter interaction regions
Bdata$PPI[which(Bdata$oeName==".")] <- "yes"  #172,263 PPI (87%) 
Bdata$distant_PEI <- rep("no",nrow(Bdata)) #promoter-enhancer interactions larger than 1Mb 
Bdata$distant_PEI[which(abs(Bdata$dist)>1000000)] <- "yes" #24,776 PEI > 1Mb (12.5%)

# Transform the coordinates from GRCh37 to GRCh38
chain <- import.chain("hg19ToHg38.over.chain")
# 1-Promoter positions
Bdata$baitChr <- paste("chr",Bdata$baitChr,sep="")
Bdata_for_granges_PR <- Bdata[,c(1,2,3)]
names(Bdata_for_granges_PR) <- c("chr","start","end")
Bdata_gr_PR <- makeGRangesFromDataFrame(unique(Bdata_for_granges_PR), keep.extra.columns = TRUE) #14,260

for(i in 1:length(Bdata_gr_PR)){
  sel <- which(Bdata$baitChr==chrom(Bdata_gr_PR)[i] & Bdata$baitStart==start(Bdata_gr_PR)[i] & Bdata$baitEnd==end(Bdata_gr_PR)[i])
  output <- as.data.frame(liftOver(Bdata_gr_PR[i], chain))
  if(nrow(output)>1){
    output2 <- output[as.character(output$seqnames)==as.character(seqnames(Bdata_gr_PR[i])),]
    Bdata$baitChr_hg18[sel] <- as.character(unique(output2$seqnames))
    Bdata$baitStart_hg18[sel] <- min(output2$start)
    Bdata$baitEnd_hg18[sel] <- max(output2$end)
  }
  if(nrow(output)==1){
    Bdata$baitChr_hg18[sel] <- as.character(output$seqnames)
    Bdata$baitStart_hg18[sel] <- output$start
    Bdata$baitEnd_hg18[sel] <- output$end
  }
}

# 2-Interacting regions
Bdata$oeChr <- paste("chr",Bdata$oeChr,sep="")
Bdata_for_granges_IR <- Bdata[,c(6,7,8)]
names(Bdata_for_granges_IR) <- c("chr","start","end")
Bdata_gr_IR <- makeGRangesFromDataFrame(unique(Bdata_for_granges_IR), keep.extra.columns = TRUE) #95,410

for(i in 1:length(Bdata_gr_IR)){
  sel <- which(Bdata$oeChr==chrom(Bdata_gr_IR)[i] & Bdata$oeStart==start(Bdata_gr_IR)[i] & Bdata$oeEnd==end(Bdata_gr_IR)[i])
  output <- as.data.frame(liftOver(Bdata_gr_IR[i], chain))
  if(nrow(output)>1){
    output2 <- output[as.character(output$seqnames)==as.character(seqnames(Bdata_gr_IR[i])),]
    Bdata$oeChr_hg18[sel] <- as.character(unique(output2$seqnames))
    Bdata$oeStart_hg18[sel] <- min(output2$start)
    Bdata$oeEnd_hg18[sel] <- max(output2$end)
  }
  if(nrow(output)==1){
    Bdata$oeChr_hg18[sel] <- as.character(output$seqnames)
    Bdata$oeStart_hg18[sel] <- output$start
    Bdata$oeEnd_hg18[sel] <- output$end
  }
}

write.table(Bdata,paste(format(Sys.Date(),"%Y%m%d"),"_PCHiC_peak_matrix_cutoff5_hg18.txt",sep=""),sep="\t",row.names = FALSE, quote=FALSE)
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> Check overlap between our regions and the promoter Biola
Bdata <- read.table("PCHiC_peak_matrix_cutoff5_hg38.txt", header=TRUE, sep="\t") #197,442 x 21

Bdata_hg18_PR <- Bdata[,c(16,17,18,5)]
colnames(Bdata_hg18_PR) <- c("chr","start","end","gene")
Bdata_hg18_PR_gr <- makeGRangesFromDataFrame(unique(Bdata_hg18_PR), keep.extra.columns = TRUE) #14,260

ocr_gene$overlap_PR_Biola <- rep("no", nrow(ocr_gene))

(promoter_overlap <- findOverlaps(peak_regulation_gr, Bdata_hg18_PR_gr, minoverlap = 6))
ocr_gene$overlap_PR_Biola[unique(promoter_overlap@from)] <- "yes"  #13,597

ocr_gene$overlap_PR_Biola_baitPeak <- rep(NA, nrow(ocr_gene))
ocr_gene$overlap_PR_Biola_baitName <- rep(NA, nrow(ocr_gene))
for(i in 1:length(promoter_overlap)){
  bb <- Bdata_hg18_PR_gr[promoter_overlap@to[i]]
  peak <- paste(chrom(bb),start(bb),end(bb),sep="_")
  gene <- bb$gene
  if(!is.na(ocr_gene$overlap_PR_Biola_baitPeak[promoter_overlap@from[i]])){
    peak <- unique(c(ocr_gene$overlap_PR_Biola_baitPeak[promoter_overlap@from[i]],peak))
    gene <- unique(c(unlist(strsplit(ocr_gene$overlap_PR_Biola_baitName[promoter_overlap@from[i]],";")),gene))
    if(length(peak)>0){peak <- paste(peak,collapse = ";")}
    if(length(gene)>0){gene <- paste(gene, collapse = ";")}
  }
  ocr_gene$overlap_PR_Biola_baitPeak[promoter_overlap@from[i]] <- peak
  ocr_gene$overlap_PR_Biola_baitName[promoter_overlap@from[i]] <- gene
}

# Are they overlapping in the same regions?
for(i in 1:nrow(ocr_gene)){
  ocr_gene$overlap_PR_Biola_ourgene_in_baitName[i] <- ocr_gene$hgnc_symbol[i]%in%unlist(strsplit(ocr_gene$overlap_PR_Biola_baitName[i],";"))                                                            
}

## Results:
# Overlapping promoter regions: 
summary(as.factor(ocr_gene$overlap_PR_Biola)) #13,597
table(ocr_gene$overlap_PR_Biola,ocr_gene$peak_type) #5,197 overlapping with enhancers and 8400 with promoters
table(ocr_gene$overlap_PR_Biola_ourgene_in_baitName,ocr_gene$peak_type) #7,651 promoters annotated to same genes and 370 enhancers annotated to same promoter gene
#-----------------------------------



##>>>>>>>>>>>>>>>>>>>>> Check overlap between our regions and the enhancer Biola
Bdata_hg18_ER <- Bdata[,c(19,20,21,5)]
colnames(Bdata_hg18_ER) <- c("chr","start","end","gene")
Bdata_hg18_ER_gr <- makeGRangesFromDataFrame(Bdata_hg18_ER, keep.extra.columns = TRUE) #197,442

ocr_gene$overlap_ER_Biola <- rep("no", nrow(ocr_gene))

(enhancer_overlap <- findOverlaps(peak_regulation_gr, Bdata_hg18_ER_gr, minoverlap = 6))
ocr_gene$overlap_ER_Biola[unique(enhancer_overlap@from)] <- "yes"  #21,476

ocr_gene$overlap_ER_Biola_oePeak <- rep(NA, nrow(ocr_gene))
ocr_gene$overlap_ER_Biola_baitName <- rep(NA, nrow(ocr_gene))
for(i in 1:length(enhancer_overlap)){
  bb <- Bdata_hg18_ER_gr[enhancer_overlap@to[i]]
  peak <- paste(chrom(bb),start(bb),end(bb),sep="_")
  gene <- bb$gene
  if(!is.na(ocr_gene$overlap_ER_Biola_oePeak[enhancer_overlap@from[i]])){
    peak <- unique(c(ocr_gene$overlap_ER_Biola_oePeak[enhancer_overlap@from[i]],peak))
    gene <- unique(c(unlist(strsplit(ocr_gene$overlap_ER_Biola_baitName[enhancer_overlap@from[i]],";")),gene))
    if(length(peak)>0){peak <- paste(peak,collapse = ";")}
    if(length(gene)>0){gene <- paste(gene, collapse = ";")}
  }
  ocr_gene$overlap_ER_Biola_oePeak[enhancer_overlap@from[i]] <- peak
  ocr_gene$overlap_ER_Biola_baitName[enhancer_overlap@from[i]] <- gene
}

# Are they overlapping in the same regions?
for(i in 1:nrow(ocr_gene)){
  ocr_gene$overlap_ER_Biola_ourgene_in_baitName[i] <- ocr_gene$hgnc_symbol[i]%in%unlist(strsplit(ocr_gene$overlap_ER_Biola_baitName[i],";"))                                                            
}

## Results:
# Overlapping promoter regions: 
summary(as.factor(ocr_gene$overlap_ER_Biola)) #21,476
table(ocr_gene$overlap_ER_Biola,ocr_gene$peak_type) #15,542 overlapping with enhancers and 5,934 with promoters
table(ocr_gene$overlap_ER_Biola_ourgene_in_baitName,ocr_gene$peak_type) #62 promoters annotated to same genes and 2,099 enhancers annotated to same promoter gene
#-----------------------------------


##>>>>>>>>>>>>>>>>>>>>>  Save comparison results 
write.table(ocr_gene,paste(format(Sys.Date(),"%Y%m%d"),"_cCREs_PCHiC_links_validated.txt",sep=""),sep="\t",row.names = FALSE, quote=FALSE)
#-----------------------------------

