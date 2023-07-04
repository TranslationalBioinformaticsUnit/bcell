#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 11_peak_annotation.sh     #
# script: 11_peak_annotation.R      #
#***********************************#

# 11. Genomic mapping/annotation of OCRs


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}


# Setting working directory
setwd(workdir)

# Library required
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#>>>>>>>>>>>>>>>>> Create the GRanges of consensus peaks

#compute the mean coverage
count_table <- read.table(count_matrix, sep="\t", header=TRUE)
mean_cov <- rowMeans(count_table)

# Read consensus peaks
peaks <- read.table(peaks_file) #105,530
colnames(peaks) <- c("chr","start","end","","")
peaksID <- paste(peaks[,1],peaks[,2], peaks[,3], sep="_")

#compute the mean coverage
count_table <- read.table(count_matrix, sep="\t", header=TRUE)
mean_cov <- rowMeans(count_table)

peaks$coverage <- mean_cov

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks[,-c(4:6)], keep.extra.columns=TRUE)


##Coverage plot
png("peaks_over_chr.png")
covplot(peak_granges, title = "ATAC-Seq Peaks over Chromosomes", weightCol="coverage")
dev.off()


# Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak_granges, windows=promoter)

# Heatmap of ChIP binding to TSS regions
png("binding_TSS.png")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
dev.off()


# Average Profile of ChIP peaks binding to TSS region
png("average_binding_TSS.png")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()


# Peak Annotation
peakAnno <- annotatePeak(peak_granges, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db", level="gene")

peak_annotation <- as.data.frame(peakAnno) #105,530 x 18

peak_annotation$annotation2 <- peak_annotation$annotation
peak_annotation$annotation2[grep("Promoter",peak_annotation$annotation)] <- "Promoter"
peak_annotation$annotation2[grep("Intron",peak_annotation$annotation)] <- "Intron"
peak_annotation$annotation2[grep("Distal Intergenic",peak_annotation$annotation)] <- "Distal Intergenic"
peak_annotation$annotation2[grep("Exon",peak_annotation$annotation)] <- "Exon"


write.csv(peak_annotation, "consensus_peaks_annotation_chipseeker_1Mb_promoter.txt")
## There are 2,082 peaks without annotation (ENSMBL ID is missing; NA)

# Visualize Genomic Annotation
png("pie_annotation.png")
plotAnnoPie(peakAnno)
dev.off()
