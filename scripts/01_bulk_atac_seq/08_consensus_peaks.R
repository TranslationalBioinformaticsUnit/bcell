#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 08_consensus_peaks.sh     #
# script: 08_consensus_peaks.R      #
#***********************************#

# 8. Consensus peaks for each cell type using GenomicRanges and rtracklayer
#https://ro-che.info/articles/2018-07-11-chip-seq-consensus


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


print(workdir)
print(filename)

#Setting working directory
setwd(workdir)

################################
# CONSENSUS PEAKS BY CELL TYPE #
################################

#https://ro-che.info/articles/2018-07-11-chip-seq-consensus
library(rtracklayer)
library(GenomicRanges)

#Get the list of .narrowPeak files in the current directory. We'll treat them all as replicates.
file_names <- read.table(filename, sep=";", header=TRUE)

#Get the general stats from BAM files (million reads per sample)
file_stats <- read.table(filestats, skip=1)
colnames(file_stats) <- c("SampleID","M_reads")

#Merge both data.frames
file_names <- merge(file_names, file_stats, by.x="SampleID", by.y="SampleID", all.x=TRUE)

#Get the different cell types
cell_types <- levels(file_names$Tissue)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

#Check the sequence deep for each cell type
aiqr <- by(file_names$M_reads, file_names$Tissue, IQR)
amedian <- by(file_names$M_reads, file_names$Tissue, median)
q25 <- by(file_names$M_reads, file_names$Tissue, FUN=function(x){quantile(x,0.25)})
q75 <- by(file_names$M_reads, file_names$Tissue, FUN=function(x){quantile(x,0.75)})
dwlim <- q25 - 1.5*aiqr
uplim <- q75 + 1.5*aiqr

dwlim <- q25 - aiqr
uplim <- q75 + aiqr

asd <- by(file_names$M_reads, file_names$Tissue, sd)
amean <- by(file_names$M_reads, file_names$Tissue, mean)
dwlim <- amean - 2*asd
uplim <- amean + 2*asd



outliers_dw <- vector()
for(i in 1:length(cell_types)){
  d <- file_names[file_names$Tissue==cell_types[i],]
  outliers_dw <- c(outliers_dw, as.character(d$SampleID[d$M_reads<dwlim[cell_types[i]]]))
}

outliers_up <- vector()
for(i in 1:length(cell_types)){
  d <- file_names[file_names$Tissue==cell_types[i],]
  outliers_up <- c(outliers_up, as.character(d$SampleID[d$M_reads>uplim[cell_types[i]]]))
}


metadata_consensus <- data.frame(SampleID=cell_types, Group=cell_types, Peaks=paste(workdir,"/consensus_",cell_types,".bed", sep=""), PeakCaller=rep("bed",length(cell_types)),ScoreCol=rep(5,length(cell_types)))
write.table(metadata_consensus,"metadata_consensus.csv", sep=",", quote=FALSE, row.names=FALSE)

pdf("peaks_exploration.pdf")
for (i in 1:length(cell_types)){

cat("Cell type: ",cell_types[i],"\n")

peak_granges <- lapply(as.list(as.character(file_names$Peaks[file_names$Tissue==cell_types[i]])), import, format = "BED", extraCols = extraCols_narrowPeak)

names(peak_granges) <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(file_names$Peaks[file_names$Tissue==cell_types[i]]),"/"),FUN=function(x){return(x[6])})),"\\."),FUN=function(x){return(x[1])}))

#Convert to a GRangesList
peak_grangeslist <- GRangesList(peak_granges)
cat("Number of GRanges objects: ",length(peak_grangeslist),"\n")


#barplot representation of the number of peaks for each case within each cell type 
barplot(unlist(lapply(peak_grangeslist, length)), xlab="Cases", ylab="Number of peaks", main=paste("Peaks in ", cell_types[i], sep=""), cex.names=0.7)

#density distribution of peak width for each sample
stats <- matrix(ncol=6, nrow=length(peak_grangeslist))
colnames(stats) <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
rownames(stats) <- names(peak_grangeslist)

for(j in 1:length(peak_grangeslist)){
  if(length(ranges(peak_grangeslist[[j]])@width)!=0){
    hist(ranges(peak_grangeslist[[j]])@width, main=paste("Peak width distribution",names(peak_grangeslist[j]),sep="\n"), xlab="Peak width (bp)")
  }
  stats[j,] <- as.numeric(summary(ranges(peak_grangeslist[[j]])@width))
}

write.table(stats, paste("width_stats_",cell_types[i],".txt", sep=""), sep="\t", dec=".", quote=FALSE)


## Start getting the consensus
#Now, find the genome regions which are covered by at least n-1 of the sets of peaks:
peak_coverage <- coverage(peak_grangeslist)

#peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. Now, using the slice function, let’s get the regions in the genome where coverage is at least n-1:
covered_ranges <- slice(peak_coverage, lower=length(peak_grangeslist)-1, rangesOnly=TRUE)

#This is a simple IRangesList object; let’s covert it to a GRanges object:
covered_granges <- GRanges(covered_ranges)
cat("************************\n",summary(covered_granges),"\n************************\n")

#Merging nearby peaks
#you might want to merge the peaks that are close to each other to avoid having many small fragments. In Bioconductor, you can do that using the min.gapwidth parameters to reduce. For example, this will merge the peaks that are separated by no more than 30 bp
red_covered_granges <- reduce(covered_granges, min.gapwidth=31)
cat("************************\n",summary(red_covered_granges),"\n************************\n")

export(red_covered_granges, paste("consensus_",cell_types[i],".bed", sep=""))
}

dev.off()


######################
# Get the whole consensus peaks for all cell types.
######################

#list of "consensus" bed files to get the whole consensus peaks
file_names_bed <- grep(".bed",dir(), value=TRUE)
#for single-cell analysis
#file_names_bed <- grep(".bed",dir(), value=TRUE)[c(5,6,11)]

peak_granges <- lapply(as.list(file_names_bed), import, format = "BED")

all_peak_grangeslist <- GRangesList(peak_granges)
#Includes a total of 181 GRanges objects, corresponding to the 18 different cell types (3 cell types removed CMP, GMP and MEP)


#Now, get the coverage
peak_coverage <- coverage(all_peak_grangeslist)

#peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. Now, using the slice function, let’s get the regions in the genome where coverage is at least in 1 sample:
covered_ranges <- slice(peak_coverage, lower=1, rangesOnly=TRUE)

#This is a simple IRangesList object; let’s covert it to a GRanges object:
covered_granges <- GRanges(covered_ranges)
cat("************************\n",summary(covered_granges),"\n************************\n")

#Merging nearby peaks
#you might want to merge the peaks that are close to each other to avoid having many small fragments. In Bioconductor, you can do that using the min.gapwidth parameters to reduce. For example, this will merge the peaks that are separated by no more than 30 bp
red_covered_granges <- reduce(covered_granges, min.gapwidth=31)
cat("************************\n",summary(red_covered_granges),"\n************************\n")


export(red_covered_granges, "consensus_cell_types.bed")
# We have a total of 83,614 intervals (peaks)

#for singel-cell comparison
#export(red_covered_granges, "consensus_3cell_types.bed")