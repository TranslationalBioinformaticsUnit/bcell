#********************************************#
# FOOTPRINT IDENTIFICATION                   #
# script: 14.3_unify_footprint.sh            #
# script: 14.4_consensus_unify_footprint.R   #
#********************************************#


# 14.4 Get the consensus footprint for each cell type


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


#Setting working directory
setwd(workdir)

################################
#  CONSENSUS FP BY CELL TYPE   #
################################

#https://ro-che.info/articles/2018-07-11-chip-seq-consensus
library(rtracklayer)
library(GenomicRanges)

#Get the list of .narrowPeak files in the current directory. We’ll treat them all as replicates.
file_names <- read.table(filename, sep=";", header=TRUE)


#Get the different cell types
cell_types <- levels(file_names$Tissue)



##>>>>>>>>>>>>>>>>>>>>>  Get the whole consensus footprinting for each cell type.

#list of "consensus" bed files to get the whole consensus peaks
file_names_bed <- grep("_unified_footprint.bed",dir(), value=TRUE) 

for(i in 1:length(cell_types)){
  samples <- file_names$SampleID[file_names$Tissue==cell_types[i]]
  file_names_bed2 <- sapply(samples, grep, x=file_names_bed, value=TRUE)
  peak_granges <- lapply(as.list(file_names_bed2), import, format = "BED")
  all_peak_grangeslist <- GRangesList(peak_granges)
  #Now, get the coverage
  peak_coverage <- coverage(all_peak_grangeslist)

  #peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. 
  #Now, using the slice function, let’s get the regions in the genome where coverage is at least in n-1 sample:
  covered_ranges <- slice(peak_coverage, lower=length(all_peak_grangeslist)-1, rangesOnly=TRUE)

  #This is a simple IRangesList object; let's covert it to a GRanges object:
  covered_granges <- GRanges(covered_ranges)
  cat("************************\n",summary(covered_granges),"\n************************\n")

  export(covered_granges, paste(cell_types[i],"consensus_unified_footprint.bed",sep="_"))
}

#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>>  Get the whole consensus footprinting for all cell type.

file_names_bed <- grep("consensus_unified_footprint.bed",dir(), value=TRUE) 

peak_granges <- lapply(as.list(file_names_bed), import, format = "BED")
all_peak_grangeslist <- GRangesList(peak_granges)
#Now, get the coverage
peak_coverage <- coverage(all_peak_grangeslist)

#peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. Now, using the slice function, let’s get the regions in the genome where coverage is at least in n-1 sample:
covered_ranges <- slice(peak_coverage, lower=1, rangesOnly=TRUE)

#This is a simple IRangesList object; let’s covert it to a GRanges object:
covered_granges <- GRanges(covered_ranges)
cat("************************\n",summary(covered_granges),"\n************************\n")

export(covered_granges, "all_cell_types_consensus_unified_footprint.bed")
#-------------------------------------
