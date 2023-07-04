#********************************************#
# FOOTPRINT IDENTIFICATION                   #
# script: 14.3_unify_footprint.sh            #
# script: 14.3_unify_footprint.R             #
#********************************************#


# 14.3 Merge footprints from HINT and Wellington outputs

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




################################
#  UNIFY FOOTPRINT BY SAMPLES  #
################################
# Footprints from HINT and Wellington algorithm

#https://ro-che.info/articles/2018-07-11-chip-seq-consensus
library(rtracklayer)
library(GenomicRanges)

#Get the list of .narrowPeak files in the current directory. We'll treat them all as replicates.
file_names <- read.table(filename, sep=";", header=TRUE)



##>>>>>>>>>>>>>>>>>>>>>  Get the union of peaks from both strategies

for(i in 1:nrow(file_names)){
  peak_granges <- lapply(list(as.character(file_names$footprint_hint[i]), 
  as.character(file_names$footprint_wellington[i])), import, format = "BED")
  all_peak_grangeslist <- GRangesList(peak_granges)
  #Now, get the coverage
  peak_coverage <- coverage(all_peak_grangeslist)

  #peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. 
  #Now, using the slice function, let's get the regions in the genome where coverage is at least in 1 of the inputs (footprint strategies):
  covered_ranges <- slice(peak_coverage, lower=1, rangesOnly=TRUE)

  #This is a simple IRangesList object; let's covert it to a GRanges object:
  covered_granges <- GRanges(covered_ranges)
  cat("************************\n",summary(covered_granges),"\n************************\n")

  export(covered_granges, paste(file_names$SampleID[i],"unified_footprint.bed",sep="_"))
}

#-------------------------------------
