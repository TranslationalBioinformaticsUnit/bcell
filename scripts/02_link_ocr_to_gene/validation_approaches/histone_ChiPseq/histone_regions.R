#***********************************************#
# LINK OCRs TO GENES                            #
# Validation approach with public ChiP-seq data #
# script: histone_regions.R                     #
#***********************************************#


# Get the histone regions from  encode ChiP data associated to CD34+ cells, naive B cells and B cells (https://www.encodeproject.org/).
# The closest data to our study population


# Required libraries
library(GenomicRanges)
library(rtracklayer)


#After download the corresponding BED files

##>>>>>>>>>>>>>>>>>>>>>  Input files
histone_files <- grep(".bed",dir(), value=TRUE)
histone_metadata <- grep(".tsv",dir(), value=TRUE)
# For example: for h3k4me1 --> experiment_report_2021_9_21_9h_42m.tsv and experiment_report_2021_9_21_9h_43m.tsv
#-----------------------------------


##>>>>>>>>>>>>>>>>>>>>>  Compute the consensus ranges within each cell type
consensus_granges <- vector(mode="list", length=length(histone_metadata))

for(i in 1:length(histone_metadata)){
  #Read metadata
  metadata <- read.table(histone_metadata[i], sep="\t", dec=",", skip=1, header=TRUE)
  target_files <- vector()
  for(j in 1:length(histone_files)){
    target_files <- c(target_files, grep(gsub(".bed","",histone_files[j]),metadata$Files))
  }
  files <- histone_files[target_files]
  
  # Read BED files in a GRanges format
  peak_granges <- vector(mode="list", length=length(files))
  
  for(k in 1:length(files)){
    data <- read.table(files[k])
    names(data)[1:4] <- c("chr","start","end","name")
    peak_granges[[k]] <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE) 
  }
  
  names(peak_granges) <- gsub(".bed","",files)
  
  peak_grangeslist <- GRangesList(peak_granges)
  cat("Number of GRanges objects: ",length(peak_grangeslist),"\n")

  # Barplot representation of the number of peaks for each case within each cell type 
  barplot(unlist(lapply(peak_grangeslist, length)), xlab="Cases", ylab="Number of peaks", main="Number of peaks per samples", cex.names=0.7)
  
  #Now, get the coverage
  peak_coverage <- coverage(peak_grangeslist)
  
  #peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. Now, using the slice function, lets get the regions in the genome where coverage is at least in n-1 sample:
  covered_ranges <- slice(peak_coverage, lower=length(files)-1, rangesOnly=TRUE)
  
  #This is a simple IRangesList object; lets covert it to a GRanges object:
  covered_granges <- GRanges(covered_ranges)
  cat("************************\n",summary(covered_granges),"\n************************\n")
  
  #Merging nearby peaks
  #you might want to merge the peaks that are close to each other to avoid having many small fragments. In Bioconductor, you can do that using the min.gapwidth parameters to reduce. For example, this will merge the peaks that are separated by no more than 30 bp
  red_covered_granges <- reduce(covered_granges, min.gapwidth=31)
  cat("************************\n",summary(red_covered_granges),"\n************************\n")
  
  consensus_granges[[i]] <- red_covered_granges
  names(consensus_granges)[i] <- unique(metadata$Biosample.term.name)
}

save(list=c("consensus_granges"), file="histone_by_cell_type.RData")

# Barplot representation 
barplot(unlist(lapply(consensus_granges, length)), xlab="Cases", ylab="Number of peaks", main="Number of peaks per samples", cex.names=0.7)
#-----------------------------------


##>>>>>>>>>>>>>>>>>>>>>  Compute the union of different consensus ranges
#Get the coverage
peak_coverage <- coverage(GRangesList(consensus_granges))

#peak_coverage is a list of numbers for every single genome position, encoded into an efficient structure Rle. Now, using the slice function, lets get the regions in the genome where coverage is at least in 1 sample:
covered_ranges <- slice(peak_coverage, lower=1, rangesOnly=TRUE)

#This is a simple IRangesList object; lets covert it to a GRanges object:
covered_granges <- GRanges(covered_ranges)
cat("************************\n",summary(covered_granges),"\n************************\n")

#Merging nearby peaks
#you might want to merge the peaks that are close to each other to avoid having many small fragments. In Bioconductor, you can do that using the min.gapwidth parameters to reduce. For example, this will merge the peaks that are separated by no more than 30 bp
red_covered_granges <- reduce(covered_granges, min.gapwidth=31)
cat("************************\n",summary(red_covered_granges),"\n************************\n")

#save results
save(list=c("red_covered_granges"), file="histone.RData")
export(red_covered_granges,con="histone.bed", format="bed")
#-----------------------------------