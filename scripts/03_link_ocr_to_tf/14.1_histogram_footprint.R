#***************************************#
# FOOTPRINT IDENTIFICATION              #
# script: 14.1_footprint_hint.sh        #
# script: 14.2_footprint_wellington.sh  #
# script: 14.1_histogram_footprint.R    #
#***************************************#


# 14.1 Histogram exploration of number of counts per footprint


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


########################################
# HISTOGRAM OF FOOTPRINTS BY CELL TYPE #
########################################


#Get the list of .narrowPeak files in the current directory. We'll treat them all as replicates.
file_names <- read.table(filename, sep=";", header=TRUE)

#Get the different cell types
cell_types <- levels(file_names$Tissue)


for(i in 1:length(cell_types)){
  list_files <- grep("footprint.bed",dir(), value=TRUE)
  samples <- file_names$SampleID[file_names$Tissue==cell_types[i]]
  list_files2 <- sapply(samples, grep, x=list_files, value=TRUE)
  stats <- vector()
  for(j in 1:length(list_files2)){
   a <- read.table(list_files2[j])
   stats[j]<-summary(a[,5])[2]
  }

  png(paste(format(Sys.Date(),"%Y%m%d"),"hist_1Q_reads_footprint_",cell_types[i],".png",sep=""))
    hist(stats)
  dev.off()
}

