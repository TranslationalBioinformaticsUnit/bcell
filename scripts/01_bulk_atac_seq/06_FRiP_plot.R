#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 06_FRiP_plot.sh           #
# script: 06_FRiP_plot.R            #
#***********************************#

# 6. Quality Control of Peak Calling - Fraction of Reads in Peaks (FRiP) PLOT

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
print(file)



#Setting working directory
setwd(workdir)

data <- read.table(file, header=FALSE)[,-c(6,8)]
colnames(data) <- c("SampleID","total_reads","reads_in_peaks","reads_in_filtered_peaks","peaks_filtered","peaks")

data$frip <- data$reads_in_peaks/data$total_reads
data$frip_filtered <- data$reads_in_filtered_peaks/data$total_reads


pdf("frip_distribution.pdf")
hist(data$frip, xlab="Fraction of reads in Peaks", main="FRiP distribution")
hist(data$frip_filtered, xlab="Fraction of reads in Filtered Peaks", main="Filtered FRiP distribution")
dev.off()

pdf("peaks_distribution.pdf")
hist(data$peaks, xlab="Number of Peaks", main="Peaks distribution")
hist(data$peaks_filtered, xlab="Number of Filtered Peaks", main="Filtered Peaks distribution")
dev.off()

write.table(data, "frip_output.txt", sep="\t", dec=".", row.names=FALSE, quote=FALSE, col.names=TRUE)