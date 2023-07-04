#***********************************#
# RNA-Seq pipeline - SINGLE-END     #
# script: 04.2_merge_count_table.sh #
# script: 04.2_merge_count_table.R  #
#***********************************#

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
setwd(indir)

Allfiles <- grep(".txt",list.files(),value=TRUE)

for (i in 1:length(Allfiles)){
  file <- Allfiles[i]
  if(i==1){
    df = read.table(file)
    colnames(df) <- c("genes", gsub("_nonempty_count_table.txt", "", file))
  }
  if(i>1){
    df_temp = read.table(file)
    colnames(df_temp) <- c("genes", gsub("_nonempty_count_table.txt", "", file))
    df <- merge(df, df_temp, by="genes", all=TRUE)
  }
}

# Setting output directory
setwd(outdir)

# Saving results
save(df, file=paste(filename,".rda",sep=""))
write.table(df, file=paste(format(Sys.time(), "%Y%m%d"),"_",filename,".txt",sep=""),sep="\t", dec=".", quote=FALSE)