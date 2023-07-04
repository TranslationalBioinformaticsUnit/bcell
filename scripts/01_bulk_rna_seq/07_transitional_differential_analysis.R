#**************************************************#
# RNA-Seq pipeline - SINGLE-END                    #
# script: 07_transitional_differential_analysis.sh #
# script: 07_transitional_differential_analysis.R  #
#**************************************************#

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


### STEP 1: LOAD DATA

#normalized data
load(rna_data)
##dimension of "v" 23,454 x 79 (just filter by number of counts)

#metadata
metadata <- read.table(metadata,sep="\t", dec=".", header=TRUE, check.names=FALSE) #213 x 23

metadata2 <- metadata[metadata$Alias%in%colnames(norm_data),]
sum(metadata2$Alias==colnames(norm_data))


metadata2$CellType <- droplevels(metadata2$CellType)
metadata2$Donor <- droplevels(metadata2$Donor)



### STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS >> VOOM+LIMMA

#Differential expression analysis
require(edgeR)

group <- as.factor(metadata2$CellType)
donor <- as.factor(metadata2$Donor)
design <- model.matrix(~0 + group + donor)

colnames(design) <- gsub("group","",colnames(design))


v_filtered <- v

fit <- lmFit(v_filtered, design)
fit2 <- eBayes(fit)


whole_stats <- topTable(fit2, n=nrow(v$E))

sum(whole_stats$adj.P.Val<0.01)
## We identify 17,059 genes DE, 97%.
## We identify 23,905 genes DE, 97%.
sum(whole_stats$adj.P.Val<0.05)
## We identify 17,277 genes DE, 98%.
## We identify 23,134 genes DE, 98%.


##Make contrasts
contr.matrix <- makeContrasts(
  HSC_CLP = HSC - CLP,
  CLP_ProB = CLP - ProB,
  ProB_PreB = ProB - PreB,
  PreB_Immature.B = PreB - Immature.B,
  Immature.B_Transitional.B = Immature.B - Transitional.B,
  Transitional.B_Naive.CD5neg = Transitional.B - Naive.CD5neg,
  Transitional.B_Naive.CD5pos = Transitional.B - Naive.CD5pos,
  Naive.CD5neg_Naive.CD5pos = Naive.CD5neg - Naive.CD5pos,
  #HSC_Transitional.B = HSC - Transitional.B,
  levels = colnames(design))

# get results
fit.contrasts <- contrasts.fit(fit, contr.matrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(v_filtered$E), sort.by="none")
sum(contrast_stats$adj.P.Val<0.01)
## We identify 9,389 genes DE, 53.4%.
## We identify 13,717 genes DE, 58.5%.

sum(contrast_stats$adj.P.Val<0.05)
## We identify 10,934 genes DE, 62.2%.
## We identify 15,813 genes DE, 67.4%.


all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i,number=nrow(v_filtered$E), sort.by="none")
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel_pos <- intersect(which(all_padj[,i]<0.05),which(all_fch[,i]>0.58))
  sig[sel_pos,i] <- 1  
  sel_neg <- intersect(which(all_padj[,i]<0.05),which(all_fch[,i]< -0.58))
  sig[sel_neg,i] <- -1
}
colnames(all_padj) <- paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) <- paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) <- paste("sig_",colnames(contr.matrix),sep="")

dea_results <- data.frame(sig,all_padj,all_fch)
rownames(dea_results) <-  rownames(contrast_stats)

sum(apply(sig,1,FUN=function(x){return(sum(x!=0))})>0) #5,940 dea / with padj<0.05 8,009 dea
#with a adj.pval<0.01 & |FC|>1.5 (logFC>0.58)
# 11,080 with a adj.pval<0.05 & |FC|>1.5 (logFC>0.58) and no CV

write.table(dea_results,paste(format(Sys.time(), "%Y%m%d"),"_dea_results.txt",sep=""),sep="\t", dec=".", quote=FALSE)



## Significant genes from binary output
dea_sig <- sig[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0,]
rownames(dea_sig) <-  rownames(contrast_stats)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
#5940 x 8
#11,080 x 8


write.table(dea_sig,paste(format(Sys.time(), "%Y%m%d"),"_dea_significant_binary.txt",sep=""),sep="\t", dec=".", quote=FALSE)

