#***************************************#
# RNA-Seq pipeline - SINGLE-END         #
# script: 06_differential_expression.sh #
# script: 06_differential_expression.R  #
#***************************************#

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

#metadata
metadata <- read.table(metadata,sep="\t", dec=".", header=TRUE, check.names=FALSE) #213 x 23
metadata2 <- metadata[metadata$Alias%in%colnames(norm_data),]
sum(metadata2$Alias==colnames(norm_data))

metadata2$CellType <- droplevels(metadata2$CellType)
metadata2$Donor <- droplevels(metadata2$Donor)
#--------------------------------------------------------

### STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS >> VOOM+LIMMA

#Differential expression analysis
require(edgeR)

group <- metadata2$CellType
donor <- metadata2$Donor
design <- model.matrix(~0 + group + donor)

colnames(design) <- gsub("group","",colnames(design))


# Filter peaks by coefficient of variation

##Set all values to positive
voom_data <- v$E+100

cv_voom <- apply(voom_data,1,FUN=function(x){a <- sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100 ; return(a)})

cat("Coefficent of varitation\n")
print(summary(cv_voom)) #1q=0.89
cat("\n")

cv_cutoff <- quantile(cv_voom, 0.25)

keep <- cv_voom>cv_cutoff
table(keep)
#FALSE  TRUE
# 5864 17590


v_filtered <- v
v_filtered$E <- v_filtered$E[keep,]
v_filtered[[3]] <- v_filtered[[3]][keep,]


write.table(v_filtered$E,"_norm_filtered_CV_voom_data.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

fit <- lmFit(v_filtered, design)
fit2 <- eBayes(fit)


whole_stats <- topTable(fit2, n=nrow(v$E))
sum(whole_stats$adj.P.Val<0.01)
## We identify 17,059 genes DE, 97%.


## Make contrasts
contr.matrix <- makeContrasts(
  HSC_CLP = CLP - HSC,
  HSC_ProB = ProB - HSC,
  HSC_PreB = PreB - HSC,
  HSC_Immature.B = Immature.B - HSC,
  HSC_Transitional.B = Transitional.B - HSC,
  HSC_Naive.CD5neg = Naive.CD5neg - HSC,
  HSC_Naive.CD5pos = Naive.CD5pos - HSC,
  CLP_ProB = ProB - CLP,
  CLP_PreB = PreB - CLP,
  CLP_Immature.B = Immature.B - CLP,
  CLP_Transitional.B = Transitional.B - CLP,
  CLP_Naive.CD5neg = Naive.CD5neg - CLP,
  CLP_Naive.CD5pos = Naive.CD5pos - CLP,
  ProB_PreB = PreB - ProB,
  ProB_Immature.B = Immature.B - ProB,
  ProB_Transitional.B = Transitional.B - ProB,
  ProB_Naive.CD5neg = Naive.CD5neg - ProB,
  ProB_Naive.CD5pos = Naive.CD5pos - ProB,
  PreB_Immature.B = Immature.B - PreB,
  PreB_Transitional.B = Transitional.B - PreB,
  PreB_Naive.CD5neg = Naive.CD5neg - PreB,
  PreB_Naive.CD5pos = Naive.CD5pos - PreB,
  Immature.B_Transitional.B = Transitional.B - Immature.B,
  Immature.B_Naive.CD5neg = Naive.CD5neg - Immature.B,
  Immature.B_Naive.CD5pos = Naive.CD5pos - Immature.B,
  Transitional.B_Naive.CD5neg = Naive.CD5neg - Transitional.B,
  Transitional.B_Naive.CD5pos = Naive.CD5pos - Transitional.B,
  Naive.CD5neg_Naive.CD5pos = Naive.CD5pos - Naive.CD5neg,
  levels = colnames(design))


## Get results
fit.contrasts <- contrasts.fit(fit, contr.matrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(v_filtered$E), sort.by="none")
sum(contrast_stats$adj.P.Val<0.01)


all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(v_filtered$E))
for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i,number=nrow(v_filtered$E), sort.by="none")
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i]<0.01),which(all_fch[,i]>0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i]<0.01),which(all_fch[,i]< -0.58))
  sig[sel2,i] <- -1  
}
colnames(all_padj) <- paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) <- paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) <- paste("sig_",colnames(contr.matrix),sep="")

dea_results <- data.frame(sig,all_padj,all_fch)
rownames(dea_results) <-  rownames(contrast_stats)

sum(rowSums(sig)!=0) 
#with a adj.pval<0.01 & |FC|>1.5 (logFC>0.58)

write.table(dea_results,paste(format(Sys.time(), "%Y%m%d"),"_dea_results_one_vs_all.txt",sep=""),sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(rownames(dea_results)[rowSums(sig)>0],"list_deg.txt",sep="\t",dec=".",quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(contrast_stats)[contrast_stats$adj.P.Val<0.01],"list_deg_anova.txt",sep="\t",dec=".",quote=FALSE, row.names=FALSE, col.names=FALSE)
