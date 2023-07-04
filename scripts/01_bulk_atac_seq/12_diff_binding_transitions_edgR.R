#**********************************************#
# ATAC-Seq pipeline - PAIR-END                 #
# script: 12_diff_binding_transitions_edgR.sh  #
# script: 12_diff_binding_transitions_edgR.R   #
#**********************************************#

# 12. Differential Binding Analysis

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
library(edgeR)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(DESeq2)


## Count table
counts <- read.table(filename, sep="\t", dec=".", header=TRUE, check.names=FALSE)

## SAMPLE TO REMOVE:
sample_to_remove <- "ImmatureB_ATAC_D215_S5"

#remove outlier samples from the analysis
counts <- counts[,!colnames(counts)%in%sample_to_remove]


# Metadata
metadata <- read.table(meta_file, sep=";", header=TRUE)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[colnames(counts),]
metadata$Tissue <- droplevels(metadata$Tissue)


# Relevant covariables
Cell <- factor(metadata$Tissue)
Donor <- factor(metadata$Factor)

### Differential binding using edgeR

# We are interested in determine differential binding in B-cell differentiation

y <- DGEList(counts=counts)

## Filtering
# low counts across all libraries provide little evidence for differential binding
# criteria: count of 5-10 in a library to be considered as a peak in that library
# Roughly speaking, the strategy keeps genes that have at least ‘min.count’ reads in a
# worthwhile number samples. More precisely, the filtering keeps genes that have count-per-million (CPM) above _k_ in _n_ samples, where _k_ is determined by ‘min.count’ and by the sample library sizes and _n_ is determined by the design matrix.
# _n_ is essentially the smallest group sample size or, more precisely, the minimum inverse leverage for any fitted value. If all the group sizes are large, then this is relaxed slightly, but with _n_ always greater than 70% of the smallest group size.

keep <- filterByExpr(y, group=Cell, min.count=10, min.total.count=15) 
table(keep)
# only 171 peaks that will be filtered out 

# keep de peaks that pass the filter and recalculate the library size.
y <- y[keep, , keep.lib.sizes=FALSE]


## Normalization
# The most obvious technical factor that affects the read counts, other than gene expression
# levels, is the sequencing depth of each RNA sample. edgeR adjusts any differential expression
# analysis for varying sequencing depths as represented by differing library sizes.
# Following Buenrostor paper, data is normalized using the TMM function in EdgeR.
y <- calcNormFactors(y, method="TMM")


# Filter peaks by coefficient of variation
tmm_data <- cpm(y)

cv_tmm <- apply(tmm_data,1,FUN=function(x){a <- sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100 ; return(a)})

cat("Coefficent of varitation\n")
print(summary(cv_tmm)) #1q=45.68
cat("\n")

cv_cutoff <- quantile(cv_tmm, 0.25)

keep <- cv_tmm>cv_cutoff
table(keep)
#FALSE  TRUE
#26340 79019

y2 <- y[keep, , keep.lib.sizes=FALSE]


## Exploration
plotMDS(y2, col=as.numeric(Cell))


## Differential binding
# Design matrix

design <- model.matrix(~0+Cell+Donor, data=y2$samples)
colnames(design) <- gsub("Cell","",colnames(design))
colnames(design) <- gsub("Donor","",colnames(design))

# Estimate dispersion
y3 <- estimateDisp(y2,design)
plotBCV(y3)


# Determining a peaks that are differently binded in any of the states
fit <- glmQLFit(y3, design)
qlf <- glmQLFTest(fit, coef=1:8)
total_diffBind <- topTags(qlf, n=nrow(cpm(y3)))
##Sig peaks
sum(total_diffBind@.Data[[1]]$FDR<0.01) #78,920 (of 79,019; 99.8%)

##List of significant peaks
stats <- total_diffBind@.Data[[1]]
sig_peaks <- rownames(stats)[total_diffBind@.Data[[1]]$FDR<0.01]


##Differential analysis by transitions
##Make contrasts
contr.matrix <- makeContrasts(
  HSC_CLP = CLP - HSC,
  CLP_ProB = proB - CLP,
  ProB_PreB = preB - proB,
  PreB_Immature.B = ImmatureB - preB,
  Immature.B_Transitional.B = Transitional_B - ImmatureB,
  Transitional.B_Naive.CD5neg = Naive_CD5neg - Transitional_B,
  Transitional.B_Naive.CD5pos = Naive_CD5pos - Transitional_B,
  Naive.CD5neg_Naive.CD5pos = Naive_CD5pos - Naive_CD5neg,
  HSC_Transitional.B = Transitional_B - HSC,
  levels = colnames(design))


qlf_contrasts_list <- vector(mode = "list", length = length(colnames(contr.matrix)))
names(qlf_contrasts_list) <- colnames(contr.matrix)

for(i in 1:length(colnames(contr.matrix))){
  contrast_output <- glmQLFTest(fit, contrast = contr.matrix[,colnames(contr.matrix)[i]])
  qlf_contrasts_list[[i]] <- topTags(contrast_output, n=nrow(contrast_output), sort.by = "none")
}



all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(y3$counts))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(y3$counts))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(y3$counts))
for(i in 1:ncol(contr.matrix)){
  target_stats <- qlf_contrasts_list[[i]]
  all_padj[,i] <- target_stats$table$FDR
  all_fch[,i] <- target_stats$table$logFC
  sel_pos <- intersect(which(all_padj[,i]<0.01),which(all_fch[,i]>0))
  sig[sel_pos,i] <- 1  
  sel_neg <- intersect(which(all_padj[,i]<0.01),which(all_fch[,i]< 0))
  sig[sel_neg,i] <- -1
}
colnames(all_padj) <- paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) <- paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) <- paste("sig_",colnames(contr.matrix),sep="")

dea_results <- data.frame(sig,all_padj,all_fch)
rownames(dea_results) <-  rownames(qlf_contrasts_list[[1]]$table)

sum(apply(sig[,-9],1,FUN=function(x){return(sum(x!=0))})>0) 
#with a adj.pval<0.01 --> 78.248 (99%)

dea_sig <- sig[apply(sig[,-9],1,FUN=function(x){return(sum(x!=0))})>0,]
rownames(dea_sig) <-  rownames(qlf_contrasts_list[[1]]$table)[apply(sig[,-9],1,FUN=function(x){return(sum(x!=0))})>0]
#78248 x 9


##Save the dea results
write.table(dea_results,paste(format(Sys.time(), "%Y%m%d"),"_dea_transitions.txt",sep=""), sep="\t", quote=FALSE)
write.table(dea_sig,paste(format(Sys.time(), "%Y%m%d"),"_sig_dea_transitions.txt",sep=""), sep="\t", quote=FALSE)
