#******************************************************#
# Archetype - TF translation                           #
# script: 18_archetype_TF_translation.R                #
#******************************************************#


# 18. Translate Archetypes to TF based on the TFBS variability, RNA coefficient of variation and correlation



##>>>>>>>>>>>>>>>>>>>>> OCR-archetype (motif) binary matrix
ocr_fp_arch_matrix <- read.table("ocr_arch_matrix.txt", sep="\t", header=TRUE, check.names = FALSE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Read correlation and variability data
result_TFBS_RNA <- read.table("archetype_TFBS_TF_expression_correlation.txt", sep="\t", dec=".", header=TRUE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Identify Archetype-TF significantly associated
#Candidate variable:
# 1- Low TFBS varibility - low CV RNA (independent of correlation result) --> NO ASSOCIATION
# 2- Low TFBS variability - high CV RNA (Poised chromatin expression changing) (independent of correlation result) --> KEEP ASSOCIATION
# 3- High TFBS variability - low CV RNA (Chromatin driving) (independent of correlation result) --> KEEP ASSOCIATION
# 4- High TFBS variability - high CV RNA - sig positive correlation --> KEEP ASSOCIATION
# 5- High TFBS variability - high CV RNA - no sig or sig negative correlation --> NO ASSOCIATION
# Candidate associations will be: candidate=2,3,4

# Cut-off for TFBS variability and RNA coefficient or variation are defined based on the 1st quantile value.
result_TFBS_RNA$candidates <- rep(5, length=nrow(result_TFBS_RNA))

sel1 <- which(result_TFBS_RNA$archetype_variability>1.6 & result_TFBS_RNA$rna_cv>0.8956 & result_TFBS_RNA$adj.pval_fdr<0.05 & result_TFBS_RNA$rho>0)
result_TFBS_RNA$candidates[sel1] <- 4

sel2 <- which(result_TFBS_RNA$archetype_variability<1.6 & result_TFBS_RNA$rna_cv<0.8956)
result_TFBS_RNA$candidates[sel2] <- 1

sel3 <- which(result_TFBS_RNA$archetype_variability<1.6 & result_TFBS_RNA$rna_cv>=0.8956)
result_TFBS_RNA$candidates[sel3] <- 2

sel4 <- which(result_TFBS_RNA$archetype_variability>=1.6 & result_TFBS_RNA$rna_cv<0.8956)
result_TFBS_RNA$candidates[sel4] <- 3


# In the case of complexes, to consider a complex significantly associated to a TFBS archetype, we want to have all the components of the complex associated
complexes <- unique(grep("\\+",result_TFBS_RNA$target_TF,value=TRUE))
for(i in 1:length(complexes)){
  aa <- result_TFBS_RNA[result_TFBS_RNA$target_TF==complexes[i],]
  if(length(unique(aa$archetype))==1){
    sel <- which(aa$candidates==1 | aa$candidates==5)
    if(length(sel)>0){
      result_TFBS_RNA$candidates[result_TFBS_RNA$target_TF==complexes[i]] <- 5
    }
  }
  if(length(unique(aa$archetype))>1){
    arch <- unique(aa$archetype)
    for(j in 1:length(arch)){
      sel <- which(aa$candidates[aa$archetype==arch[j]]==1 | aa$candidates[aa$archetype==arch[j]]==5)
      if(length(sel)>0){
        result_TFBS_RNA$candidates[result_TFBS_RNA$target_TF==complexes[i] & result_TFBS_RNA$archetype==arch[j]] <- 5
      }
    }
  }
}
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Candidate pairs for ARCH-TF translation
result_TFBS_RNA2 <- na.omit(result_TFBS_RNA[result_TFBS_RNA$candidates%in%c(2,3,4),]) #384 x 11
write.table(result_TFBS_RNA2,paste(format(Sys.Date(),"%Y%m%d"),"_archetype_TF_translation.txt",sep=""),sep="\t",dec=",",row.names=FALSE, col.names=TRUE, quote=FALSE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> OCR-TF binary matrix based on translation strategy
tfs <- unique(result_TFBS_RNA2$target_TF) #285

#TF matrix
ocr_tf_matrix <- matrix(nrow=nrow(ocr_fp_arch_matrix), ncol=length(tfs)) #105,530 x 285
rownames(ocr_tf_matrix) <- rownames(ocr_fp_arch_matrix)
colnames(ocr_tf_matrix) <- tfs

for(i in 1:length(tfs)){
  archetypes <- as.character(unique(result_TFBS_RNA2$archetype[result_TFBS_RNA2$target_TF==tfs[i]]))
  if(length(archetypes)>1){
    tmp <- rowSums(ocr_fp_arch_matrix[,archetypes])
    tmp[tmp>1] <- 1
    ocr_tf_matrix[,i] <- tmp
  }
  if(length(archetypes)==1){
    ocr_tf_matrix[,i] <- ocr_fp_arch_matrix[,archetypes]
  }
}

write.table(ocr_tf_matrix,paste(format(Sys.Date(),"%Y%m%d"),"_ocr_tf_matrix.txt",sep=""),sep="\t",dec=",",row.names=TRUE, col.names=TRUE, quote=FALSE)
#-------------------------------------
