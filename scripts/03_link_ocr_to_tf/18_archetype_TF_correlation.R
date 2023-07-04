#******************************************************#
# Archetype - TF correlation                           #
# script: 18_archetype_TF_correlation.R                #
#******************************************************#


# 18. Assess the correlation between archetype TFBS and associated TF expression



# Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVAR)
library(TFBSTools)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pals)
library(BuenColors)
library(circlize)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load internal function
source("chromVAR_internal_functions.R")



##>>>>>>>>>>>>>>>>>>>>> OCR-archetype (motif) binary matrix
ocr_fp_arch_matrix <- read.table("ocr_arch_matrix.txt", sep="\t", header=TRUE, check.names = FALSE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Aggregated archetype score data
# Read the TFBS deviation scores for each archetype
archetype_TFBS_deviations <- read.table("archetype_TFBS_deviations.txt", sep="\t", header=TRUE, dec=".", check.names = FALSE)
archetype_TFBS_score <- read.table("archetype_TFBS_score.txt", sep="\t", header=TRUE, dec=".", check.names = FALSE)

# Read the TFBS variability for each archetype
archetype_TFBS_variability <- read.table("archetype_TFBS_variability.txt", sep="\t", header=TRUE, dec=".")
hist(archetype_TFBS_variability$variability)
quantile(archetype_TFBS_variability$variability,0.25) #cut-off: 1.6

# Get archetype/motif annotation
TF_family <- read.table("vierstra_motif_annotations.txt", header=TRUE, sep="\t")
TF_cluster <- read.table("vierstra_motif_annotations_archetypes_clusters.txt", header=TRUE, sep="\t")
TF_cluster <- TF_cluster[order(TF_cluster$Cluster_ID),]
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> ATAC counts and metadata
# counts
counts_file <- "count_data_atac.txt"
my_counts_matrix <- as.matrix(read.delim(counts_file, check.names = FALSE))
#remove outlier sample
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix)=="ImmatureB_ATAC_D215_S5")]

# metadata
atac_metadata <- read.table("metadata_atac.txt", sep="\t", header=TRUE)
rownames(atac_metadata) <- atac_metadata$SampleID
atac_metadata2 <- atac_metadata[colnames(my_counts_matrix),]
#order cell type labels
atac_metadata2$Tissue <- factor(as.character(atac_metadata2$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> RNA data and metadata
# Read the RNA expression of TF genes
rna <- read.table("voom_normalized_filtered_data_rna.txt", sep="\t", header=TRUE, dec=".", check.names = FALSE)
rna_metadata <- read.table("metadata_rna.txt", sep="\t", header=TRUE, dec=".")
rownames(rna_metadata) <- rna_metadata$Sample.Name
rna <- rna[,as.character(rna_metadata$Alias)] 
colnames(rna) <- rownames(rna_metadata) 

# Get gene annotation
gene_anno <- read.table("ensmbl_gene_annotation.txt", sep="\t", header=TRUE, dec=".")
rownames(gene_anno) <- gene_anno$ensembl_gene_id
rna_gene_anno <- gene_anno[rownames(rna),]
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Common samples between RNA and ATAC
rownames(atac_metadata2) <- paste(atac_metadata2$Factor,atac_metadata2$Tissue,sep="_")
colnames(archetype_TFBS_score) <- rownames(atac_metadata2)
colnames(archetype_TFBS_deviations) <- rownames(atac_metadata2)

ids_atac_rna <- intersect(colnames(archetype_TFBS_deviations), colnames(rna)) 
length(ids_atac_rna) #59 share samples

rna2 <- as.matrix(rna[,ids_atac_rna]) 
rownames(rna2) <- rna_gene_anno$hgnc_symbol

archetype_TFBS_score2 <- as.matrix(archetype_TFBS_score[,ids_atac_rna])
archetype_TFBS_deviations2 <- as.matrix(archetype_TFBS_deviations[,ids_atac_rna])

atac_metadata3 <- atac_metadata2[ids_atac_rna,]
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Correlation between Archetype TFBS and the RNA expression of TFs from Archetype family
arch_tf_correlation <- vector()

bcell_colors <- jdb_color_map(c("MPP","LMPP","CD8","CD4","GMP-B","Ery","mono","MEP"))
names(bcell_colors) <- levels(atac_metadata3$Tissue)


pdf(paste(format(Sys.Date(),"%Y%m%d"),"_archetype_TFBS_TF_expression_correlation.pdf",sep=""))

col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))

for(i in 1:nrow(TF_cluster)){ 
  sel_arch <- colnames(ocr_fp_arch_matrix)[which(colnames(ocr_fp_arch_matrix)==TF_cluster$Name[i])]
  if(length(sel_arch)==1){
    motif_id_gene0 <- unique(TF_family$TF_Gene_ID[TF_family$Cluster_ID==i])
    tf_tmp_matrix <- matrix(NA, ncol=2+ncol(rna2), nrow=length(motif_id_gene0))
    colnames(tf_tmp_matrix) <- c("pval","rho",colnames(rna2))
    rownames(tf_tmp_matrix) <- paste(sel_arch,motif_id_gene0, sep="_")

    ha = HeatmapAnnotation(stage = atac_metadata3$Tissue, 
                           TFBS=anno_barplot(archetype_TFBS_score2[sel_arch,]),
                           col = list(stage = bcell_colors))
    
    for(j in 1:nrow(tf_tmp_matrix)){
      sel_rna_gene <- which(rna_gene_anno$hgnc_symbol%in%motif_id_gene0[j])
      if(length(sel_rna_gene)==1){
        tf_tmp_matrix[j,3:ncol(tf_tmp_matrix)] <- rna2[sel_rna_gene,]
        cor_output <- cor.test(rna2[sel_rna_gene,],archetype_TFBS_score2[sel_arch,], method="spearman")
        tf_tmp_matrix[j,1] <- cor_output$p.value
        tf_tmp_matrix[j,2] <- cor_output$estimate
      }
    }
    
    cat("Cluster ",i,">> num motifs: ",length(motif_id_gene0)," - num exprs.genes: ",sum(!is.na(tf_tmp_matrix[,1])),"\n")
    
    complex <- grep("\\+",motif_id_gene0)
    if(length(complex)!=0){
      cat("TF complexes in the cluster\n")
      complex_to_remove <- vector()
      for(k in 1:length(complex)){
        individual_gene_complex <- unlist(strsplit(motif_id_gene0[complex[k]],"\\+"))
        sel_rna_gene <- vector()
        for(h in 1:length(individual_gene_complex)){
          sel_rna_gene <- c(sel_rna_gene, which(rna_gene_anno$hgnc_symbol==individual_gene_complex[h]))
        }
        if(length(individual_gene_complex)==length(sel_rna_gene)){
          complex_info <- matrix(NA, ncol=2+ncol(rna2), nrow=length(individual_gene_complex))
          colnames(complex_info) <- c("pval","rho",colnames(rna2))
          rownames(complex_info) <- paste(sel_arch,"_",motif_id_gene0[complex[k]]," (",individual_gene_complex,")", sep="")
          complex_to_remove <- c(complex_to_remove,complex[k])
          for(w in 1:length(sel_rna_gene)){
            complex_info[w,3:ncol(tf_tmp_matrix)] <- rna2[sel_rna_gene[w],]
            cor_output <- cor.test(rna2[sel_rna_gene[w],],archetype_TFBS_score2[sel_arch,], method="spearman")
            complex_info[w,1] <- cor_output$p.value
            complex_info[w,2] <- cor_output$estimate
          }
          tf_tmp_matrix <- rbind(tf_tmp_matrix,complex_info)
        }
      }
      if(length(complex_to_remove)!=0){
        tf_tmp_matrix <- tf_tmp_matrix[-complex_to_remove,]
      }
    }
    
    arch_tf_correlation <- rbind(arch_tf_correlation,tf_tmp_matrix)
    
    if(!is.matrix(tf_tmp_matrix[,3:ncol(tf_tmp_matrix)])){
      tf_tmp_matrix <- rbind(tf_tmp_matrix,tf_tmp_matrix)
    }
    data_to_plot <- t(scale(t(tf_tmp_matrix[,3:ncol(tf_tmp_matrix)]), center = TRUE, scale = TRUE))
    data_to_plot[is.na(data_to_plot)] <- 0
    h_cell_size <- 0.47*nrow(data_to_plot)
    rha = rowAnnotation(rho = tf_tmp_matrix[,2],
                        col = list(rho = col_fun1))     
    draw(Heatmap(data_to_plot, top_annotation = ha, right_annotation = rha,
                 show_row_names = TRUE, show_column_names =FALSE, 
                 column_names_gp = gpar(fontsize = 6),
                 heatmap_legend_param = list(title = "gene expression "),
                 column_title = paste(sel_arch," (clust_",i,") - archetype expression",sep=""), 
                 width = unit(6, "cm"), height = unit(h_cell_size, "cm"),
                 row_names_gp = gpar(fontsize = 8),
                 column_order = order(atac_metadata3$Tissue)))
    
  }
}
dev.off()
 


result_TFBS_RNA <- data.frame(ID=rownames(arch_tf_correlation),pval=arch_tf_correlation[,1],rho=arch_tf_correlation[,2])


# Correct for multiple testing
result_TFBS_RNA$adj.pval_bonferroni <- p.adjust(result_TFBS_RNA$pval, method="bonferroni")
result_TFBS_RNA$adj.pval_fdr <- p.adjust(result_TFBS_RNA$pval, method="fdr")
result_TFBS_RNA$archetype <- unlist(lapply(strsplit(result_TFBS_RNA$ID,"_"),FUN=function(x){return(x[1])}))
result_TFBS_RNA$TF <- unlist(lapply(strsplit(result_TFBS_RNA$ID,"_"),FUN=function(x){return(x[2])}))


# Include TFBS variability
result_TFBS_RNA$archetype_variability <- rep(NA, length=nrow(result_TFBS_RNA))

arch <- unique(result_TFBS_RNA$archetype)
for(i in 1:length(arch)){
  result_TFBS_RNA$archetype_variability[result_TFBS_RNA$archetype==arch[i]] <- archetype_TFBS_variability$variability[archetype_TFBS_variability$name==arch[i]]
}


# Include RNA deviation
## Before compute CV need to have all variables positive
rna_pos <- rna+100
rna_cv <- apply(rna_pos,1,FUN=function(x){return(sd(x)/mean(x))})*100
hist(rna_cv)
quantile(rna_cv,0.25) #cut-off: 0.8956

result_TFBS_RNA$rna_cv <- rep(NA, length=nrow(result_TFBS_RNA))
result_TFBS_RNA$target_TF <- rep(NA, length=nrow(result_TFBS_RNA))
tf <- unique(result_TFBS_RNA$TF)
for(i in 1:length(tf)){
  target_tf <- tf[i]
  target_tf2 <- tf[i]
  if(length(grep("\\+",target_tf))>0){
    target_tf <- gsub(")","",unlist(strsplit(tf[i],"\\("))[2])
    target_tf2 <- gsub(")","",unlist(strsplit(tf[i],"\\("))[1])
  }
  result_TFBS_RNA$target_TF[result_TFBS_RNA$TF==tf[i]] <- target_tf2
  if(!is.na(target_tf)){
    gene_id <- rna_gene_anno$ensembl_gene_id[rna_gene_anno$hgnc_symbol==target_tf]
    if(length(gene_id)>0){
      result_TFBS_RNA$rna_cv[result_TFBS_RNA$TF==tf[i]] <- rna_cv[gene_id]
    }
  }
}


# Save results
write.table(result_TFBS_RNA,paste(format(Sys.Date(),"%Y%m%d"),"_correlation_ARCH_TFBS_RNA_TF.txt",sep=""),sep="\t",dec=".",row.names=FALSE, col.names=TRUE, quote=FALSE)
#-------------------------------------
