#****************************************#
# OCR-TF analysis - with chromVAR tool   #
# script: 21_TF_OCR_exploration.R        #
#****************************************#


# 21. OCR - TF exploration


# Load required libraries
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pals)
library(BuenColors)
library(circlize)
library(dplyr)


##>>>>>>>>>>>>>>>>>>>>> Inputs required

# 1- Normalized OCR (TMM) and metadata
atac_norm <- read.table("tmm_normalized_filtered_data_atac.txt", sep="\t", dec=".", header=TRUE, check.names = FALSE)
atac_metadata <- read.table("metadata_atac.txt", sep="\t", header=TRUE)

atac_norm <- atac_norm[,atac_metadata$SampleID]

atac_metadata$Tissue <- factor(as.character(atac_metadata$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))
rownames(atac_metadata) <- paste(atac_metadata$Factor,atac_metadata$Tissue,sep="_")
colnames(atac_norm) <- rownames(atac_metadata)


# 2- Normalized RNA data
rna <- read.table("voom_normalized_filtered_data_rna.txt", sep="\t", header=TRUE, check.names = FALSE)
rna_metadata <- read.table("metadata_rna.txt", sep="\t", header=TRUE)
rna_metadata$CellType <- factor(as.character(rna_metadata$CellType),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))
rna <- rna[,as.character(rna_metadata$Alias)] #23,454 x 79
rownames(rna_metadata) <- paste(rna_metadata$Donor,rna_metadata$CellType,sep="_")
colnames(rna) <- rownames(rna_metadata)
  
anno <- read.table("ensmbl_gene_annotation.txt", sep="\t", header=TRUE, check.names = FALSE)
anno_rna <- anno[anno$ensembl_gene_id %in% rownames(rna),]
rownames(anno_rna) <- anno_rna$ensembl_gene_id
anno_rna <- anno_rna[rownames(rna),]


# 3- Common samples ATAC-RNA
ids_atac_rna <- intersect(colnames(atac_norm), colnames(rna)) 
length(ids_atac_rna) #59 share samples

rna2 <- as.matrix(rna[,ids_atac_rna])
rownames(rna2) <- anno_rna$hgnc_symbol

atac_norm2 <- as.matrix(atac_norm[,ids_atac_rna])


# 4 - OCR-TF (motif) binary matrix
ocr_tf_matrix <- read.table("ocr_tf_matrix.txt", sep="\t", header=TRUE, check.names = FALSE)


# 5 - Differenital Analysis clustering
cluster_ocr <- read.table("clustering_of_OCRs_based_on_atacseq_dba_transitions.txt")
cluster_rna <- read.table("clustering_of_genes_based_on_rnaseq_dea_transitions.txt")
#-------------------------------------



# >>>>>>>>>>>> Heatmap of OCRs were target TFs binds (Figure 3E)

#Visualze select OCRs for specific TFs motifs

bcell_colors <- jdb_color_map(c("MPP","LMPP","CD8","CD4","GMP-B","Ery","mono","MEP"))
names(bcell_colors) <- levels(atac_metadata$Tissue)

target_TF_list <- c("TCF3","GATA2","GATA1","PAX5","HOXA9","CBFB","EBF1","NFATC2","SPIB")

for(i in 1:length(target_TF_list)){
  target_TF <- target_TF_list[i]
  
  target_ocrs <- rownames(ocr_tf_matrix)[ocr_tf_matrix[target_TF]!=0]
  target_ocrs_de <- target_ocrs[target_ocrs%in%rownames(cluster_ocr)]
  cat(target_TF,": ",length(target_ocrs_de),"target ocrs\n")

  data_to_plot <- log2(atac_norm2[target_ocrs_de,]+1)
  dim(data_to_plot)
  
  ha = HeatmapAnnotation(rna= anno_lines(rna2[target_TF,], smooth = TRUE, pch=16, gp = gpar(col = "red")),
                         stage = atac_metadata[ids_atac_rna,]$Tissue,
                         col = list(stage=bcell_colors))

  cluster_col <- as.character(jdb_palette("Darjeeling"))
  names(cluster_col) <- c(1:5)
  
  rha = rowAnnotation(clust = cluster_ocr[target_ocrs_de,]$clust_following_transitions,
                      col=list(clust=cluster_col))
  
  aa <- table(cluster_ocr[target_ocrs_de,]$clust_following_transitions)
  data_for_piechart <- data.frame(cluster=as.factor(names(aa)),
                                  value=as.numeric(aa))
  

  # Compute the position of labels
  library(dplyr)
  data_for_piechart <- data_for_piechart %>% 
    arrange(desc(cluster)) %>%
    mutate(prop = value / sum(data_for_piechart$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_OCRs_with_",target_TF,"_motif.pdf",sep=""))
  # Heatmap
  draw(Heatmap(data_to_plot, top_annotation = ha, right_annotation = rha,
               show_row_names = FALSE, show_column_names =FALSE, 
               col=as.character(jdb_palette("algae_earth"))[5:9],
               name="log2(ATAC)", row_names_gp = gpar(fontsize = 6),
               column_title = paste("OCRs with ",target_TF,"motif",sep=""),
               column_order = order(atac_metadata[ids_atac_rna,]$Tissue)))
  
  # Piechart
  print(ggplot(data_for_piechart, aes(x="", y=prop, fill=cluster)) +
          geom_bar(stat="identity", width=1, color="white") +
          coord_polar("y", start=0) +
          theme_void() +
          scale_fill_manual(values=cluster_col) +
          geom_text(aes(y = ypos, label = paste(round(data_for_piechart$prop),"%",sep="")), color = "white", size=6) +
          theme(legend.position="none"))

  dev.off()
}
#-------------------------------------



# >>>>>>>>>>>> Jaccard index between clustering of ATAC-Seq and RNA-Seq profiles (TFs) (Figure 3F)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


##RNA clusters

cl1 <- "HSC"
cl2 <- c("HSC","CLP")
cl3 <- c("HSC","CLP","proB")
cl4 <-  c("HSC","CLP","proB","preB")
cl5 <-  c("HSC","CLP","proB","ImmatureB")
cl6 <-  c("HSC","CLP","proB","ImmatureB", "TransitionalB")
cl7 <- c("NaiveCD5pos","NaiveCD5neg")
cl8 <- "NaiveCD5pos"
cl9 <- c("TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl10 <- c("ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl11 <- c("preB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl12 <- c("proB","preB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl13 <- c("CLP","proB","preB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl14 <- c("HSC","CLP","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl15 <- c("HSC","CLP","proB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl16 <- c("CLP","proB","preB","ImmatureB")
cl17 <- c("CLP","proB","preB")
cl18 <- "CLP"
cl19 <- c("proB","preB")
cl20 <- "proB"
cl21 <- "preB"
cl22 <- c("proB","preB","ImmatureB")

rna_clusters <- list(cl1,cl2,cl3,cl4,cl5,cl6,cl7,cl8,cl9,cl10,cl11,cl12,cl13,cl14,cl15,cl16,cl17,cl18,cl19,cl20,cl21,cl22)

#ATAC clusters

cl1a <- c("HSC","CLP")
cl2a <- c("HSC","CLP","proB")
cl3a <- c("HSC","CLP","proB","preB")
cl4a <- c("preB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")
cl5a <- c("CLP","proB","preB","ImmatureB","TransitionalB","NaiveCD5pos","NaiveCD5neg")

atac_clusters <- list(cl1a,cl2a,cl3a,cl4a,cl5a)


jaccard_clusters <- matrix(ncol=22,nrow=5)
colnames(jaccard_clusters) <- c(1:22)
rownames(jaccard_clusters) <- c(1:5)

for(i in 1:length(rna_clusters)){
  for(j in 1:length(atac_clusters)){
    jaccard_clusters[j,i] <- jaccard(rna_clusters[[i]],atac_clusters[[j]])
  }
}

write.table(jaccard_clusters,paste(format(Sys.Date(),"%Y%m%d"),"_jaccard_between_RNA_ATAC_clusters.txt",sep=""),
            sep="\t",dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)


pdf(paste(format(Sys.Date(),"%Y%m%d"),"_jaccard_clusters.pdf",sep=""))
Heatmap(t(jaccard_clusters), cluster_columns = FALSE, cluster_rows = FALSE,
        col=as.character(jdb_palette("solar_rojos"))[-1],
        rect_gp = gpar(col = "darkslategray", lwd = 0.3),
        column_names_side = "top")
dev.off()
#-------------------------------------



# >>>>>>>>>>>> Heatmap of OCRs were target TFs binds that can be relevant sites of actions (Figure 3G)

target_TF_list <- c("TCF3","GATA2","GATA1","PAX5","HOXA9","CBFB","EBF1","NFATC2","SPIB")

for(i in 1:length(target_TF_list)){
  target_TF <- target_TF_list[i]
  
  target_ocrs <- rownames(ocr_tf_matrix)[ocr_tf_matrix[target_TF]!=0]
  
  TF_cluster <- cluster_rna$subclust_following_transitions[cluster_rna$gene==anno_rna$ensembl_gene_id[anno_rna$hgnc_symbol==target_TF]]
  if(length(TF_cluster)==0){
    cat(target_TF, " without DE gene\n")
    next
  }
  OCR_cluster <- jaccard_clusters[,TF_cluster]
  OCR_cluster <- which(OCR_cluster>=0.5 | OCR_cluster==max(OCR_cluster)) #Decision based on jaccard index
  
  target_ocrs_de <- target_ocrs[target_ocrs%in%rownames(cluster_ocr[cluster_ocr$clust_following_transitions%in%OCR_cluster,])]
  cat(target_TF,": ",length(target_ocrs_de),"target ocrs\n")

  data_to_plot <- log2(atac_norm2[target_ocrs_de,]+1)
  dim(data_to_plot)
  
  ha = HeatmapAnnotation(rna= anno_lines(rna2[target_TF,], smooth = TRUE, pch=16, gp = gpar(col = "red")),
                         stage = atac_metadata[ids_atac_rna,]$Tissue,
                         col = list(stage=bcell_colors))
  
  
  cluster_col <- as.character(jdb_palette("Darjeeling"))
  names(cluster_col) <- c(1:5)
  
  rha = rowAnnotation(clust = cluster_ocr[target_ocrs_de,]$clust_following_transitions,
                      col=list(clust=cluster_col))
  
  aa <- table(cluster_ocr[target_ocrs_de,]$clust_following_transitions)
  data_for_piechart <- data.frame(cluster=as.factor(names(aa)),
                                  value=as.numeric(aa))
  
  # Compute the position of labels
  data_for_piechart <- data_for_piechart %>% 
    arrange(desc(cluster)) %>%
    mutate(prop = value / sum(data_for_piechart$value) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_candidate_OCRs_with_",target_TF,"_motif.pdf",sep=""))
  # Heatmap
  draw(Heatmap(data_to_plot, top_annotation = ha, right_annotation = rha,
               show_row_names = FALSE, show_column_names =FALSE, 
               col=as.character(jdb_palette("algae_earth"))[5:9],
               name="log2(ATAC)", row_names_gp = gpar(fontsize = 6),
               column_title = paste("OCRs with ",target_TF,"motif",sep=""),
               column_order = order(atac_metadata[ids_atac_rna,]$Tissue)))
  
  # Piechart
  print(ggplot(data_for_piechart, aes(x="", y=prop, fill=cluster)) +
          geom_bar(stat="identity", width=1, color="white") +
          coord_polar("y", start=0) +
          theme_void() +
          scale_fill_manual(values=cluster_col) +
          geom_text(aes(y = ypos, label = paste(round(data_for_piechart$prop),"%",sep="")), color = "white", size=6) +
          theme(legend.position="none"))

  dev.off()
}
#-------------------------------------
