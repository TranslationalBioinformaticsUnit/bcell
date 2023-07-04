#****************************************************#
# ATAC-Seq pipeline - PAIR-END                       #
# script: 13_diff_binding_transitions_clustering.sh  #
# script: 13_diff_binding_transitions_clustering.R   #
#****************************************************#

# 13. Differential Binding Analysis Clustering

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


# Required libraries
library(mclust)
library(edgeR)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

set.seed(12345)


##>>>>>>>>>>>>>>>>>>>>> Input data

## significant transitions from DBA 
dba_sig <- read.table(sig_peaks, sep="\t", dec=".", header=TRUE, check.names=FALSE)

# Metadata
metadata <- read.table(meta_file, sep=";", header=TRUE)
rownames(metadata) <- as.character(metadata$SampleID)
metadata <- metadata[colnames(counts),]
metadata$Tissue <- droplevels(metadata$Tissue)


##>>>>>>>>>>>>>>>>>>>>> Clustering

data_scaled <- dba_sig[,-9]
# Clustering with mclust

# BIC <- mclustBIC(data_scaled)
mod1 <- Mclust(data_scaled, G = 9, modelName = "EII")

summary(as.factor(mod1$classification))

clust <- data.frame(gene=names(mod1$classification), cluster=mod1$classification)

# save cluster classification
write.table(clust,paste(format(Sys.time(), "%Y%m%d"),"_cluster_of_binary_genes_based_on_atacseq.txt",sep=""), sep="\t", quote=FALSE)


##>>>>>>>>>>>>>>>>>>>>> Heatmap

#1. Order data by cell type
data_scaled <- dba_sig
Cell <- colnames(data_scaled)

#2. Order data by cluster
cluster_order <- order(clust[,2])

#3. Ordered data
heatmap_matrix <- data_scaled[cluster_order,]

#4. Heatmap
my_palette1 <- c(brewer.pal(11, "Spectral")[c(1:11)])

ha1 <- HeatmapAnnotation(
  CellType = Cell, 
  col = list(CellType = c("sig_HSC_CLP"=my_palette1[1], 
                          "sig_CLP_ProB"=my_palette1[2], 
                          "sig_ProB_PreB"=my_palette1[3],
                          "sig_PreB_Immature.B"=my_palette1[4], 
                          "sig_Immature.B_Transitional.B"=my_palette1[5], 
                          "sig_Transitional.B_Naive.CD5neg"=my_palette1[6], 
                          "sig_Transitional.B_Naive.CD5pos"=my_palette1[7], 
                          "sig_Naive.CD5neg_Naive.CD5pos"=my_palette1[8], 
                          "sig_HSC_Transitional.B"=my_palette1[9])),
  show_annotation_name = TRUE
)


png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_cluster_binary_peaks.png",sep=""), width=1000, height=1000)
draw(Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(-1,0,1), c("pink","gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clust[cluster_order,2],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()



##>>>>>>>>>>>>>>>>>>>>> Sub-Clustering

## Get first clustering
data_scaled <- dba_sig[-9]
# Clustering with mclust

# Identify optimal number of groups within cluster
for(i in 5:9){
  BIC <- mclustBIC(data_scaled[rownames(clust)[clust$cluster==i],])
  pdf(paste(format(Sys.time(), "%Y%m%d"),"_BIC_plot_C",i,".pdf",sep=""))
    plot(BIC)
  dev.off()
}

subcluster <- vector(mode = "list", length =9)
names(subcluster) <- paste("C",1:9,sep="")

for(i in 1:max(clust$cluster)){
  dd <- dba_sig[clust$cluster==i,-9]
  dd_plot <- dba_sig[clust$cluster==i,]

  if(i%in%c(3:8)){
  mod1 <- Mclust(dd, G = 9, modelName = "EII")

  summary(as.factor(mod1$classification))

  subclust <- data.frame(gene=names(mod1$classification), cluster=mod1$classification)
  subcluster[[i]] <- mod1

##>>>>>>>>>>>>>>>>>>>>> Heatmap

#1. Order data by cell type
Cell <- colnames(dd_plot)

#2. Order data by cluster
cluster_order <- order(subclust[,2])

#3. Ordered data
heatmap_matrix <- dd_plot[cluster_order,]

#4. Heatmap
my_palette1 <- c(brewer.pal(11, "Spectral")[c(1:11)])

ha1 <- HeatmapAnnotation(
  CellType = Cell, 
  col = list(CellType = c("sig_HSC_CLP"=my_palette1[1], 
                          "sig_CLP_ProB"=my_palette1[2], 
                          "sig_ProB_PreB"=my_palette1[3],
                          "sig_PreB_Immature.B"=my_palette1[4], 
                          "sig_Immature.B_Transitional.B"=my_palette1[5], 
                          "sig_Transitional.B_Naive.CD5neg"=my_palette1[6], 
                          "sig_Transitional.B_Naive.CD5pos"=my_palette1[7], 
                          "sig_Naive.CD5neg_Naive.CD5pos"=my_palette1[8], 
                          "sig_HSC_Transitional.B"=my_palette1[9])),
  show_annotation_name = TRUE
)

png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_subcluster_binary_peaks_C",i,"_G9.png",sep=""), width=1000, height=1000)
draw(Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(-1,0,1), c("pink","gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=subclust[cluster_order,2],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()
}

}

saveRDS(subcluster,file=paste0(format(Sys.time(), "%Y%m%d"),"_subclustering.rds"))

# Generate a dataframe with all sublcusters
gene <- vector()
subclusters <- vector()
cluster <- vector()

subcluster[[1]]$classification <- rep(1, sum(clust$cluster==1))
names(subcluster[[1]]$classification) <- rownames(dba_sig)[clust$cluster==1]

subcluster[[2]]$classification <- rep(1, sum(clust$cluster==2))
names(subcluster[[2]]$classification) <- rownames(dba_sig)[clust$cluster==2]

subcluster[[9]]$classification <- rep(1, sum(clust$cluster==9))
names(subcluster[[9]]$classification) <- rownames(dba_sig)[clust$cluster==9]


for(i in 1:length(subcluster)){
    gene <- c(gene, names(subcluster[[i]]$classification))
    subclusters <- c(subclusters, subcluster[[i]]$classification)
    cluster <- c(cluster,rep(i,length(subcluster[[i]]$classification)))
}

subcluster_id <- paste0(cluster,"_",subclusters)

clustering_df <- data.frame(gene,cluster,subclusters,subcluster_id)


### Get the subcluster expression patterns

clust_id <- vector()
HSC <- vector()
CLP <- vector()
ProB <- vector()
PreB <- vector()
Immature.B <- vector()
Transitional.B <- vector()
Naive.CD5neg <- vector()
Naive.CD5pos <- vector()

for(i in 1:9){
  sub_clust_info <- read.table(paste("20221108_pattern_subcluster",i,".txt", sep=""),
   header=FALSE)
  clust_id <- c(clust_id, paste0(i,"_",1:9)) 
  HSC <- c(HSC,sub_clust_info[,1])
  CLP <- c(CLP,sub_clust_info[,2])
  ProB <- c(ProB,sub_clust_info[,3])
  PreB <- c(PreB,sub_clust_info[,4])
  Immature.B <- c(Immature.B,sub_clust_info[,5])
  Transitional.B <- c(Transitional.B,sub_clust_info[,6])
  Naive.CD5neg <- c(Naive.CD5neg,sub_clust_info[,7])
  Naive.CD5pos <- c(Naive.CD5pos,sub_clust_info[,8])
}

cluster_pattern_df <- data.frame(clust_id, HSC=(HSC-1)*-1, CLP=(CLP-1)*-1, ProB=(ProB-1)*-1, PreB=(PreB-1)*-1, Immature.B=(Immature.B-1)*-1, Transitional.B=(Transitional.B-1)*-1, Naive.CD5neg=(Naive.CD5neg-1)*-1, Naive.CD5pos=(Naive.CD5pos-1)*-1)
cluster_pattern_df$pattern_id <- rep(NA, nrow(cluster_pattern_df))
##Unique patterns
unique_patterns <- unique(cluster_pattern_df[,2:9]) #29 in total

#High in HSC:
sel1 <- which(unique_patterns[,1]==1)
unique_patterns_ord_high <- order(unique_patterns[sel1,2],
                    unique_patterns[sel1,3],
                    unique_patterns[sel1,4],
                    unique_patterns[sel1,5],
                    unique_patterns[sel1,6],
                    unique_patterns[sel1,7],
                    unique_patterns[sel1,8],
                    decreasing=FALSE)
#Low in HSC:
sel2 <- which(unique_patterns[,1]==0)
unique_patterns_ord_low <- order(unique_patterns[sel2,2],
                    unique_patterns[sel2,3],
                    unique_patterns[sel2,4],
                    unique_patterns[sel2,5],
                    unique_patterns[sel2,6],
                    unique_patterns[sel2,7],
                    unique_patterns[sel2,8],
                    decreasing=FALSE)
                    
unique_patterns_ord <- unique_patterns[c(sel1[unique_patterns_ord_high],sel2[unique_patterns_ord_low]),]

for(i in 1:nrow(unique_patterns_ord)){
  ss <- apply(cluster_pattern_df[,2:9],1,FUN=function(x){
    return(sum(x==unique_patterns_ord[i,]))})
  sel <- which(ss==8)
  cluster_pattern_df$pattern_id[sel] <- i 
}


clustering_df2 <- merge(clustering_df,cluster_pattern_df, by.x="subcluster_id",by.y="clust_id", sort=FALSE)
rownames(clustering_df2) <- clustering_df2$gene

saveRDS(clustering_df2,file=paste0(format(Sys.time(), "%Y%m%d"),"_clustering_patterns_after_subclustering.rds"))



##>>>>>>>>>>>>>>>>>>>>> Heatmap
cluster_order <- order(clustering_df2$pattern_id)

#1. Ordered data
heatmap_matrix <- clustering_df2[cluster_order,c(5:10,12,11)]

#2. Heatmap
my_palette1 <- c(brewer.pal(11, "Spectral")[c(1:11)])

ha1 <- HeatmapAnnotation(
  CellType = colnames(heatmap_matrix), 
  col = list(CellType = c("HSC"=my_palette1[1], 
                          "CLP"=my_palette1[2], 
                          "ProB"=my_palette1[3],
                          "PreB"=my_palette1[4], 
                          "Immature.B"=my_palette1[5], 
                          "Transitional.B"=my_palette1[6], 
                          "Naive.CD5pos"=my_palette1[7], 
                          "Naive.CD5neg"=my_palette1[8])),
  show_annotation_name = TRUE
)


png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_cluster_binary_pattern_peaks.png",sep=""), width=1000, height=1500)
draw(Heatmap(heatmap_matrix, name = "Peaks vs Cell Type",
        col = colorRamp2(c(0,1), c("gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clustering_df2$pattern_id[cluster_order],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()
