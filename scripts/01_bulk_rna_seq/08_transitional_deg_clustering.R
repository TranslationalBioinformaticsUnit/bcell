#*******************************************#
# RNA-Seq pipeline - SINGLE-END             #
# script: 08_transitional_deg_clustering.sh #
# script: 08_transitional_deg_clustering.R  #
#*******************************************#


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
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

set.seed(12345)



##>>>>>>>>>>>>>>>>>>>>> Input data
#transitional dea output 
dea_sig <- read.table(trans_dea_sig,sep="\t", dec=".", header=TRUE)
dea_sig_plot <- read.table(trans_dea_sig_plot,sep="\t", dec=".", header=TRUE)


##>>>>>>>>>>>>>>>>>>>>> Clustering

data_scaled <- dea_sig
# Clustering with mclust

BIC <- mclustBIC(data_scaled)
pdf(paste(format(Sys.time(), "%Y%m%d"),"_BIC_plot.pdf",sep=""))
  plot(BIC)
dev.off()

mod1 <- Mclust(data_scaled, G = 9, modelName = "EII")

summary(as.factor(mod1$classification))

clust <- data.frame(gene=names(mod1$classification), cluster=mod1$classification)

# save cluster classification
write.table(clust,paste(format(Sys.time(), "%Y%m%d"),"_cluster_of_binary_genes_based_on_rnaseq.txt",sep=""),sep="\t", dec=".", quote=FALSE)



##>>>>>>>>>>>>>>>>>>>>> Heatmap

#1. Order data by cell type
data_scaled <- dea_sig_plot
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

png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_cluster_binary_genes.png",sep=""), width=1000, height=1000)
Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(-1,0,1), c("pink","gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clust[cluster_order,2],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()




##>>>>>>>>>>>>>>>>>>>>> Sub-Clustering


## Clustering
data_scaled <- dea_sig
# Clustering with mclust

# Identify optimal number of groups within cluster
for(i in 1:9){
  BIC <- mclustBIC(data_scaled[rownames(clust)[clust$cluster==i],])
  pdf(paste(format(Sys.time(), "%Y%m%d"),"_BIC_plot_C",i,".pdf",sep=""))
    plot(BIC)
  dev.off()
}

g <- c(3,4,4,3,4,3,3,8,3)

subcluster <- vector(mode = "list", length =9)
names(subcluster) <- paste("C",1:9,sep="")

for(i in 1:max(clust$cluster)){
  dd <- data_scaled[clust$cluster==i,]
  dd_plot <- dea_sig[clust$cluster==i,]
  
# BIC <- mclustBIC(dd)
mod1 <- Mclust(dd, G = g[i], modelName = "EII")

summary(as.factor(mod1$classification))

subclust <- data.frame(gene=names(mod1$classification), cluster=mod1$classification)
subcluster[[i]] <- mod1


##>>>>>>>>>>>>>>>>>>>>> Heatmap

#1. Order data by cell type?
Cell <- colnames(dd_plot)

#2. Order data by cluster
cluster_order <- order(subclust[,2])

#3. Ordered data
heatmap_matrix <- dd_plot[cluster_order,]

#4. Heatmap
png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_subcluster_binary_genes_cluster",i,".png",sep=""), width=1000, height=1000)
draw(Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(-1,0,1), c("pink","gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=subclust[cluster_order,2],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

}


### Unify the subclusters of each cluster in a unique list and heatmap

cont <- 0
clust$subclusters <- clust$cluster
for(i in 1:length(subcluster)){
clust$subclusters[clust$cluster==i] <- cont + subcluster[[i]]$classification
cont <- max(cont + subcluster[[i]]$classification)
cat(cont,"\n")
}


##>>>>>>>>>>>>>>>>>>>>> Heatmap
data_scaled <- dea_sig

#1. Order data by cell type
Cell <- colnames(data_scaled)

#2. Order data by cluster
cluster_order <- order(clust$subclusters)

#3. Ordered data
heatmap_matrix <- data_scaled[cluster_order,]

#4. Heatmap
png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_subcluster_binary_genes.png",sep=""), width=1000, height=1000)
Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(-1,0,1), c("pink","gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clust[cluster_order,3],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()


# save cluster classification
write.table(clust,paste(format(Sys.time(), "%Y%m%d"),"_subcluster_of_genes_based_on_binary_rnaseq.txt",sep=""), sep="\t", quote=FALSE)



##>>>>>>>>>>>>>>>>>>>>> Identifying patterns

heatmap_matrix2 <- dea_sig[,-9]

patterns <- data.frame(c1=c(0,0,0,0,0,0,1,1),
                       c2=c(0,0,0,0,0,1,1,1),
                       c3=c(0,0,0,0,0,0,1,0),
                       c4=c(1,1,1,1,0,0,0,0),
                       c5=c(1,1,1,1,0,0,0,0),
                       c6=c(1,1,1,1,0,0,0,0),
                       c7=c(0,0,0,1,0,0,0,0),
                       c8=c(1,0,0,0,0,0,0,0),
                       c9=c(1,1,1,0,0,0,0,0),
                       c10=c(1,1,1,0,0,0,0,0),
                       c11=c(1,1,0,0,0,0,0,0),
                       c12=c(0,1,0,0,0,0,0,0),
                       c13=c(1,1,0,0,0,0,0,0),
                       c14=c(1,1,0,0,1,1,1,1),
                       c15=c(0,0,0,0,1,1,1,1),
                       c16=c(0,0,1,1,1,1,1,1),
                       c17=c(0,0,0,0,1,1,1,1),
                       c18=c(0,0,0,1,1,1,1,1),
                       c19=c(1,1,1,0,0,0,0,0),
                       c20=c(1,1,1,0,0,0,0,0),
                       c21=c(1,1,1,0,1,1,1,1),
                       c22=c(0,0,1,1,1,1,1,1),
                       c23=c(0,0,1,1,0,0,0,0),
                       c24=c(0,0,1,0,0,0,0,0),
                       c25=c(1,1,1,1,1,0,0,0),
                       c26=c(0,0,0,1,1,1,1,1),
                       c27=c(0,0,0,1,1,1,1,1),
                       c28=c(0,0,1,1,1,1,1,1),
                       c29=c(0,0,0,1,1,1,1,1),
                       c30=c(1,1,1,1,1,1,0,0),
                       c31=c(0,0,1,1,1,1,1,1),
                       c32=c(0,0,1,1,1,0,0,0),
                       c33=c(0,1,1,1,1,0,0,0),
                       c34=c(0,1,1,1,0,0,0,0),
                       c35=c(0,1,1,1,1,1,1,1))

for(k in 1:max(clust$subclusters)){
  sel <- which(clust$subclusters==k)
  for(i in 1:length(sel)){
    heatmap_matrix2[sel[i],] <- patterns[,k]
  }
}

colnames(heatmap_matrix2) <- c("HSC","CLP","proB","preB","ImmatureB","TransitionalB","NaiveCD5neg","NaiveCD5pos")


##>>>>>>>>>>>>>>>>>>>>> Heatmap

cluster_order <- order(clust$subclusters)

#1. Ordered data
heatmap_matrix <- heatmap_matrix2[cluster_order,]


#2. Heatmap
my_palette1 <- c(brewer.pal(11, "Spectral")[c(1:11)])

ha1 <- HeatmapAnnotation(
  CellType = colnames(heatmap_matrix), 
  col = list(CellType = c("HSC"=my_palette1[1], 
                          "CLP"=my_palette1[2], 
                          "proB"=my_palette1[3],
                          "preB"=my_palette1[4], 
                          "ImmatureB"=my_palette1[5], 
                          "TransitionalB"=my_palette1[6], 
                          "NaiveCD5neg"=my_palette1[7], 
                          "NaiveCD5pos"=my_palette1[8])),
  show_annotation_name = TRUE
)


png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_subcluster_binary_pattern_genes.png",sep=""), width=1000, height=1000)
draw(Heatmap(heatmap_matrix, name = "Genes vs Cell Type",
        col = colorRamp2(c(0,1), c("gray", "darkblue")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clust[cluster_order,3],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()




##>>>>>>>>>>>>>>>>>>>>> Tunning heatmap - to show the clusters according to cell transisiton

clust$subclust_following_transitions <- clust$subclusters

clust$subclust_following_transitions[clust$subclusters==8] <- 1
clust$subclust_following_transitions[clust$subclusters==11] <- 2
clust$subclust_following_transitions[clust$subclusters==13] <- 2
clust$subclust_following_transitions[clust$subclusters==9] <- 3
clust$subclust_following_transitions[clust$subclusters==10] <- 3
clust$subclust_following_transitions[clust$subclusters==19] <- 3
clust$subclust_following_transitions[clust$subclusters==20] <- 3
clust$subclust_following_transitions[clust$subclusters==4] <- 4
clust$subclust_following_transitions[clust$subclusters==5] <- 4
clust$subclust_following_transitions[clust$subclusters==6] <- 4
clust$subclust_following_transitions[clust$subclusters==25] <- 5
clust$subclust_following_transitions[clust$subclusters==30] <- 6
clust$subclust_following_transitions[clust$subclusters==1] <- 7
clust$subclust_following_transitions[clust$subclusters==3] <- 8
clust$subclust_following_transitions[clust$subclusters==2] <- 9
clust$subclust_following_transitions[clust$subclusters==15] <- 10
clust$subclust_following_transitions[clust$subclusters==17] <- 10
clust$subclust_following_transitions[clust$subclusters==18] <- 11
clust$subclust_following_transitions[clust$subclusters==26] <- 11
clust$subclust_following_transitions[clust$subclusters==27] <- 11
clust$subclust_following_transitions[clust$subclusters==29] <- 11
clust$subclust_following_transitions[clust$subclusters==16] <- 12
clust$subclust_following_transitions[clust$subclusters==22] <- 12
clust$subclust_following_transitions[clust$subclusters==28] <- 12
clust$subclust_following_transitions[clust$subclusters==31] <- 12
clust$subclust_following_transitions[clust$subclusters==35] <- 13
clust$subclust_following_transitions[clust$subclusters==14] <- 14
clust$subclust_following_transitions[clust$subclusters==21] <- 15
clust$subclust_following_transitions[clust$subclusters==33] <- 16
clust$subclust_following_transitions[clust$subclusters==34] <- 17
clust$subclust_following_transitions[clust$subclusters==12] <- 18
clust$subclust_following_transitions[clust$subclusters==23] <- 19
clust$subclust_following_transitions[clust$subclusters==24] <- 20
clust$subclust_following_transitions[clust$subclusters==7] <- 21
clust$subclust_following_transitions[clust$subclusters==32] <- 22

cluster_order <- order(clust$subclust_following_transitions)

heatmap_matrix <- heatmap_matrix2[cluster_order,]

png(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_subcluster_binary_pattern_genes_by_transitions.png",sep=""), width=1000, height=1000)
Heatmap(heatmap_matrix, name = "Peaks vs Cell Type",
        col = colorRamp2(c(0, 1), c("darkgreen","darkorange")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = FALSE, row_split=clust$subclust_following_transitions[cluster_order],
        cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()


# save cluster classification
write.table(clust,paste(format(Sys.time(), "%Y%m%d"),"_subcluster_of_genes_based_on_binary pattern_rnaseq_by_transitions.txt",sep=""), sep="\t", quote=FALSE)
