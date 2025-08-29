#***************************************************#
# GRN generation                                    #
# script: 03_grn_properties_igraph.R                #
#***************************************************#

# GRN generation


# Load required libraries
library(reshape2)
library(igraph)
library(GenomicRanges)


##>>>>>>>>>>>>>>>>>>>>> Inputs required

# TF-OCR-gene associations
tf_ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated_extended.rds"))

# Differential Analysis clustering
cluster_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.rds"))
#-------------------------------------


# >>>>>>>>>>>> B-cell GRN properties (igraph) --------------------
edges <- data.frame(from=tf_gene$tf_ensembl, to=tf_gene$gene_ensembl, width=tf_gene$gene_tf_rho)
nodes <- data.frame(id=unique(c(edges$from, edges$to)), label=unique(c(edges$from, edges$to)))

## GRN 
net <- graph.data.frame(edges, nodes, directed=T)

#Centrality measures
V(net)$degree <- igraph::degree(net)                        # Degree centrality
V(net)$indegree <- igraph::degree(net, mode="in")           
V(net)$outdegree <- igraph::degree(net, mode="out")           
V(net)$eig <- evcent(net)$vector                    # Eigenvector centrality
V(net)$hubs <- hub.score(net)$vector                # "Hub" centrality
V(net)$authorities <- authority.score(net)$vector   # "Authority" centrality
V(net)$closeness <- closeness(net)                  # Closeness centrality
V(net)$betweenness <- betweenness(net)              # Vertex betweenness centrality

save("net",file=paste(format(Sys.Date(),"%Y%m%d"),"_bcell_grn.RData",sep=""))

# To get adjacent matrix:
#adj_matrix <- as_adjacency_matrix(net)
#-------------------------------------

# >>>>>>>>>>>> B-cell GRN per cell type: properties (igraph) ---------------------------
cell_type <- c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5neg","Naive_CD5pos")
list_grn <- vector(mode="list", length=length(cell_type))

for(i in 1:length(cell_type)){
  sel <- which(colnames(tf_ocr_gene_extended)==paste0(cell_type[i],"_interaction"))
  aa <- tf_ocr_gene_extended[tf_ocr_gene_extended[,sel]=="yes",]
  exp_genes <- rownames(cluster_rna[cluster_rna[,cell_type[i]]==1,])
  sel2 <- which(aa$tf_ensembl%in%exp_genes & aa$gene.x%in%exp_genes)
  aa2 <- aa[sel2,]
  
  edges <- unique(data.frame(from=aa2$tf_ensembl, 
                             to=aa2$gene.x, 
                             weight=abs(aa2$gene_tf_rho))) 
  nodes <- data.frame(id=unique(c(edges$from, edges$to)), label=unique(c(edges$from, edges$to)))
  
  ## GRN 
  net <- graph.data.frame(edges, nodes, directed=T)
  
  #Centrality measures
  V(net)$degree <- igraph::degree(net)                        # Degree centrality
  V(net)$indegree <- igraph::degree(net, mode="in")           
  V(net)$outdegree <- igraph::degree(net, mode="out")           
  V(net)$eig <- evcent(net)$vector                    # Eigenvector centrality
  V(net)$hubs <- hub.score(net)$vector                # "Hub" centrality
  V(net)$authorities <- authority.score(net)$vector   # "Authority" centrality
  V(net)$closeness <- closeness(net)                  # Closeness centrality
  V(net)$betweenness <- betweenness(net)              # Vertex betweenness centrality
  
  list_grn[[i]] <- net
}
names(list_grn) <- cell_type

saveRDS(list_grn_properties,"bcell_grn_by_cell_type.rds")

# To get adjacent matrix:
#adj_matrix <- as_adjacency_matrix(net)
#-------------------------------------
