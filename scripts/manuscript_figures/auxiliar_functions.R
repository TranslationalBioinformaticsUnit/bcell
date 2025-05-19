#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its 
## implications in the pathogenesis of B-cell acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2024
## Code: Auxiliar function
#***************************************************


# Plot TreeDot for GSEA comparision results -----------------------------------
plot_treedot = function(comp_pair_term, top_paths=5, clust_num=3) {
  # libraries
  library(ggplot2)
  library(tidyverse)
  library(ggdendro)
  library(cowplot)
  library(clusterProfiler)
  library(ggtree)
  library(patchwork) 
  library(org.Hs.eg.db)
  library(RColorBrewer)
  
  
  # Estruture Data
  comp_pair_term_fort <- fortify(comp_pair_term, showCategory = top_paths, includeAll = TRUE, split = NULL)
  comp_pair_term_fort$Cluster <- sub("\n.*", "", comp_pair_term_fort$Cluster)
  comp_pair_term_fort$geneID <- comp_pair_term_fort$core_enrichment
  
  # Merge Clusters
  merge_compareClusterResult <- function(yy) {
    yy_union <- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
      ids <- unique(unlist(strsplit(x$geneID, "/")))
      cnt <- length(ids)
      list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })
    
    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))
    
    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
  }
  merged_ggData <- merge_compareClusterResult(comp_pair_term_fort)
  
  # Prepare data for IDs Cluster
  prepare_pie_category <- function(enrichDf, pie = "equal") {
    pie <- match.arg(pie, c("equal", "count", "Count"))
    if (pie == "count") pie <- "Count"
    
    pie_data <- enrichDf[,c("Cluster", "Description", "Count")]
    pie_data[,"Description"] <- as.character(pie_data[,"Description"])
    prepare_pie_data(pie_data, pie = pie)
  }
  prepare_pie_data <- function(pie_data, pie = "equal",type = "category") {
    if(type == "category"){
      ID_unique <- unique(pie_data[,2])
    } else {
      ID_unique <- unique(pie_data[,3])
    }
    
    Cluster_unique <- unique(pie_data[,1])
    ID_Cluster_mat <- matrix(0, nrow = length(ID_unique), ncol = length(Cluster_unique))
    rownames(ID_Cluster_mat) <- ID_unique
    colnames(ID_Cluster_mat) <- Cluster_unique
    ID_Cluster_mat <- as.data.frame(ID_Cluster_mat, stringAsFactors = FALSE)
    if(pie == "Count") {
      for(i in seq_len(nrow(pie_data))) {
        ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- pie_data[i,3]
      }
      for(kk in seq_len(ncol(ID_Cluster_mat))) {
        ID_Cluster_mat[,kk] <- as.numeric(ID_Cluster_mat[,kk])
      }
      return(ID_Cluster_mat)
    }
    for(i in seq_len(nrow(pie_data))) {
      if(type == "category"){
        ID_Cluster_mat[pie_data[i,2],pie_data[i,1]] <- 1
      } else {
        ID_Cluster_mat[pie_data[i,3],pie_data[i,1]] <- 1
      }
      
    }
    return(ID_Cluster_mat)
  }
  ID_Cluster_mat <- prepare_pie_category(comp_pair_term_fort,pie = "equal")
  
  # Hierarchical Clustering 
  fill_termsim <- function(x, keep) {
    termsim <- x@termsim[keep, keep]
    termsim[which(is.na(termsim))] <- 0
    termsim2 <- termsim + t(termsim)
    for ( i in seq_len(nrow(termsim2)))
      termsim2[i, i] <- 1
    return(termsim2)
  }
  termsim2 <- fill_termsim(comp_pair_term, rownames(ID_Cluster_mat))
  hc <- stats::hclust(stats::as.dist(1- termsim2),method = "ward.D")
  
  # Info of Clusters wanted 
  clus <- stats::cutree(hc, clust_num)
  
  # ORA for clusters
  keywords <- vector()
  for (i in 1:clust_num) {
    paths_of_clust <- clus[clus == i]
    
    cluster_paths_info <- comp_pair_term_fort[comp_pair_term_fort$Description %in% names(paths_of_clust) ,]
    genes2check <-  unique(unlist(str_split(cluster_paths_info$geneID, "/")))
    
    if (comp_pair_term@fun == "gseGO"){
      ora <- enrichGO(genes2check, OrgDb= "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "BP", pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, 
                      minGSSize = 10,maxGSSize = 500, readable = FALSE, pool = FALSE)
      
    } else if (comp_pair_term@fun %in% c("gseKEGG", "gseDO")) {
      ora <- enrichGO(genes2check, OrgDb= "org.Hs.eg.db", keyType= "ENTREZID",
                      ont = "BP", pvalueCutoff = 1, pAdjustMethod= "BH", qvalueCutoff = 1, 
                      minGSSize = 10,maxGSSize = 500, readable = FALSE, pool = FALSE)
      
    } else {print("object@keytype is neither: SYMBOL or kegg")}
    
    keywords[i] <- ora@result$Description[1]
  }
  
  
  # Plot dendogram
  g <- split(names(clus), clus)
  p <- ggtree(hc, size=1.2)
  clades <- sapply(g, function(n) MRCA(p, n))
  p <- groupClade(p, clades, group_name='SubTree_ORA') + aes(color=SubTree_ORA)
  mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 6))
  ggtree_plot <- p + scale_color_manual(values=mycolors, breaks=1:clust_num,labels = keywords)+
    theme(legend.position='right',legend.justification = c(0,1.5))
  ggtree_plot_noLegend <- ggtree_plot + theme(legend.position = "none")
  
  
  # Prepare for dotplot
  comp_pair_term_fort$log_p.adjust <- -log10(comp_pair_term_fort$p.adjust)
  comp_pair_term_fort$Description = factor(comp_pair_term_fort$Description, levels = hc$labels[hc$order])
  comp_pair_term_fort$Cluster <- factor(comp_pair_term_fort$Cluster,levels=levels(comp_pair_term@compareClusterResult$Cluster))
  # Plot dotplot
  dotplot <- comp_pair_term_fort %>% 
    ggplot(aes(x=Cluster, y =Description , color = NES, size = log_p.adjust)) + 
    geom_point() +
    scale_y_discrete(position = "right")+
    scale_color_gradient2(low="blue4",mid="white", high="red")+
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    guides(size=guide_legend(title="-log10(p.adj)"))+
    theme(axis.ticks = element_blank(),legend.position = "right",legend.justification = c(0,0))
  dotplot_noLegend <- dotplot + theme(legend.position = "none")
  
  # Legends 
  legend_tree <- get_legend(ggtree_plot)
  legend_dot <- get_legend(dotplot)
  
  # Merge Plot
  combine_plot <- plot_grid(ggtree_plot_noLegend,NULL, dotplot_noLegend, nrow = 1, rel_widths = c(0.3,-0.05,2), align = 'h')
  combine_legend <- plot_grid(legend_dot,NULL,legend_tree, ncol=1, rel_heights = c(1,-0.5,1))
  big_plot <- plot_grid(combine_plot,combine_legend,NULL, nrow = 1, rel_widths = c(1,0.1, 0.16))
  return(big_plot)
  
}
# ----------------------------------------------------------------------------- #


# Human - Mouse conversion ----------------------------------------------------
# Basic function to convert mouse to human gene names
convertMouseGeneList = function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  attributes = listAttributes(mouse)
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows = TRUE)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}
# ----------------------------------------------------------------------------- #


# Regulon specificity score (RSS) ---------------------------------------------
# Calculates Regulon specificity score (RSS) from binary regulon activity.
# Code obtained from: https://rdrr.io/github/FloWuenne/scFunctions/src/R/calculate_rrs.R
calculate_rss = function(metadata,binary_regulons,cell_type_column){
  
  cell_types <- unique(metadata[,cell_type_column])
  regulons <- rownames(binary_regulons)
  
  
  jsd_matrix_ct <- data.frame("regulon" = c(),
                              "cell_type" = c(),
                              "jsd" = c())
  
  
  cell_type_counter <- 0
  for(ct in unique(cell_types)) {
    
    cell_type_counter <- cell_type_counter + 1
    print(paste("Processing cell type:",cell_type_counter,ct,sep=" "))
    
    for(regulon_no in 1:length(regulons)) {
      
      regulon <- regulons[regulon_no]
      
      regulon_vec <- binary_regulons[regulon,]
      
      regulon_vec_sum <- sum(regulon_vec)
      
      ## Check that there are cells with binary activity > 0 for this regulon
      if(regulon_vec_sum > 0){
        
        #progress(regulon_no)
        
        regulon_norm <- regulon_vec/regulon_vec_sum
        
        genotype_vec <- metadata[colnames(binary_regulons),]
        genotype_vec <- genotype_vec %>%
          mutate("cell_class" = if_else(get(cell_type_column) == ct,1,0))
        
        genotype_vec <- genotype_vec$cell_class
        genotype_norm <- genotype_vec/sum(genotype_vec)
        
        dist_df <- rbind(regulon_norm,genotype_norm)
        
        ## Calculate the Jensen-Shannon divergence
        jsd_divergence <- suppressMessages(philentropy::JSD(dist_df))
        
        ## Calculate Jensen-Shannon distance
        rss <- 1-sqrt(jsd_divergence)
        
        regulon_jsd <- data.frame("regulon" = regulon,
                                  "cell_type" = ct,
                                  "RSS" = rss[1])
        
        jsd_matrix_ct <- rbind(jsd_matrix_ct,regulon_jsd)
        
      }else if(regulon_vec_sum == 0){
        print(paste("Filtered out:",regulon,". No cells with binary activity > 0 identified. Please check your threshold for this regulon!",sep=""))
      }
    }
    
  }
  
  jsd_matrix_ct <- jsd_matrix_ct %>%
    arrange(desc(RSS))
  
  return(jsd_matrix_ct)
}
# ----------------------------------------------------------------------------- #


# MDS representation from GRN pearson correlation -----------------------------
pearson_mds = function(GRN, csGRN, name, pca_colors){
  library(ComplexHeatmap)
  library(igraph)
  edges_GRN <- as_data_frame(GRN, what = "edges")
  rownames(edges_GRN) <- paste0(edges_GRN$from,"_",edges_GRN$to)
  output <- matrix(0, ncol=length(names(csGRN)), nrow=nrow(edges_GRN),
                   dimnames = list(rownames(edges_GRN),names(csGRN)))
  
  for(i in names(csGRN)){
    edges_csGRN_df <- as_data_frame(csGRN[[i]], what = "edges")
    edges_csGRN <- paste0(edges_csGRN_df$from,"_",edges_csGRN_df$to)
    output[which(rownames(edges_GRN)%in%edges_csGRN),i] <- 1
  }  
  
  cor_output <- cor(output, method="pearson")
  data_to_plot <- cor_output

  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_",name,"_correlation_heatmap_and_mds.pdf",sep=""))
  bcell_colors <- pca_colors
  
  ha = HeatmapAnnotation(stage = colnames(output),
                         col = list(stage=bcell_colors))
  
  rha = rowAnnotation(stage = colnames(output),
                      col = list(stage=bcell_colors))
  
  draw(Heatmap(data_to_plot, column_title=name, 
               show_row_names = TRUE, show_column_names = TRUE,
               cluster_rows = FALSE, cluster_columns = FALSE,
               top_annotation = ha, left_annotation = rha))
  
  draw(Heatmap(data_to_plot, column_title=name, 
               show_row_names = TRUE, show_column_names = TRUE,
               cluster_rows = TRUE, cluster_columns = TRUE,
               top_annotation = ha, left_annotation = rha))
  
    fit <- cmdscale(1-data_to_plot,eig=TRUE, k=2) # k is the number of dim
  
  # plot solution
  mplot <- data.frame(PC1=fit$points[,1], PC2=fit$points[,2], Cell_type=names(csGRN))
  print(ggplot(mplot, aes(x = PC1, y = PC2, color = Cell_type)) +
    geom_point(size = 3, show.legend=FALSE) + 
    coord_fixed(ratio=1) + 
    scale_color_manual(values = bcell_colors) + 
    directlabels::geom_dl(aes(label = Cell_type), method = list("smart.grid", cex=1)) + 
    labs(title="MDS",x = "Coordinate 1", y = "Coordinate 2") +
    theme_classic())

    dev.off()
}
# ----------------------------------------------------------------------------- #


# Over representation analysis ------------------------------------------------
ORA=function(list_query,list_terms, x_lab="", y_lab=""){
  
  library(viridis)
  library(enrichplot)
  library(clusterProfiler)
  
  # Term to gene data frame
  c1 <- vector()
  c2 <- vector()
  for(i in 1:length(list_terms)){
    c1 <- c(c1,rep(names(list_terms)[i],length(list_terms[[i]])))
    c2 <- c(c2,list_terms[[i]])
  }
  Gx_GeneSet_x2 <- data.frame(term=c1,gene=c2)
  
  
  #ORA for each cell type
  ORA_Padj <- vector()
  ORA_GeneRatio <- vector()
  ORA_GR <- vector()
  ORA_Bg <- vector()
  
  Give_GR <- function(x) {
    if (!is.na(x)){
      separate <- strsplit(x, "/")
      division <- as.numeric(separate[[1]][1]) / as.numeric(separate[[1]][2])
      return(division)
    } else {
      return(NA)
    }
  }
  
  for(i in 1:length(list_query)){
    target_genes <- unique(list_query[[i]])
    
    ora_output <- enricher(gene = target_genes, 
                           pvalueCutoff = 1, 
                           pAdjustMethod = "BH",
                           minGSSize = 10, 
                           maxGSSize = 10000, 
                           qvalueCutoff = 1, 
                           TERM2GENE = Gx_GeneSet_x2, 
                           TERM2NAME = NA)
    ora_output_ord <- ora_output@result[names(list_terms),]
    ORA_Padj <- c(ORA_Padj, log10(ora_output_ord$p.adjust)*-1)
    ORA_GeneRatio <- c(ORA_GeneRatio,sapply(ora_output_ord$GeneRatio,Give_GR))
    ORA_GR <-  c(ORA_GR, ora_output_ord$GeneRatio)
    ORA_Bg <-  c(ORA_Bg, ora_output_ord$BgRatio)
  }
  
  ORA_Padj[ORA_Padj>10] <- 10
  
  data_to_plot <- data.frame(Cell=factor(rep(names(list_query),each=length(list_terms)),levels = names(list_query)), 
                             Group=factor(rep(names(list_terms), length(list_query)), levels=names(list_terms)),
                             Gene_ratio = ORA_GR,
                             Bg_ratio = ORA_Bg,
                             GR=as.numeric(ORA_GeneRatio),
                             Padj=ORA_Padj)
  
  data_to_plot <- data_to_plot[!data_to_plot$Cell=="global",]
  
  myPalette <- colorRampPalette(viridis(12))
  
  p <- ggplot(data_to_plot, aes(x=Cell, y=Group)) +
    geom_count(mapping=aes(color=GR, size=Padj))+
    scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))+
    scale_size(range = c(0, 10), limits = c(0,10))+
    labs(color='Gene Ratio', size="-log(p.adj)")+
    xlab(x_lab) +
    ylab(y_lab) +
    scale_y_discrete(limits=rev)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust=-0.05),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_discrete(position = "top")
  
  return(list(data_to_plot,p))
}  
ORA2=function(list_query,list_terms, x_lab="", y_lab="", y_elements_to_highlight=NULL, col_rect_highlight="gray", order="both"){
  
  library(viridis)
  library(enrichplot)
  library(clusterProfiler)
  
  # Term to gene data frame
  c1 <- vector()
  c2 <- vector()
  for(i in 1:length(list_terms)){
    c1 <- c(c1,rep(names(list_terms)[i],length(list_terms[[i]])))
    c2 <- c(c2,list_terms[[i]])
  }
  Gx_GeneSet_x2 <- data.frame(term=c1,gene=c2)
  
  
  #ORA for each cell type
  ORA_Padj <- vector()
  ORA_GeneRatio <- vector()
  ORA_GR <- vector()
  ORA_Bg <- vector()
  
  Give_GR <- function(x) {
    if (!is.na(x)){
      separate <- strsplit(x, "/")
      division <- as.numeric(separate[[1]][1]) / as.numeric(separate[[1]][2])
      return(division)
    } else {
      return(NA)
    }
  }
  
  for(i in 1:length(list_query)){
    target_genes <- unique(list_query[[i]])
    
    ora_output <- enricher(gene = target_genes, 
                           pvalueCutoff = 1, 
                           pAdjustMethod = "BH",
                           minGSSize = 10, 
                           maxGSSize = 10000, 
                           qvalueCutoff = 1, 
                           TERM2GENE = Gx_GeneSet_x2, 
                           TERM2NAME = NA)
    ora_output_ord <- ora_output@result[names(list_terms),]
    ORA_Padj <- c(ORA_Padj, log10(ora_output_ord$p.adjust)*-1)
    ORA_GeneRatio <- c(ORA_GeneRatio,sapply(ora_output_ord$GeneRatio,Give_GR))
    ORA_GR <-  c(ORA_GR, ora_output_ord$GeneRatio)
    ORA_Bg <-  c(ORA_Bg, ora_output_ord$BgRatio)
  }
  
  ORA_Padj[ORA_Padj>10] <- 10
  
  data_to_plot <- data.frame(Cell=factor(rep(names(list_query),each=length(list_terms)),levels = names(list_query)), 
                             Group=factor(rep(names(list_terms), length(list_query)), levels=names(list_terms)),
                             Gene_ratio = ORA_GR,
                             Bg_ratio = ORA_Bg,
                             GR=as.numeric(ORA_GeneRatio),
                             Padj=ORA_Padj)
  
  data_to_plot <- data_to_plot[!data_to_plot$Cell=="global",]
  
  if(order=="GR"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  
  if(order=="pval"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  if(order=="both"){
    dd1 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd1[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd1[is.na(dd1)] <- 0
    
    dd2 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd2[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd2[is.na(dd2)] <- 0
    
    dd <- dd2*dd1
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  myPalette <- colorRampPalette(viridis(12))
  
  p <- ggplot(data_to_plot, aes(x=Cell, y=Group)) +
    geom_count(mapping=aes(color=GR, size=Padj))+
    scale_colour_gradientn(colours = myPalette(100))+
    scale_size(range = c(0, 5), limits = c(0,10))+
    labs(color='Gene Ratio', size="-log(p.adj)")+
    xlab(x_lab) +
    ylab(y_lab) +
    scale_y_discrete(limits=rev)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = -45, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_discrete(position = "top")
  
  p2 <- NULL
  if(!is.null(y_elements_to_highlight)){
    sel <- which(rev(levels(p$data$Group))%in%y_elements_to_highlight)
    p2 <- p
    for(i in 1:length(sel)){
      p2 <- p2 + annotate("rect",xmin = 0, xmax = 9, ymin = sel[i]-0.5, ymax = sel[i]+0.5, fill = col_rect_highlight, alpha = 0.2)
    }
  }
  return(list(data_to_plot,p,p2))
}  

ORA2_rev=function(list_query,list_terms, x_lab="", y_lab="", y_elements_to_highlight=NULL, col_rect_highlight="gray", order="both"){
  
  library(viridis)
  library(enrichplot)
  library(clusterProfiler)
  
  # Term to gene data frame
  c1 <- vector()
  c2 <- vector()
  for(i in 1:length(list_terms)){
    c1 <- c(c1,rep(names(list_terms)[i],length(list_terms[[i]])))
    c2 <- c(c2,list_terms[[i]])
  }
  Gx_GeneSet_x2 <- data.frame(term=c1,gene=c2)
  
  
  #ORA for each cell type
  ORA_Padj <- vector()
  ORA_GeneRatio <- vector()
  ORA_GR <- vector()
  ORA_Bg <- vector()
  
  Give_GR <- function(x) {
    if (!is.na(x)){
      separate <- strsplit(x, "/")
      division <- as.numeric(separate[[1]][1]) / as.numeric(separate[[1]][2])
      return(division)
    } else {
      return(NA)
    }
  }
  
  for(i in 1:length(list_query)){
    target_genes <- unique(list_query[[i]])
    
    ora_output <- enricher(gene = target_genes, 
                           pvalueCutoff = 1, 
                           pAdjustMethod = "BH",
                           minGSSize = 10, 
                           maxGSSize = 10000, 
                           qvalueCutoff = 1, 
                           TERM2GENE = Gx_GeneSet_x2, 
                           TERM2NAME = NA)
    ora_output_ord <- ora_output@result[names(list_terms),]
    ORA_Padj <- c(ORA_Padj, log10(ora_output_ord$p.adjust)*-1)
    ORA_GeneRatio <- c(ORA_GeneRatio,sapply(ora_output_ord$GeneRatio,Give_GR))
    ORA_GR <-  c(ORA_GR, ora_output_ord$GeneRatio)
    ORA_Bg <-  c(ORA_Bg, ora_output_ord$BgRatio)
  }
  
  ORA_Padj[ORA_Padj>10] <- 10
  
  data_to_plot <- data.frame(Cell=factor(rep(names(list_query),each=length(list_terms)),levels = names(list_query)), 
                             Group=factor(rep(names(list_terms), length(list_query)), levels=names(list_terms)),
                             Gene_ratio = ORA_GR,
                             Bg_ratio = ORA_Bg,
                             GR=as.numeric(ORA_GeneRatio),
                             Padj=ORA_Padj)
  
  data_to_plot <- data_to_plot[!data_to_plot$Cell=="global",]
  
  if(order=="GR"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  
  if(order=="pval"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  if(order=="both"){
    dd1 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd1[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd1[is.na(dd1)] <- 0
    
    dd2 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd2[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd2[is.na(dd2)] <- 0
    
    dd <- dd2*dd1
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  myPalette <- colorRampPalette(viridis(12))
  
  data_to_plot2 <- data_to_plot
  data_to_plot2$Padj[data_to_plot2$Padj<1.3] <- 0
  
  aa <- by(data_to_plot2$Padj,data_to_plot2$Group,sum, na.rm=TRUE)
  TF_to_keep <- names(aa)[aa>0]
  data_to_plot3 <- data_to_plot2[data_to_plot2$Group%in%TF_to_keep,]
  data_to_plot3$Group <- droplevels(data_to_plot3$Group)
  
  p <- ggplot(data_to_plot3, aes(x=Cell, y=Group)) +
    geom_count(mapping=aes(color=GR, size=Padj))+
    scale_colour_gradientn(colours = myPalette(100))+
    scale_size(range = c(0, 5), limits = c(1,10))+
    labs(color='Gene Ratio', size="-log(p.adj)")+
    xlab(x_lab) +
    ylab(y_lab) +
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -45, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_discrete(position = "top")
  
  p2 <- NULL
  if(!is.null(y_elements_to_highlight)){
    sel <- which(rev(levels(p$data$Group))%in%y_elements_to_highlight)
    p2 <- p
    for(i in 1:length(sel)){
      p2 <- p2 + annotate("rect",xmin = 0, xmax = 9, ymin = sel[i]-0.5, ymax = sel[i]+0.5, fill = col_rect_highlight, alpha = 0.2)
    }
  }
  return(list(data_to_plot,p,p2))
}  


ORA2_rev2=function(list_query,list_terms, x_lab="", y_lab="", TF_order, y_elements_to_highlight=NULL, col_rect_highlight="gray", order="both"){
  
  library(viridis)
  library(enrichplot)
  library(clusterProfiler)
  
  # Term to gene data frame
  c1 <- vector()
  c2 <- vector()
  for(i in 1:length(list_terms)){
    c1 <- c(c1,rep(names(list_terms)[i],length(list_terms[[i]])))
    c2 <- c(c2,list_terms[[i]])
  }
  Gx_GeneSet_x2 <- data.frame(term=c1,gene=c2)
  
  
  #ORA for each cell type
  ORA_Padj <- vector()
  ORA_GeneRatio <- vector()
  ORA_GR <- vector()
  ORA_Bg <- vector()
  
  Give_GR <- function(x) {
    if (!is.na(x)){
      separate <- strsplit(x, "/")
      division <- as.numeric(separate[[1]][1]) / as.numeric(separate[[1]][2])
      return(division)
    } else {
      return(NA)
    }
  }
  
  for(i in 1:length(list_query)){
    target_genes <- unique(list_query[[i]])
    
    ora_output <- enricher(gene = target_genes, 
                           pvalueCutoff = 1, 
                           pAdjustMethod = "BH",
                           minGSSize = 10, 
                           maxGSSize = 10000, 
                           qvalueCutoff = 1, 
                           TERM2GENE = Gx_GeneSet_x2, 
                           TERM2NAME = NA)
    ora_output_ord <- ora_output@result[names(list_terms),]
    ORA_Padj <- c(ORA_Padj, log10(ora_output_ord$p.adjust)*-1)
    ORA_GeneRatio <- c(ORA_GeneRatio,sapply(ora_output_ord$GeneRatio,Give_GR))
    ORA_GR <-  c(ORA_GR, ora_output_ord$GeneRatio)
    ORA_Bg <-  c(ORA_Bg, ora_output_ord$BgRatio)
  }
  
  ORA_Padj[ORA_Padj>10] <- 10
  
  data_to_plot <- data.frame(Cell=factor(rep(names(list_query),each=length(list_terms)),levels = names(list_query)), 
                             Group=factor(rep(names(list_terms), length(list_query)), levels=names(list_terms)),
                             Gene_ratio = ORA_GR,
                             Bg_ratio = ORA_Bg,
                             GR=as.numeric(ORA_GeneRatio),
                             Padj=ORA_Padj)
  
  data_to_plot <- data_to_plot[!data_to_plot$Cell=="global",]
  
  if(order=="GR"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  
  if(order=="pval"){
    dd <- matrix(NA,
                 ncol=length(levels(data_to_plot$Cell)),
                 nrow=length(levels(data_to_plot$Group)),
                 dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd[is.na(dd)] <- 0
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  if(order=="both"){
    dd1 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd1[,i] <- data_to_plot$Padj[data_to_plot$Cell==i]
    }
    
    dd1[is.na(dd1)] <- 0
    
    dd2 <- matrix(NA,
                  ncol=length(levels(data_to_plot$Cell)),
                  nrow=length(levels(data_to_plot$Group)),
                  dimnames = list(levels(data_to_plot$Group),levels(data_to_plot$Cell)))
    
    for(i in levels(data_to_plot$Cell)){
      dd2[,i] <- data_to_plot$GR[data_to_plot$Cell==i]
    }
    
    dd2[is.na(dd2)] <- 0
    
    dd <- dd2*dd1
    rr <- hclust(dist(dd))
    
    data_to_plot$Group <- factor(data_to_plot$Group, levels=rownames(dd)[rr$order])
  }
  
  data_to_plot$Group <- factor(data_to_plot$Group, levels=TF_order)
  
  myPalette <- colorRampPalette(viridis(12))
  
  data_to_plot2 <- data_to_plot
  data_to_plot2$Padj[data_to_plot2$Padj<1.3] <- 0
  
  aa <- by(data_to_plot2$Padj,data_to_plot2$Group,sum, na.rm=TRUE)
  TF_to_keep <- names(aa)[aa>0]
  data_to_plot3 <- data_to_plot2[data_to_plot2$Group%in%TF_to_keep,]
  data_to_plot3$Group <- droplevels(data_to_plot3$Group)
  
  p <- ggplot(data_to_plot3, aes(x=Cell, y=Group)) +
    geom_count(mapping=aes(color=GR, size=Padj))+
    scale_colour_gradientn(colours = myPalette(100))+
    scale_size(range = c(0, 5), limits = c(1,10))+
    labs(color='Gene Ratio', size="-log(p.adj)")+
    xlab(x_lab) +
    ylab(y_lab) +
    scale_y_discrete(limits=rev)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = -45, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_discrete(position = "top")
  
  p2 <- NULL
  if(!is.null(y_elements_to_highlight)){
    sel <- which(rev(levels(p$data$Group))%in%y_elements_to_highlight)
    p2 <- p
    for(i in 1:length(sel)){
      p2 <- p2 + annotate("rect",xmin = 0, xmax = 9, ymin = sel[i]-0.5, ymax = sel[i]+0.5, fill = col_rect_highlight, alpha = 0.2)
    }
  }
  return(list(data_to_plot,p,p2))
}  


# ----------------------------------------------------------------------------- #


# Peak overlap  ---------------------------------------------------------------
my_enrichPeakOverlap=function(gr1,gr2, n_ocrs){
  cat("....Starting enrichment\n \n")
  a <-findOverlaps(gr1,gr2)
  f11 <- length(a@from)
  f01 <- a@nLnode-f11
  f10 <- a@nRnode-f11
  f00 <- n_ocrs - sum(f11,f01,f10)
  tt <- matrix(c(f00,f01,f10,f11),ncol=2, byrow=TRUE)
  colnames(tt) <- c(0,1)
  rownames(tt) <- c(0,1)
  print(tt)
  cat("subject (row; gr2) x query (col; gr1)\n")
  cat("\n.... Fisher test\n")
  ftt <- fisher.test(tt)
  print(ftt)
  cat("....End enrichment\n")
  return(c(pval=ftt$p.value, estimate=ftt$estimate, f11=f11, f01=f01, f10=f10, f00=f00, Lquery=a@nLnode, Lsubj=a@nRnode))
}
# ----------------------------------------------------------------------------- #
