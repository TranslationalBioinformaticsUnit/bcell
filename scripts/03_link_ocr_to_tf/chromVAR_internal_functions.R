#******************************************************#
# Aggregated archetype score - with chromVAR           #
# script: 17_archetype_chromVAR.R                      #
#******************************************************#


##>>>>>>>>>>>>>>>>>>>>> chromVAR_assay function
# Run the full analysis of chromVAR
chromVAR_assay = function(ocr_granges, ocr_counts, ocr_metadata, ocr_motif_anno, assay_name, target_tfs){
  #Arguments:
  #ocr_granges -> granges with OCRs to explore
  #ocr_counts -> counts for OCRs provided in ocr_granges
  #ocr_metadata -> metadata of OCRs provided in ocr_granges
  #ocr_motif_anno -> binary matrix with OCR-Motif information
  #assay_name -> name of analysis done
  #target_tfs -> tf name or general name to highlight in heatmap plot and ridgeplot
  
  #Generate summarized experiment object
  fragment_counts <- SummarizedExperiment(assays = 
                                            list(counts = ocr_counts[names(ocr_granges),]),
                                          rowRanges = ocr_granges,
                                          colData = ocr_metadata)
  #Adding GC content
  fragment_counts2 <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  #Filtering peaks
  counts_filtered <- filterPeaks(fragment_counts2, non_overlapping = TRUE)
  
  anno_ix <- getAnnotations(ocr_motif_anno[rownames(counts_filtered),], 
                            rowRanges = rowRanges(counts_filtered))
  
  # computing deviations
  dev <- computeDeviations(object = counts_filtered, 
                           annotations = anno_ix)
  save(list=c("dev"),file=paste(format(Sys.Date(),"%Y%m%d"),"_",assay_name,"_dev.RData",sep=""))

  # variability
  variability <- computeVariability(dev)
  save(list=c("variability"),file=paste(format(Sys.Date(),"%Y%m%d"),"_",assay_name,"_variability.RData",sep=""))
  cat("Top variable motifs:\n")
  print(variability[order(variability$variability, decreasing = TRUE)[1:10],])
  plotVariability(variability, use_plotly = FALSE)
  
  # clustering samples
  sample_cor <- getSampleCorrelation(dev)
  #top annotation
  bcell_colors <- jdb_color_map(c("MPP","LMPP","CD8","CD4","GMP-B","Ery","mono","MEP"))
  names(bcell_colors) <- levels(ocr_metadata$Tissue)
  ha = HeatmapAnnotation(stage = ocr_metadata$Tissue,
                         col = list(stage =bcell_colors))
  #plot color
  range_color <- as.character(jdb_palette("solar_flare",9))
  f1 <- colorRamp2(seq(-1, 1, length = 9), range_color)
  
  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_",assay_name,"_sample_correlation.pdf",sep=""))
  draw(Heatmap(sample_cor, top_annotation = ha, 
          show_row_names = FALSE, show_column_names =FALSE, 
          row_names_gp = gpar(fontsize = 3), col=f1, 
          name="rho",
          column_title = "Correlation for bias corrected deviations"))
  dev.off()
  
  # Get deviations and scores
  #bias corrected deviations in accessibility
  TF_deviations <- deviations(dev) 
  colnames(TF_deviations) <- gsub("-",".",colnames(TF_deviations))
  #deviationScores are the Z-scores for each bias corrected deviations
  TF_scores <- deviationScores(dev)
  colnames(TF_scores) <- gsub("-",".",colnames(TF_scores))
  
  # Graphical exploration
  #heatmap
  ha = HeatmapAnnotation(stage = ocr_metadata$Tissue,
                         col = list(stage =bcell_colors))
  
  target_TF <- target_tfs
  target_TF_position <- vector()
  for(i in 1:length(target_TF)){
    target_TF_position <- c(target_TF_position,grep(target_TF[i],rownames(TF_scores)))
  }
  
  rha = rowAnnotation(foo = anno_mark(at = target_TF_position, 
                                      labels = rownames(TF_scores)[target_TF_position],
                                      labels_gp = gpar(fontsize = 6)))
  
  # heatmap color 
  range_color <- as.character(jdb_palette("brewer_celsius",9))
  f1 <- colorRamp2(seq(-1*min(abs(min(TF_scores)),max(TF_scores)), min(abs(min(TF_scores)),max(TF_scores)), length = 9), range_color)
  
  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_",assay_name,"_TFBS_heatmap.pdf",sep=""))
  draw(Heatmap(TF_scores, top_annotation = ha, right_annotation = rha, 
          show_row_names = FALSE, show_column_names =FALSE, 
          row_names_gp = gpar(fontsize = 3), col=f1, 
          name="Z-score",
          column_title = "TFBS accessibility score"))
  
  draw(Heatmap(TF_scores, top_annotation = ha, 
          show_row_names = TRUE, show_column_names =FALSE, 
          row_names_gp = gpar(fontsize = 1), col=f1, 
          name="Z-score",
          column_title = "TFBS accessibility score"))
  
  draw(Heatmap(TF_scores, top_annotation = ha, right_annotation = rha,
          show_row_names = FALSE, show_column_names =FALSE, 
          column_names_gp = gpar(fontsize = 6), col=f1, 
          name="Z-score",
          column_title = "TFBS accessibility score",
          column_order = order(ocr_metadata$Tissue)))
  dev.off()
  
  #ridgeplot for specific TFs
  library(ggridges)
  library(hrbrthemes)
  library(viridis)
  
  group <- factor(as.character(ocr_metadata$Tissue),
                  levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))
  
  target_TF <- target_tfs
  
  pdf(paste(format(Sys.Date(),"%Y%m%d"),"_",assay_name,"_TFBS_ridgeplot.pdf",sep=""))
  for(i in 1:length(target_TF)){
    target_TF_position <- grep(target_TF[i],rownames(TF_scores))
    if(length(target_TF_position)>0){
      for(j in 1:length(target_TF_position)){
        M <- data.frame(cell_type=group, score=TF_scores[target_TF_position[j],])
        TF_name <- rownames(TF_scores)[target_TF_position[j]]
        print(ggplot(M, aes(x = score, y = cell_type, fill= ..x..)) +
                geom_density_ridges_gradient(scale=2) +
                scale_fill_gradientn(colors = jdb_palette("brewer_celsius")) +
                labs(title = paste(TF_name," - TFBS accessibility score",sep="")) +
                theme_ridges())
      }
    }
  }
  dev.off()
  
}
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> TFBS_RNA_heatmap function
## function to plot the heatmaps 
TFBS_RNA_heatmap = function(TF_gene_association, atac_metadata, rna_metadata, color_group, atac_group, rna_group, tfbs_data, rna_data, add_raw_rna, tfbs_variability, rna_sd, target_tfs, assay_name){
  #correlation color
  col_fun1 = colorRamp2(seq(-1, 1, length=9), jdb_palette("solar_flare",9,"discrete"))
  #variability color
  col_fun2 = colorRamp2(seq(min(tfbs_variability$variability),max(tfbs_variability$variability), length = 9), jdb_palette("brewer_violet",9,"discrete"))
  #sd color
  col_fun3 = colorRamp2(seq(min(rna_sd, na.rm=TRUE),max(rna_sd, na.rm=TRUE), length = 9), jdb_palette("brewer_marine",9,"discrete"))
  
  ##TFBS heatmap
  # Column annotation
  ha = HeatmapAnnotation(stage = atac_metadata[,atac_group],
                         col = list(stage =color_group), show_annotation_name = FALSE)
  # Row annotation
  rha = rowAnnotation(tfbs_var = tfbs_variability[rownames(tfbs_data),"variability"],
                      rna_cv = rna_sd,
                      rho = TF_gene_association[rownames(tfbs_data),"spearman_rho"],
                      col = list(tfbs_var=col_fun2, rna_cv=col_fun3, rho=col_fun1))
  
  #Split positive vs negative correlation
  split <- TF_gene_association$spearman_rho
  split[TF_gene_association$spearman_rho>0] <- 1
  split[TF_gene_association$spearman_rho<0] <- 2
  
  # Color 
  f1 <- colorRamp2(seq(-1*min(abs(min(tfbs_data)),max(tfbs_data)), min(abs(min(tfbs_data)),max(tfbs_data)), length = 9), jdb_palette("brewer_celsius",9,"discrete"))
  # Plot
  ht_tfbs <- Heatmap(tfbs_data, top_annotation = ha, right_annotation = rha,
                     show_row_names = FALSE, show_column_names =FALSE, 
                     column_names_gp = gpar(fontsize = 6),
                     row_names_gp = gpar(fontsize = 3),
                     col=f1, 
                     row_split = split,
                     name="Z-score TFBS",
                     column_title = "TFBS score",
                     column_title_gp = gpar(fontsize = 9),
                     column_order = order(atac_metadata[,atac_group]))
  
  ##RNA heatmap
  # Column annotation
  ha = HeatmapAnnotation(stage = rna_metadata[,rna_group],
                         col = list(stage =color_group), show_annotation_name = FALSE)
  # Row annotation
  target_TF <- target_tfs
  target_TF_position <- vector()
  for(i in 1:length(target_TF)){
    target_TF_position <- c(target_TF_position,which(rownames(rna_data)==target_TF[i]))
  }
  rha = rowAnnotation(foo = anno_mark(at = target_TF_position, 
                                      labels = target_TF,
                                      labels_gp = gpar(fontsize = 7)))
  
  if(add_raw_rna==FALSE){
    #Plot
    M <- t(scale(t(rna_data),center=TRUE, scale=TRUE))
    ht_rna <- Heatmap(M, top_annotation = ha, right_annotation = rha, 
                      show_row_names = TRUE, show_column_names =FALSE, 
                      row_names_gp = gpar(fontsize = 4), #col=f1, 
                      name="Z-score RNA",
                      row_split = split,
                      row_names_side = "right",
                      column_title = "RNA",
                      column_title_gp = gpar(fontsize = 9),
                      column_order = order(rna_metadata[,rna_group]))
    
    # Joint plot
    pdf(paste(format(Sys.Date(),"%Y%m%d"),"_tfbs_and_rna_heatmap_",assay_name,".pdf",sep=""))
    draw(ht_tfbs + ht_rna)
    dev.off()
  }
  
  if(add_raw_rna==TRUE){
    #Plot
    M <- t(scale(t(rna_data),center=TRUE, scale=TRUE))
    ht_rna <- Heatmap(M, top_annotation = ha,  
                      show_row_names = FALSE, show_column_names =FALSE, 
                      row_names_gp = gpar(fontsize = 4), #col=f1, 
                      name="Z-score RNA",
                      row_split = split,
                      row_names_side = "right",
                      column_title = "RNA",
                      column_title_gp = gpar(fontsize = 9),
                      column_order = order(rna_metadata[,rna_group]))
    
    # Joint plot
    #Plot
    M <- rna_data
    ht_raw_rna <- Heatmap(M, top_annotation = ha, right_annotation = rha, 
                          show_row_names = TRUE, show_column_names =FALSE, 
                          row_names_gp = gpar(fontsize = 4), 
                          col=jdb_palette("flame_macaw"), 
                          name="RNA level",
                          row_split = split,
                          row_names_side = "right",
                          column_title = "",
                          column_order = order(rna_metadata[,rna_group]))
    
    pdf(paste(format(Sys.Date(),"%Y%m%d"),"_tfbs_and_rna_heatmap_",assay_name,".pdf",sep=""))
    draw(ht_tfbs + ht_rna + ht_raw_rna)
    dev.off()
  }
}
#-------------------------------------