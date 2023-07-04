
bw_plot=function(list_data, peaks, gene, ylims, start, end, chr){
  genomeAxis <- GenomeAxisTrack(name="MyAxis") 
  customFromTxDb <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, start=start, end=end, chromosome=chr, name = "Gene Region") 
  atrack <- AnnotationTrack(peaks[which(seqnames(peaks)==chr),], start=start, end=end, name=gene)
  plotTracks(c(atrack,list_data,customFromTxDb,genomeAxis), 
             from=start,to=end, 
             chromosome=chr,type="l", ylim=ylims, main=gene)
}         


bw_plot2=function(list_data, peaks, gene, ylims, start, end, chr,h3k4me1, h3k4me3, h3k27ac){
  genomeAxis <- GenomeAxisTrack(name="MyAxis") 
  customFromTxDb <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, start=start, end=end, chromosome=chr, name = "Gene Region") 
  atrack <- AnnotationTrack(peaks[which(seqnames(peaks)==chr),], start=start, end=end, name=gene)
  atrack1 <- AnnotationTrack(h3k4me1[which(seqnames(h3k4me1)==chr),], start=start, end=end, name="H3K4me1")
  atrack2 <- AnnotationTrack(h3k4me3[which(seqnames(h3k4me3)==chr),], start=start, end=end, name="H3k4me3")
  atrack3 <- AnnotationTrack(h3k27ac[which(seqnames(h3k27ac)==chr),], start=start, end=end, name="H3k27ac")
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  plotTracks(c(itrack,atrack,list_data,customFromTxDb,atrack2,atrack1,atrack3,genomeAxis), 
             from=start,to=end, 
             chromosome=chr,type="l", ylim=ylims, main=gene)
}         


bw_plot3=function(list_data, peaks, gene, ylims, start, end, chr,h3k4me1, h3k4me3, h3k27ac, snp){
  genomeAxis <- GenomeAxisTrack(name="MyAxis") 
  customFromTxDb <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, start=start, end=end, chromosome=chr, name = "Gene Region") 
  atrack <- AnnotationTrack(peaks[which(seqnames(peaks)==chr),], start=start, end=end, name=gene)
  atrack1 <- AnnotationTrack(h3k4me1[which(seqnames(h3k4me1)==chr),], start=start, end=end, name="H3K4me1")
  atrack2 <- AnnotationTrack(h3k4me3[which(seqnames(h3k4me3)==chr),], start=start, end=end, name="H3k4me3")
  atrack3 <- AnnotationTrack(h3k27ac[which(seqnames(h3k27ac)==chr),], start=start, end=end, name="H3k27ac")
  atrack_snp <- AnnotationTrack(snp[which(seqnames(snp)==chr),], start=start, end=end, name=gene)
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  plotTracks(c(itrack,atrack_snp,atrack,list_data,customFromTxDb,atrack2,atrack1,atrack3,genomeAxis), 
             from=start,to=end, 
             chromosome=chr,type="l", ylim=ylims, main=gene)
}         


ORA=function(list_query,list_terms){
  
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
    xlab("Bcell type") +
    ylab("Leukemia Subgroup") +
    scale_y_discrete(limits=rev)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust=-0.05),
          plot.margin = unit(c(1,1,1,1), "cm"))+
    scale_x_discrete(position = "top")
  
  return(list(data_to_plot,p))
}  


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


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}