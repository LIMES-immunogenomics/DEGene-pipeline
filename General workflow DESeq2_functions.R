list.of.packages <- c("tidyverse", "reshape2","dplyr","stringr","pheatmap","colourlovers",'VennDiagram',"ggbeeswarm","scales","magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
list.of.bioc.packages<- c("rhdf5","DESeq2","IHW","clusterProfiler","GSEABase","DOSE","AnnotationDbi","sva","tximport","limma","geneplotter","genefilter","org.Mm.eg.db","org.Hs.eg.db","biomaRt")
new.packages.bioc <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
source("https://bioconductor.org/biocLite.R")
if(length(new.packages.bioc)>0) biocLite(new.packages.bioc,suppressUpdates=TRUE)
lapply(c(list.of.packages,list.of.bioc.packages), require, character.only = TRUE)

df <- data.frame(x=c(48:129))


filter_counts <- function(count_mat=count_matrix, 
                          keytype_opt="ENSEMBL",
                          change_names=TRUE,
                          organism_option=org.Mm.eg.db,
                          threshold=10){
  list_genes <- list()
  
  list_genes_use <- rowSums(count_mat)>=threshold
  
  count_matrix <- count_mat[list_genes_use,]
  
  if(change_names==TRUE){
    rownames(count_matrix) <- mapIds(x=organism_option,
                                     keys=rownames(count_matrix),
                                     column="SYMBOL",
                                     keytype=keytype_opt,
                                     multiVals="first")
    
    count_matrix <- as.data.frame(count_matrix[!is.na(rownames(count_matrix)),])
    count_matrix<-aggregate(.~rownames(count_matrix), data=count_matrix, sum)
    rownames(count_matrix) <- count_matrix[,1]
    count_matrix[,1]<-NULL
    count_matrix[1:ncol(count_matrix)]<- lapply(count_matrix, as.integer)
    count_matrix
  }
  else{
    count_matrix[1:ncol(count_matrix)]<- lapply(count_matrix, as.integer)
    count_matrix
  }
}

Var_stab<- function(deseq_count_matrix=dsfm, blind_param=T, count_mat=count_matrix){
  if(ncol(count_mat)<40){
    rld<-rlog(deseq_count_matrix, blind=blind_param)
  }
  else{
    vst <- vst(deseq_count_matrix, blind = blind_param)
  }
}
Checking_normalized_counts <- function(raw_data=count_matrix, 
                                       normalized_source=dds, 
                                       scaling_opt="log2",
                                       anno_colour=anno,
                                       condition=condition){
  if (is.null(anno_colour)==F){
  if(scaling_opt=="log2"){ 
    dat <- stack(as.data.frame(log2(raw_data+1)))
    colnames(dat) <- c("Reads", "Sample")
    dat$new <- annotation[[paste0(condition)]][dat$Sample]
    legend_title=paste0(condition)
    raw_counts_log2_boxplots <- ggplot(dat, aes(x=Sample,y=Reads, fill=new))+
      geom_boxplot()+
      labs(title="Unnormalized reads")+
      ylab("Reads (log2)")+
      scale_fill_manual(values=anno_colour,name=legend_title)+
      theme(panel.background = element_rect(fill=NA, color = "black"),
            plot.background = element_blank(),
            aspect.ratio = 0.5,
            legend.background = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = NA, color = "black"),
            axis.text.x = element_text(angle = 90, vjust=0.5),
            axis.line = element_line(colour = "black", size = 1),
            axis.text = element_text(colour = "black", size = 10, face = "bold"),
            axis.ticks = element_line(color = "black",size = 1),
            axis.title = element_text(colour = "black", size = 10, face = "bold"),
            axis.title.x = element_blank(),
            legend.key = element_rect(fill = "white"))
    
    
    normalized_data <- stack(as.data.frame(log2(counts(normalized_source, normalized=T)+1)))
    normalized_data$values <- as.integer(normalized_data$values)
    colnames(normalized_data) <- c("Reads", "Sample")
    normalized_data$new <- annotation[[paste0(condition)]][normalized_data$Sample]
    legend_title=paste0(condition)
    normalized_counts_log2_boxplots <- ggplot(normalized_data, aes(x=Sample,y=Reads, fill=new))+
      geom_boxplot()+
      labs(title="Normalized reads")+
      scale_fill_manual(values=anno_colour,name=legend_title)+
      ylab("Reads (log2)")+
      theme(panel.background = element_rect(fill=NA, color = "black"),
            plot.background = element_blank(),
            aspect.ratio = 0.5,
            legend.background = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = NA, color = "black"),
            axis.text.x = element_text(angle = 90, vjust=0.5),
            axis.line = element_line(colour = "black", size = 1),
            axis.text = element_text(colour = "black", size = 10, face = "bold"),
            axis.ticks = element_line(color = "black",size = 1),
            axis.title.x = element_blank(),
            axis.title = element_text(colour = "black", size = 10, face = "bold"),
            legend.key = element_rect(fill = "white"))
    
    multiplot(raw_counts_log2_boxplots, normalized_counts_log2_boxplots, cols=2)
  }
  else{
    # using untransformed data to plot
    dat <- stack(as.data.frame(raw_data+1))
    colnames(dat) <- c("Reads", "Sample")
    dat$new <- annotation[[paste0(condition)]][dat$Sample]
    legend_title=paste0(condition)
    raw_counts_boxplot <- ggplot(dat, aes(x=Sample,y=Reads, fill=new))+
      geom_boxplot()+
      scale_y_log10()+
      labs(title="Unnormalized reads")+
      scale_fill_manual(values=anno_colour,name=legend_title)+
      theme(panel.background = element_rect(fill=NA, color = "black"),
            plot.background = element_blank(),
            aspect.ratio = 0.5,
            legend.background = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = NA, color = "black"),
            axis.text.x = element_text(angle = 90, vjust=0.5),
            axis.line = element_line(colour = "black", size = 1),
            axis.text = element_text(colour = "black", size = 10, face = "bold"),
            axis.ticks = element_line(color = "black",size = 1),
            axis.title = element_text(colour = "black", size = 10, face = "bold"),
            axis.title.x = element_blank(),
            legend.key = element_rect(fill = "white"))
    
    
    normalized_data <- stack(as.data.frame(counts(normalized_source, normalized=T)+1))
    normalized_data$values <- as.integer(normalized_data$values)
    colnames(normalized_data) <- c("Reads", "Sample")
    normalized_data$new <- annotation[[paste0(condition)]][normalized_data$Sample]
    legend_title=paste0(condition)
    normalized_counts_boxplot <- ggplot(normalized_data, aes(x=Sample,y=Reads, fill=new))+
      geom_boxplot()+
      scale_y_log10()+
      labs(title="Normalized reads")+
      scale_fill_manual(values=anno_colour,name=legend_title)+
      theme(panel.background = element_rect(fill=NA, color = "black"),
            plot.background = element_blank(),
            aspect.ratio = 0.5,
            legend.background = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            strip.background = element_rect(fill = NA, color = "black"),
            axis.text.x = element_text(angle = 90, vjust=0.5),
            axis.line = element_line(colour = "black", size = 1),
            axis.text = element_text(colour = "black", size = 10, face = "bold"),
            axis.ticks = element_line(color = "black",size = 1),
            axis.title.x = element_blank(),
            axis.title = element_text(colour = "black", size = 10, face = "bold"),
            legend.key = element_rect(fill = "white"))
    
    multiplot(raw_counts_boxplot, normalized_counts_boxplot, cols=2)
  }}
  #option if no colors given
  else{
    if(scaling_opt=="log2"){ 
      dat <- stack(as.data.frame(log2(raw_data+1)))
      colnames(dat) <- c("Reads", "Sample")
      dat$new <- annotation[[paste0(condition)]][dat$Sample]
      legend_title=paste0(condition)
      raw_counts_log2_boxplots <- ggplot(dat, aes(x=Sample,y=Reads, fill=new))+
        geom_boxplot()+
        labs(title="Unnormalized reads")+
        ylab("Reads (log2)")+
        scale_fill_brewer(name=legend_title)+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              plot.background = element_blank(),
              aspect.ratio = 0.5,
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.text.x = element_text(angle = 90, vjust=0.5),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              axis.title.x = element_blank(),
              legend.key = element_rect(fill = "white"))
      
      
      normalized_data <- stack(as.data.frame(log2(counts(normalized_source, normalized=T)+1)))
      normalized_data$values <- as.integer(normalized_data$values)
      colnames(normalized_data) <- c("Reads", "Sample")
      normalized_data$new <- annotation[[paste0(condition)]][normalized_data$Sample]
      legend_title=paste0(condition)
      normalized_counts_log2_boxplots <- ggplot(normalized_data, aes(x=Sample,y=Reads, fill=new))+
        geom_boxplot()+
        labs(title="Normalized reads")+
        scale_fill_brewer(name=legend_title)+
        ylab("Reads (log2)")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              plot.background = element_blank(),
              aspect.ratio = 0.5,
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.text.x = element_text(angle = 90, vjust=0.5),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title.x = element_blank(),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white"))
      
      multiplot(raw_counts_log2_boxplots, normalized_counts_log2_boxplots, cols=2)
    }
    else{
      # using untransformed data to plot
      dat <- stack(as.data.frame(raw_data+1))
      colnames(dat) <- c("Reads", "Sample")
      dat$new <- annotation[[paste0(condition)]][dat$Sample]
      legend_title=paste0(condition)
      raw_counts_boxplot <- ggplot(dat, aes(x=Sample,y=Reads, fill=new))+
        geom_boxplot()+
        scale_y_log10()+
        labs(title="Unnormalized reads")+
        scale_fill_brewer(name=legend_title)+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              plot.background = element_blank(),
              aspect.ratio = 0.5,
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.text.x = element_text(angle = 90, vjust=0.5),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              axis.title.x = element_blank(),
              legend.key = element_rect(fill = "white"))
      
      
      normalized_data <- stack(as.data.frame(counts(normalized_source, normalized=T)+1))
      normalized_data$values <- as.integer(normalized_data$values)
      colnames(normalized_data) <- c("Reads", "Sample")
      normalized_data$new <- annotation[[paste0(condition)]][normalized_data$Sample]
      legend_title=paste0(condition)
      normalized_counts_boxplot <- ggplot(normalized_data, aes(x=Sample,y=Reads, fill=new))+
        geom_boxplot()+
        scale_y_log10()+
        labs(title="Normalized reads")+
        scale_fill_brewer(name=legend_title)+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              plot.background = element_blank(),
              aspect.ratio = 0.5,
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.text.x = element_text(angle = 90, vjust=0.5),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title.x = element_blank(),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white"))
      
      multiplot(raw_counts_boxplot, normalized_counts_boxplot, cols=2)
    }
  }}
sample_cor<- function(rld, limit_set = c(0.9,1),midpoint_set = 0.95,condition="Treatment"){
  # melt correlatiopn matrix
  cor_mat<-cor(rld, method = "pearson")
  
  melt_correlation_matrix<-melt(as.matrix(cor_mat))
  annotation[[paste0(condition)]] <- as.character(annotation[[paste0(condition)]])
  
  melt_correlation_matrix$Var1 <- paste( melt_correlation_matrix$Var1, annotation[as.character(melt_correlation_matrix$Var1),][[paste0(condition)]], sep = " - " )
  melt_correlation_matrix$Var2 <- paste( melt_correlation_matrix$Var2, annotation[as.character(melt_correlation_matrix$Var2),][[paste0(condition)]], sep = " - " )
  #ggplot heatmap
  ggplot(melt_correlation_matrix, aes(Var2, Var1, fill = value)) + 
    geom_tile(color = "white") + 
    scale_fill_gradient2(limit = limit_set, 
                         midpoint = midpoint_set,low = "blue", high = "red", mid = "white",
                         name = "Pearson\nCorrelation") + 
    theme_minimal() + 
    labs(title="Pearson correlation between all conditions (ungrouped)")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          aspect.ratio = 1,
          axis.title.x = element_blank(), axis.title.y = element_blank())
}
sorted_pearson_correaltion<- function(rld,midpoint_set = 0.9995,limit_set = c(0.999,1),condition="Treatment"){
  cor_mat<-cor(rld, method = "pearson")
  #Reorder the correlation matrix
  dd <- as.dist((1-cor_mat)/2)
  hc <- hclust(dd)
  cormat <-cor_mat[hc$order, hc$order]
  # Melt the correlation matrix
  melted_cormat <- melt(cormat, na.rm = TRUE)
  annotation[[paste0(condition)]] <- as.character(annotation[[paste0(condition)]])
  melted_cormat$Var1 <- paste( melted_cormat$Var1, annotation[as.character(melted_cormat$Var1),][[paste0(condition)]], sep = " - " )
  melted_cormat$Var2 <- paste( melted_cormat$Var2, annotation[as.character(melted_cormat$Var2),][[paste0(condition)]], sep = " - " )
  melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = unique(as.character(melted_cormat$Var1)))
  melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = unique(as.character(melted_cormat$Var2)))
  ggplot(data=melted_cormat, aes(x=Var2, y=Var1, fill = value))+
    geom_tile(color = "white")+
    labs(title="Pearson correlation of all samples")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = midpoint_set,oob=squish, limit = limit_set, space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          aspect.ratio = 1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust=0.5))+
    coord_fixed()
} 
PCA_rawdata <-function(object=rld, color_obj="condition", 
                       ntop=500, PC_1=1, PC_2=2,
                       count_mat=count_matrix,
                       shape_opt="NULL",
                       continuous=F,
                       colour_gradient=c("red","green","blue")){
  if(ntop=="all"){
    pca <- prcomp(t(log10(as.data.frame(count_mat)+1))) 
  }
  else{
    # calculate the variance for each gene
    rv <- rowVars(assay(object))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(log10(as.data.frame(count_mat)+1)[select,]))
  }
  
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(color_obj %in% names(colData(object)))) {
    stop("the argument 'color_obj' should specify columns of colData(dds)")
  }
  
  color_obj.df <- as.data.frame(colData(object)[, color_obj, drop=FALSE])
  
  # add the color_obj factors together to create a new grouping factor
  group <- if (length(color_obj) > 1) {
    factor(apply( color_obj.df, 1, paste, collapse=":"))
  } 
  else {
    colData(object)[[color_obj]]
  }
  
  # assembly the data for the plot
  pcaData <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, color_obj.df, name=colnames(object))
  colnames(pcaData)[1] <- paste("PC_1")
  colnames(pcaData)[2] <- paste("PC_2")
  attr(pcaData, "percentVar") <- percentVar[c(PC_1,PC_2)]
  pcaData$condition <- pcaData[[paste0(color_obj)]]
  # transform variance to percent
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  color_title<-paste0(color_obj)
  if (continuous==F){
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
        geom_point(size =3) +
        scale_color_discrete(name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete(name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete(name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }
    }}
  else{
    color_obj<-as.integer(annotation[[paste0(color_obj)]])
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =3) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}
  }
  pca_summary <- plot(pca, main="Variance explained by components", type="l")
  pca <- list()
  pca <- list(pca_plot,pca_summary)
}

PCA_DEseq2 <- function(object=rld, color_obj="condition", 
                       ntop=500, PC_1=2, PC_2=3,
                       shape_opt = "NULL",
                       continuous=F,
                       colour_gradient=c("red","green","blue")){
  if(ntop=="all"){
    pca <- prcomp(t(assay(object))) 
  }
  else{
    # calculate the variance for each gene
    rv <- rowVars(assay(object))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select,]))
  }
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(color_obj %in% names(colData(object)))) {
    stop("the argument 'color_obj' should specify columns of colData(dds)")
  }
  
  color_obj.df <- as.data.frame(colData(object)[, color_obj, drop=FALSE])
  
  # add the color_obj factors together to create a new grouping factor
  group <- if (length(color_obj) > 1) {
    factor(apply( color_obj.df, 1, paste, collapse=":"))
  } 
  else {
    colData(object)[[color_obj]]
  }
  
  # assembly the data for the plot
  pcaData <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, color_obj.df, name=colnames(object))
  colnames(pcaData)[1] <- paste("PC_1")
  colnames(pcaData)[2] <- paste("PC_2")
  attr(pcaData, "percentVar") <- percentVar[c(PC_1,PC_2)]
  pcaData$condition <- pcaData[[paste0(color_obj)]]
  # transform variance to percent
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  color_title<-paste0(color_obj)
  if(continuous ==F){
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
        geom_point(size =3) +
        scale_color_discrete()+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}}
  else{
    if(shape_opt=="NULL"){
      color_obj<-as.integer(annotation[[paste0(color_obj)]])
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =3) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        color_obj<-as.integer(annotation[[paste0(color_obj)]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        color_obj<-as.integer(annotation[[paste0(color_obj)]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}
  }
  pca_summary <- plot(pca, main="Variance explained by components", type="l")
  pca <- list()
  pca <- list(pca_plot,pca_summary)
}

Limma_batch <- function(rld_obj=rld,
                        object=rld_df, 
                        batch_obj=annotation$Sex,
                        batch2_obj=NULL,
                        shape_opt="Sex",
                        ntop=500, PC_1=1, PC_2=2,
                        color_obj="Treatment",
                        continuous=F,
                        colour_gradient=c("red","green","blue")){
  
  removedbatch_rld <- removeBatchEffect(x=object, batch=batch_obj, batch2 = batch2_obj)
  if(ntop=="all"){
    pca <- prcomp(t(removedbatch_rld)) 
  }
  else{
    # select the ntop genes by variance
    rv <- rowVars(assay(rld_obj))
    
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    pca <- prcomp(t(removedbatch_rld[select,]))
  }
  
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  color_title<-paste0(color_obj)
  
  color_obj.df <- as.data.frame(annotation[paste0(color_obj)])
  group <- annotation[paste0(color_obj)]
  pcaData <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, color_obj.df, name=colnames(object))
  colnames(pcaData)[1] <- paste("PC_1")
  colnames(pcaData)[2] <- paste("PC_2")
  attr(pcaData, "percentVar") <- percentVar[c(PC_1,PC_2)]
  pcaData$condition <- pcaData[[paste0(color_obj)]]
  
  # transform variance to percent
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #plot PCA
  if(continuous ==F){
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
        geom_point(size =3) +
        scale_color_discrete()+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}}
  else{
    if(shape_opt=="NULL"){
      color_obj<-as.integer(annotation[[paste0(color_obj)]])
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =3) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        color_obj<-as.integer(annotation[[paste0(color_obj)]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        color_obj<-as.integer(annotation[[paste0(color_obj)]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}
  }
  pca_plot
}
Limma_batch_sva <- function(rld_obj=rld,
                            object=rld_df, 
                            batch_obj=annotation$Sex,
                            shape_opt="Sex",
                            ntop=500, PC_1=1, PC_2=2,
                            color_obj="Treatment",
                            continuous=F,
                            colour_gradient=c("red","blue")){
  
  removedbatch_rld <- removeBatchEffect(x=as.matrix(assay(rld)),covariates = batch_obj,design = model.matrix(~annotation$merged))
  if(ntop=="all"){
    pca <- prcomp(t(removedbatch_rld)) 
  }
  else{
    # select the ntop genes by variance
    rv <- rowVars(assay(rld_obj))
    
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    pca <- prcomp(t(removedbatch_rld[select,]))
  }
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  color_title<-paste0(color_obj)
  
  color_obj.df <- as.data.frame(annotation[paste0(color_obj)])
  group <- annotation[paste0(color_obj)]
  pcaData <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, color_obj.df, name=colnames(object))
  colnames(pcaData)[1] <- paste("PC_1")
  colnames(pcaData)[2] <- paste("PC_2")
  attr(pcaData, "percentVar") <- percentVar[c(PC_1,PC_2)]
  pcaData$condition <- pcaData[[paste0(color_obj)]]
  
  # transform variance to percent
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  #plot PCA
  if(continuous ==F){
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
        geom_point(size =3) +
        scale_color_discrete()+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition, shape=new)) +
          geom_point(size =3) +
          scale_color_discrete()+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}}
  else{
    if(shape_opt=="NULL"){
      color_obj<-as.integer(annotation[[color_obj]])
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =3) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.background = element_rect(fill=NA, color = "black"),
              aspect.ratio = 1,
              plot.background = element_blank(),
              legend.background = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5),
              strip.background = element_rect(fill = NA, color = "black"),
              axis.line = element_line(colour = "black", size = 1),
              axis.text = element_text(colour = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black",size = 1),
              axis.title = element_text(colour = "black", size = 10, face = "bold"),
              legend.key = element_rect(fill = "white", colour = "black"))
    }
    else{
      if (length(levels(annotation[[paste0(shape_opt)]]))<6){
        color_obj<-as.integer(annotation[[color_obj]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
      }
      else{
        color_obj<-as.integer(annotation[[color_obj]])
        pcaData["new"]<- as.character(annotation[[paste0(shape_opt)]])
        legend_title <- paste0(shape_opt)
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj, shape=new)) +
          geom_point(size =3) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.background = element_rect(fill=NA, color = "black"),
                aspect.ratio = 1,
                plot.background = element_blank(),
                legend.background = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0.5),
                strip.background = element_rect(fill = NA, color = "black"),
                axis.line = element_line(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 10, face = "bold"),
                axis.ticks = element_line(color = "black",size = 1),
                axis.title = element_text(colour = "black", size = 10, face = "bold"),
                legend.key = element_rect(fill = "white", colour = "black"))
        
      }}
  }
  pca_plot
}

setClass(Class = "DESeq2_analysis_object",
         slots = c(parameters="list", results="list", DE_genes="list", Number_DE_genes="data.frame"))
Dea_analysis <- function(annotation_file=annotation,
                         count_matrix_file=count_matrix,
                         IHW_option=F,
                         alpha_option=0.05, 
                         lfc_Threshold=0, 
                         control=list("CON_GM","CON_WM"), 
                         condition="merged",
                         design_variable=design){
  DE_objects <- list()
  
  for(j in control){
    cond_list <- list()
    list_conditions<-list()
    list_DE_genes <- list()
    list_DE_sum <- list()
    list_test_all <- list()
    list_DE_genes_names <- list()
    control_list<-list()
    #creat DE_object
    DE_object <- new(Class = "DESeq2_analysis_object")
    # Define parameters
    
    # relevel condition column to have control as place number one in the levels
    anno_DE_obj<-annotation_file
    relevel(anno_DE_obj[[condition]], ref=j)
    # create new DESeqDataSet
    dsfm <- DESeqDataSetFromMatrix(countData = count_matrix_file,
                                   colData = anno_DE_obj,
                                   design = design_variable )
    # Run DESeq function
    dds <- DESeq(dsfm)
    
    DE_object@parameters <- list(dds,IHW_option, alpha_option, lfc_Threshold, j, condition)
    
    #Define the levels that should be compared against the control
    cond_list <- levels(anno_DE_obj[[condition]])
    cond_list<-cond_list[cond_list != j]
  # Define parameters
  # Create results table as data.frame
  
  for (i in cond_list){
    if (IHW_option==T) {
      res_deseq_lfc <- results(dds,contrast = c(condition, i, j),
                               lfcThreshold = lfc_Threshold,
                               alpha = alpha_option,
                               filterFun = ihw,
                               altHypothesis = "greaterAbs")
      res_deseq_lfc <- lfcShrink(dds, contrast = c(condition, i, j),
                                 res=res_deseq_lfc)
      res_deseq_lfc <- as.data.frame(res_deseq_lfc)
      list_conditions[[paste(i)]] <- assign(  paste(i), res_deseq_lfc )
    }
    if (IHW_option==F) {
      res_deseq_lfc <- results(dds,contrast = c(condition, i, j),
                               lfcThreshold = lfc_Threshold,
                               alpha = alpha_option,
                               altHypothesis = "greaterAbs")
      res_deseq_lfc <- lfcShrink(dds, contrast = c(condition, i, j),
                                 res=res_deseq_lfc)
      res_deseq_lfc <- as.data.frame(res_deseq_lfc)
      list_conditions[[paste(i)]] <- assign(  paste(i), res_deseq_lfc )
    } 
  DE_object@results <- list_conditions
    list_DE_genes <- list(rownames(list_conditions[[paste(i)]][!is.na(list_conditions[[paste(i)]]$padj)&
                                                                 list_conditions[[paste(i)]]$padj<alpha_option&
                                                                 list_conditions[[paste(i)]]$log2FoldChange>lfc_Threshold,]),
                          rownames(list_conditions[[paste(i)]][!is.na(list_conditions[[paste(i)]]$padj)&
                                                                 list_conditions[[paste(i)]]$padj<alpha_option&
                                                                 list_conditions[[paste(i)]]$log2FoldChange<(-lfc_Threshold),]))
    names(list_DE_genes) = c(paste("up-regulated genes"), 
                             paste("down-regulated genes"))
    list_DE_genes_names[[paste(i)]] <- assign(  paste(i), list_DE_genes )
    list_DE_sum <- list(nrow(list_conditions[[paste(i)]][!is.na(list_conditions[[paste(i)]]$padj)&
                                                           list_conditions[[paste(i)]]$padj<alpha_option&
                                                           list_conditions[[paste(i)]]$log2FoldChange>lfc_Threshold,]) ,
                        nrow(list_conditions[[paste(i)]][!is.na(list_conditions[[paste(i)]]$padj)&
                                                           list_conditions[[paste(i)]]$padj<alpha_option&
                                                           list_conditions[[paste(i)]]$log2FoldChange<(-lfc_Threshold),]))
    df_DE_genes <- data.frame(matrix(unlist(list_DE_sum), nrow=2, byrow=T),stringsAsFactors=FALSE, 
                              row.names = c(paste("up-regulated padj. <", alpha_option, sep = "_"), 
                                            paste("down-regulated padj. <", alpha_option, sep = "_")))
    list_test_all[[paste(i)]] <- assign(  paste(i), df_DE_genes )
  }
  df_DE_genes <- data.frame(matrix(unlist(list_DE_sum), nrow=2, byrow=T),stringsAsFactors=FALSE, 
                            row.names = c(paste("up-regulated padj. <", alpha_option, sep = " "), 
                                          paste("down-regulated padj. <", alpha_option, sep = " ")))
  DE_object@DE_genes <- list_DE_genes_names
  df <- do.call("cbind", list_test_all)
  colnames(df) <- levels(anno_DE_obj[[condition]])[-1]
  DE_object@Number_DE_genes <- df
  DE_objects[[paste(j)]] <- assign(  paste(j), DE_object)}
  return(DE_objects)
}

multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
df_function <- function(results_obj=DE_object_01@results){
  list_df<- list()
  for (j in names(results_obj)){
    res_df <- results_obj[[j]]
    res_df <- assign(  paste(j, sep = ""),res_df)
    list_df[[paste(j)]] <- res_df
  }
  list_df
}
MA_function <- function(condition=NULL,df_names=df_names_CON_GM, DE_obj=DE_object, y_lim=c(-10,10)){
  list_store <- list()
  if(!is.null(condition)){ multi_plot =F} else {multi_plot=T}
  plot_MA <- function(condition, lfc_shrinkage=T){
    label <- as.character(condition)
    
    main_name = paste(label," vs ", DE_obj@parameters[[5]])
    ylab_name = "MAP log2 fold change"
    
    # get shrunken lfc object as data frame
    res_deseq_lfc_df <- as.data.frame(list_df[[condition]])
    
    res_deseq_lfc_df$significant <- res_deseq_lfc_df$padj < DE_obj@parameters[[3]]
    # res_deseq_lfc_df$significant_lFC <- (res_deseq_lfc_df$padj < 0.1 & res_deseq_lfc_df$log2FoldChange>1) 
    
    # plot shrunken lfc vs mean 
    ggplot(res_deseq_lfc_df, aes(x = baseMean, y = log2FoldChange)) + 
      geom_point(aes(color = padj<DE_obj@parameters[[3]]), size = 1.5, alpha = 0.5) + 
      geom_hline(yintercept = 0, color = "slategray", size=1) + 
      ylim(y_lim) + ylab(ylab_name) + 
      theme_bw() +
      labs(title = main_name, color="Significance") +
      scale_x_log10("Mean expression",breaks=c(1,10, 100,1000,10000,100000)) +
      scale_colour_manual(labels = c(paste0("padj > ",DE_obj@parameters[[3]]), paste0("padj < ",DE_obj@parameters[[3]])),
                          values = c("#2B2D42","#6E0B21", "#7A7D7F"))
  }
  if (multi_plot==F){
    plot_MA(condition=condition)
    
  } 
  else {
    for (i in df_names){
      list_store[[paste("MA_",i,"plot_lfc_shrinkage",sep = "_")]] <- plot_MA(condition = i)
    }
    multiplot(plotlist = list_store,cols = 2)
  }
}

plot_pval <- function(df_names=df_names_CON_GM,condition=NULL,ylim_obj=c(0,10000), DE_obj=DE_object){
  list_pval <- list()
  if(!is.null(condition)){ multi_plot =F} else {multi_plot=T}
  histo <- function(condition){
    label <- as.character(condition)
    condition <- as.data.frame(list_df[[condition]])
    ggplot(condition, aes(x = pvalue)) + 
      geom_histogram(binwidth = 0.025, boundary = 0)+
      theme_bw() + labs(title = paste(label,"vs", DE_obj@parameters[[5]],sep=" "))+
      ylab("Frequency") + ylim(ylim_obj)
  }
  if (multi_plot==F){
    histo_gram <- histo(condition=condition)
    return(histo_gram)
  }
  if (multi_plot==T){
    for (i in df_names){
      list_pval[[paste("pval",i,sep = "_")]] <- histo(condition = i)
    }
    multiplot(plotlist = list_pval,cols = 2)
  }
}
DE_genes_plot <- function(DE_genes_df = DE_object@Number_DE_genes, DE_obj=DE_object){
  # melt
  melt_DE_genes_df <- melt(as.matrix(DE_genes_df))
  # add column up down regulated and conditions
  melt_DE_genes_df$regulation <- ifelse(grepl("up",melt_DE_genes_df$Var1, ignore.case = T), "upregulated", "downregulated")
  colnames(melt_DE_genes_df)[2] <- "condition"
  melt_DE_genes_df$condition <- factor(melt_DE_genes_df$condition)
  Geom_bar_DE <- ggplot(melt_DE_genes_df , 
                        aes(Var1, value)) + 
    geom_bar(stat = "identity", aes(fill = regulation)) +
    theme_bw() + facet_grid(.~condition)+
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_blank())+
    scale_fill_manual(values = c("#B1BBCF", "#FC9D9A"))+
    guides(fill = guide_legend(title = "DE genes"))
  
  ##Add gene numbers to the bars
  Geom_bar_DE + geom_text(aes(label = value), size = 3,
                          hjust = 0.5, vjust = 3, position="stack")+
    labs(title="Number of DE genes", 
         subtitle= c(paste0("padj < ", DE_obj@parameters[[3]], 
                            ", Log2FC threshold: ", DE_obj@parameters[[4]])))+
    ylab("Number of genes")
}


Volcano_plot <- function(input_file=DE_object, condition="CpG", x_limit=c(-10,10), y_limit=c(NA,NA)){
  df <- as.data.frame(input_file@results[[condition]])
  df$declogp <- -log10(df$padj)
  
  ggplot(df, aes(x=df$log2FoldChange, y=df$declogp))+
    geom_point(colour=ifelse(df$log2FoldChange< -input_file@parameters[[4]]&df$declogp>10^-(input_file@parameters[[3]]) | 
                               df$log2FoldChange>input_file@parameters[[4]]&df$declogp>10^-(input_file@parameters[[3]]),"red","grey"))+
    scale_y_discrete(expression('-log'['10']*' adjusted p-value'), limits=y_limit)+
    scale_x_discrete(expression('Fold change log'['2']), limits=x_limit)+
    coord_cartesian(xlim = x_limit)+
    labs(title=paste(condition, "VS. ",input_file@parameters[[5]]))+
    theme(panel.background = element_rect(fill=NA, color = "black"),
          plot.background = element_blank(),
          legend.background = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          strip.background = element_rect(fill = NA, color = "black"),
          axis.line = element_line(colour = "black", size = 1),
          axis.text = element_text(colour = "black", size = 10, face = "bold"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.title = element_text(colour = "black", size = 10, face = "bold"))
}
Volcano_plot_display_gene <- function(input_file=DE_object, condition="CpG", 
                                      gene="Tnf",size_def=3, x_limit=c(-10,10),
                                      y_limit=c(NA,NA)){
  df <- as.data.frame(input_file@results[[condition]])
  df$declogp <- -log(df$padj)
  
  ggplot(df, aes(x=df$log2FoldChange, y=df$declogp))+
    geom_point(colour=ifelse(df$log2FoldChange< -input_file@parameters[[4]]&df$declogp>10^-(input_file@parameters[[3]]) | 
                               df$log2FoldChange>input_file@parameters[[4]]&df$declogp>10^-(input_file@parameters[[3]]),"red","grey"))+
    scale_y_discrete(expression('-log'['10']*' p-value'), limits = y_limit)+
    scale_x_discrete(expression('Fold change log'['2']),limits = x_limit)+
    geom_text(data=df[gene,], mapping=aes(x=df[gene,][,"log2FoldChange"], 
                                          y=df[gene,][,"declogp"], 
                                          label=gene),
              size=size_def, hjust=0.5, vjust=1)+
    geom_point(data=df[gene,], mapping=aes(x=df[gene,][,"log2FoldChange"], 
                                           y=df[gene,][,"declogp"]), 
               colour="blue")+
    
    labs(title=paste(condition, "VS. ",input_file@parameters[[6]]))+
    theme(panel.background = element_rect(fill=NA, color = "black"),
          plot.background = element_blank(),
          legend.background = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          strip.background = element_rect(fill = NA, color = "black"),
          axis.line = element_line(colour = "black", size = 1),
          axis.text = element_text(colour = "black", size = 10, face = "bold"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.title = element_text(colour = "black", size = 10, face = "bold"))
}
cluster_Top_genes_output <- function(object=rld, dds_obj=dds, anno_color,
                                     first_annotation="condition",
                                     second_annotation="Sex",
                                     heatmap_title, gene_list=all_DE_genes,
                                     ntop=1000){
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  df_ <- as.data.frame(colData(dds)[,c(first_annotation,second_annotation)])
  colnames(df_) <- c(first_annotation,second_annotation)
  row.names(df_) <- colnames(object)
  
  map <- pheatmap(assay(object)[select,],
                  cluster_rows=TRUE, cluster_cols=T,
                  clustering_distance_rows="correlation",
                  clustering_distance_cols="correlation",
                  show_rownames=F, show_colnames = TRUE,
                  scale = "row",
                  main = heatmap_title,
                  annotation_col=df_,
                  annotation_colors = anno_color) 
  name_output <<- assay(object)[map$tree_row$order,]
}
Cluster_genelist_output <-function(object=rld, dds_obj=dds, anno_color,
                                   first_annotation="condition",
                                   second_annotation="Sex",
                                   heatmap_title, gene_list=all_DE_genes,
                                   display_row=F){
  
  df_ <- as.data.frame(colData(dds)[,c(first_annotation,second_annotation)])
  colnames(df_) <- c(first_annotation,second_annotation)
  row.names(df_) <- colnames(object)
  
  map <- pheatmap(assay(object)[gene_list,],
                  cluster_rows=TRUE, cluster_cols=T,
                  clustering_distance_rows="correlation",
                  clustering_distance_cols="correlation",
                  show_rownames=display_row, show_colnames = TRUE,
                  scale = "row",
                  main = heatmap_title,
                  annotation_col=df_,
                  annotation_colors =anno_color ) 
  Clustering_output <<- assay(object)[map$tree_row$order,]
}




plot_enrichGO <- function(downreg_genes_list, 
                          upreg_genes_list, 
                          downreg_title, upreg_title, 
                          nr_of_BPs_to_plot=20,
                          minGSSize_no = 10,
                          maxGSSize_no = 500,
                          keytype_opt="ENSEMBL",
                          pvalueCutoff_opt=0.05,
                          ontology="BP",
                          organism="org.Mm.eg.db"){
  # make first plot for upreg genes
  # GO Enrichment Analysis
  enrichGO_upreg <- enrichGO(upreg_genes_list, 
                             OrgDb = organism,
                             ont = ontology, 
                             keyType = keytype_opt, 
                             minGSSize = minGSSize_no,
                             maxGSSize = maxGSSize_no,
                             pvalueCutoff = pvalueCutoff_opt)
  # make dataframe
  upreg_df <- as.data.frame(enrichGO_upreg)
  #order by pvalue
  upreg_ordered <- upreg_df[order(upreg_df$pvalue),]
  # select nr of biological processes(BPs) to be ploted
  upreg_nr_selected_BPs <- upreg_ordered[1:nr_of_BPs_to_plot,]
  # remove all columns exept for Description, pvalue and Count
  # upreg_nr_selected_BPs[,c(1,3:4,6:8)] <- list(NULL)
  
  # order Description column according to pvalues for proper plotting
  upreg_nr_selected_BPs$Description <- factor(upreg_nr_selected_BPs$Description, 
                                              levels = c(upreg_nr_selected_BPs[order(upreg_nr_selected_BPs$pvalue, 
                                                                                     decreasing = T), "Description"]))
  levels(upreg_nr_selected_BPs$Description) <- gsub("(.{25,}?)\\s", "\\1\n", levels(upreg_nr_selected_BPs$Description))  
  
  # ggplot upreg genes
  gigi_up <- ggplot(upreg_nr_selected_BPs, aes(Description))+
    geom_bar(aes(y=Count, fill = pvalue),
             stat = "identity")  + 
    scale_x_discrete(position = "top") + 
    coord_flip() + 
    scale_fill_gradient(low="#6E0B21",
                        high="#B1BBCF",
                        guide = guide_colorbar(reverse=T)) +
    theme(plot.title = element_text(hjust = 1), 
          aspect.ratio = 2.5,
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          axis.text.x = element_text(size = 10),
          panel.background = element_rect(fill = "white", colour = "grey85")) +
    labs(title=upreg_title) + 
    ylab("Count") + 
    xlab(label = "")
  
  # make second plot for downreg genes
  # GO Enrichment Analysis
  enrichGO_downreg <- enrichGO(downreg_genes_list, 
                               OrgDb = organism,
                               ont = ontology, 
                               keyType = keytype_opt,
                               minGSSize = minGSSize_no,
                               maxGSSize = maxGSSize_no,
                               pvalueCutoff = pvalueCutoff_opt)
  # make dataframe
  downreg_df <- as.data.frame(enrichGO_downreg)
  #order by pvalue
  downreg_ordered <- downreg_df[order(downreg_df$pvalue),]
  # select nr of biological processes(BPs) to be ploted
  downreg_nr_selected_BPs <- downreg_ordered[1:nr_of_BPs_to_plot,]
  # remove all columns exept for Description, pvalue and Count
  # downreg_nr_selected_BPs[,c(1,3:4,6:8)] <- list(NULL)
  
  # order Description column according to pvalues for proper plotting
  downreg_nr_selected_BPs$Description <- factor(downreg_nr_selected_BPs$Description, 
                                                levels = c(downreg_nr_selected_BPs[order(downreg_nr_selected_BPs$pvalue, 
                                                                                         decreasing = T), "Description"]))
  
  levels(downreg_nr_selected_BPs$Description) <- gsub("(.{25,}?)\\s", "\\1\n", levels(downreg_nr_selected_BPs$Description))  
  
  # ggplot downreg genes
  gigi_down <- ggplot(downreg_nr_selected_BPs, aes(Description)) + 
    geom_bar(aes(y=-Count, fill=pvalue),
             stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low="#53777A", high="#C8C8A9",
                        guide = guide_colorbar(reverse=T)) +
    theme(legend.justification=c(0,0), 
          aspect.ratio = 2.5,
          legend.position=c(0,0),
          axis.text.x = element_text(size = 10),
          panel.background = element_rect(fill = "white", colour = "grey85")) +
    labs(title=downreg_title) +
    ylab("Count") +
    xlab(label = "")
  
  
  # use the multiplot function to plot gigi_upreg and gigi_downreg together
  multiplot(gigi_down, gigi_up, cols = 2)
}
Plot_pathway_cluster <- function(gene_list=all_upreg_genes,
                                 organism="org.Mm.eg.db",
                                 ontology="CC",
                                 keytype_opt="SYMBOL",
                                 pvalue_cutoff=0.05,
                                 pathway="anchoring junction",
                                 first_annotation="condition",
                                 second_annotation="Sex",
                                 anno_color=list(condition=anno,Sex=anno_batch),
                                 object=rld,dds_obj=dds,
                                 display_genes=T){
  
  GO_CC_all_DE_up <- as.data.frame(enrichGO(gene_list, 
                                            OrgDb = organism,
                                            keyType=keytype_opt,
                                            ont=ontology,
                                            pvalueCutoff=pvalue_cutoff))
  
  # Make a genes_list out of them# 
  
  Genes_focal_adhesion <- filter(GO_CC_all_DE_up, Description==pathway) %>% dplyr::select(geneID)
  
  Genes_focal_adhesion <- Genes_focal_adhesion$geneID
  
  Genes_focal_adhesion <-  gsub("/", " ", Genes_focal_adhesion)
  Genes_focal_adhesion
  
  Genes_focal_adhesion <- unlist(strsplit(Genes_focal_adhesion, " "))
  class(Genes_focal_adhesion)
  
  df_ <- as.data.frame(colData(dds_obj)[,c(first_annotation,second_annotation)])
  colnames(df_) <- c(first_annotation,second_annotation)
  row.names(df_) <- colnames(object)
  
  pheatmap(assay(object)[Genes_focal_adhesion,],
           cluster_rows=TRUE, cluster_cols=T,
           show_rownames=display_genes, show_colnames = TRUE,
           scale = "row",
           main = pathway,
           annotation_col=df_,
           annotation_colors = anno_color) 
}
plot_KEGG <- function(upreg_genes_list, downreg_genes_list,
                      upreg_title="LPS_up",downreg_title="LPS_down",
                      organism_KEGG = "mmu", keyType = "SYMBOL",
                      minGSSize_no = 10,
                      maxGSSize_no = 500,
                      pvalueCutoff_opt = 0.05, 
                      nr_of_BPs_to_plot = 20,
                      organism_opt=org.Mm.eg.db){
  
  # Transform keytypes of DE genes into ncbi_gene_IDs
  
  upreg_genes_list <- mapIds(x=organism_opt,
                             keys=upreg_genes_list,
                             column="ENTREZID",
                             keytype=keyType,
                             multiVals="first")
  upreg_genes_list <- as.character(upreg_genes_list)
  
  downreg_genes_list <- mapIds(x=organism_opt,
                               keys=downreg_genes_list,
                               column="ENTREZID",
                               keytype=keyType,
                               multiVals="first")
  downreg_genes_list <- as.character(downreg_genes_list)
  
  # make first plot for upreg genes
  # GO Enrichment Analysis
  enrichGO_upreg <- enrichKEGG(gene=upreg_genes_list, 
                               organism = organism_KEGG, keyType = "ncbi-geneid", 
                               minGSSize = minGSSize_no,
                               maxGSSize = maxGSSize_no,
                               pvalueCutoff = pvalueCutoff_opt)
  # make dataframe
  upreg_df <- as.data.frame(enrichGO_upreg)
  #order by pvalue
  upreg_ordered <- upreg_df[order(upreg_df$pvalue),]
  # select nr of biological processes(BPs) to be ploted
  upreg_nr_selected_BPs <- upreg_ordered[1:nr_of_BPs_to_plot,]
  # remove all columns exept for Description, pvalue and Count
  # upreg_nr_selected_BPs[,c(1,3:4,6:8)] <- list(NULL)
  
  # order Description column according to pvalues for proper plotting
  upreg_nr_selected_BPs$Description <- factor(upreg_nr_selected_BPs$Description, 
                                              levels = c(upreg_nr_selected_BPs[order(upreg_nr_selected_BPs$pvalue, 
                                                                                     decreasing = T), "Description"]))
  levels(upreg_nr_selected_BPs$Description) <- gsub("(.{25,}?)\\s", "\\1\n", levels(upreg_nr_selected_BPs$Description))  
  
  # ggplot upreg genes
  gigi_up <- ggplot(upreg_nr_selected_BPs, aes(Description))+
    geom_bar(aes(y=Count, fill = pvalue),
             stat = "identity")  + 
    scale_x_discrete(position = "top") + 
    coord_flip() + 
    scale_fill_gradient(low="#6E0B21",
                        high="#B1BBCF",
                        guide = guide_colorbar(reverse=T)) +
    theme(plot.title = element_text(hjust = 1), 
          aspect.ratio = 2.5,
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          axis.text.x = element_text(size = 10),
          panel.background = element_rect(fill = "white", colour = "grey85")) +
    labs(title=upreg_title) + 
    ylab("Count") + 
    xlab(label = "")
  
  # make second plot for downreg genes
  # GO Enrichment Analysis
  enrichGO_downreg <- enrichKEGG(gene=downreg_genes_list, 
                                 organism = organism_KEGG, keyType = "ncbi-geneid", 
                                 minGSSize = minGSSize_no,
                                 maxGSSize = maxGSSize_no,
                                 pvalueCutoff = pvalueCutoff_opt)
  # make dataframe
  downreg_df <- as.data.frame(enrichGO_downreg)
  #order by pvalue
  downreg_ordered <- downreg_df[order(downreg_df$pvalue),]
  # select nr of biological processes(BPs) to be ploted
  downreg_nr_selected_BPs <- downreg_ordered[1:nr_of_BPs_to_plot,]
  # remove all columns exept for Description, pvalue and Count
  # downreg_nr_selected_BPs[,c(1,3:4,6:8)] <- list(NULL)
  
  # order Description column according to pvalues for proper plotting
  downreg_nr_selected_BPs$Description <- factor(downreg_nr_selected_BPs$Description, 
                                                levels = c(downreg_nr_selected_BPs[order(downreg_nr_selected_BPs$pvalue, 
                                                                                         decreasing = T), "Description"]))
  
  levels(downreg_nr_selected_BPs$Description) <- gsub("(.{25,}?)\\s", "\\1\n", levels(downreg_nr_selected_BPs$Description))  
  
  # ggplot downreg genes
  gigi_down <- ggplot(downreg_nr_selected_BPs, aes(Description)) + 
    geom_bar(aes(y=-Count, fill=pvalue),
             stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low="#53777A", high="#C8C8A9",
                        guide = guide_colorbar(reverse=T)) +
    theme(legend.justification=c(0,0), 
          aspect.ratio = 2.5,
          legend.position=c(0,0),
          axis.text.x = element_text(size = 10),
          panel.background = element_rect(fill = "white", colour = "grey85")) +
    labs(title=downreg_title) +
    ylab("Count") +
    xlab(label = "")
  
  
  # use the multiplot function to plot gigi_upreg and gigi_downreg together
  multiplot(gigi_down, gigi_up, cols = 2)
}
plot_single_gene <- function(dds_object=dds, gene_symbol="Tnf", 
                             condition="condition", pc_cond=F, anno_colour=anno,
                             order=c("suspension_wt",
                                     "suspension_Cyth2_KO",
                                     "fibronectin_Cyth2_KO",
                                     "fibronectin_wt",
                                     "gelatin_wt",
                                     "gelatin_Cyth2_KO"),shape_opt="Treatment") {
  geneCounts_lfc <- plotCounts(dds_object, gene = gene_symbol, 
                               color_obj = condition, pc=pc_cond,
                               returnData = TRUE)
  geneCounts_lfc$condition <- annotation[[paste0(condition)]]
  geneCounts_lfc$condition <- factor(geneCounts_lfc$condition, levels =order )
  geneCounts_lfc$sign <- annotation[[paste0(shape_opt)]]
  legend_shape<-paste0(shape_opt)
  if (shape_opt=="NULL"){
    ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition)) +
      scale_y_discrete(limits = c(0,max(geneCounts_lfc$count))) +  
      scale_color_manual(values=anno_colour)+
      geom_beeswarm(cex = 3, na.rm=T)+
      ylab("Normalized counts")+
      labs(title=paste0(gene_symbol),colour=condition)+
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))}
  else{
    ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition, shape=sign)) +
      scale_y_discrete(limits = c(0,max(geneCounts_lfc$count))) +  
      scale_color_manual(values=anno_colour)+
      scale_shape(name=legend_shape)+
      geom_beeswarm(cex = 3, na.rm=T)+
      ylab("Normalized counts")+
      labs(title=paste0(gene_symbol),colour=condition)+
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))
  }
}

Compare_FC <- function(input_file=DE_object@results, condition1="CpG", condition2="LPS_IFNg"){
  df <- as.data.frame(input_file)
  df_cond1 <- df[[condition1]]
  df_cond2 <- df[[condition2]]
  ggplot(df, aes(x=df_cond1, y=df_cond2))+
    geom_point(colour="blue")+
    scale_x_discrete(expression(paste('log'['2']*'fold change', condition1, 'VS control')))+
    scale_y_discrete(expression(paste('log'['2']*'fold change', condition2, 'VS control')))+
    coord_cartesian()+
    labs(title=paste("FC", condition1, "VS", condition2, sep=" "))+
    theme(panel.background = element_rect(fill=NA, color = "black"),
          plot.background = element_blank(),
          legend.background = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          strip.background = element_rect(fill = NA, color = "black"),
          axis.line = element_line(colour = "black", size = 1),
          axis.text = element_text(colour = "black", size = 10, face = "bold"),
          axis.ticks = element_line(color = "black",size = 1),
          axis.title = element_text(colour = "black", size = 10, face = "bold"))
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
DE_genes_Top <- function(ntop=500, condi="NULL", input_file=DE_object$CON_GM, condition=annotation$condition){
  if (condi=="NULL"){
    cond_list <- list()
    list_top_DE <- list()
    cond_list <- names(input_file@results)
    for (i in cond_list){
      select <- order(input_file@results[[i]]$padj, decreasing = F, na.last=T )[seq_len(min(ntop, length(input_file@results[[i]]$padj)))]
      Top500_DE_genes <- rownames(input_file@results[[i]][select,])
      list_top_DE[[paste(i)]] <- assign(paste(i),Top500_DE_genes )}
    list_top_DE<-do.call(c, list_top_DE)
    list_top_DE <- unique(list_top_DE)
  }
  else {
    select <- order(input_file@results[[condi]]$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(DE_object@results[[condi]]$padj)))]
    
    Top500_DE_genes <- input_file@results[[condi]][select,]
  }
}
DE_genes_Top_up <- function(ntop=500, condi="NULL",input_file=DE_object$CON_GM, condition="merged"){
  if (condi=="NULL"){
    cond_list <- list()
    list_top_DE <- list()
    cond_list <- names(input_file@results)
    for (i in cond_list){
      df<-input_file@results[[i]][input_file@results[[i]]$log2FoldChange>0,]
      select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
      
      Top500_DE_genes <- rownames(df[select,])
      list_top_DE[[paste(i)]] <- assign(paste(i),Top500_DE_genes )}
    list_top_DE<-do.call(c, list_top_DE)
    list_top_DE <- unique(list_top_DE)
  }
  else {
    df<-input_file@results[[condi]][input_file@results[[condi]]$log2FoldChange>0,]
    select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
    
    Top500_DE_genes <- rownames(df[select,])
  }
}

DE_genes_Top_down <- function(ntop=10, input_file=DE_object$CON_GM,condition="CON_GM", condi="NULL"){
  if (condi=="NULL"){
    cond_list <- list()
    list_top_DE <- list()
    cond_list <- names(input_file@results)
    for (i in cond_list){
      df<-input_file@results[[i]][input_file@results[[i]]$log2FoldChange<0,]
      select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
      
      Top500_DE_genes <- rownames(df[select,])
      list_top_DE[[paste(i)]] <- assign(paste(i),Top500_DE_genes )}
    list_top_DE<-do.call(c, list_top_DE)
    list_top_DE <- unique(list_top_DE)
  }
  else {
    df<-input_file@results[[condi]][input_file@results[[condi]]$log2FoldChange<0,]
    select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
    
    Top500_DE_genes <- rownames(df[select,])
  }
}
