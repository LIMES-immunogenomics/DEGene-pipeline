list.of.packages <- c("tidyverse", "reshape2","dplyr","stringr","pheatmap","colourlovers",'VennDiagram',
                      "ggbeeswarm","scales","magrittr", "tidyr","bit", "data.table","RSQLite","Rcpp", "rJava",
                      "devtools", "matrixStats", "xlsx", "gridBase", "survival", "rlang", "MASS", "ggforce", "tibble",
                      "stringi","gtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
library(r2excel)
list.of.bioc.packages<- c("rhdf5","DESeq2","IHW","clusterProfiler","GSEABase","AnnotationDbi",
                          "sva","tximport","limma","geneplotter","genefilter","org.Mm.eg.db","org.Hs.eg.db",
                          "biomaRt", "ReactomePA", "grid","pcaGoPromoter","pcaGoPromoter.Hs.hg19",
                          "pcaGoPromoter.Mm.mm9", "limma", "DOSE")
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
    count_matrix<-aggregate(.~rownames(count_matrix), data=count_matrix, FUN=sum)
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
  if(ncol(count_mat)<30){
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
sample_cor<- function(rld, limit_set = c(0.9,1),midpoint_set = 0.95,condition="Treatment",annotation=annotation){
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
                       anno_colour=anno_merged,
                       ntop=500, PC_1=1, PC_2=2,
                       count_mat=count_matrix,
                       shape_opt="NULL",
                       continuous=F,
                       colour_gradient=c("red","green","blue"),
                       point_size=3){
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
    if(is.null(anno_colour)){
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_discrete(name=color_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_discrete(name=color_title)+
            scale_shape(name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_discrete(name=color_title)+
            scale_shape_manual(values=df$x,name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_manual(values=anno_colour, name=color_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
            scale_shape(name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
            scale_shape_manual(values=df$x,name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
      }}}
  else{
    color_obj<-as.integer(annotation[[paste0(color_obj)]])
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
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
          geom_point(size =point_size) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
          geom_point(size =point_size) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
                       ntop=500, PC_1=2, PC_2=2,
                       anno_colour=anno_merged,
                       shape_opt = "NULL",
                       continuous=F,
                       colour_gradient=c("red","green","blue"),
                       point_size=3){
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
    if(is.null(anno_colour)){
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_discrete(name=color_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_discrete(name=color_title)+
            scale_shape(name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_discrete(name=color_title)+
            scale_shape_manual(values=df$x,name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_manual(values=anno_colour, name=color_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
            scale_shape(name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
            scale_shape_manual(values=df$x,name=legend_title)+
            xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
            ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
            coord_fixed()+
            labs(title="PCA")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
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
      }}}
  else{
    color_obj<-as.integer(annotation[[paste0(color_obj)]])
    if(shape_opt=="NULL"){
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = colour_gradient,name=color_title)+
        xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
        ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
        coord_fixed()+
        labs(title="PCA")+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
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
          geom_point(size =point_size) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape(name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
          geom_point(size =point_size) +
          scale_color_gradientn(colours = colour_gradient,name=color_title)+
          scale_shape_manual(values=df$x,name=legend_title)+
          xlab(paste0("PC ",PC_1, ": ", percentVar[1], "% variance")) +
          ylab(paste0("PC ",PC_2,": ", percentVar[2], "% variance")) +
          coord_fixed()+
          labs(title="PCA")+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
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
                        column_interest=annotation$merged,
                        batch_obj=annotation$Sex,
                        batch2_obj=NULL,
                        shape_opt="Sex",
                        ntop=500, PC_1=1, PC_2=2,
                        color_obj="Treatment",
                        anno_colour=anno_merged,
                        continuous=F,
                        colour_gradient=c("red","green","blue"),
                        point_size=3){
  
  removedbatch_rld <- removeBatchEffect(x=object, batch=batch_obj, batch2 = batch2_obj, model=model.matrix(~paste0(column_interest)))
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
    if(is.null(anno_colour)){
      
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
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
            geom_point(size =point_size) +
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
            geom_point(size =point_size) +
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
          
        }}
    }
    else{
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_manual(values=anno_colour, name=color_title)+
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
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
  }
  if(continuous==T){
    if(shape_opt=="NULL"){
      color_obj<-as.integer(annotation[[paste0(color_obj)]])
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =point_size) +
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
          geom_point(size =point_size) +
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
          geom_point(size =point_size) +
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
                            condition=annotation$merged,
                            batch_obj=annotation$Sex,
                            shape_opt="Sex",
                            ntop=500, PC_1=1, PC_2=2,
                            color_obj="Treatment",
                            anno_colour=anno_merged,
                            continuous=F,
                            colour_gradient=c("red","blue"),
                            point_size=3){
  
  removedbatch_rld <- removeBatchEffect(x=as.matrix(assay(rld_obj)),covariates = batch_obj,design = model.matrix(~condition))
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
    if(is.null(anno_colour)){
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
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
            geom_point(size =point_size) +
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
            geom_point(size =point_size) +
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
          
        }# end of length
      }# end of shape
    }
    else{
      if(shape_opt=="NULL"){
        pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=condition)) +
          geom_point(size =point_size) +
          scale_color_manual(values=anno_colour, name=color_title)+
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
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
            geom_point(size =point_size) +
            scale_color_manual(values=anno_colour, name=color_title)+
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
          
        }# end of length
      }# end of shape
    }
  }#continuous ends here
  else{
    if(shape_opt=="NULL"){
      color_obj<-as.integer(annotation[[color_obj]])
      pca_plot <- ggplot(pcaData, aes(x = PC_1, y = PC_2,colour=color_obj)) +
        geom_point(size =point_size) +
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
          geom_point(size =point_size) +
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
          geom_point(size =point_size) +
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
                         design_variable=design,
                         fdr_correction="BH"){
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
    
    dds <- dds[which(mcols(dds)$betaConv),]
    
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
        res_deseq_lfc$FC <- logratio2foldchange(res_deseq_lfc$log2FoldChange,base=2)
        list_conditions[[paste(i)]] <- assign(  paste(i), res_deseq_lfc )
      }
      if (IHW_option==F) {
        res_deseq_lfc <- results(dds,contrast = c(condition, i, j),
                                 lfcThreshold = lfc_Threshold,
                                 alpha = alpha_option,
                                 independentFiltering = T,
                                 altHypothesis = "greaterAbs",
                                 pAdjustMethod=fdr_correction)
        res_deseq_lfc <- lfcShrink(dds, contrast = c(condition, i, j),
                                   res=res_deseq_lfc)
        res_deseq_lfc <- as.data.frame(res_deseq_lfc)
        res_deseq_lfc$FC <- logratio2foldchange(res_deseq_lfc$log2FoldChange,base=2)
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
      df_DE_genes <- data.frame(matrix(unlist(list_DE_sum), nrow=2, byrow=T),stringsAsFactors=FALSE, 
                                row.names = c(paste("up-regulated padj. <", alpha_option, sep = " "), 
                                              paste("down-regulated padj. <", alpha_option, sep = " ")))
      DE_object@DE_genes <- list_DE_genes_names}
    df <- do.call("cbind", list_test_all)
    colnames(df) <- cond_list
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
MA_function <- function(condition=NULL,
                        df_names=df_names_CON_GM, 
                        DE_obj=DE_object, 
                        y_lim=c(-10,10),
                        list_df=data.statistics){
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

plot_pval <- function(df_names=df_names_CON_GM,
                      condition=NULL,
                      ylim_obj=c(0,10000), 
                      DE_obj=DE_object,
                      list_df=list_df){
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
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values = c("#B1BBCF", "#FC9D9A"))+
    guides(fill = guide_legend(title = "DE genes"))
  
  ##Add gene numbers to the bars
  Geom_bar_DE + geom_text(aes(label = value), size = 3,
                          hjust = 0.5, vjust = 3, position="stack")+
    labs(title=c(paste0("Number of DE genes compared against ", DE_obj@parameters[[5]])), 
         subtitle= c(paste0("padj < ", DE_obj@parameters[[3]], 
                            ", Log2FC threshold: ", DE_obj@parameters[[4]])))+
    ylab("Number of genes")
}

Volcano_plot <- function(input_file=DE_object, condition="CpG", x_limit=c(-10,10), y_limit=c(NA,NA)){
  df <- as.data.frame(input_file@results[[condition]])
  
  ggplot(df, aes(x=df$log2FoldChange, y=log10(df$padj)))+
    geom_point(colour=ifelse(df$log2FoldChange< -input_file@parameters[[4]]&df$padj<input_file@parameters[[3]] | 
                               df$log2FoldChange>input_file@parameters[[4]]&df$padj<input_file@parameters[[3]],"red","grey"))+
   
    scale_x_discrete(expression('Fold change log'['2']), limits=x_limit)+
    scale_y_reverse(expression('-log'['10']*' adjusted p-value'), limits=y_limit)+
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
  
  ggplot(df, aes(x=df$log2FoldChange, y=log10(df$padj)))+
    geom_point(colour=ifelse(df$log2FoldChange< -input_file@parameters[[4]]&df$padj<input_file@parameters[[3]] | 
                               df$log2FoldChange>input_file@parameters[[4]]&df$padj<input_file@parameters[[3]],"red","grey"))+
    scale_y_reverse(expression('-log'['10']*' adjusted p-value'), limits=y_limit)+
    scale_x_discrete(expression('Fold change log'['2']),limits = x_limit)+
    geom_text(data=df[gene,], mapping=aes(x=df[gene,][,"log2FoldChange"], 
                                          y=log10(df[gene,][,"padj"]), 
                                          label=gene),
              size=size_def, hjust=0.5, vjust=1)+
    geom_point(data=df[gene,], mapping=aes(x=df[gene,][,"log2FoldChange"], 
                                           y=log10(df[gene,][,"padj"])), 
               colour="blue")+
    
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
removeBatchEffect_function <- function(x=rld_df,
                                       batch = annotation$`Customer ID`,
                                       model = model.matrix(~annotation$merged)){
  if(is.numeric(batch[,1])==T){
    removeBatchEffect(x=as.matrix(x),
                      covariates = batch,
                      design = model)}
  else{
    removeBatchEffect(x=x,
                      batch = batch,
                      design = model)
  }
}
cluster_Top_genes_output <- function(object=rld_df, dds_obj=dds, anno_color,
                                     first_annotation="condition",
                                     second_annotation="Sex",
                                     heatmap_title, gene_list=all_DE_genes,
                                     ntop=1000){
  # calculate the variance for each gene
  rv <- rowVars(object)
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  if(is.null(second_annotation)==F){
    df_ <- as.data.frame(colData(dds_obj)[,c(first_annotation,second_annotation)])
    colnames(df_) <- c(first_annotation,second_annotation)
    row.names(df_) <- colnames(object)
    
    map <- pheatmap(object[select,],
                    cluster_rows=TRUE, cluster_cols=T,
                    clustering_distance_rows="correlation",
                    clustering_distance_cols="correlation",
                    show_rownames=F, show_colnames = TRUE,
                    scale = "row",
                    main = heatmap_title,
                    annotation_col=df_,
                    annotation_colors = anno_color) 
    name_output <<- object[map$tree_row$order,]
  }
  else{
    df_ <- as.data.frame(colData(dds_obj)[,first_annotation])
    colnames(df_) <- first_annotation
    row.names(df_) <- colnames(object)
    
    map <- pheatmap(object[select,],
                    cluster_rows=TRUE, cluster_cols=T,
                    clustering_distance_rows="correlation",
                    clustering_distance_cols="correlation",
                    show_rownames=F, show_colnames = TRUE,
                    scale = "row",
                    main = heatmap_title,
                    annotation_col=df_,
                    annotation_colors = anno_color) 
    name_output <<- object[map$tree_row$order,]
  }
}
Cluster_genelist_output <-function(object=rld_df, dds_obj=dds, anno_color,
                                   first_annotation="condition",
                                   second_annotation="Sex",
                                   heatmap_title, gene_list=all_DE_genes,
                                   display_row=F){
  
  if(is.null(second_annotation)==F){
    df_ <- as.data.frame(colData(dds_obj)[,c(first_annotation,second_annotation)])
    colnames(df_) <- c(first_annotation,second_annotation)
    row.names(df_) <- colnames(object)
    
    map <- pheatmap(object[gene_list,],
                    cluster_rows=TRUE, cluster_cols=T,
                    clustering_distance_rows="correlation",
                    clustering_distance_cols="correlation",
                    show_rownames=display_row, show_colnames = TRUE,
                    scale = "row",
                    main = heatmap_title,
                    annotation_col=df_,
                    annotation_colors =anno_color ) 
    Clustering_output <<- object[map$tree_row$order,]
  }
  else{
    df_ <- as.data.frame(colData(dds_obj)[,first_annotation])
    colnames(df_) <- c(first_annotation)
    row.names(df_) <- colnames(object)
    
    map <- pheatmap(object[gene_list,],
                    cluster_rows=TRUE, cluster_cols=T,
                    clustering_distance_rows="correlation",
                    clustering_distance_cols="correlation",
                    show_rownames=display_row, show_colnames = TRUE,
                    scale = "row",
                    main = heatmap_title,
                    annotation_col=df_,
                    annotation_colors =anno_color ) 
    Clustering_output <<- object[map$tree_row$order,] 
  }
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
plot_single_gene <-function(input=norm_anno, gene_symbol="TNF", 
                            condition="Genotype_Age", anno_colour=Genotype_Age_col,
                            shape_opt="age") {
  input<-as.data.frame(input)
  geneCounts_lfc <- as.data.frame(t(input[rownames(input)=="TNF",]))
  geneCounts_lfc$condition <- annotation[[condition]]
  colnames(geneCounts_lfc)<-c("count","condition")
  if(is.null(anno_colour)==F){
    if (is.null(shape_opt)==T){
      ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition)) +
        scale_y_continuous(expand=c(0.05,0.25)) +  
        scale_color_manual(values=anno_colour)+
        geom_beeswarm(cex = 3, na.rm=T)+
        ylab("Batch-correct rlog transformed counts")+
        labs(title=paste0(gene_symbol),colour=condition)+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5))+
        expand_limits(y=0)+
        geom_boxplot(width=.5,alpha=0)}
    else{
      geneCounts_lfc$sign <- annotation[[paste0(shape_opt)]]
      legend_shape<-paste0(shape_opt)
      ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition, shape=sign)) +
        scale_y_continuous(expand=c(0.05,0.25)) +   
        scale_color_manual(values=anno_colour)+
        scale_shape(name=legend_shape)+
        geom_beeswarm(cex = 3, na.rm=T)+
        ylab("Batch-correct rlog transformed counts")+
        labs(title=paste0(gene_symbol),colour=condition)+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5))+
        expand_limits(y=0)+
        geom_boxplot(width=.5,alpha=0)}}
  else{
    if (is.null(shape_opt)==T){
      ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition)) +
        scale_y_continuous(expand=c(0.05,0.25)) +  
        scale_color_brewer(palette = "Spectral")+
        geom_beeswarm(cex = 3, na.rm=T)+
        ylab("Batch-correct rlog transformed counts")+
        labs(title=paste0(gene_symbol),colour=condition)+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5))+
        expand_limits(y=0)+
        geom_boxplot(width=.5,alpha=0)}
    else{
      geneCounts_lfc$sign <- annotation[[paste0(shape_opt)]]
      legend_shape<-paste0(shape_opt)
      ggplot(geneCounts_lfc, aes(x = condition, y = count, colour=condition, shape=sign)) +
        scale_y_continuous(expand=c(0.05,0.25)) +  
        scale_color_brewer(palette = "Spectral")+
        scale_shape(name=legend_shape)+
        geom_beeswarm(cex = 3, na.rm=T)+
        ylab("Batch-correct rlog transformed counts")+
        labs(title=paste0(gene_symbol),colour=condition)+
        theme_bw()+
        theme(plot.title = element_text(hjust=0.5))+
        expand_limits(y=0)+
        geom_boxplot(width=.5,alpha=0)
    }
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

DE_genes_Top_up <- function(ntop=500, condi="NULL",input_file=DE_object$CON_GM, condition="merged"){
  if (condi=="NULL"){
    cond_list <- list()
    list_top_DE <- list()
    cond_list <- names(input_file@results)
    for (i in cond_list){
      df<-input_file@results[[i]][input_file@results[[i]]$log2FoldChange>input_file@parameters[[4]],]
      select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
      
      Top500_DE_genes <- rownames(df[select,])
      list_top_DE[[paste(i)]] <- assign(paste(i),Top500_DE_genes )}
    list_top_DE<-do.call(c, list_top_DE)
    list_top_DE <- unique(list_top_DE)
  }
  else {
    df<-input_file@results[[condi]][input_file@results[[condi]]$log2FoldChange>input_file@parameters[[4]],]
    select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
    
    Top500_DE_genes <- rownames(df[select,])
  }
}

DE_genes_Top_down <- function(ntop=500, condi="NULL",input_file=DE_object$CON_GM, condition="merged"){
  if (condi=="NULL"){
    cond_list <- list()
    list_top_DE <- list()
    cond_list <- names(input_file@results)
    for (i in cond_list){
      df<-input_file@results[[i]][input_file@results[[i]]$log2FoldChange<(-input_file@parameters[[4]]),]
      select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
      
      Top500_DE_genes <- rownames(df[select,])
      list_top_DE[[paste(i)]] <- assign(paste(i),Top500_DE_genes )}
    list_top_DE<-do.call(c, list_top_DE)
    list_top_DE <- unique(list_top_DE)
  }
  
  
  else {
    df<-input_file@results[[condi]][input_file@results[[condi]]$log2FoldChange<(-input_file@parameters[[4]]),]
    select <- order(df$padj ,decreasing = F, na.last=T )[seq_len(min(ntop, length(df$padj)))]
    
    Top500_DE_genes <- rownames(df[select,])
  }
}
functional_prediction <- function(input_file= dds,
                                  cond= condition,
                                  orga_type="human",
                                  DE_gene_list= list('OAvsCtrl_UP' = DE_object$` 24h_ctrl`@DE_genes$` 24h_OA`$`up-regulated genes`,
                                                     'OAvsCtrl_DOWN' = DE_object$` 24h_ctrl`@DE_genes$` 24h_OA`$`down-regulated genes`),
                                  gmtfile=file("h.all.v5.2.symbols.gmt")){
  normalized_counts <- counts(input_file, normalized=T)
  # create normdata file
  t <- t(normalized_counts)
  a <-as.data.frame(annotation[[cond]], rownames=rownames(annotation))
  rownames(a)<-rownames(annotation)
  colnames(a) <- c("merged")
  normdata <- merge(a, t, by = "row.names", all = TRUE)
  normdata <- data.frame(normdata[,-1],row.names = normdata[,1])
  #read in the expression data and save genes names as universe - to be used as background for GOEA
  universe <- colnames(normdata)
  
  ### define DE gene lists####
  #convert gene lists to dataframe
  if(length(DE_gene_list)==1){
    cluster_genes_original <- as.data.frame(DE_gene_list)
    l <- as.data.frame(gsub(".", "-", cluster_genes_original[[1]], fixed = TRUE))
    colnames(l) <- colnames(cluster_genes_original)
    cluster_genes_original <- l
    #class(cluster_genes_original)
    
  }
  if(length(DE_gene_list)>1){
    library(qpcR)
    cluster_genes_original <- as.data.frame(do.call(qpcR:::cbind.na,DE_gene_list))
    colnames <- colnames(cluster_genes_original)
    cluster_genes_original <- data.frame(lapply(cluster_genes_original, function(x) {gsub(".", "-",x, fixed=TRUE) }))
    
    #write.table(cluster_genes_original, "cluster_genes_original.txt", sep="\t", na="")
    #class(cluster_genes_original)
    #cluster_genes_original <- as.data.frame(t(as.data.frame.list(DE_gene_list)))
    # will also save a txt of your DE genes
    #write.table(cluster_genes_original, "cluster_genes_original.txt", sep="\t", na="")
    #cluster_genes_original <- read.delim("cluster_genes_original.txt", stringsAsFactor = FALSE, na.strings = "", check.names = FALSE, header = T)
    #cluster_genes_original <- cluster_genes_original[ , order(names(cluster_genes_original))]
  }
  
  # define where results are saved
  plotPath = file.path(getwd(), "clusterprofiler");
  dir.create(file.path(getwd(), "clusterprofiler"), showWarnings = FALSE)
  cluster_folder <- "clusterprofiler"
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  
  if(orga_type == "mouse"){
    universe_mouse_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = universe, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    universe_mouse_human <- universe_mouse_human[,2]
    universe_Entrez_mouse = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    universe_Entrez_mouse = unlist(universe_Entrez_mouse[2],use.names = FALSE)
    universe_Entrez_mouse_human = bitr(universe_mouse_human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez_mouse_human = unlist(universe_Entrez_mouse_human[2],use.names = FALSE)
  }else{
    universe_Entrez = bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez = unlist(universe_Entrez[2],use.names = FALSE)
  }
  
  #gmtfile <- file("h.all.v5.2.symbols.gmt")
  c1_hallmark_genes <- read.gmt(gmtfile)
  
  
  
  list_of_entrez <- list()
  
  
  
  ##
  for(id_it in 1:ncol(cluster_genes_original)){
    
    list_of_genes <- list(cluster_genes_original[,id_it])
    
    #list_of_genes <- cluster_genes_original[,1]
    na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
    list_of_genes <- na.omit.list(list_of_genes)
    #list_of_genes <- na.omit(list_of_genes)
    #length <- length(list_of_genes)
    #list_of_genes <- list(list_of_genes)
    #class(list_of_genes)
    
    #if mouse, mouse symbols are translated to human, because some gene sets works only with human
    if(orga_type == "mouse"){
      universe_orignal <- universe
      genes_cluster = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unlist(list_of_genes), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      cluster_genes <- genes_cluster[,2]
      list_of_genes_mouse_human <- list(cluster_genes)
      entrez_de = bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
      module_entrez_mouse <- unlist(entrez_de[2],use.names = FALSE)
      entrez_de = bitr(unlist(list_of_genes_mouse_human), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez_mouse_human <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez_mouse
      
    }else{
      entrez_de = bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez
    }
    
    wb <- createWorkbook(type="xlsx")
    #wb <- createWorkbook()
    
    
    pdf(paste(cluster_folder,"/ClusterProfiler_",colnames(cluster_genes_original)[id_it],".pdf",sep=""), onefile=FALSE,  width = 15, height = 15)
    
    font_size <- 8
    
    
    #Create figure window and layout
    plot.new()
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 2)))
    
    ####Hallmark gene sets####
    #start
    #egmt <- enricher(module_entrez, TERM2GENE=c1_hallmark, universe = universe_genes, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    
    if(orga_type == "mouse"){
      egmt <- enricher(unlist(list_of_genes_mouse_human), TERM2GENE=c1_hallmark_genes, universe = universe_mouse_human, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    }else{
      egmt <- enricher(unlist(list_of_genes), TERM2GENE=c1_hallmark_genes, universe = universe, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    }
    
    if(!is.null(egmt)){
      hallmark_plot <- dotplot(egmt, font.size = font_size, title = "Hallmark enrichment")
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
      print(hallmark_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Hallmark")
      #writeData(wb, 1, egmt@result)
      
      sheet <- xlsx::createSheet(wb, sheetName = "Hallmark")
      xlsx.addTable(wb, sheet, egmt@result)
    }
    #ende
    
    #####KEGG#####
    #start
    
    if(orga_type == "mouse"){
      KEGG <- enrichKEGG(gene = module_entrez_mouse, organism = 'mmu', pvalueCutoff=0.05,universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_Kegg<-NULL
      Enriched_Kegg_obj<-NULL
      reg_Mm=AnnotationDbi::select(org.Mm.eg.db,as.character(unlist(list_of_genes)),"ENTREZID","SYMBOL",multiVals="first")
      
      
      if(!is.null(KEGG))
      {
        if(nrow(summary(KEGG))>0)
        {
          df_kk<-as.data.frame(summary(KEGG))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_Kegg<-df_kk
          Enriched_Kegg_obj<-KEGG
        }
        else
        {
          Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_Kegg_obj<-NULL
        }
        
      }
      
    }else{
      KEGG <- enrichKEGG(gene = module_entrez, organism = 'hsa', pvalueCutoff=0.05,universe = universe_Entrez, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_Kegg<-NULL
      Enriched_Kegg_obj<-NULL
      reg_Hs=AnnotationDbi::select(org.Hs.eg.db,as.character(unlist(list_of_genes)),"ENTREZID","SYMBOL",multiVals="first")
      
      
      if(!is.null(KEGG))
      {
        if(nrow(summary(KEGG))>0)
        {
          df_kk<-as.data.frame(summary(KEGG))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_Kegg<-df_kk
          Enriched_Kegg_obj<-KEGG
        }
        else
        {
          Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_Kegg_obj<-NULL
        }
        
      }
      
    }
    
    
    if(!is.null(KEGG)){
      KEGG_plot <- dotplot(KEGG, font.size = font_size, title = "KEGG enrichment")
      pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
      print(KEGG_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "KEGG")
      #writeData(wb, 2, Enriched_Kegg)
      sheet <- xlsx::createSheet(wb, sheetName = "KEGG")
      xlsx.addTable(wb, sheet, Enriched_Kegg)
    }
    #ende
    
    #####Reactome#####
    #start
    if(orga_type == "mouse"){
      ReacTome <- enrichPathway(gene=module_entrez_mouse, organism = "mouse", pvalueCutoff=0.05, universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      
      Enriched_ReacTome<-NULL
      Enriched_ReacTome_obj<-NULL
      
      
      if(!is.null(ReacTome))
      {
        if(nrow(summary(ReacTome))>0)
        {
          df_kk<-as.data.frame(summary(ReacTome))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_ReacTome<-df_kk
          Enriched_ReacTome_obj<-ReacTome@result
        }
        else
        {
          Enriched_ReacTome<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_ReacTome_obj<-NULL
        }
        
      }
      
    }else{
      ReacTome <- enrichPathway(gene=module_entrez, organism = "human", pvalueCutoff=0.05, universe = universe_Entrez, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_ReacTome<-NULL
      Enriched_ReacTome_obj<-NULL
      
      
      if(!is.null(ReacTome))
      {
        if(nrow(summary(ReacTome))>0)
        {
          df_kk<-as.data.frame(summary(ReacTome))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_ReacTome<-df_kk
          Enriched_ReacTome_obj<-ReacTome@result
        }
        else
        {
          Enriched_ReacTome<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_ReacTome_obj<-NULL
        }
        
      }
    }
    
    
    if(!is.null(Enriched_ReacTome)){
      
      #ReacTome$Description <- factor(ReacTome$Description, levels = rev(unique(ReacTome$Description)))
      
      # tryCatch(ReacTome_plot <- dotplot(ReacTome, font.size = font_size, title = "Reactome enrichment"))
      #ReacTome_plot <- dotplot(KEGG, font.size = font_size, title = "KEGG enrichment")
      ReacTome_plot <- dotplot(ReacTome, font.size = font_size, title = "Reactome enrichment")
      pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
      print(ReacTome_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Reactome")
      #writeData(wb, 3, Enriched_ReacTome)
      sheet <- xlsx::createSheet(wb, sheetName = "Reactome")
      xlsx.addTable(wb, sheet, Enriched_ReacTome)
    }
    #ende
    
    ####enrichDO####
    #start
    if(orga_type == "mouse"){
      enrichDO <- enrichDO(gene= module_entrez_mouse_human,
                           ont           = "DO", 
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "none",
                           universe      = universe_Entrez_mouse_human, 
                           qvalueCutoff  = 1.0,
                           readable      = FALSE)
      
      Enriched_enrichDO<-NULL
      Enriched_enrichDO_obj<-NULL
      
      
      if(!is.null(enrichDO))
      {
        if(nrow(summary(enrichDO))>0)
        {
          df_kk<-as.data.frame(summary(enrichDO))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_enrichDO<-df_kk
          Enriched_enrichDO_obj<-enrichDO@result
        }
        else
        {
          Enriched_enrichDO<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_enrichDO_obj<-NULL
        }
        
      }
      
      
    }else{
      enrichDO <- enrichDO(gene= module_entrez,
                           ont           = "DO", 
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "none",
                           universe      = universe_Entrez, 
                           qvalueCutoff  = 1.0,
                           readable      = FALSE)
      
      Enriched_enrichDO<-NULL
      Enriched_enrichDO_obj<-NULL
      
      
      if(!is.null(enrichDO))
      {
        if(nrow(summary(enrichDO))>0)
        {
          df_kk<-as.data.frame(summary(enrichDO))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_enrichDO<-df_kk
          Enriched_enrichDO_obj<-enrichDO@result
        }
        else
        {
          Enriched_enrichDO<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_enrichDO_obj<-NULL
        }
        
      }
      
    }
    
    
    
    if(!is.null(Enriched_enrichDO)){
      enrichDO_plot <- dotplot(enrichDO, font.size = font_size, title = "Disease enrichment")
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
      print(enrichDO_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Disease")
      #writeData(wb, 4, enrichDO@result)
      sheet <- xlsx::createSheet(wb, sheetName = "Disease")
      xlsx.addTable(wb, sheet, Enriched_enrichDO)
    }
    #ende
    
    ####enrichGO_dotplot####
    #start
    if(orga_type == "mouse"){
      enrichGO <- enrichGO(gene = module_entrez_mouse,
                           universe = universe_Entrez_mouse,
                           OrgDb = org.Mm.eg.db,
                           # keytype = 'ENTREZID',
                           ont = "BP",
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 1,
                           readable      = T)
      
    }else{
      enrichGO <- enrichGO(gene = module_entrez,
                           universe = universe_Entrez,
                           OrgDb = org.Hs.eg.db,
                           #keytype = "ENTREZID",
                           ont = "BP",
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = T)
      
    }
    
    
    if(!is.null(enrichGO)){
      # BAUSTELLE: shorten the labels on the y axis
      enrichGO_plot <- dotplot(enrichGO, font.size = font_size, title = "GO enrichment",showCategory=20)
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
      print(enrichGO_plot, new = FALSE)
      popViewport()
      
      #addWorksheet(wb, "GO")
      #writeData(wb, 5, enrichGO@result)
      sheet <- xlsx::createSheet(wb, sheetName = "GO")
      xlsx.addTable(wb, sheet, enrichGO@result)
    }
    #ende
    
    #enrichGO_enrichment_map
    #start
    
    # pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    # par(fig = gridFIG(), new = TRUE)
    # enrichMap(enrichGO, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.8, n=10)
    # popViewport()
    #}
    #ende
    # 
    ####TF prediction####
    #plotFunction<-function(x=TFs_matrix_all_plot.melt.mean){p<- ggplot2:::ggplot(x,aes(x=variable,y=mean,fill=merged)) +
    #  geom_bar(stat="identity",position = "dodge") + facet_wrap(~variable,scales = "free")
    #  print(p)
    #}
    
    #start
    if(orga_type == "mouse"){
      
      r <- as.vector(list_of_genes[[1]])
      
      TFtable <- primo(as.vector(list_of_genes[[1]]), inputType = "geneSymbol", org = "Mm")
      #class(list_of_genes)
      
      
      TFs <- str_split(TFtable$overRepresented[,4], "_")
      
      TFs <- do.call(rbind, TFs)[,1]
      
      TFs <- str_split(TFs, "::")
      
      TFs <- do.call(rbind, TFs)
      
      if(grepl(":", TFs)&&!is.null(TFs)){
        
        TFs <- c(TFs[,1],TFs[,2])
        
      }
      
      TFs <-unique(TFs)
      
      TFs <- str_to_title(TFs, locale = "")
      
      
      
    }else{
      
      #TFtable <- primo(unlist(list_of_genes), inputType = "geneSymbol", org = "Hs")
      TFtable <- primo(as.vector(list_of_genes[[1]]), inputType = "geneSymbol", org = "Hs")
      
      
      TFs <- str_split(TFtable$overRepresented[,4], "_")
      
      TFs <- do.call(rbind, TFs)[,1]
      
      TFs <- str_split(TFs, "::")
      
      TFs <- do.call(rbind, TFs)
      
      if(grepl(":", TFs)&&!is.null(TFs)){
        
        TFs <- c(TFs[,1],TFs[,2])
        
      }
      
      TFs <-unique(TFs)
      
      TFs <- str_to_upper(TFs, locale = "")
      
    }  
    
    
    
    
    
    if(!is.null(TFs)){
      
      
      
      #addWorksheet(wb, "TFoverrepresented")
      
      #writeData(wb, 6, TFtable$overRepresented)
      
      
      
      
      
      
      
      if(length(TFs) != 0L){
        
        
        
        sheet <- xlsx::createSheet(wb, sheetName = "TFoverrepresented")
        
        xlsx.addTable(wb, sheet, TFtable$overRepresented)
        
        
        
        TFS_intersect <- intersect(TFs,colnames(normdata))
        
        
        
        if(length(TFS_intersect) != 0L){
          
          
          
          TFs_matrix_all_plot <- normdata[,c("merged",TFS_intersect)]
          
          
          
          TFs_matrix_all_plot.melt <- melt(TFs_matrix_all_plot, id="merged")
          
          # calculate means
          
          require(dplyr)
          library(dplyr)
          TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt %>% dplyr::group_by(merged,variable) %>% dplyr::summarise(mean=mean(value))
          # # get the mean over merged group (90 to 30 samples)
          # TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt %>% group_by(merged,variable)
          # # commas into dots
          # TFs_matrix_all_plot.melt.mean$value <- as.numeric(sub(",", ".", TFs_matrix_all_plot.melt.mean$value, fixed = TRUE))
          # # reorder the results according to group and not TF (!?)
          # TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt.mean %>%  dplyr::summarise(mean=mean(value)) 
          
          
          condis <- unique(normdata$merged)
          
          TFs_matrix_all_plot.melt.mean$merged <- factor(TFs_matrix_all_plot.melt.mean$merged,levels = condis)
          
          
          
          
          
          
          ####plotFunction####
          plotTheThing<-function(){
            
            p<- ggplot(TFs_matrix_all_plot.melt.mean,aes(x=variable,y=mean,fill=merged)) +
              
              geom_bar(stat="identity",position = "dodge") +
              
              facet_wrap(~variable,scales = "free")
            
            print(p)
            
          }
          
          
          # here you put the plot into the excel sheet
          r2excel::xlsx.addPlot(wb, sheet, plotTheThing)
          
        }
        
        
        
        TFS_intersect <- intersect(TFs,colnames(normdata))
        
        
        
        if(length(TFS_intersect) != 0L){
          
          
          
          
          
          TFs_matrix <- normdata[,c("merged",head(TFS_intersect,4))]
          
          
          
          TFs_matrix.melt <- melt(TFs_matrix, id="merged")        # calculate means
          
          #require(dplyr)
          
          
          
          TFs_matrix.melt %>% group_by(merged,variable) -> TFs_matrix.melt.mean
          
          TFs_matrix.melt.mean$value<-as.numeric(sub(",", ".", TFs_matrix.melt.mean$value, fixed = TRUE))
          
          TFs_matrix.melt.mean %>% summarise(mean=mean(value)) -> TFs_matrix.melt.mean
          
          condis <- unique(normdata$merged)
          
          #TFs_matrix.melt.mean$merged <- factor(TFs_matrix.melt.mean$merged,levels = condis)
          
          
          # for pdf
          TF_plot <- ggplot(TFs_matrix.melt.mean,aes(x=variable,y=mean,fill=merged)) +
            
            geom_bar(stat="identity",position = "dodge") +
            
            facet_wrap(~variable,scales = "free")
          
          
          
          
          # embed in pdf
          pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
          
          print(TF_plot, newpage = FALSE)
          
          popViewport()
          
        }
        
        
      }
      
    }
    
    xlsx::saveWorkbook(wb, paste(cluster_folder,"/ClusterProfiler_",colnames(cluster_genes_original)[id_it],".xlsx",sep=""))
    
    dev.off()
    
    
    
    
    
  }
}
