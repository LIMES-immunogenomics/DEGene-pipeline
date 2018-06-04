##### Final general workflow DESeq2  #### 

###### General workflow analysing RNA-seq data using DESeq2 ######
setwd("E:/Nico/GitHub-Deseq2")

##### Install & load required packages ####
##### Installation
source("DESeq2 Edition_2_functions.R")


##### Annotation file #####

# Annotation file -> samples in rows & parameters in columns,
# Annotation file should have conditions called "conditions" that you want to compare,
# Make sure to set rownames as character, so they are not sorted by DESeq2

annotation <- read.table("anno_wo11-087.txt",
                         header = T,
                         sep = "\t",
                         check.names = F,
                         row.names = 1)

annotation[] <- lapply(annotation, factor)
str(annotation)

##### Define the colours of your conditions #####
anno <- c(CON="#B1BBCF",
          MS="#757575")

anno_batch <- c(GM="#0288D1",WM="#F44336")

anno_merged<-c(CON_GM="#B1BBCF", CON_WM="#757575", MS_GM="#0288D1",  MS_WM="#F44336" )

#################################################################################
##### Kallisto reads (jump this step if you have a Symbol Gene count table) ##### 
#################################################################################

# read in file with transcripts and genes
tx2gene = read.csv("E:/Nico/Shiny App/MS project/tx2genes_human.csv")

# If you do not have this, you should use the latest bioMart transcript lists
# dataset: mmusculus_gene_ensembl, hsapiens_gene_ensembl,rnorvegicus_gene_ensembl
# or if you don't know the name of the species search the following list

# mart = useMart("ensembl")
#listDatasets(mart)

#ensembl <- useMart(biomart="ensembl", dataset='hsapiens_gene_ensembl', host = "dec2017.archive.ensembl.org")
#bm <- getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id"), mart = ensembl)
#tid_gid <- data.frame(transcript_id=paste(bm[,1], bm[,2], sep = "."), gene_id=bm[,3])

#txi.kallisto <- tximport(files, type="kallisto", tx2gene=tid_gid)
#counts <- txi.kallisto$counts


# Define path where the Kallisto files are stored
basepath = "C:/Users/Nico/Documents/Master theesis/DESeq2/DESeq2/kallisto_output"
samples = dir(path=basepath, full.names=FALSE, no..=TRUE)
files <- file.path(basepath, samples, "abundance.h5")
names(files) <- samples

# Perform the distribution of transcripts to genes
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
# Extract the table with the counts as Ensemble Gene IDs only
counts <- txi.kallisto$counts

# Filter out genes with <10 counts in all conditions & 
# transform Ensembl Gene IDs to Gene Symbols
count_matrix <- filter_counts(count_mat=counts, 
                              keytype_opt="ENSEMBL",
                              change_names=T,
                              organism_option = org.Hs.eg.db,
                              threshold=10)

#################################################################################
#### Symbol Gene count table (use this if you already have a count table) ####
#################################################################################

#load the count table and call it count_matrix

count_matrix <- read.table("reads_wo11-087.txt",
                           header = T,
                           sep = "\t",
                           check.names = F,
                           row.names = 1)

# Filter out genes with <10 counts in all conditions
count_matrix <- filter_counts(count_mat=count_matrix, 
                              keytype_opt="SYMBOL",
                              change_names=F,
                              organism_option = org.Hs.eg.db,
                              threshold=10)



#################################################################################
###### Preparing DESeqDataSetFromMatrix #####
#################################################################################

# design formula: if you have a known batch effect, 
# you can enter that parameter from the annotation table in the first place of the design formula, 
# the last place in the design formula is the parameter over which the DE genes are later calculated


dsfm <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = annotation,
                               design = ~ merged )


# Run DESeq function
dds <- DESeq(dsfm)
dds

design <- design(dds)

# In case you need them for a different analysis
normalized_counts <- counts(dds, normalized=T)

#### Quality control ####
# Samples distribution before and after normalization 
# using log2 data to plot
# The +1 pseudocount step is neccessary so the value that fall below 0 due to the 
# multiplication with the sij, type in anything in scaling option and you will get 
# untransformed counts
# if you did not define colors in the beginning just type in NULL for anno_color
# enter the condition of the samples you want to be displayed in different colors

Checking_normalized_counts(raw_data=count_matrix, 
                           normalized_source=dds, 
                           scaling_opt="NULL",
                           anno_colour=NULL,
                           condition="merged")



##### Visual exploration of sample relationship ####

##### Variance stabilization ####

# Perform rlog transformation if n<30 samples 
# and variance stabilizing transformations (VST) if n > 30 samples,
# when applying this for visually exploring the data, blind should be set TRUE,
# if using the rld for later downstream analysis 
# (e.g. displaying DE genes in a heatmap) it should be set FALSE

rld<-Var_stab(dsfm,blind_param = T)
rld_df<-as.data.frame(assay(rld))


##### Sample-to-Sample Pearson correlation coefficient mattrix #####

### The limits will be set automatically, 
### but you need to edit the midpoint after seing how the scale is distributed

sample_cor(rld=rld_df, limit_set = c(NA,NA),midpoint_set = 0.99, condition="merged")

##### Sort samples by correlation #####

# Create a ggheatmap

sorted_pearson_correaltion(rld=rld_df,midpoint_set= 0.99,
                           limit_set = c(NA,NA), condition="merged")


##### PCA #####
# Since PCA can be slightly problematic with high dimensional data, 
# we first select only the 500 genes showing the highest variance. 
# The principal components to be ploted can be defined by PC_1 & PC_2
# if you have a continuous variable for the colour, change continuous to T

# PCA of raw, unnormalized counts
pca_raw <- PCA_rawdata(object=rld,
                       ntop=500, PC_1=1, PC_2=2,
                       color_obj="PMD",
                       count_mat=count_matrix,
                       shape_opt="type",
                       continuous=T,
                       colour_gradient=c("red","white","blue"),
                       point_size=3)
pca_raw[1]

# PCA of normalized, variance stabilized counts, you can either set shape_opt="NULL" or 
# a column of annotation that you would like to plot as shape

pca <- PCA_DEseq2(object=rld, 
                  ntop=500, PC_1=1, PC_2=2,
                  color_obj="PMD",
                  shape_opt = "merged",
                  continuous=T,
                  colour_gradient=c("red","yellow","blue"),
                  point_size=3)
pca[1]

##################################################################
##### Display the compensation for known batch effect #####
##################################################################

# Display compensation for batch effect using Limma, if no second batch name it "NULL"

Limma_batch(rld_obj=rld,
            object=rld_df, 
            batch_obj=annotation$`Customer ID`,
            batch2_obj=NULL,
            shape_opt="NULL",
            ntop=500, PC_1=1, PC_2=2,
            color_obj="merged",
            continuous=F,
            colour_gradient=c("blue","green","red"),
            point_size=3)


###############################################################################
###### Calculate surrogate variables in case you have unknown batches ######
###############################################################################

# Enter your variable to compare (your biological question) 
# behind the ~ & SVA defines a matrix for your conditions
mod  <- model.matrix(~ merged, annotation)
mod0 <- model.matrix(~   1, colData(dds))

# Calculate the marix corresponding to the hidden batch (surrogate variables),
# first run it with n.sv=NULL to see how many surrogate variables 
# the algorithm recommands you, if the results are not as you expected them to be
# you can try out different amounts of SVAs
# n.sv will define how many itterations SVA will perform

# find number of sva
n.sv = num.sv(dat,mod,method="leek")
svseq_NULL<- sva(assay(rld), mod, mod0, n.sv = 7)

#svaobj <-svaseq(assay(rld), mod, mod0, n.sv=7)

par(mfrow = c(2, 2), mar = c(3,5,3,1))
for (i in 1:7) {
  stripchart(svseq_NULL$sv[, i] ~ dds$merged, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

# Add surrogate variables to annotation table to reanalyse the batch effect via Limma

annotation$SV1 <- svseq_NULL$sv[,1]
annotation$SV2 <- svseq_NULL$sv[,2]
annotation$SV3 <- svseq_NULL$sv[,3]
annotation$SV4 <- svseq_NULL$sv[,4]
annotation$SV5 <- svseq_NULL$sv[,5]
annotation$SV6 <- svseq_NULL$sv[,6]
annotation$SV7 <- svseq_NULL$sv[,7]




##### Check on a PCA what influence the SVA has on the clustering of the samples ####

Limma_batch_sva(rld_obj=rld,
                object=rld_df, 
                condition=annotation$merged,
                batch_obj=annotation[c("SV1","SV2","SV3","SV4","SV5","SV6","SV7")],
                shape_opt="NULL",
                ntop=500, PC_1=1, PC_2=2,
                color_obj="merged",
                continuous=F,
                colour_gradient=c("red","green","blue"),
                point_size=3)

##### SVA into design formula ####

dds <- dds
dds$SV1 <- svseq_NULL$sv[,1]
dds$SV2 <- svseq_NULL$sv[,2]
dds$SV3 <- svseq_NULL$sv[,3]
dds$SV4 <- svseq_NULL$sv[,4]
dds$SV5 <- svseq_NULL$sv[,5]
dds$SV6 <- svseq_NULL$sv[,6]
dds$SV7 <- svseq_NULL$sv[,7]

###############################################################################
##### Redesigning the design matrix (if a batch was found) #####
###############################################################################


design(dds) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + merged
design <- design(dds)

dds %<>% DESeq

# In case you get the following error:
# "rows did not converge in beta, labelled in mcols(object)$betaConv. 
# Use larger maxit argument with nbinomWaldTest"
# run the following code

#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds, maxit=1000)

# If still some rows cannot be converged, you can just delet them
# as they are probably very low expressed genes: 
dds <- dds[which(mcols(dds)$betaConv),]




##### Prepare results table
# By default: FDR threshold: alpha=0.05 (5% chance that DE gene is not true); 
# IHW=T; lfcThreshold=0.58 corresponding to genes with more than 1.5 fold-change and 
# an adjusted p value below 0.05 are defined as significant DEG; 
# alternative options: independent filtering=TRUE if you set IHW_option=F; 
# lfcThreshold=0 (no log2 fold change threshold) 
# pAdjustMethod="BH"; extrem outliers will be marked as "NA"; 
# condition: the name of the variable from annotation that contains 
# control: here you can enter all the conditions that you want to compare the other groups against
# (or just one if you like...)


DE_object <- Dea_analysis(annotation_file=annotation,
                          count_matrix_file=count_matrix,
                          IHW_option=F,
                          alpha_option=0.05, 
                          lfc_Threshold=0, 
                          control=list("CON_GM"), 
                          condition="merged",
                          design_variable=design)

names(DE_object$CON_GM@results)
##### Get shrunken lfc object as data frame ####

list_df_CON_GM <- df_function(results_obj=DE_object$CON_GM@results)

##### MA-plot #####

df_names_CON_GM <- names(list_df)

### Either enter a name of a specific condition e.g.: fibronectin_df from the list_df here 
### or just type MA_function(condi=NULL) for multiplot of all conditions, 
# missing p values (NA) will not be displayed

MA_function(condition=NULL,
df_names=df_names_CON_GM, 
DE_obj=DE_object$CON_GM, 
y_lim=c(-10,10),
list_df=data.statistics)


#####p value distribution #####

# plot p values histograms 
# either set condi=NULL for multiplot or 
# set condition as one of the conditions from df_names;
# the statistical test is designed so that the p values 
# of 1 correspond to genes in the span between +LFC threshold 
# and -LFC threshold, thus non-DEG

plot_pval(df_names=df_names_CON_GM,
condition=NULL,
ylim_obj=c(0,10000), 
DE_obj=DE_object$CON_GM,
          list_df=data.statistics)

##### Plot up- & downregulated DE genes for all conditions ####

DE_genes_plot(DE_object$CON_GM@Number_DE_genes, DE_obj=DE_object$CON_GM)


##### How does the multiplot function work? #####

# Either plot a list and just define the amounts of columns
# multiplot(plotlist = list_volcano,cols = 2)

# Plot a list of plots with the order of the plot indicated in the matrix, 
# each number representing the position of the plot in the list,
# and define the amount of columns and rows behind and whether you want the plots entered by row 
# or by column (the order how they will be filled up)
# multiplot(plotlist = list_volcano,layout=matrix(c(2,1,3,4,5,6),3, 2, byrow = F))
#layout=matrix(c(2,1,3,4,5,6),3, 2, byrow = T)
#layout=matrix(c(5,2,3,4,1,6),3, 2, byrow = TRUE)


#### Volcano plot ####

# Enter the names of the conditions that you want to compare
Volcano_plot(input_file = DE_object$CON_GM, condition = "CON_WM", x_limit=c(-5,5), y_limit = c(NA,NA))

list_volcano <- list(Volcano_plot(input_file=DE_object$CON_GM, condition="CON_WM", x_limit=c(-5,5)), 
                     Volcano_plot(input_file=DE_object$CON_GM, condition="MS_GM", x_limit=c(-5,5)),
                     Volcano_plot(input_file=DE_object$CON_GM, condition="MS_WM", x_limit=c(-5,5)))

multiplot(plotlist = list_volcano,cols = 3)

# Volcano plot displaying name of a gene or a gene list

# Define and plot top genes for all conditions

Top_10_CON_WM_up <- DE_genes_Top_up(ntop=10, input_file=DE_object$CON_GM,condition="CON_GM", condi="CON_WM")
Top_10_CON_WM_up

Top_10_CON_WM_down <- DE_genes_Top_down(ntop=10, input_file=DE_object$CON_GM,condition="CON_GM", condi="NULL")

Top<-DE_genes_Top(ntop=50, condi="NULL", input_file=DE_object$CON_GM, condition=annotation$condition)
Top_10_CON_WM <- union(Top_10_CON_WM_up,Top_10_CON_WM_down)

Top_all <- DE_genes_Top(ntop=50, condi="NULL", input_file=DE_object$CON_GM, condition=annotation$condition)
# This plot will display the Volcano plot with a certain gene highlighted,
# condition equals the title of the plot
Volcano_plot_display_gene(input_file=DE_object$CON_GM, condition="CON_WM", 
                          gene=Top_10_CON_WM_up,size=5, 
                          x_limit=c(-10,10), y_limit = c(NA,NA))
Volcano_plot_display_gene(input_file=DE_object, condition="AB_H1N1", 
                          gene=Top_10_H1N1,size=2, 
                          x_limit=c(-10,10), y_limit = c(NA,NA))

# This plot displays the Top 10 varying genes (lowest p value for all genes) on the Volcano plot

#list_volcano <- list(Volcano_plot_display_gene(input_file=DE_object, 
#                                               condition="suspension_Cyth2_KO",Top_10_suspension_Cyth2_KO), 
#                     Volcano_plot_display_gene(input_file=DE_object, 
#                                               condition="fibronectin_wt",Top_10_fibronectin_wt),
#                     Volcano_plot_display_gene(input_file=DE_object, 
#                                               condition="fibronectin_Cyth2_KO",Top_10_fibronectin_Cyth2_KO), 
#                     Volcano_plot_display_gene(input_file=DE_object, 
#                                               condition="gelatin_wt",Top_10_gelatin_wt), 
#                     Volcano_plot_display_gene(input_file=DE_object, 
#                                               condition="gelatin_Cyth2_KO",Top_10_gelatin_Cyth2_KO))

#multiplot(plotlist = list_volcano,cols = 2)


##### Compar FCs ####



##### Venn Diagram #####
# Enter the list of DE_genes of each condition
# Define the colours in col and fill according to the anno vector defined at the beginning
venn.diagram(list(DE_object$CON_GM@DE_genes$CON_WM$`up-regulated genes`,
                  DE_object$CON_GM@DE_genes$MS_GM$`up-regulated genes`,
                  DE_object$CON_GM@DE_genes$MS_WM$`up-regulated genes`),
             filename = "Venn Diagram of upregulated genes.tiff",
             category.names = c("CON_WM", "MS_GM",
                                "MS_WM"),
             col = anno[-c(1)],
             fill = anno[-c(1)],
             alpha = 0.5,
             main = "Venn Diagram of upregulated genes")

venn.diagram(list(DE_object$CON_GM@DE_genes$CON_WM$`down-regulated genes`,
                  DE_object$CON_GM@DE_genes$MS_GM$`down-regulated genes`,
                  DE_object$CON_GM@DE_genes$MS_WM$`down-regulated genes`),
             filename = "Venn Diagram of downregulated genes.tiff",
             category.names = c("CON_WM", "MS_GM",
                                "MS_WM"),
             col = anno[-c(1)],
             fill = anno[-c(1)],
             alpha = 0.5,
             main = "Venn Diagram of downregulated genes")


##############################################
##### Not gone over yet #####
##############################################


#### Defining all DE genes ####

# Dinfe the DEG list that you later want to compare in a GO analysis
# union:union of the DEG lists, every gene just once; detdiff: DEG that differ between two lists;
# intersec: genes that overlap in two gene lists

# For H1N1
AB_H1N1_up <- DE_object@DE_genes$AB_H1N1$`up-regulated genes`
AB_H1N1_down <- DE_object@DE_genes$AB_H1N1$`down-regulated genes`
all_DE_genes <- union(AB_H1N1_up,AB_H1N1_down)
#

fibronectin_wt_genes_great <- DE_object@DE_genes$fibronectin_wt$`up-regulated genes`
fibronectin_Cyth2_KO_genes_great <- DE_object@DE_genes$`fibronectin_Cyth2_KO`$`up-regulated genes`
gelatin_wt_genes_great <- DE_object@DE_genes$gelatin_wt$`up-regulated genes`
gelatin_Cyth2_KO_genes_great <- DE_object@DE_genes$`gelatin_Cyth2_KO`$`up-regulated genes`

fibronectin_wt_genes_less <- DE_object@DE_genes$fibronectin_wt$`down-regulated genes`
fibronectin_Cyth2_KO_genes_less <- DE_object@DE_genes$`fibronectin_Cyth2_KO`$`down-regulated genes`
gelatin_wt_genes_less <- DE_object@DE_genes$gelatin_wt$`down-regulated genes`
gelatin_Cyth2_KO_genes_less <- DE_object@DE_genes$`gelatin_Cyth2_KO`$`down-regulated genes`

all_upreg_genes <- union(DE_object$CON_GM@DE_genes$CON_WM$`up-regulated genes`,
                         union(DE_object$CON_GM@DE_genes$MS_GM$`up-regulated genes`,DE_object$CON_GM@DE_genes$MS_WM$`up-regulated genes`))

all_downreg_genes <- union(DE_object$CON_GM@DE_genes$CON_WM$`down-regulated genes`,
                           union(DE_object$CON_GM@DE_genes$MS_GM$`down-regulated genes`,DE_object$CON_GM@DE_genes$MS_WM$`down-regulated genes`))

all_DE_genes<- union(all_downreg_genes, all_upreg_genes)

inter_fibro_wt_gel_wt_great <- intersect(fibronectin_wt_genes_great,gelatin_wt_genes_great)
inter_fibro_wt_gel_wt_less <- intersect(fibronectin_wt_genes_less, gelatin_wt_genes_less)

dif_fib_gel_ko_up <- setdiff(fibronectin_Cyth2_KO_genes_great, gelatin_Cyth2_KO_genes_great)
dif_fib_gel_ko_down <- setdiff(fibronectin_Cyth2_KO_genes_less, gelatin_Cyth2_KO_genes_less)

unique_fib_wt_great <- setdiff(fibronectin_wt_genes_great,union(fibronectin_Cyth2_KO_genes_great, union(gelatin_wt_genes_great,gelatin_Cyth2_KO_genes_great)))
unique_fib_wt_less <- setdiff(fibronectin_wt_genes_less,union(fibronectin_Cyth2_KO_genes_less, union(gelatin_wt_genes_less,gelatin_Cyth2_KO_genes_less)))

unique_fib_ko_great <- setdiff(fibronectin_Cyth2_KO_genes_great,union(fibronectin_wt_genes_great, union(gelatin_wt_genes_great,gelatin_Cyth2_KO_genes_great)))
unique_fib_ko_less <- setdiff(fibronectin_Cyth2_KO_genes_less,union(fibronectin_wt_genes_less, union(gelatin_wt_genes_less,gelatin_Cyth2_KO_genes_less)))

unique_gel_wt_great <- setdiff(gelatin_wt_genes_great,union(fibronectin_Cyth2_KO_genes_great, union(fibronectin_wt_genes_great,gelatin_Cyth2_KO_genes_great)))
unique_gel_wt_less <- setdiff(gelatin_wt_genes_less,union(fibronectin_Cyth2_KO_genes_less, union(fibronectin_wt_genes_less,gelatin_Cyth2_KO_genes_less)))

unique_gel_ko_great <- setdiff(gelatin_Cyth2_KO_genes_great,union(fibronectin_Cyth2_KO_genes_great, union(gelatin_wt_genes_great,fibronectin_wt_genes_great)))
unique_gel_ko_less <- setdiff(gelatin_Cyth2_KO_genes_less,union(fibronectin_Cyth2_KO_genes_less, union(gelatin_wt_genes_less,fibronectin_wt_genes_less)))

spec_gel_wt_up <- setdiff(gelatin_wt_genes_great, inter_fibro_wt_gel_wt_great)
spec_gel_wt_down <- setdiff(gelatin_wt_genes_less, inter_fibro_wt_gel_wt_less)

spec_fib_wt_up <- setdiff(fibronectin_wt_genes_great, inter_fibro_wt_gel_wt_great)
spec_fib_wt_down <- setdiff(fibronectin_wt_genes_less, inter_fibro_wt_gel_wt_less)

inter_fibronectin_wt_less <- setdiff(intersect(fibronectin_Cyth2_KO_genes_less,fibronectin_wt_genes_less),gelatin_wt_genes_less)
inter_gelatin_wt_less <- setdiff(intersect(fibronectin_Cyth2_KO_genes_less, gelatin_wt_genes_less),fibronectin_wt_genes_less)
inter_fibronectin_wt_gelatin_wt_less <- setdiff(intersect(fibronectin_wt_genes_less, gelatin_wt_genes_less),fibronectin_Cyth2_KO_genes_less)
inter_fibronectin_wt_LPS_gelatin_wt_less <- intersect(intersect(fibronectin_Cyth2_KO_genes_less,fibronectin_wt_genes_less),
                                                      gelatin_wt_genes_less)
inter_fibronectin_wt_LPS_gelatin_wt_less <- intersect(intersect(fibronectin_Cyth2_KO_genes_less,fibronectin_wt_genes_less),
                                                      gelatin_wt_genes_less)

spec_fib_wt_comp_ko_up <- setdiff(fibronectin_wt_genes_great, intersect(fibronectin_wt_genes_great,fibronectin_Cyth2_KO_genes_great))
spec_fib_wt_comp_ko_down <- setdiff(fibronectin_wt_genes_less, intersect(fibronectin_wt_genes_less, fibronectin_Cyth2_KO_genes_less))

spec_fib_ko_comp_wt_up <- setdiff(fibronectin_Cyth2_KO_genes_great, intersect(fibronectin_wt_genes_great, fibronectin_Cyth2_KO_genes_great))
spec_fib_ko_comp_wt_down <- setdiff(fibronectin_Cyth2_KO_genes_less, intersect(fibronectin_Cyth2_KO_genes_less, fibronectin_wt_genes_less))

spec_gel_wt_comp_ko_up <- setdiff(gelatin_wt_genes_great, intersect(gelatin_wt_genes_great, gelatin_Cyth2_KO_genes_great))
spec_gel_wt_comp_ko_down <- setdiff(gelatin_wt_genes_less, intersect(gelatin_wt_genes_less, gelatin_Cyth2_KO_genes_less))

spec_gel_ko_comp_wt_up <- setdiff(gelatin_Cyth2_KO_genes_great, intersect(gelatin_wt_genes_great, gelatin_Cyth2_KO_genes_great))
spec_gel_ko_comp_wt_down <- setdiff(gelatin_Cyth2_KO_genes_less, intersect(gelatin_Cyth2_KO_genes_less, gelatin_wt_genes_less))


##### Cluster genes function, as arg data to use and list of DE genes #####

# Plots all DE genes of all conditions &
# Output the genes table of the clustering of the DE genes

# Perform rlog transformation if n<30 samples 
# and variance stabilizing transformations (VST) if n > 30 samples,
# when applying this for visually exploring the data, blind should be set TRUE,
# if using the rld for later downstream analysis 
# (e.g. displaying DE genes in a heatmap) it should be set FALSE

rld<-Var_stab(dsfm,blind_param = F)
rld_df<-as.data.frame(assay(rld))

# create a data frame with the batch-corrected values of the expression
batch_corrected_rld <- removeBatchEffect(x=rld_df,batch = annotation$`Customer ID`,design = model.matrix(~annotation$merged))


# Will display and cluster a gene list and its samples
# if you have no second annotation, just add second_annotation=NULL
# either add rld as input to display the normalized, rlog-transformed data
# or add batch-corrected_rld as input to display the batch-corrected rlog-transformed data
# change second_annotation to NULL if you do not want the second batch do be displayed

Clustering_all_DEgenes_output <- Cluster_genelist_output(object=rld_df, dds_obj=DE_object$CON_GM@parameters[[1]], 
                                                         heatmap_title="All DE genes",
                                                         first_annotation="status",
                                                         second_annotation=NULL,
                                                         anno_color=list(status=anno,type=anno_batch), 
                                                         gene_list=all_DE_genes,
                                                         display_row=F)

# Only cluster the top 1000 varying gene
Cluster_top1000 <- cluster_Top_genes_output(object=rld_df, dds_obj=DE_object$CON_GM@parameters[[1]], 
                                            heatmap_title="Top 1000 genes", 
                                            ntop=1000,
                                            first_annotation="type",
                                            second_annotation=NULL,
                                            anno_color=list(type=anno_batch,status=anno), 
                                            gene_list=all_DE_genes)



##### GO enrichment plot #####
# clusterProfiler implements methods to analyze and visualize 
# functional profiles of gene and gene clusters, by default this will display the
# biological pathway "BP", but can be edited to "CC" for cellular compartment
# or "MF" for molecular function, make sure you change the organism to "org.Mm.eg.db" for mouse

# Plot all DE genes
GO_all <- plot_enrichGO(upreg_genes_list=AB_H1N1_up,
                        downreg_genes_list=AB_H1N1_down, 
                        downreg_title="GO BP downregulated all DEG", 
                        upreg_title="GO BP upregulated all DEG", 
                        nr_of_BPs_to_plot=20,
                        minGSSize_no = 10,
                        maxGSSize_no = 500,
                        keytype_opt="SYMBOL",
                        pvalueCutoff_opt=0.05,
                        ontology="BP",
                        organism="org.Hs.eg.db")

##### Return GO analysis as table #####

GO_all_up <- as.data.frame(enrichGO(AB_H1N1_up, 
                                    OrgDb = "org.Hs.eg.db",
                                    keyType="SYMBOL",
                                    ont="BP",
                                    pvalueCutoff=0.05))

##### Plot genes of GO pathways ###
# If you have no second annotation, enter NULL WITHOUT ''!
Plot_pathway_cluster(gene_list=AB_H1N1_up,
                     organism="org.Hs.eg.db",
                     ontology="BP",
                     keytype_opt="SYMBOL",
                     pvalue_cutoff=0.05,
                     pathway="defense response to virus",
                     first_annotation="Treatment",
                     second_annotation="Type",
                     anno_color=list(Treatment=anno,Type=anno_batch),
                     object=rld,dds_obj=DE_object$CON_GM@parameters[[1]])


##### KEGG enrichment plot #####
# organism_KEGG= "mmu" for mouse & "hsa" for human
KEGG_LPS <- plot_KEGG(upreg_genes_list=AB_H1N1_up,
                      downreg_genes_list=AB_H1N1_down, 
                      downreg_title="KEGG downregulated all DEG", 
                      upreg_title="KEGG upregulated all DEG", 
                      organism_KEGG = "hsa", keyType = "SYMBOL",
                      minGSSize_no = 10,
                      maxGSSize_no = 500,
                      pvalueCutoff_opt = 0.05, 
                      nr_of_BPs_to_plot = 20,
                      organism_opt = org.Hs.eg.db)


##### Plot single genes ####

#plots a single gene indicating the normalized counts for all conditions
# This duntion will plot the normalized count
# you can set anno_colour to NULL to make it chose the colours itself with color_brewer

plot_single_gene(dds_object=DE_object$CON_GM@parameters[[1]], gene_symbol="IL6", 
                 condition="merged", pc_cond=F, anno_colour=NULL, 
                 order=c("CON_GM", "CON_WM", "MS_GM",  "MS_WM"),shape_opt="type")

# This function plots the batch-corrected counts
plot_batch_corrected_counts(batch_rld=batch_corrected_rld, gene_symbol="IL6", 
                            condition="merged", anno_colour=NULL,
                            order=c("CON_GM", "CON_WM", "MS_GM",  "MS_WM"),shape_opt="type")

##### Creating tables of GO terms #####

Table_GO <- cbind(head(GO_CpG_up[,"Description"], 10),
                  head(GO_LPS_up[,"Description"], 10),
                  head(GO_PolyIC_up[,"Description"], 10),
                  head(GO_CpG_down[,"Description"], 10),
                  head(GO_LPS_down[,"Description"], 10),
                  head(GO_PolyIC_down[,"Description"], 10))

colnames(Table_GO)<- c("CpG","LPS+IFNg","PolyI:C","CpG","LPS+IFNg","PolyI:C")

write.csv(Table_GO, "GO pathways.csv")

Table_GO_up <- cbind(head(GO_CpG_up[,"Description"], 10),
                     head(GO_LPS_up[,"Description"], 10),
                     head(GO_PolyIC_up[,"Description"], 10))

colnames(Table_GO_up)<- c("CpG","LPS+IFNg","PolyI:C")


##### Genes of interest #####

GOI <- c("Cnr1", "Cnr2","Dagla", "Daglb","Napepld","Faah","Mgll","AbhD6","AbhD12",  "Cd40", "Icam1", "Il6", "Tnf", "Ccl2")

# Add FC values to the DE_object for all conditions
DE_object@results$CpG$FC <- 2^(DE_object@results$CpG$log2FoldChange)
DE_object@results$LPS_IFNg$FC <- 2^(DE_object@results$LPS_IFNg$log2FoldChange)
DE_object@results$PolyIC$FC <- 2^(DE_object@results$PolyIC$log2FoldChange)

# Extract the FC of the GOI from the different conditions


write.csv(GOI_FC_table, "FC of GOI.csv")

# Look at the corresponding p values
list_df$CpG_df[GOI,]$padj
list_df$LPS_IFNg_df[GOI,]$padj
list_df$PolyIC_df[GOI,]$padj

# Create a table with the FC of all conditions for those genes
GOI_padj_table <- cbind(list_df$CpG_df[GOI,]$padj,list_df$LPS_IFNg_df[GOI,]$padj,list_df$PolyIC_df[GOI,]$padj)
colnames(GOI_padj_table) <- c("CpG","LPS_IFNg", "PolyI:C")
rownames(GOI_padj_table) <- GOI

write.csv(GOI_padj_table, "Adjusted p values of GOI.csv")
