# DEGene-pipeline

This pipeline includes the General workflow to find differentially expressed genes with DESeq starting from Kallisto files or Gene count tables including helpful plots.

#08.06.2018:
1) Fold change is now also included in the results table
2) Excel files with the DEgenes and their pvalues and Fold-changes can now be created automatically.
3) Creating an excel file with the normalized counts

# 07-06-2018: 
1) Now you can either add a colour to the PCA under anno_colour matching a colour to each factor in colour_obj, or not by settin anno_colour=NULL
2) Adjusted the Limma_batch command so you can add your column of interest here so the variance explained by this variable never gets lost.
