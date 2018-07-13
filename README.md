# DEGene-pipeline

This pipeline includes the General workflow to find differentially expressed genes with DESeq starting from Kallisto files or Gene count tables including helpful plots.

# 13.07.2018
1) DEA_analysis: added Fold-change as column in results of DE_object
2) fixed batch_corrected_table when using SVA
3) Added boxplots to single gene plotting 
4) PCA looks prettier now :D

# 22.06.2018_
1) Changed Volcano plot because it was displaying the wrong genes in red!
2) Changed the PCA of the normalized counts so you can now also plot it without a shape_option
3) Changed the Top_10_genes function as it calculated the wrong genes! (did not take fold-change threshold into account)

# 19.06.2018:
Changed the n.sv command

# 11.06.2018:
1) Dea_analysis: now the DEgene_number data-frame displays the right names.
2) DE_genes_plot now has the comparison of the conditions in the heading.

# 08.06.2018:
1) Fold change is now also included in the results table
2) Excel files with the DEgenes and their pvalues and Fold-changes can now be created automatically.
3) Creating an excel file with the normalized counts

# 07-06-2018: 
1) Now you can either add a colour to the PCA under anno_colour matching a colour to each factor in colour_obj, or not by settin anno_colour=NULL
2) Adjusted the Limma_batch command so you can add your column of interest here so the variance explained by this variable never gets lost.
