###################################
### Working with ComplexHeatmap ###
###################################

###library and WD
library(ComplexHeatmap)
library(rstudioapi)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tibble)
setwd(dirname(getActiveDocumentContext()$path)) 

### sourcing of in-house functions
source("./functions/dendo_info_extractor.R")

#### Import the data
## ANOVA prots this information is used just to remove some rows that we do not want to work with in the in the matrix
anova_dat <- read.table("./raw_data/Grouped_Heatmap/TopTable_Clean.tsv", sep = "\t", dec = ".", header = T)
anova_dat <- anova_dat %>% 
  filter(adj.P.Val < 0.05)
anova_prots <- anova_dat$Accession

## This is used to remove rows
my_patients <- read.table("./raw_data/Grouped_Heatmap/patients_meta_pearson_removed_out_k4.tsv", sep = "\t", dec = ".", header = T)
my_patients <- my_patients$patient

## raw_expression data
raw_data <- openxlsx::read.xlsx("./raw_data/Grouped_Heatmap/raw_data_ITACAT.xlsx",sheet = 1)
raw_data <- raw_data[raw_data$Accession %in% anova_prots,]

## expression matrix - formating the expression matrix
expression_mat <- openxlsx::read.xlsx("./raw_data/Grouped_Heatmap/raw_data_ITACAT.xlsx",sheet = 1)
expression_mat <- expression_mat[expression_mat$Accession %in% anova_prots,]
expression_mat <- expression_mat[,-c(1:3)]
rownames(expression_mat) <- raw_data$Accession
expression_mat <- expression_mat[,colnames(expression_mat) %in% my_patients]


#### Annotation 
## annotation for the patients - meta data limma
meta_data <- read.table("./raw_data/Grouped_Heatmap/patients_meta_pearson_removed_out_k4.tsv", sep = "\t", dec = ".", header = T)
rownames(meta_data) <- meta_data$patient
meta_data <- meta_data[,c(2,3)]

## generate the annotation df
group_annot <- meta_data
rownames(group_annot) <- group_annot$patient
group_annot <- group_annot[,2, drop = F]
colnames(group_annot) <- c("patient_group")

## color annotation
colours <- list(
  "patient_group" = c("group1" = "#ff00ff", "group2" = "#ffff00", 
                      "group3" = "#00ff99", "group4" = "#ff8100"))

## Generate the row annotation
row_annot <- HeatmapAnnotation(df = group_annot,
                               col = colours,
                               which = "row")

#### Reshape data for the heatmap
## prepare the matrix
numeric_values <- expression_mat[,1:ncol(expression_mat)]
rownames(numeric_values) <- raw_data$Accession
## save cols and rows
cols <- colnames(numeric_values)
rows <- rownames(numeric_values)

## make the numeric matrix a matrix
numeric_values <- as.matrix(numeric_values)
numeric_values <- matrix(as.numeric(numeric_values),
                         ncol = ncol(numeric_values))

## columns and rows for the numeric matrix
colnames(numeric_values) <- cols
rownames(numeric_values) <- rows

# prepare numeric data
numeric_values <- t(numeric_values)
numeric_values <- scale(numeric_values)

# Force numeric values to have the same ordering as group_annot dataframe
numeric_values <- numeric_values[rownames(group_annot),]

#### Grouping of columns
## First we have to import a dataframe that will need three columns.
# A column with the name of the proteins/genes.
# A column with the group of proteins/genes it belongs.
# A column with with the name of the Heatmap this group belongs.
### Import the data
proteins_to_individual_heatmap <- openxlsx::read.xlsx("./raw_data/Grouped_Heatmap/proteins_to_plot.xlsx", sheet = 2)
# Remove those entries with No heatmap assigned
proteins_to_individual_heatmap <- proteins_to_individual_heatmap[!is.na(proteins_to_individual_heatmap$Which_heatmap),]

## Change the colnames of numeric values
# Numeric values is annotated as Protein Group and we want to reffer ourselves to Genes, so we do this change.
x <- colnames(numeric_values)
colnames(numeric_values) <- anova_dat$Gene.names[match(x, anova_dat$Accession)]

## Subset the protein grouping
# Get the interesting useful columns
groupings_of_prots <- proteins_to_individual_heatmap %>% 
  subset(select = c(Gene.names, group_of_prot, Which_heatmap))
# If the proteins are not grouped but we want to add them to the plot
# the Which Heatmap column won't be NA but the group_of_prot column will be 
# NA, the group of prot will be used to the grouping calculation, thus, we paste
# the name of the heatmap and the gene to create an unique value and the mean of its'
# intensity will be calculated with only its value and it won't change
for (i in 1:length(rownames(groupings_of_prots))){
  if (is.na(groupings_of_prots$group_of_prot[i])){
    groupings_of_prots$group_of_prot[i] <- paste0(groupings_of_prots$Which_heatmap[i],
                                                  groupings_of_prots$Gene.names[i])}
}

## Get the different proteins that will be plotted in distincts heatmaps
# Proteasome Heatmap proteins
groupings_of_prots_proteasome <- groupings_of_prots[groupings_of_prots$Which_heatmap == 
                                                      "Proteasome",]
# Heme Metabolism proteins
groupings_of_prots_heme <- groupings_of_prots[groupings_of_prots$Which_heatmap == 
                                                "Heme_metabolism",]

## Subset the different proteins and create the different subsets of "numeric_values"
# proteasome
proteins_proteasome <- proteins_to_individual_heatmap$Gene.names[
  proteins_to_individual_heatmap$Which_heatmap == "Proteasome"]
numeric_proteasome <- numeric_values[,colnames(numeric_values) %in% proteins_proteasome]
# heme
proteins_heme <- proteins_to_individual_heatmap$Gene.names[
  proteins_to_individual_heatmap$Which_heatmap == "Heme_metabolism"]
numeric_heme <- numeric_values[,colnames(numeric_values) %in% proteins_heme]

## Grouping the proteins
# proteasome
gene_groups <- groupings_of_prots_proteasome
expression_matrix <- numeric_proteasome

result_matrix <- sapply(unique(gene_groups$group_of_prot), function(group) {
  group_genes <- gene_groups$Gene.names[gene_groups$group_of_prot == group]
  group_expression <- expression_matrix[, grepl(paste0("^", group_genes, "$", collapse = "|"), colnames(expression_matrix)), drop = FALSE]
  group_mean <- rowMeans(group_expression, na.rm = TRUE)
  group_mean[is.na(group_mean)] <- group_expression[which(is.na(group_mean), arr.ind = TRUE)]
  return(group_mean)
})

numeric_proteasome <- result_matrix
remove(result_matrix)
numeric_proteasome <- as.data.frame(numeric_proteasome)
numeric_proteasome <- numeric_proteasome %>% 
  subset(select = c(Activator,inhibitor,RPNs,RPTs,
                    alpha,beta,Ubiquitin_E1, Ubiquitin_E2, Ubiquitin_E3))
colnames(numeric_proteasome)[2] <- "Inhibitor"
numeric_proteasome <- as.matrix(numeric_proteasome)

# heme
gene_groups <- groupings_of_prots_heme
expression_matrix <- numeric_heme

result_matrix <- sapply(unique(gene_groups$group_of_prot), function(group) {
  group_genes <- gene_groups$Gene.names[gene_groups$group_of_prot == group]
  group_expression <- expression_matrix[, grepl(paste0("^", group_genes, "$", collapse = "|"), colnames(expression_matrix)), drop = FALSE]
  group_mean <- rowMeans(group_expression, na.rm = TRUE)
  group_mean[is.na(group_mean)] <- group_expression[which(is.na(group_mean), arr.ind = TRUE)]
  return(group_mean)
})

for (i in 1:length(colnames(result_matrix))){
  colnames(result_matrix)[i] <- gsub("Heme_metabolism","",colnames(result_matrix)[i])
}
numeric_heme <- result_matrix
numeric_heme <- as.data.frame(numeric_heme)
numeric_heme <- numeric_heme[,order(colnames(numeric_heme))]
numeric_heme <- numeric_heme %>% 
  relocate(c(Hemoglobin, Spectrins), .before = ADD1)
numeric_heme <- as.matrix(numeric_heme)

#### Heatmaps
## splits
split_reorder <- factor(group_annot$patient_group, levels = c("group1","group2","group3","group4"))

# proteasome
hmap_proteasome <- Heatmap(numeric_proteasome, 
                           cluster_row_slices = FALSE, 
                           column_names_rot = 90, column_title = "Proteasome and Ubiquitins",
                           show_row_dend = F, show_column_dend = F, column_title_rot = 30,
                           column_names_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 12),
                           row_split = split_reorder, cluster_columns = F)

hmap_proteasome_drawn <- draw(hmap_proteasome)

# heme metabolism
hmap_heme <- Heatmap(numeric_heme, 
                     cluster_row_slices = FALSE, 
                     row_split = split_reorder, column_names_rot = 90, column_title = "Heme metabolism",
                     show_row_dend = F, show_column_dend = F, show_heatmap_legend = FALSE, column_title_rot = 30,
                     column_names_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 12),
                     right_annotation = row_annot, cluster_columns = F)
hmap_heme_drawn <- draw(hmap_heme)

### HMAP list
list_of_hm <- hmap_proteasome + hmap_heme
draw(list_of_hm)
