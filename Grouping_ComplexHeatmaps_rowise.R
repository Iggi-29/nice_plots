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

#### Grouping of rows
## First we have to import a dataframe that will need two columns.
# A column with the name of the patients accessions
# A column with the group of patients it belongs.
### Import the data
patients_grouping <- as.data.frame(readr::read_tsv("./raw_data/Grouped_Heatmap/groups_of_patients.tsv"))

## Subset the protein grouping
# Get the interesting useful columns
patients_grouping <- patients_grouping %>% 
  subset(select = c(patient, grouping))
rownames(patients_grouping) <- patients_grouping$patient

## Grouping the patients
expression_matrix <- numeric_values

result_matrix <- sapply(unique(patients_grouping$grouping), function(group) {
  groups_patients <- patients_grouping$patient[patients_grouping$grouping == group]
  group_expression <- expression_matrix[grepl(paste0("^", groups_patients, "$", collapse = "|"), rownames(expression_matrix)), 
                                                , drop = FALSE]
  group_mean <- colMeans(group_expression, na.rm = TRUE)
  group_mean[is.na(group_mean)] <- group_expression[which(is.na(group_mean), arr.ind = TRUE)]
  return(group_mean)
})

numeric_values <- result_matrix
remove(result_matrix)
numeric_values <- as.data.frame(numeric_values)
numeric_values <- as.matrix(numeric_values)
numeric_values <- t(numeric_values)
###numeric_values <- scale(numeric_values)

### subset numeric_values
proteins_to_plot <- openxlsx::read.xlsx("./raw_data/Grouped_Heatmap/proteins_to_plot.xlsx", sheet = 1)
numeric_values <- numeric_values[,
                                 colnames(numeric_values) %in% 
                                   proteins_to_plot$Accession[!is.na(proteins_to_plot$Which_heatmap)]]

colnames(numeric_values) <- proteins_to_plot$Gene.names[match(colnames(numeric_values),
                                                              proteins_to_plot$Accession)]

### Annotation
## dataframe
group_annot <- as.data.frame(tibble::tibble(
  pats = c(rownames(numeric_values)),
  Groups = c(rownames(numeric_values))
))
group_annot <- group_annot[,-1, drop = F]

## colors
colours <- list("Groups" = c("rest" = "yellow",
                            "C110" = "purple"))

## annotation formation
row_annot <- HeatmapAnnotation(df = group_annot,
                               col = colours,
                               which = "row")

## splits
split_reorder <- factor(group_annot$Groups, levels = c("C110","rest"))

### Do the heatmap
hmap_rowise <- Heatmap(numeric_values, 
                     cluster_row_slices = FALSE, 
                     row_split = split_reorder,
                     show_row_dend = F, show_column_dend = F, 
                     show_heatmap_legend = T,
                     column_names_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 12),
                     right_annotation = row_annot, column_title = "Rowwise grouped heatmap", name = "log2\n intensity",
                     cluster_columns = F, show_row_names = T)
hmap_rowise_drawn <- draw(hmap_rowise)
