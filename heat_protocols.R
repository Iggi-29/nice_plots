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
### Expression matrix
expression_matrix <- as.data.frame(readr::read_tsv("./raw_data/expression_matrix_noExos.AA.Pellet.tsv"))

## Reshape the expression matrix a bit, this is not entirelly necessary but it can make our graphs a bit cooler
# Create the column Gene.names_1
expression_matrix <- expression_matrix %>%
  group_by(`Protein Group`) %>%  
  mutate(Gene.names_1 = unlist(strsplit(split = "\\;", Gene.names),use.names = F)[1]) %>% 
  ungroup() %>% 
  relocate(Gene.names_1, .after = Gene.names)

# Now change the column names of the expression matrix
colnames(expression_matrix) <- gsub(pattern = "Exos\\.", replacement = "", colnames(expression_matrix))
colnames(expression_matrix) <- gsub(pattern = "\\.pellet", replacement = "", colnames(expression_matrix))
colnames(expression_matrix) <- gsub(pattern = "\\.\\.Pellet", replacement = "", colnames(expression_matrix))
colnames(expression_matrix) <- gsub(pattern = "\\.Pellet", replacement = "", colnames(expression_matrix))
colnames(expression_matrix) <- gsub(pattern = "\\.[Mm]edian", replacement = "", colnames(expression_matrix))

# Filter for the repeated Gene.names_1 and for NAs
expression_matrix <- expression_matrix %>% 
  filter(!duplicated(Gene.names_1)) %>% 
  filter(Gene.names != "")

### Meta data file
meta_data <- as.data.frame(readr::read_tsv("./raw_data/meta_data.tsv"))

### Reshape the data to median-normalization and use this normalized data to work on
samples <- colnames(expression_matrix)[7:ncol(expression_matrix)]

long_format <- expression_matrix %>%
  tidyr::pivot_longer(cols = all_of(samples),
                      names_to = "sample_name",
                      values_to = "intens")

# median of all the experiment
median_all <- median(long_format$intens)
# median of each group
long_format <- long_format %>%
  group_by(sample_name) %>%
  mutate(MED = median(intens))
# do the normalization
long_format$normalized_intensity <- (long_format$intens - long_format$MED) + median_all 
# check for the normalization w a ggpplot
ggplot(long_format, mapping = aes(x = sample_name, y = normalized_intensity))+
  geom_boxplot(data = long_format, mapping = aes(x = sample_name, y = normalized_intensity, fill = sample_name))+
  theme_bw()+
  xlab("Sample") +
  ylab("log2 Detection Intensity")+
  ggtitle("Intensity of detection (Normalized intensities)")+
  theme(legend.position= "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#### Numeric values creation
# retrive the expression matrix and save it as numeric_values
# numeric values will be the matrix we will be able to use ComplexHeatmap with
expression_matrix <- reshape2::dcast(long_format, Gene.names_1 ~ sample_name, value.var="normalized_intensity")
rownames(expression_matrix) <- expression_matrix$Gene.names_1
expression_matrix <- expression_matrix[,-c(1)] 
numeric_values <- expression_matrix[,1:ncol(expression_matrix)]

## reshape numeric_values and scale it
# save cols and rows
cols <- colnames(numeric_values)
rows <- rownames(numeric_values)
# make the numeric matrix a matrix
numeric_values <- as.matrix(numeric_values)
# this will make it numeric
numeric_values <- matrix(as.numeric(numeric_values),
                         ncol = ncol(numeric_values))
# columns and rows for the numeric matrix
colnames(numeric_values) <- cols
rownames(numeric_values) <- rows
# scale the numeric values data
numeric_values <- t(numeric_values)
numeric_values <- scale(numeric_values)

########### NUMERIC VALUES IS OK TO GO!! ##################################
### Now it is time to go for group annotation,                          ###
### for this, we will need the meta_data that we have imported earlier  ###
###########################################################################

#### Annotation for the sample groups
# group annot creation, this will be the dataframe to annotate the rows of the Heatmap
group_annot <- meta_data
rownames(group_annot) <- group_annot$samples
group_annot <- group_annot[,-1, drop = F]
colnames(group_annot) <- "Group"
group_annot <- group_annot %>% 
  arrange(Group)

# color as a named list different entries of the columns of group_annot
colors = list("Group" = c("group1" = "#009E73",
                          "another_sample_1" = "red",
                          "another_sample_2" = "green",
                          "post" = "orange", "prev" = "blue"))
## Do the annotation as a ComplexHeatmap HeatmapAnnotation oject
row_annot <- ComplexHeatmap::HeatmapAnnotation(df = group_annot,
                                               col = colors,
                                               which = "row") 
# !! which, indicates if, in the matrix of the Heatmap items to annotate or columns are in rows

## Map the rownames of the numeric values to the ones in the group_annot
numeric_values <- numeric_values[rownames(group_annot),]

###################################
### WORKING WITH COMPLEXHEATMAP ###
###################################

### Very basic ComplexHeatmap
ssGSEA_standard_heatmap <- ComplexHeatmap::Heatmap(numeric_values,                            ## The matrix to plot
                                                   show_column_names = F, show_row_names = T, ## Show rownames and do not show colnames 
                                                   row_names_gp = gpar(fontsize = 10),        ## Rownames graphic options
                                                   column_names_gp = gpar(fontsize = 6),      ## Colnames graphic options
                                                   cluster_rows = F, cluster_columns = F,     ## Do not cluster rows and colums, the order of the matrix will prevail
                                                   right_annotation = row_annot,              ## The annotation of the 
                                                  ####column_km = 4, row_km = 5,            #### BORRAR AIXO !!
                                                   name = "Zscored protein \n expression",    ## This is the name that will appear in the legend
                                                   column_title = "Integrated protein expression \n profile",) ## Column title acts like if it was a general title
ssGSEA_standard_heatmap <- draw(ssGSEA_standard_heatmap, padding = unit(c(25, ## Down padding
                                                                          5, ## Right padding
                                                                          12.5,## Top padding
                                                                          5), ## Left padding
                                                                        "mm")) ## Padding units

### Correlation-based Heatmap
# This heatmap uses a correlation index line the Pearson, but can use any other to do the work
## Things to define before
# Do the correlation matrix
pearson.cor <- (expression_matrix) %>% 
  cor(use = "pairwise.complete.obs", method = "pearson")
pearson.cor <- as.matrix(pearson.cor)
# Set the color scheme
color.scheme <- rev(brewer.pal(11,"RdBu"))
# Heatmap
pearson_hm <- ComplexHeatmap::pheatmap(pearson.cor,                             ## The correlation 
                                       color = color.scheme,                    ## Color scheme     
                                       cluster_rows = T, cluster_cols = T,      ## Do the clustering
                                       cutree_cols = 3, cutree_rows = 3,        ## Cut the tree to make n groups
                                       border_color = NA,                     ## Border of the squares, hint, put it as NA
                                       show_rownames = F, show_colnames = T,    ## Show or not rownames
                                       name = "k6",                             ## Name on the legend
                                       main = "Pearson Correlation")            ## Title of the plot
pearson_hm <- draw(pearson_hm)
# https://jokergoo.github.io/ComplexHeatmap-reference/book/integrate-with-other-packages.html 
# At the 10.1.2 Translation section there is the translation between functions of pheatmap and Heatmap of ComplexHeatmap

### Extract information from a dendogram !
## Supose you have done a heatmap on the correlation, you do 3 sample groups and you want to extract the info from this analysis
## We will use this function to extract this data source this from this repo <3
# It also works with Heatmaps type heatmaps
groupping_info <- dendo_info_extractor(heatm = pearson_hm,
                           side = "column",
                           mat = pearson.cor, 
                           names_of_groups = c("group1","group2","group3"), grouped_element = "pacient")

#### Reshape the ComplexHeatmap
numeric_values <- numeric_values[rownames(group_annot),
                                 proteins_data$protein]
split1 <- factor(group_annot$Cluster, levels = c("Cluster1","Cluster2","Cluster3",
                                                 "s093","s097"))
split2 <- factor(c(proteins_data$cluster), levels = c("group1","group2","group3",
                                                      "group4","group5","group6"))

hm <-Heatmap(numeric_values, show_column_names = FALSE,
             right_annotation = row_annot,
             row_labels = rownames(numeric_values), 
             row_names_gp = gpar(fontsize = 8),
             column_split = split2, cluster_column_slices = F,
             row_split = split1, cluster_row_slices = F,
             name = "Zscored \n protein abundance")
hm <- draw(hm)
