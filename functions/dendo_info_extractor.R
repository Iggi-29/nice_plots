############################
### dendo_info_extractor ###
############################

# heatm = a heatmap done w ComplexHeatmap, it has to be drawn !
# mat = the matrix that creates the map
# names_of_groups = names of the different groups
# grouped_element = the elements you are grouping: "proteins","patients" etc.

dendo_info_extractor <- function(heatm,side,
                                 mat,
                                 names_of_groups,
                                 grouped_element = "elements") {
  ### Reset the heatmap variable name, this is taken from the code that originated the matrix, and may be changed in the future
  hm <- heatm
  ## Pick the dendogram info depending wether if we retrive info from the columns or the rows 
  if (side == "column"){
    col_dend <- column_dend(hm)
    col_list <- column_order(hm)} 
  if (side == "row"){
    col_dend <- row_dend(hm)
    col_list <- row_order(hm)}
  # Add the group names
  names(col_list) <- names_of_groups
  # Print the message of the numbers
  for (i in 1:length(col_list)){
    print(paste(names(col_list)[i], " has ", length(col_list[[i]]), grouped_element))}
  # First dataframe
  element_data <- as.data.frame(tibble(
    element_number = unlist(col_dend),              ## The index of the cols or rows that we want to map 
    ##cluster = rep("Not done",
    ##              times = length(unlist(col_dend))), ## Initialize the cluster column
    element = rep("Not done",
                  times = length(unlist(col_dend))))) ## Initialize the element column
  ## Map the elements to the matrix
  for (element_ in 1:length(rownames(element_data))){
    if (side == "column"){
      element_data$element[element_] <- colnames(mat)[element_data$element_number[element_]]}
    if (side == "row"){
      element_data$element[element_] <- rownames(mat)[element_data$element_number[element_]]}
  }
  element_data$element_number <- as.character(element_data$element_number)
  ### Cluster info
  # Transform col_list to a dataframe
  col_list_dat <- plyr::ldply (col_list, data.frame)
  
  colnames(col_list_dat)[1] <- "group"
  colnames(col_list_dat)[2] <- "element_number"
  col_list_dat$element_number <- as.character(col_list_dat$element_number)
  
  element_data <- merge(element_data,col_list_dat, 
                        by = "element_number")
  
  ### Map the cluster name
  ###for (k in 1:length(rownames(element_data))){
  ###  element_data$cluster[k] <- col_list_dat$group[k][col_list_dat$element_number == element_data$element_number]}
  
  element_data <- element_data[,c(2,3)]
  colnames(element_data)[1] <- grouped_element
  return(element_data)
}