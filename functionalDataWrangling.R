# install and load packages
if (!requireNamespace("ggalt", quietly = TRUE))
  install.packages("ggalt")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("ggfortify", quietly = TRUE))
  install.packages("ggfortify")
if (!requireNamespace("plotly", quietly = TRUE))
  install.packages("plotly")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("PCAtools", quietly = TRUE))
  install.packages("PCAtools")
if (!requireNamespace("docstring", quietly = TRUE))
  install.packages("docstring")
if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"vegan" %in% installed.packages()) install.packages("vegan")
if (!"stringr" %in% installed.packages()) install.packages("stringr")
if (!"readxl" %in% installed.packages()) install.packages("readxl")

if (!"corrplot" %in% installed.packages()) install.packages("corrplot")

library("docstring")
library("dplyr")
library("ggplot2")

##### File parsing
# Read and prepare FIA pos/neg table
read_fia_table <- function(table_path, sheet = "pos", fix_names = FALSE, sort_table = TRUE){
  feature_table <- readxl::read_excel(path = table_path, sheet = "pos", col_names = TRUE)
  # Retain only necessary columns
  fia_df <- cbind(feature_table[, 2], feature_table[, 5:ncol(feature_table)])
  # Transpose table
  fia_df_t <- t(fia_df)
  # Set column names as metabolites names
  colnames(fia_df_t) <- fia_df_t[1,]
  # Remove metnames column
  fia_df_t <- fia_df_t[2:nrow(fia_df_t),]
  # Reconvert to dataframe
  fia_df_t <- as.data.frame(fia_df_t)
  # Transform to numeric all columns
  fia_df_t[,1:ncol(fia_df_t)] <- sapply(fia_df_t[,1:ncol(fia_df_t)],as.numeric)
  # Check the data types of columns
  #print(sapply(fia_df_t, class))
  if (fix_names) {
    colnames(fia_df_t) <- make.names(colnames(fia_df_t), unique=TRUE)
  }
  
  if (isTRUE(sort_table)) {
    fia_df_t <- fia_df_t[order(row.names(fia_df_t)), ] # sort my row names (sample names)
  }
  
  return(fia_df_t)
}

read_ft_1 <- function(path, sort_by_names = FALSE){
  # Read Hitschickers guide style export feature table
  ft <- read.csv(path, header = TRUE, row.names = 1) #read csv table
  if(isTRUE(sort_by_names)){
    ft <- ft[order(row.names(ft)), ] # sort my row names (sample names)
  }
  return(ft)
}

read_metadata <- function(path, sort_table = FALSE){
  md <- read.csv(path, row.names = 1)
  if(isTRUE(sort_table)){
    md <- md[order(row.names(md)), ] # sort my row names (sample names)
  }
  return(md)
}

read_metadata_xls <- function(path, sheet, sort_table = FALSE){
  md <- readxl::read_excel(path, sheet = sheet)
  
  md <- as.data.frame(md)
  row.names(md) <- md[,1]
  md <- md[2:ncol(md)]
  
  if(isTRUE(sort_table)){
    md <- md[order(row.names(md)), ] # sort my row names (sample names)
  }
  return(md)
}

# Remove highly variable metabolites
filter_by_error <- function(feature_table, metadata_table, grouping_var = NULL, error_threshold = 25){
  # Step 1: Join grouping column in metadata to feature table.
  if (!all.equal(rownames(feature_table), rownames(metadata_table))) {
    print("Sample names in feature table and metadatable are not identical")
    return()
  }else{
    print("Sample names in feature table and metadatable are identical :)")
  }
  
  feature_table <- cbind(metadata_table[grouping_var], feature_table)
  
  # Step 2: Calculate the error for each variable per type
  errors <- feature_table %>%
    group_by(ATTRIBUTE_Sample) %>%
    summarise(across(where(is.numeric), ~ (sd(.) / mean(.)) * 100, .names = "error_{col}"))
  
  # View the errors dataframe
  #print(head(errors))
  
  # Step 3: Average the error for each variable for all types
  avg_errors <- errors %>%
    summarise(across(starts_with("error_"), mean, na.rm = TRUE))
  
  # View the average errors dataframe
  #print(head(avg_errors))
  
  variables_to_keep <- names(avg_errors)[avg_errors <= error_threshold]
  variables_to_keep <- gsub("error_", "", variables_to_keep)
  
  # Keep the metadata columns and the variables with error below the threshold
  df_filtered <- feature_table %>%
    select(all_of(variables_to_keep))
  
  # View the filtered dataframe
  return(df_filtered)
}

##### Graphing functions
get_palette <- function(nColors = 50){
  #' @title Get random colors for plotting.
  #'
  #' @description Returns a vector of "nCOlors" number random of colors to work in ggplot2. Max colors to ask for is 60.
  #' @param nCOlors integer. The number of colors to return.
  #' @value A vector of nColors.
  #' @details
  #' The list of colors is:
  #' "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2",
  #' "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3",
  #' "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
  #' "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
  #' "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
  #' "lavenderblush4", "lawngreen", "indianred1", "lightblue1", "honeydew4",
  #' "hotpink", "#e3ae78", "#a23f3f", "#290f76", "#ce7e00",
  #' "#386857", "#738564", "#e89d56", "#cd541d", "#1a3a46",
  #' "#9C4A1A", "#ffe599", "#583E26", "#A78B71", "#F7C815",
  #' "#EC9704", "#4B1E19", "firebrick2", "#C8D2D1", "#14471E",
  #' "#6279B8", "#DA6A00", "#C0587E", "#FC8B5E", "#FEF4C0",
  #' "#EA592A", "khaki3", "lavenderblush3", "indianred4", "lightblue",
  #' "honeydew1", "hotpink4", "ivory3", "#49516F", "#502F4C",
  #' "#A8C686", "#669BBC", "#29335C", "#E4572E", "#F3A712",
  #' "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD
  colors_vec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2",
                  "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3",
                  "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
                  "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
                  "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
                  "lavenderblush4", "lawngreen", "indianred1", "lightblue1", "honeydew4",
                  "hotpink", "#e3ae78", "#a23f3f", "#290f76", "#ce7e00",
                  "#386857", "#738564", "#e89d56", "#cd541d", "#1a3a46",
                  "#9C4A1A", "#ffe599", "#583E26", "#A78B71", "#F7C815",
                  "#EC9704", "#4B1E19", "firebrick2", "#C8D2D1", "#14471E",
                  "#6279B8", "#DA6A00", "#C0587E", "#FC8B5E", "#FEF4C0",
                  "#EA592A", "khaki3", "lavenderblush3", "indianred4", "lightblue",
                  "honeydew1", "hotpink4", "ivory3", "#49516F", "#502F4C",
                  "#A8C686", "#669BBC", "#29335C", "#E4572E", "#F3A712",
                  "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD")
  
  #set.seed(1)
  
  return(colors_vec[sample(1:length(colors_vec), size = nColors)])
}

# Do PCA plot, prev. fia_pca
ft_pca_1 <- function(feature_table, metadata_table, grouping_col, encircle = FALSE){
  # transposing feature table
  ft_t <- t(feature_table)
  if (isTRUE(all.equal(colnames(ft_t),row.names(metadata)))) {
    print("Sample names in feature table and metadatable are identical :)")
    
    fia_pca <- PCAtools::pca(ft_t, scale = TRUE, metadata = metadata_table, transposed = FALSE)
    
    p2 <- PCAtools::biplot(fia_pca, showLoadings = TRUE, ntopLoadings = 2, lab = NULL, colby = grouping_col,
                           legendPosition = "right", axisLabSize = 8, legendLabSize = 8, legendIconSize = 2, pointSize = 1.5,
                           colkey = get_palette(nColors = 60), encircle = encircle)
    
    p2
  }else{
    print("Sample names in feature table and metadatable are not identical")
    #print(all.equal(colnames(ft_t),metadata_table$Sample))
  }
  
}

ft_pca_2 <- function(feature_table, metadata_table, grouping_col = NULL, p_shape = NULL, dist_method = "euclidean"){
  # Compute distance matrix according to dist_method
  dist_matrix <- vegan::vegdist(feature_table, method = dist_method)
  
  # Perform PCoA
  pcoa_results <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  
  pcoa_df <- as.data.frame(pcoa_results$points)
  colnames(pcoa_df) <- c("PC1", "PC2")  # Rename axes
  pcoa_df$Sample <- rownames(feature_table)  # Add sample names to pcoa_df
  
  metadata_table$Sample <- rownames(metadata_table)  # Add sample names to pcoa_df
  
  #print(pcoa_df$Sample == row.names(metadata_table))
  if (identical(pcoa_df$Sample,row.names(metadata_table))) {
    print("Sample names in feature table and metadatable are identical :)")
    
  } else{
    print("Sample names in feature table and metadatable are not identical")
    return()
  }
  
  pcoa_df <- left_join(pcoa_df, metadata_table, by = c("Sample" = "Sample"))
  colour_palette <- get_palette(nColors = 20)
  #print(colour_palette)
  
  if (is.null(grouping_col) && is.null(p_shape)) {
    ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2")
  }else if(!is.null(grouping_col) && is.null(p_shape)){
    print("Plotting with grouping variable")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }else if(!is.null(grouping_col) && !is.null(p_shape)){
    print("Plotting with two grouping variables")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]], shape = .data[[p_shape]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }
  
}

############ Extracts data that matches a list of compounds froma  feature table.
extract_features_comparison <- function(feature_table, sig_features, columns_to_preserve = NULL){
  # Create list of columns to extract
  if (!is.null(columns_to_preserve)) {
    columns_to_preserve <- c(columns_to_preserve, sig_features)
  }else{
    columns_to_preserve <- sig_features
  }
  
  #print(columns_to_preserve)
  
  temp_df <- dplyr::select(feature_table, any_of(columns_to_preserve))
  #print(head(temp_df))
  #temp_df <- temp_df[2:length(temp_df)]
  ###temp_df_t <- as.data.frame(t(temp_df))
  temp_df_t <- as.data.frame((temp_df))
  
  #colnames(temp_df_t) <- temp_df$ATTRIBUTE_GROUPS
  
  ###colnames(temp_df_t) <- make.unique(colnames(temp_df_t), sep = "_")
  
  #temp_df_t <- temp_df_t[2:nrow(temp_df_t),]
  
  # Convert to numeric
  i = colnames(temp_df_t)
  temp_df_t[ , i] <- apply(temp_df_t[ , i], 2,  # Specify own function within apply
                           function(x) as.numeric(as.character(x)))
  
  #temp_df_t$Metabolite <- row.names(temp_df_t)
  return((temp_df_t))
}

##########################################################
normalize_by_od <- function(feature_table, metadata_table, samples_group_to_exclude = NULL, grouping_variable = NULL, select_by = "all", od_column = "od"){
  #' @title Normalize by OD
  #' @description
  #' Normalizes a feature table using OD measurements.
  #' 
  #' @param feature_table A dataframe containing samples in rows and features in columns.
  #' @param metadata_table A dataframe containing samples in columns and variables in rows. One row should be the OD values for all samples in feature_table (default = "OD").
  #' @param od_column String. Name of column containing OD values for all samples. By default is named "OD".
  #' @param select_by String. String indicating the way to select columns to normalize ("all", "range", "pattern").
  #' @param select_info Vector/string. If other than "all" selected as select_by, this is used to select columns to normalize.
  #' @return Dataframe where rows are samples and columns are features. Features are normalized by OD values in metadata_table
  
  #print(row.names(feature_table))
  #print(metadata_table$sample)
  #print(isTRUE(all.equal(row.names(feature_table), metadata_table$sample)))
  
  if (!isTRUE(all.equal(row.names(feature_table), row.names(metadata_table)))) {
    print("Sample names in feature table and metadatable are not identical")
    return()
  }else{
    print("Sample names in feature table and metadatable are identical :)")
  }
  
  #feature_table <- cbind(metadata_table[od_column], feature_table)
  
  #print(head(feature_table))
  
  if (!is.null(grouping_variable) && !is.null(samples_group_to_exclude)) {
    # To do, get the sample names that are part of the grouping variable
    samples_exclude <- rownames(feature_table)[metadata_table[[grouping_variable]] %in% samples_group_to_exclude]
    #print(samples_exclude)
    #filtered_df <- feature_table[!metadata_table$syncom %in% samples_group_to_exclude, ]
    #filtered_metadata <- metadata_table[!metadata_table$syncom %in% samples_group_to_exclude, ]
    filtered_df <- feature_table[!row.names(metadata_table) %in% samples_exclude, ]
    filtered_metadata <- metadata_table[!row.names(metadata_table) %in% samples_exclude, ]
    od_values <- as.vector(filtered_metadata$od)
  }else{
    od_values <- as.vector(metadata_table$od)
  }
  
  print(od_values)
  print(nrow(filtered_df) == length(od_values))
  
  #print("Select by all runing")
  df_normalized <- filtered_df / od_values # to divide based on names list or indices
  #return(select(df_normalized, c(2:length(df_normalized))))
  df_normalized <- rbind(df_normalized, feature_table[row.names(metadata_table) %in% samples_exclude, ])
  return(df_normalized)
}

#######################
log2_convert <- function(metabolites_table){
  ###metabolites_table_log2 <- metabolites_table %>% dplyr::select(-Metabolite) %>% log2()
  metabolites_table_log2 <- metabolites_table %>% log2()
  #metabolites_table_log2 <- metabolites_table %>% log2()
  #rownames(metabolites_table_log2) <- metabolites_table$Metabolite
  #metabolites_table_log2 <- tibble::add_column(metabolites_table_log2, "Metabolite" = metabolites_table$Metabolite, .before = 1)
  return(metabolites_table_log2)
}

#######################
normalize_table_to_treatment <- function(feature_table, metadata_table, grouping_variable = NULL, samples_group_to_norm = NULL){
  #metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite) %>% t()
  ###metabolites_table_norm <- metabolites_table %>% dplyr::select(-Metabolite)
  feature_table_norm <- as.data.frame(feature_table)
  
  
  samples_norm <- rownames(feature_table)[metadata_table[[grouping_variable]] %in% samples_group_to_norm]
  print(samples_norm)
  # Get the means of the treatment used to normalize the rest of treatments
  norm_treatment_means <- colMeans(feature_table[metadata_table[[grouping_variable]] %in% samples_group_to_norm, ])
  ###norm_treatment_means <- rowMeans(dplyr::select(metabolites_table, samples_norm))
  print(norm_treatment_means)
  # Get the sds of the treatment used to normalize the rest of treatments
  #norm_treatment_stds <- apply(dplyr::select(metabolites_table, starts_with(treatment_to_normilize)), MARGIN = 1, FUN = sd)
  ###norm_treatment_stds <- apply(dplyr::select(metabolites_table, samples_norm), MARGIN = 1, FUN = sd)
  #norm_treatment_stds <- apply(metabolites_table[metadata_table$syncom %in% samples_norm, ], MARGIN = 2, FUN = sd)
  #print(norm_treatment_stds)
  
  for (feature in 1:ncol(feature_table)) {
    #print(metabolite)
    for (cSample in 1:nrow(feature_table)) {
      #print(metabolites_table_norm[metabolite, cSample]*1000)
      #metabolites_table_norm[metabolite, cSample] <- (metabolites_table_norm[metabolite, cSample]-norm_treatment_means[metabolite])/norm_treatment_stds[metabolite]
      feature_table_norm[cSample, feature] <- (feature_table_norm[cSample, feature]-norm_treatment_means[feature])
    }
  }
  
  #print(head(metabolites_table_norm))
  #rownames(metabolites_table_norm) <- metabolites_table$Metabolite
  ###metabolites_table_norm <- tibble::add_column(metabolites_table_norm, "Metabolite" = metabolites_table$Metabolite, .before = 1)
  return(as.data.frame(feature_table_norm))
}

#######################
graph_metabolites <- function(feature_table, y1 = -10, y2= 10, dotsize = 0.5, binwidth = 0.5, ylab = "Log fold change"){
  # Data shape
  feature_table_t <- as.data.frame(t(feature_table))
  #print(head(feature_table_t))
  #print(colnames(feature_table))
  feature_table_t <- tibble::add_column(feature_table_t, "Metabolite" = colnames(feature_table), .before = 1)
  #feature_table_t$Metabolite <- colnames(feature_table)
  print(head(feature_table_t))
  
  metabolite_table_g <- tidyr::gather(data = feature_table_t, key = "Sample", value = "Change", colnames(feature_table_t[,2:ncol(feature_table_t)]))
  
  print(head(metabolite_table_g))
  
  metabolite_table_g$sample_type <- sapply(strsplit(metabolite_table_g$Sample, ".", fixed = TRUE), "[", 1)
  
  dotplot <- ggplot(data = metabolite_table_g, aes(x = Metabolite, y = Change, fill = sample_type))
  dotplot + geom_dotplot(binaxis = 'y', dotsize = dotsize, stackdir='center', position=position_dodge(0.8), binwidth = binwidth) +
    labs(y= ylab, x = "Metabolites") +
    ylim(y1, y2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
}

scale_0_1 <- function(df, scale_by = "variable") {
  if (!scale_by %in% c("variable", "sample")) {
    stop("scale_by must be either 'variable' or 'sample'")
  }
  
  # Scale by columns (variables)
  if (scale_by == "variable") {
    df_scaled <- as.data.frame(apply(df, 2, function(x) (x - min(x)) / (max(x) - min(x))))
  } 
  # Scale by rows (samples)
  else if (scale_by == "sample") {
    df_scaled <- as.data.frame(t(apply(df, 1, function(x) (x - min(x)) / (max(x) - min(x)))))
  }
  
  # Preserve row and column names
  rownames(df_scaled) <- rownames(df)
  colnames(df_scaled) <- colnames(df)
  
  return(df_scaled)
}

###################################################################################################

get_replicate_means <- function(replicate_data, id_col_name = NULL, id = NULL, replicates, subreplicates = 2, standard = TRUE){
  #' @title Calculate the means of treatments .
  #'
  #' @description Returns .
  #' @param replicate_data dataframe. .
  #' @value A DF .
  #' @details
  # Takes .
  
  # create means df, first col is compound names
  means <- replicate_data[1]
  # second col is id name if provided
  if (!is.null(id) && !is.null(id_col_name)) {
    means[id_col_name] <- c(id)
  }
  
  # getting the mean of the areas of the standard for these samples.
  if (standard) {
    std_means <- c()
    for(row in 1:nrow(replicate_data)){
      # Initializing sum and actual summed count
      temp_sum = 0
      rep_count = 0
      # For loop iterating the number of subreplicates
      for (replicate_count in 0:(subreplicates-1)) {
        # current replicate beggins in col 2 where std is. SUBSETTING
        current_replicate = replicate_data[[row,2+replicate_count]]
        # if na or null in data lets not sum these and not get their mean.
        if (!is.null(current_replicate) && !is.na(current_replicate)) {
          # adding current replicate to sum and increasing summed replicate count
          temp_sum = temp_sum + current_replicate
          rep_count = rep_count + 1
        }
      }
      # adding current meand to the final mean vector
      std_means <- c(std_means, temp_sum/rep_count)
    }
    # assigning means vector to means df
    means["std"] <- std_means
    # if std existed starting point for next step is increased by subreplicates
    starting_point = 2 + subreplicates
  }else{
    starting_point = 2
  }
  end_point = 1 + (replicates*subreplicates)
  # if std existed ending point for next step is increased by subreplicates
  if (standard) {
    end_point = end_point + subreplicates
  }
  
  # getting the means for each reaplicate for these samples.
  for(col in seq(starting_point, end_point, by = subreplicates)){
    current_means <- c()
    for(row in 1:nrow(replicate_data)){
      # Initiliazing sum and actual summed count
      temp_sum = 0
      rep_count = 0
      # For loop iterating the number of subreplicates
      for (replicate_count in 0:(subreplicates-1)) {
        # take current replicate where col. SUBSETTING.
        current_replicate = replicate_data[[row,col + replicate_count]]
        if (!is.null(current_replicate) && !is.na(current_replicate)) {
          # adding current replicate to sum and increasing summed replicate count
          temp_sum = temp_sum + current_replicate
          rep_count = rep_count + 1
        }
      }
      # adding current mean to the final mean vector
      current_means <- c(current_means, temp_sum/rep_count)
    }
    # assigning means vector to means df
    means[stringr::str_sub(colnames(replicate_data[col]), 1, 3)] <- current_means
  }
  return(means)
}

###################################################################################################
stability_per_period <- function(replicate_data, period_length, time_unit = "day") {
  #' Stability per period
  #'
  #' Calculates the mean stability index mean of a functional parameter of replicate data for a period of time.
  #'
  #' @param replicate_data 
  #' @param period_length
  #' @param time_unit
  #'
  #' @return Dataframe where rows are replicates and columns are "time_unit"(s) (days, etc). Each entry corresponds to the stability index of the period including the current time and previous "period_length" time units.
  #' @export
  #'
  #' @examples

  # Creating empty dataframe where stability measures for each replicate will be stored.
  # Rows are replicates.
  # Columns are "time_unit"(s) (days, etc).
  # Each entry corresponds to the stability index of the period including the current time and previous "period_length" time units.
  results_df <- data.frame(row.names = colnames(replicate_data_biogas))
  # Iterating over each time in which functional parameter was measured.
  for (day in seq(from = period_length, to = nrow(replicate_data), by= 1)) {
    # Empty vector that becomes the column containing stability indices for each replicate for the current time.
    temp_col <- c()
    # starting day is the first day included in stability calculation
    starting_day <- day-period_length+1
    #print(starting_day)
    # starting day is the last day included in stability calculation
    end_day <- day
    #print(end_day)
    # Subsetting the data included in the stability indices calculation.
    temp_data <- replicate_data[starting_day:end_day, ]
    #print(temp_data)
    # Mean of the functional parameter for each replicate for the current period.
    means <- sapply(temp_data, mean)
    # standard deviation of the functional parameter for each replicate for the current period.
    std_dev <- sapply(temp_data, sd)
    # Stability indices of the functional parameter for each replicate for the current period.
    period_stability <- 1 -(std_dev/means)
    #print(period_stability)
    # append the stability data for current period
    results_df[paste(paste(time_unit, " "), day)] <- period_stability
  }
  return(results_df)
}

feature_table_heatmap <- function(ft1, ft2){
  ft1_t <- t(ft1)
  ft1_t <- ft1_t[order(row.names(ft1_t)), ] # Ordering by row names
  
  
  ft2_t <- t(ft2)
  ft2_t <- ft2_t[order(row.names(ft2_t)), ] # Ordering by row names
  
  
  ft1Xft2 <- cor(ft1_t, ft2_t, method = "pearson")
  #heatmap(ft1Xft2)
  corrplot::corrplot(ft1Xft2, , tl.cex=0.5)
}


feature_table_heatmap_w_sig <- function(ft1, ft2){
  ft1_t <- t(ft1)
  ft1_t <- ft1_t[order(row.names(ft1_t)), ] # Ordering by row names
  
  
  ft2_t <- t(ft2)
  ft2_t <- ft2_t[order(row.names(ft2_t)), ] # Ordering by row names
  
  
  ft1Xft2 <- cor(ft1_t, ft2_t, method = "spearman")
  
  ft1Xft2_df <- as.data.frame(ft1Xft2)
  
  ft1Xft2_df <- cbind(rownames(ft1Xft2_df), data.frame(ft1Xft2_df, row.names=NULL))
  colnames(ft1Xft2_df)[1] <- "species"
  
  # Gathering data by species
  adXm_g <- ft1Xft2_df %>%
    pivot_longer(cols = where(is.numeric), names_to = "compound", values_to = "correlation")
  
  #print(adXm_g)
  
  adxm.pval <- as.data.frame(Hmisc::rcorr(ft1_t, ft2_t, type = "spearman")$P)
  
  if (identical(ft1, ft2)) {
    adxm.pval <- adxm.pval[1:(nrow(adxm.pval)/2), 1:(ncol(adxm.pval)/2)]
    # Converting rownames to column 1
    adxm.pval <- cbind(rownames(ft1), data.frame(adxm.pval, row.names=NULL))
  }else{
    adxm.pval <- adxm.pval[1:(nrow(ft1)), (nrow(ft1)+1):ncol(adxm.pval)]
    # Converting rownames to column 1
    adxm.pval <- cbind(rownames(adxm.pval), data.frame(adxm.pval, row.names=NULL))
  }
  colnames(adxm.pval)[1] <- "species"

  # Gathering data by species
  adXm_pvs_g <- adxm.pval %>%
    pivot_longer(cols = where(is.numeric), names_to = "compound", values_to = "p_value")
  
  # store pvalues
  pvals <- adXm_pvs_g$p_value
  
  # Remove NAs
  pvals[is.na(pvals)] <- 1
  
  # Create column of significance labels
  adXm_g["pvals"] <- pvals
  
  # Create column of significance labels
  adXm_g["stars"] <- cut(pvals, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  # Plotting
  ggplot(aes(x=compound, y=species, fill=correlation), data=adXm_g) +
    geom_tile() +
    scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") +
    geom_text(aes(label=stars), color="black", size=2) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
  
}
