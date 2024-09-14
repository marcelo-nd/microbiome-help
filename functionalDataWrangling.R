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

get_palette <- function(nColors = 50){
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

###################################################################################################

get_replicate_means <- function(replicate_data, id_col_name = NULL, id = NULL, replicates, subreplicates = 2, standard = TRUE){
  # Takes the mean of "subreplicates" number of observations for "replicates" number of replicates
  # Receives optional ID name (id_col_name), ID for group (id), replicates number and sureplicates to teake the mean.
  # Returns df
  
  # Installing required packages
  if (!"stringr" %in% installed.packages()) install.packages("stringr")
  
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
stability_per_period <- function(replicate_data, period_length, time_unit = "day") {
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

# Read and prepare FIA pos/neg table
read_fia_table <- function(table_path, sheet = "pos"){
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
  colnames(fia_df_t) <- make.names(colnames(fia_df_t), unique=TRUE)
  return(fia_df_t)
}

##########################################################
normalize_by_od <- function(feature_table, metadata_table, od_column = "OD", select_by = "all", select_info = NULL){
  if (!all.equal(rownames(feature_table),metadata_table$Sample)) {
    print("Sample names in feature table and metadatable are not identical")
    return()
  }else{
    print("Sample names in feature table and metadatable are identical :)")
  }
  
  feature_table <- cbind(metadata_table[od_column], feature_table)

  if (select_by == "all") {
    df_normalized <- feature_table %>%
      mutate(across(where(is.numeric), ~ . / OD)) # to divide based on names list or indices
    return(select(df_normalized, c(2:length(df_normalized))))
  }
  else if (select_by == "range") {
    #print(c((select_info[1]+1):(select_info[2]+1)))
    #test_df <- select(dataframe, c((select_info[1]+1):(select_info[2])))
    #print(select(dataframe, c((select_info[1]+1):(select_info[2]+1))))
    df_normalized <- feature_table %>%
      mutate(across(all_of(c((select_info[1]+1):(select_info[2]))), ~ . / OD)) # to divide based on names list or indices
      #mutate(across(starts_with("var"), ~ . / div_value)) %>% # to divide variables by pattern in name
      #select(-div_value)  # Remove the division value column if not needed
      # Return the normalized dataframe
    return(select(df_normalized, c(2:length(df_normalized))))
  }
  else if (select_by == "numeric") {
    df_normalized <- feature_table %>%
      mutate(across(where(is.numeric), ~ . / div_value)) %>% # to divide all numeric variables
      # Return the normalized dataframe
      return(df_normalized)
  }
  else if (select_by == "pattern") {
    df_normalized <- feature_table %>%
      mutate(across(starts_with(select_info), ~ . / div_value)) %>% # to divide variables by pattern in name e.g. "var"
      # Return the normalized dataframe
      return(df_normalized)
  }
  else{
    print("select_by value not valid")
  }
}

# Remove highly variable metabolites
filter_by_error <- function(feature_table, metadata_table, grouping_var = NULL, error_threshold = 50){
  # Step 1: Join grouping column in metadata to feature table.
  if (!all.equal(rownames(feature_table),metadata_table$Sample)) {
    print("Sample names in feature table and metadatable are not identical")
    return()
  }else{
    print("Sample names in feature table and metadatable are identical :)")
  }

  feature_table <- cbind(metadata_table[grouping_var], feature_table)
  
  # Step 2: Calculate the error for each variable per type
  errors <- feature_table %>%
    group_by(SynCom) %>%
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

# Do PCA plot
fia_pca <- function(feature_table, metadata_table, grouping_col, encircle = FALSE){
  # transposing feature table
  ft_t <- t(feature_table)
  #print(sapply(ft_t, class))
  # Step 1: Join grouping column in metadata to feature table.
  if (!all.equal(colnames(ft_t),metadata_table$Sample)) {
    print("Sample names in feature table and metadatable are not identical")
    return()
  }else{
    print("Sample names in feature table and metadatable are identical :)")
  }
  metadata_table2 <- as.data.frame(metadata_table)
  rownames(metadata_table2) <- metadata_table2$Sample
  #print(sapply(metadata_table2, class))
  
  fia_pca <- PCAtools::pca(ft_t, scale = TRUE, metadata =metadata_table2, transposed = FALSE)

  p2 <- PCAtools::biplot(fia_pca, showLoadings = TRUE, ntopLoadings = 2, lab = NULL, colby = grouping_col,
                         legendPosition = "right", axisLabSize = 8, legendLabSize = 8, legendIconSize = 2, pointSize = 1.5,
                         colkey = get_palette(nColors = 60), encircle = encircle)
  
  #p2 <- PCAtools::biplot(fia_pca, showLoadings = TRUE, ntopLoadings=2)
  p2
}


