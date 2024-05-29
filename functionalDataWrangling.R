

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

get_qnty_from__std <- function(raw_qnts, id_col = NULL, replicates, std_values){
  # Calcula el valor de una cantidad con base en un std.
  
  # create temp df, first col is compound names
  convtd_qnts <- raw_qnts[1]
  
  # if id_col was provided transfer it to convtd_qnts DF.
  if (!is.null(id_col)) {
    convtd_qnts[colnames(raw_qnts[id_col])] <- raw_qnts[id_col]
    # set starting point accordingly if id_col exists in raw_qnts DF
    starting_point <- 4
  }else{
    starting_point <-  3
  }
  
  ending_point <- starting_point + (replicates - 1)
  
  # Calculate quantities from stds
  # for each replicate from staring point to ending point
  for(col in starting_point:ending_point){
    # create temporary column
    temp_col <- c()
    # for all compound rows in DF
    for(row in 1:nrow(raw_qnts)){
      # append value for each calculated compound.
      # Calculation is a rule of three between compound equivalencies,  standard values. SUBSETTING.
      temp_col <- c(temp_col, c(std_values[[row, 2]] * raw_qnts[[row, col]] / raw_qnts[[row, 3]]))
    }
    # Adding temporary column to cnvt_qnts DF
    convtd_qnts[colnames(raw_qnts[col])] <- temp_col
  }
  return(convtd_qnts)
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

###################################################################################################

# Normalize by OD of SynComs

normalize_by_od <- function(dataframe, select_by = "range", select_info){
  if (select_by == "range" || select_by == "names") {
    df_normalized <- dataframe %>%
      mutate(across(all_of(select_info), ~ . / OD)) %>% # to divide based on names list or indices
      #mutate(across(starts_with("var"), ~ . / div_value)) %>% # to divide variables by pattern in name
      #select(-div_value)  # Remove the division value column if not needed
      # Return the normalized dataframe
      return(df_normalized)
  }
  else if (select_by == "numeric") {
    df_normalized <- dataframe %>%
      mutate(across(where(is.numeric), ~ . / div_value)) %>% # to divide all numeric variables
      # Return the normalized dataframe
      return(df_normalized)
  }
  else if (select_by == "pattern") {
    df_normalized <- dataframe %>%
      mutate(across(starts_with(select_info), ~ . / div_value)) %>% # to divide variables by pattern in name e.g. "var"
      # Return the normalized dataframe
      return(df_normalized)
  }
  else{
    print("select_by value not valid")
  }
}

###################################################################################################
# Remove highly variable metabolites

filter_by_error <- function(dataframe, error_threshold = 50){
  
  # Step 1: Calculate the error for each variable per type
  errors <- dataframe %>%
    group_by(SynCom) %>%
    summarise(across(where(is.numeric), ~ (sd(.) / mean(.)) * 100, .names = "error_{col}"))
  
  # View the errors dataframe
  #print(head(errors))
  
  # Step 2: Average the error for each variable for all types
  avg_errors <- errors %>%
    summarise(across(starts_with("error_"), mean, na.rm = TRUE))
  
  # View the average errors dataframe
  #print(head(avg_errors))
  
  variables_to_keep <- names(avg_errors)[avg_errors <= error_threshold]
  variables_to_keep <- gsub("error_", "", variables_to_keep)
  
  # Keep the metadata columns and the variables with error below the threshold
  df_filtered <- dataframe %>%
    select(all_of(c("Sample", "SynCom", "Time", "OD", variables_to_keep)))
  
  # View the filtered dataframe
  return(df_filtered)
}