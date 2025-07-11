# Sort otu table in barcodes numeration
sort_nanopore_table_by_barcodes <- function(df, new_names = NULL){
  cn <- colnames(df) # store column names
  sorted_names <- cn[order(nchar(cn), cn)] # order columns names
  df_sorted <- df[, sorted_names] # order data frame using colnames
  if (!is.null(new_names) && ncol(df_sorted == length(new_names))) {
    colnames(df_sorted) <- new_names
  }
  return(df_sorted)
}

###################################################################################################

filter_otus_by_counts_nas <- function(otu_table, min_count, percentage){
    return(otu_table[which(rowMeans(! is.na(otu_table)) >= percentage), ])
}

filter_otus_by_counts_col_percent <- function(otu_table, min_count, percentage){
    return(otu_table[which(rowMeans(otu_table >= min_count) >= percentage), ])
}

filter_otus_by_counts_col_counts <- function(otu_table, min_count, col_number){
  if (ncol(otu_table) > 1) {
    return(otu_table[which(rowSums(otu_table >= min_count) >= col_number), ])
  }else{
    return(otu_table)
  }
}

filter_species_by_col_counts <- function(otu_table, min_count, col_number) {
  if (ncol(otu_table) == 1) { # 1 column for species and 1 for data
    # Filter rows directly for the single sample column
    filtered_df <- otu_table[otu_table[1] >= min_count, ,drop = FALSE]
  } else {
    # Check rows that meet the criteria: at least `m` columns with values >= `n
    filtered_df <- otu_table[rowSums(otu_table >= min_count) >= col_number, ,drop = FALSE]
  }
  return(filtered_df)
}

sub_otutable <- function(otu_table, sample_indices, sample_names){
  # Select samples defined by user
  sub_otut <- dplyr::select(otu_table, any_of(sample_indices))
  # Remove rows that have only 0s
  sub_otut <- sub_otut[rowSums(sub_otut != 0) > 0, ]
  # if user provided sample names, apply them
  if (!is.null(sample_names) && length(sample_indices) == length(sample_names)) {
    colnames(sub_otut) <- sample_names
  }
  return(sub_otut)
}

# To Do
join_otu_tables <- function(otu_table, otu_table2){
  
}

prep_otutable <- function(otu_table, rename_cols, set_row_names, col_pattern){
  otu_table <- as.data.frame(otu_table)
  if (rename_cols) {
    colnames(otu_table) <- paste(colnames(otu_table), col_pattern, sep = "_")
  }
  if (set_row_names) {
    row.names(otu_table) <- otu_table[,1]
    otu_table <- otu_table[2:ncol(otu_table_v1)]
  }
  return(otu_table)
}

### function for renaming samples (columns) in a table with a prefix given by user + the current sample name
rename_sample_names_prefix <- function(feature_tables, prefixes){
  new_tables <- list()
  #print(feature_tables)
  for (table in seq(from = 1, to = length(feature_tables), by=1)) {
    rownames_saved <- rownames(feature_tables[[table]])
    print(rownames_saved)
    #print(colnames(feature_tables[[table]]))
    colnames(feature_tables[[table]]) <- paste0(prefixes[table], colnames(feature_tables[[table]]))
    rownames(feature_tables[[table]]) <- rownames_saved
    print(head(feature_tables[[table]]))
    new_tables <- append(new_tables, list(feature_tables[[table]]), 0)
    
  }
  return(new_tables)
}

# to do replace NAs with 0s
# code on sequencing bar plots nasal syncom script:
#syncomall2 <- syncomall %>% replace(is.na(.), 0)

transform_feature_table <- function(feature_table, transform_method){
  if (transform_method == "zscale") {
    # Z-Scaling
    df_transformed <- as.data.frame(scale(feature_table))
  } else if (transform_method == "min_max"){
    df_transformed <- feature_table
    normalize = function(x) (x- min(x))/(max(x) - min(x))
    cols <- sapply(df_transformed, is.numeric)
    df_transformed[cols] <- lapply(df_transformed[cols], normalize)
  }else if (transform_method == "rel_abundance"){
    # Relative abundance
    df_transformed <- sweep(feature_table, 2, colSums(feature_table), FUN = "/")
  } else{
    "Transform method not valid"
  }
  return(df_transformed)
}

cumulative_half_addition <- function(df) {
  # Get numeric part of the dataframe (excluding species names)
  abundance_matrix <- as.matrix(df[, -1])
  
  df2 <- as.matrix(df[, -1])
  
  # Apply calculations for each column (sample)
  for (col in 1:ncol(abundance_matrix)) {
    for (row in nrow(abundance_matrix):1) { # Iterate from bottom to top
      if (row == nrow(abundance_matrix)) {
        # Last row: itself divided by 2
        if (df2[row, col] == 0) {
          abundance_matrix[row, col] <- 0
        } else{
          abundance_matrix[row, col] <- df2[row, col] / 2
        }
      } else { # Other rows: cumulative sum of rows below + itself/2
        if (df2[row, col] == 0) {
          abundance_matrix[row, col] <- 0
        } else {
          abundance_matrix[row, col] <- sum(df2[(row + 1):nrow(df2), col]) +
            (df2[row, col] / 2)
        }
      }
    }
  }
  
  # Replace values back into the dataframe
  df[, -1] <- abundance_matrix
  return(df)
}

calculate_relative_abundance <- function(df) {
  species <- rownames(df)
  # Calculate relative abundance
  relative_abundance <- sweep(df, 2, colSums(df), "/")
  # Combine species names back with the relative abundance data
  rownames(relative_abundance) <- species
  # Return the result as a dataframe
  return(as.data.frame(relative_abundance))
}

order_samples_by_clustering <- function(feature_table){
  # Takes feature_table and returns the list of samples ordered according to the clustering algorithm
  df_otu <- feature_table %>% rownames_to_column(var = "Species")
  
  df_t <- as.matrix(t(df_otu[, -1]))  # Exclude the "Species" column after moving it to row names
  
  # Perform hierarchical clustering
  d <- dist(df_t, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  # Get the order of samples based on clustering
  ordered_samples_cluster <- colnames(df_otu)[-1][hc$order]  # Remove "Species" again
  
  return(ordered_samples_cluster)
}
