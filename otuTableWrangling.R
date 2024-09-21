###################################################################################################
    
filter_otus_by_counts_nas <- function(otu_table, min_count, percentage){
    return(otu_table[which(rowMeans(! is.na(otu_table)) >= percentage), ])
}
    
###################################################################################################

filter_otus_by_counts_col_percent <- function(otu_table, min_count, percentage){
    return(otu_table[which(rowMeans(otu_table >= min_count) >= percentage), ])
}
    
###################################################################################################

filter_otus_by_counts_col_counts <- function(otu_table, min_count, col_number){
  if (ncol(otu_table) > 1) {
    return(otu_table[which(rowSums(otu_table >= min_count) >= col_number), ])
  }else{
    return(otu_table)
  }
}

# to do replace NAs with 0s
# code on sequencing bar plots nasal syncom script

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
