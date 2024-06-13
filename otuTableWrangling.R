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

#to do replace NAs with 0s
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
