
# Funtion to remove features by name
remove_feature_by_prefix <- function(df, patterns) {
  # Create a single regex pattern that matches any of the species names at the start
  combined_pattern <- paste0("^(", paste(patterns, collapse = "|"), ")")
  
  # Filter the dataframe: keep rows that do NOT match the pattern
  df_filtered <- df[!grepl(combined_pattern, rownames(df)), ]
  
  return(df_filtered)
}

##### Table Filtering
filter_features_by_col_counts <- function(feature_table, min_count, col_number){
  if (ncol(feature_table) > 1) {
    return(feature_table[which(rowSums(feature_table >= min_count) >= col_number), ])
  }
  else if(ncol(feature_table) == 1){
    ft <- feature_table[feature_table >= min_count, ,drop=FALSE]
    return(ft)
  }
  else{
    print("Dataframe has no columns")
  }
}

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

# Orders the samples of a feature table using ward clustering
order_samples_by_clustering <- function(feature_table, diss_method = "euclidean", clust_method = "ward.D2"){
  # Takes feature_table and returns the list of samples ordered according to the clustering algorithm
  df_otu <- feature_table %>% rownames_to_column(var = "Species")
  
  df_t <- as.matrix(t(df_otu[, -1]))  # Exclude the "Species" column after moving it to row names
  
  # Perform hierarchical clustering
  d <- dist(df_t, method = diss_method)
  hc <- hclust(d, method = clust_method)
  
  # Get the order of samples based on clustering
  ordered_samples_cluster <- colnames(df_otu)[-1][hc$order]  # Remove "Species" again
  
  return(ordered_samples_cluster)
}

##### Function to convert a OTU table to a strain-level-table
# It takes the otu table at species level and a second dataframe including strain-level data.
# First dataframe should be a dataframe containing species level data.
# Second dataframe should be a dataframe containing
# the strain inoculation data in the following format:
merge_abundance_by_strain <- function(df1, df2) {
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)
  # Extract species names from df1
  species_names_df1 <- rownames(df1)
  #print(species_names_df1)
  
  # Extract species names and strain names from df2
  strain_names_df2 <- df2[, 1]  # Full strain names (including strain number)
  #print(strain_names)
  
  # Create an empty matrix to store the new abundance data
  new_abundance_matrix <- matrix(0, nrow = nrow(df2), ncol = ncol(df1))
  rownames(new_abundance_matrix) <- strain_names_df2
  colnames(new_abundance_matrix) <- colnames(df1)
  
  #print(head(new_abundance_matrix))
  #print(nrow(new_abundance_matrix))
  #print(ncol(new_abundance_matrix))
  
  samples <- colnames(new_abundance_matrix)
  # Iterate over each sample of the new DF
  for (i in seq_along(samples)) {
    #print(samples[i])
    # Get the SC to which it belongs to.
    current_sc <- strsplit(x = samples[i], split = "_")[[1]][1]
    print(current_sc)
    # get the list of strains inoculated in that sample.
    inoc_strains_per_sample <- get_inoculated_strains(df2 = df2, sample_name = current_sc)
    print(inoc_strains_per_sample)
    
    for (x in seq_along(inoc_strains_per_sample)) {
      strain_name <- inoc_strains_per_sample[x]
      print(strain_name)
      # get the index where the data is going to be inserted. The index is the same row as in the df2
      index_strain_df2 <- which(strain_names_df2 == strain_name) # this is also the same in the new df
      #print(index_strain_df2)
      # get the name of the species.
      species_name <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", strain_name)  # Remove strain number, keeping species
      print(species_name)
      if (species_name %in% species_names_df1) {
        index_species_df1 <- which(species_names_df1 == species_name)
        #print(index_species_df1)
        #print(species_names_df1[index_species_df1])
        # get the actual data, that corresponds to the species in df1
        current_abundance <- df1[index_species_df1, i]
        #print(current_abundance)
        # paste the data
        new_abundance_matrix[index_strain_df2, i] <- current_abundance
      }
    }
  }
  return(as.data.frame(new_abundance_matrix))
}

get_inoculated_strains <- function(df2, sample_name) {
  # Select the column corresponding to the sample
  sample_column <- df2[[sample_name]]
  
  # Get row indices where the value is 1 (inoculated strains)
  inoculated_indices <- which(sample_column == 1)
  
  # Extract the strain names based on the indices
  inoculated_strains <- df2[inoculated_indices, 1]  # First column contains strain names
  
  return(inoculated_strains)
}

# This function takes a dataframe where the rownames are strain level OTUs/ASVs in the form:
# Genera species strain data. The two first words are used a the Species names that are numbered then as:
# Genera species 1; Genera species 2; Genera species 3
strain_name2strain_number <- function(df){
  # Extract only the "Genus species" part
  species_names <- sub(" \\S+$", "", rownames(df))  
  
  # Create a numeric ID for each strain within the same species
  species_ids <- ave(species_names, species_names, FUN = function(x) seq_along(x))
  
  # Create new rownames with species + strain ID
  new_rownames <- paste(species_names, species_ids)
  
  # Assign new rownames to the dataframe
  rownames(df) <- new_rownames
  
  # Print the updated dataframe
  #print(df)
  return(df)
}

# Sort otu table in barcodes numeration
sort_nanopore_table_by_barcodes <- function(df, new_names = NULL){
  cn <- colnames(df) # store column names
  sorted_names <- cn[order(nchar(cn), cn)] # order columns names
  df_sorted <- df[, sorted_names] # order data frame using colnames
  if (!is.null(new_names) && ncol(df_sorted) == length(new_names)) {
    colnames(df_sorted) <- new_names
  }
  return(df_sorted)
}

# Function to set selected species/sample combinations to zero
zero_out_species_in_samples <- function(df, species_name, sample_names) {
  # Safety check: does the species exist?
  if (!(species_name %in% rownames(df))) {
    stop(paste("Species", species_name, "not found in rownames"))
  }
  
  # Safety check: do all samples exist?
  if (!all(sample_names %in% colnames(df))) {
    missing_samples <- sample_names[!sample_names %in% colnames(df)]
    stop(paste("Samples not found in dataframe:", paste(missing_samples, collapse = ", ")))
  }
  
  # Set the selected cells to zero
  df[species_name, sample_names] <- 0
  
  return(df)
}

# Read Hitchhikers guide style export feature table
read_ft <- function(path, sort_by_names = FALSE, p_sep = ","){
  ft <- read.csv2(path, header = TRUE, row.names = 1, sep = p_sep, dec = ".") #read csv table
  if(isTRUE(sort_by_names)){
    ft <- ft[order(row.names(ft)), ] # sort my row names (sample names)
  }
  
  rownames(ft) <- gsub("\\.mzML$", "", rownames(ft))
  #col_names <- colnames(ft)
  #ft <- sapply(ft, as.numeric)
  #ft <- as.data.frame(ft)
  #colnames(ft) <- col_names
  return(t(ft))
}

read_metadata <- function(path, sort_table = FALSE){
  md <- read.csv(path, row.names = 1)
  if(isTRUE(sort_table)){
    md <- md[order(row.names(md)), ] # sort my row names (sample names)
  }
  return(md)
}
