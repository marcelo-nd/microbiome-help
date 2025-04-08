
get_inoculated_strains <- function(df2, sample_name) {
  # Select the column corresponding to the sample
  sample_column <- df2[[sample_name]]
  
  # Get row indices where the value is 1 (inoculated strains)
  inoculated_indices <- which(sample_column == 1)
  
  # Extract the strain names based on the indices
  inoculated_strains <- df2[inoculated_indices, 1]  # First column contains strain names
  
  return(inoculated_strains)
}

##### Function to convert a OTU table to a strain-level-table
# It takes the otu table at species level and a second dataframe including strain-level data.
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
    #print(current_sc)
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
