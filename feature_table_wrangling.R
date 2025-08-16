###################################################################################################
##### Install packages
if (!"collections" %in% installed.packages()) install.packages("collections")
if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"tidyr" %in% installed.packages()) install.packages("tidyr")
if (!"tibble" %in% installed.packages()) install.packages("tibble")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!"phyloseq" %in% installed.packages()) BiocManager::install("phyloseq", update = FALSE)

###################################################################################################

# This function takes a "tax_rank" string that correspond to a taxonomic rank in Greengenes format.
# Returns a list of strings which represent the columns in the tax_table of a biom file 
# that have to be joined to get the taxonomy assignment of each AVS/OTU as a string.
# If a not valid tax_rank is provided it returns an error.
get_colNames_per_rank <- function(tax_rank){
  colNames = NULL
  switch(tax_rank,
         Strain = {
           colNames = c("Genus", "Species", "Strain")
         },
         Species = {
           # Species level
           colNames = c("Genus", "Species")
         },
         Genus = {
           # Genus level
           colNames = c("Genus")
         },
         Family = {
           # Family
           colNames = c("Family")
         },
         Order = {
           # Order
           colNames = c("Order")}
  )
  if (!is.null(colNames)){
    return(colNames)
  }else{
    stop("Please choose a valid taxonomy rank!", call. = FALSE)
  }
}

###################################################################################################

# This function takes a character vector containing the result of splitting a taxonomy vector in the greenegenes format.
# It returns a named vector where each field is a taxonomic rank for the passed taxonomy entry.
# The taxonomic ranks are the same as in the greengenes taxonomy format but include a "Strain" rank.
# This function is used by phyloseq's "import_biom" function to parse taxonomy.
# import_biom splits taxonomy vectors automatically when they are the in the greengenes format.
parse_taxonomy_strain <- function(char.vec){
  # Remove the greengenes taxonomy rank id prefix.
  named.char.vec <- substring(char.vec, first = 4)
  # Set the names for each rank.
  names(named.char.vec) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  return(named.char.vec)
}

###################################################################################################

# This function takes a biom object and extracts it's "tax_table" and the "otu table".
# Then cbinds both dataframes to obtain a dataframe where the 1st column is the taxonomy and the 
extract_table <- function(biom_object, tax_rank, col_names){
  # Agglomerate tax_table by the chosen tax_rank
  biom_object <- phyloseq::tax_glom(biom_object, taxrank = tax_rank, NArm=TRUE)
  # cbind tax_table and otu_table
  feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(phyloseq::tax_table(biom_object)),
                                                    taxonomy,
                                                    all_of(col_names),
                                                    sep = "_"), "taxonomy"),
                         data.frame(phyloseq::otu_table(biom_object)))
}

###################################################################################################

# This function takes a feature table. First it make all the values in the "taxonomy" column unique.
# Then it makes the "taxonomy" column the rownames of the table.
# If "order_table" is TRUE it orders the table by ASVs/OTUs abundance.
clean_table <- function(feature_table, order_table){
  # Get valid (unique) names for all ASVs/OTUs.
  feature_table["taxonomy"] <- make.unique(feature_table$taxonomy, sep = "_")
  # Set taxonomy column as rownames
  feature_table <- tibble::column_to_rownames(tibble::remove_rownames(feature_table), var = "taxonomy")
  if (order_table) {
    # Order by abundances mean, from higher to lower.
    feature_table <- feature_table[order(rowMeans(feature_table), decreasing = TRUE),]
  }else{
    return(feature_table)
  }
}

# This function takes a "biom_path" to a biom file with an otu_table and a tax_table,
# a string "tax_rank" which indicates the level of analyses, and bool "order_table".
# tax_rank parameter must be a value of greengenes ranks format; if not an error is returned.
# The ASVs/OTUs in the biom file are agglomerated by the "tax_rank" provided
# "order_table" indicates if the table should be ordered by larger to smaller values of rowMeans.
# Generally, ASV/OTU tables from QIIME2 are already ordered by row sums.
# This function returns a dataframe where rows are the ASVs and the columns are samples,
# "rownames" are ASVs taxonomy at the selected rank, and "colnames" are samples names.
# Taxonomy is dereplicated so that no row has the same name (which is not allowed in R dataframes).
# The output format is useful for using in other packages like vegan and to generate plots like barplots and heatmaps.
load_biom_as_table <- function(biom_path, tax_rank = "Species", strain_taxonomy = FALSE, order_table = FALSE){
  
  unite_colNames <- get_colNames_per_rank(tax_rank)
  
  if(strain_taxonomy) {
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_strain)
  }else{
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=phyloseq::parse_taxonomy_greengenes)
  }
  
  extracted_feature_table <- extract_table(biom_object, tax_rank, unite_colNames)
  
  return(clean_table(extracted_feature_table, order_table = order_table))
}

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

flag_samples_by_abundance <- function(feature_df, metadata_df, feature_name, percentage_threshold) {
  # 1. Check if feature_df and metadata_df have the same samples in the same order
  if (identical(colnames(feature_df), rownames(metadata_df))) {
    print("Both dataframes have the same sample names")
  }else{
    print("Both dataframes have different sample names")
    return()
  }
  
  # 2. Check that the feature_name exists in feature_df
  if (!feature_name %in% rownames(feature_df)) {
    stop("The feature name is not found in the feature dataframe.")
  }
  
  # 3. Get the vector of abundances for that feature
  feature_abundances <- feature_df[feature_name, ]
  
  # 4. Calculate percentage per sample
  # First calculate total abundance per sample
  total_abundance <- colSums(feature_df)
  
  percentage_abundance <- (feature_abundances / total_abundance) * 100
  
  # 5. Create logical vector indicating whether each sample is above threshold
  above_threshold <- percentage_abundance >= percentage_threshold
  
  # 6. Make a copy of the metadata and add the new column
  metadata_copy <- metadata_df
  new_colname <- paste0("ATTRIBUTE_above_", percentage_threshold, "pct_", gsub(" ", "_", feature_name))
  metadata_copy[[new_colname]] <- above_threshold[1,]
  
  # Return modified metadata
  return(metadata_copy)
}

# Takes an feature table (OTUs) and removes the strain information from the species NOT in the passed vector of species
merge_non_target_strains <- function(df, target_species) {
  # Extract species names (first two words) from rownames
  species_names <- sapply(strsplit(rownames(df), " "), function(x) paste(x[1:2], collapse = " "))
  #print(species_names)
  # Identify which rows belong to target or non-target species
  is_target <- species_names %in% target_species
  #print(is_target)
  # Separate target and non-target
  target_df <- df[is_target, , drop = FALSE]
  #print(target_df)
  non_target_df <- df[!is_target, , drop = FALSE]
  #print(non_target_df)
  non_target_species <- species_names[!is_target]
  #print(non_target_species)
  # 
  # # Aggregate non-target strains by species
  if (nrow(non_target_df) > 0) {
    aggregated <- aggregate(non_target_df, by = list(Species = non_target_species), FUN = sum)
    # Set the species name as rownames and remove the Group column
    # >>> THIS IS WHERE YOU ADD " 1" TO THE SPECIES NAMES <<<
    rownames(aggregated) <- paste(aggregated$Species, "1")
    #rownames(aggregated) <- aggregated$Species
    aggregated$Species <- NULL
  } else {
    aggregated <- NULL
  }
  print(aggregated)
  # Combine target and aggregated non-target dataframes
  result <- rbind(target_df, aggregated)
  
  return(result)
}

filter_low_abundance <- function(rel_abundance, threshold = 0.01) {
  # rel_abundance: species x samples matrix/dataframe of relative abundances
  # threshold: minimum relative abundance (e.g., 0.01 = 1%)
  
  # Compute maximum abundance of each species across samples
  species_max <- apply(rel_abundance, 1, max)
  
  # Keep only species with max abundance >= threshold
  filtered_df <- rel_abundance[species_max >= threshold, ]
  
  message(paste("Filtered from", nrow(rel_abundance), 
                "to", nrow(filtered_df), "species"))
  
  return(filtered_df)
}
