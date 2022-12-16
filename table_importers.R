###################################################################################################
##### Install packages
if (!"collections" %in% installed.packages()) install.packages("collections")
if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"tidyr" %in% installed.packages()) install.packages("tidyr")
if (!"tibble" %in% installed.packages()) install.packages("tibble")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!"phyloseq" %in% installed.packages()) BiocManager::install("phyloseq", update = FALSE)

###################################################################################################

# This function takes a character vector containing the result of splitting a taxonomy vector.
# It returns a named vector where each field is a taxonomic rank for the passed taxonomy entry.
# The taxonomic ranks are the same as in the greengenes taxonomy format but include a "Strain" rank.
# This function is used by phyloseq's "import_biom" function to parse taxonomy
# import_biom splits taxonomy vectors atuomatically when in greengenes' format.
parse_taxonomy_strain <- function(char.vec){
  # Remove the greengenes taxonomy rank id prefix.
  named.char.vec <- substring(char.vec, first = 4)
  # Set the names for each rank.
  names(named.char.vec) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  return(named.char.vec)
}

###################################################################################################
# This function takes a vector of strings. Each string is the taxonomic id for a feature (OTU/ASV) in an ASV/OTU table.
# Dereplication in this context means to make every string in the vector unique.
# Strings (taxonomic ids) with repeated occurrences are enumerated (a sequential number n is added to the end of the string).
# The function requires the "collections" package to store and count the string entries.
dereplicate_taxonomy <- function(tax_vector, first_numerated = FALSE){
  
  features_taxonomy_list_dereplicated <- c() # list to store final names of bacterial features
  features_counts <- collections::dict() # dictionary to help store names of bacterial features and their counts.
  
  # for each feature in the taxonomy column
  for(taxonomic_id in tax_vector){
    # if feature's taxonomic id is already in the dictionary (it has already been encountered before).
    if(features_counts$has(taxonomic_id)){
      # set the feature's count +1 in the "features_counts" dictionary.
      features_counts$set(taxonomic_id, features_counts$get(taxonomic_id) + 1)
      # Add the taxonomic id with it's respective count to the features_taxonomy_list_dereplicated.
      features_taxonomy_list_dereplicated <- c(features_taxonomy_list_dereplicated, c( paste(taxonomic_id, features_counts$get(taxonomic_id), sep = "_")))
    }else{
      # if species has not been found before, add species to dictionary.
      features_counts$set(taxonomic_id, 1)
      # Add species to features_taxonomy_list_dereplicated.
      if (first_numerated == TRUE) {
        features_taxonomy_list_dereplicated <- c(features_taxonomy_list_dereplicated, c( paste(taxonomic_id, "1", sep = "_") ))
      }else{
        features_taxonomy_list_dereplicated <- c(features_taxonomy_list_dereplicated, c(taxonomic_id))
      }
    }
  }
  return(features_taxonomy_list_dereplicated)
}

# This function takes a "biom_path" to a biom file with an otu_table and a tax_table, a string "tax_rank",
# and bool "order_table".
# tax_rank parameter must be a value of greengenes ranks format.
# "order_table" indicates if the table should be ordered by larger to smaller values of rowMeans.
# Generally, ASV/OTU tables from QIIME2 are already ordered.
# The ASVs in the biom file are agglomerated by this rank.
# This function returns a dataframe where rows are the ASVs and the columns are samples.
# "rownames" are ASVs taxonomy at the selected rank. "colnames" are samples names.
# Taxonomy is dereplicated so that no row has the same name (which is not allowed).
# The output format is useful for using in other packages like vegan an to generate plots like barplots and heatmaps.

load_biom <- function(biom_path, tax_rank = "Species", order_table = FALSE){
  unite_colNames = get_colNames_per_rank(tax_rank)
  
  if (!is.null(unite_colNames)) {
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_strain)
    extracted_feature_table <- extract_table(biom_object, tax_rank)
  }else{
    print("Please choose a valid taxonomy rank!")
    return()
  }
  
  return(clean_table(extracted_feature_table, order_table = order_table))
  }

extract_table <- function(biom_object, tax_rank){
  unite_colNames = get_colNames_per_rank(tax_rank)
  # Agglomerate tax_table by the chosen tax_rank
  biom_object <- phyloseq::tax_glom(biom_object, taxrank = tax_rank, NArm=TRUE)
  # cbind tax_table and otu_table
  feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(phyloseq::tax_table(biom_object)),
                                                    taxonomy,
                                                    all_of(unite_colNames),
                                                    sep = "_"), "taxonomy"),
                         data.frame(phyloseq::otu_table(biom_object)))
}

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
         Family = {
           # Order
           colNames = c("Order")}
  )
  if (!is.null(colNames)){
    return(colNames)
  }else{
    return(NULL)
  }
}

clean_table <- function(feature_table, order_table){
  # Dereplicate taxonomy
  feature_table["taxonomy"] <- dereplicate_taxonomy(feature_table$taxonomy, first_numerated = FALSE)
  # Set taxonomy column as rownames
  feature_table <- tibble::column_to_rownames(tibble::remove_rownames(feature_table), var = "taxonomy")
  
  if (order_table) {
    # Order by abundances mean, from higher to lower.
    feature_table <- extracted_feature_table[order(rowMeans(extracted_feature_table), decreasing = TRUE),]
  }else{
    return(feature_table)
  }
}

###################################################################################################

get_otu_table_dada <- function(biom_file, starting_col, level = "Species"){
  if(!requireNamespace("BiocManager")){
    install.packages("BiocManager")
  }
  
  if(!requireNamespace("phyloseq")){
    BiocManager::install("phyloseq")
  }
  
  library(phyloseq)
  
  if(!requireNamespace("metagMisc")){
    devtools::install_github("vmikk/metagMisc")
  }
  
  biom_otu_tax <- phyloseq::import_biom(biom_file, parseFunction = phyloseq::parse_taxonomy_greengenes)
  #rank_names(biom_otu_tax)
  tax_sp <- phyloseq::tax_glom(biom_otu_tax, taxrank="Species")
  tax_sp_df <- metagMisc::phyloseq_to_df(tax_sp, addtax = T, addtot = F, addmaxrank = T, sorting = "abundance")
  tax_sp_df["taxonomy"] <- dereplicate_taxonomy(extract_tax_gg_from_biom(table = tax_sp_df))
  out_df <- tax_sp_df[1:nrow(tax_sp_df), starting_col:(ncol(tax_sp_df)-1)]
  #print(tax_sp_df$taxonomy)
  rownames(out_df) <- tax_sp_df$taxonomy
  return(out_df)
}

###################################################################################################

read_qiime_otu_table <- function(table_path, level = "Species"){
  if (!"readr" %in% installed.packages()) install.packages("readr")
  if (!"collections" %in% installed.packages()) install.packages("collections")
  if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
  if (!"tibble" %in% installed.packages()) install.packages("tibble")
  
  # read otu_table "as is"
  otu_table <- readr::read_delim(table_path, skip = 1 ,delim = "\t")
  
  # getting a vector of parsed taxonomy
  if (level == "Species") {
    tax_col <- apply(otu_table["taxonomy"], 1, greengenes_parser_sp_level)
  }
  else if (level == "Strain") {
    tax_col <- apply(otu_table["taxonomy"], 1, greengenes_parser_strain_level)
  }
  
  
  # renaming species in taxonomy.
  # To avoid duplicates, if OTUs with the same names are found more than once "sp + number" is added to name.
  species_list <- c() # list to store final names of OTUs
  species_counts <- collections::dict() # dictionary to help store names of OTUs and their counts
  
  # for each "species" row in taxonomy column
  for(species in tax_col){
    # if species is already in dictionary, has already been found before
    if(species_counts$has(species)){
      # set the species count +1 and add the species including it's count to the species list.
      species_counts$set(species, species_counts$get(species) + 1)
      species_list <- c(species_list, c( paste(species, species_counts$get(species), sep = "_")))
    }else{
      # if species has not been found before, add species to dictionary, and add species to list.
      species_counts$set(species, 1)
      species_list <- c(species_list, c( paste(species, "1", sep = "_") ))
    }
  }
  
  # assigning taxonomy to column parsed_taxonomy
  otu_table["parsed_taxonomy"] <- species_list
  
  # moving tax column to the first column
  otu_table <- cbind(otu_table[, ncol(otu_table)], otu_table[1:nrow(otu_table), 2:(ncol(otu_table)-2)])
  # renaming tax to taxonomy. rename() is a dplyr function.
  otu_table <- dplyr::rename(otu_table, taxonomy = parsed_taxonomy)
  # setting row names and dropping rownames column
  otu_table <- tibble::column_to_rownames(otu_table, var = "taxonomy")
  return(otu_table)
}

###################################################################################################

read_qiime_otu_table2 <- function(table_path){
  if (!"readr" %in% installed.packages()) install.packages("readr")
  if (!"collections" %in% installed.packages()) install.packages("collections")
  if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
  if (!"tibble" %in% installed.packages()) install.packages("tibble")
  
  # Read otu_table "as is"
  user_table <- readr::read_delim(table_path, skip = 1 ,delim = "\t")
  
  # Getting a vector of parsed taxonomy
  
  tax_col <- apply(user_table["taxonomy"], 1, greengenes_parser_strain_level)
  
  # Assigning taxonomy to column taxonomy
  user_table["taxonomy"] <- tax_col
  
  # Moving tax column to the first column
  user_table <- cbind(user_table[, ncol(user_table)], user_table[1:nrow(user_table), 2:(ncol(user_table)-1)])
  
  # Collapse all the ASVs with the same taxonomy assignation
  user_table <- user_table %>%
    group_by(taxonomy) %>%
    summarise_all(sum)
  
  # Setting row names and dropping rownames column
  user_table <- tibble::column_to_rownames(user_table, var = "taxonomy")
  return(user_table)
}