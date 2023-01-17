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
    biom_object <- phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_greengenes)
  }
  
  extracted_feature_table <- extract_table(biom_object, tax_rank, unite_colNames)
  
  return(clean_table(extracted_feature_table, order_table = order_table))
}
