parse_taxonomy_strain <- function(char.vec){
  named.char.vec <- substring(char.vec, first = 4)
  names(named.char.vec) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
  return(named.char.vec)
}

load_biom <- function(path, tax_rank = "Species", dereplicate = FALSE){
  biom <- phyloseq::import_biom(biom_path, parseFunction=parse_taxonomy_strain)
  biom <- phyloseq::tax_glom(biom, taxrank = tax_rank, NArm=TRUE)
  # cbind tax_table and otu_table
  if (tax_rank == "Strain") {
    feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(tax_table(biom)),
                                                      taxonomy,
                                                      c("Genus", "Species", "Strain"),
                                                      sep = "_"), "taxonomy"),
                           data.frame(otu_table(biom)))
  }else if (tax_rank == "Species") {
    feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(tax_table(biom)),
                                                      taxonomy,
                                                      c("Genus", "Species"),
                                                      sep = "_"), "taxonomy"),
                           data.frame(otu_table(biom)))
  }else if (tax_rank == "Family" || tax_rank == "Order"){
    feature_table <- cbind(dplyr::select(tidyr::unite(data.frame(tax_table(biom)),
                                                      taxonomy,
                                                      tax_rank,
                                                      sep = "_"), "taxonomy"),
                           data.frame(otu_table(biom)))
  }
  
  # Dereplicate taxonomy
  return(tibble::column_to_rownames(tibble::remove_rownames(feature_table), var = "taxonomy"))
}


tax_table_test <- data.frame(tax_table(tax_test))

tax_table_test

test_str <- strsplit("k__Bacteria; p__Actinobacteria; c__Actinomycetia; o__Propionibacteriales; f__Propionibacteriaceae; g__Cutibacterium; s__acnes; n__ NR_113028.1", split  = "; ")[[1]]
test_str

parse_taxonomy_strain(test_str)

typeof(parse_taxonomy_default(test_str))

tax_table(test_biom)








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