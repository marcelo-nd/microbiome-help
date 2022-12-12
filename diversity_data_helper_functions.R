extract_tax_gg_from_biom <- function(table, level = "Species"){
  tax_list <- c()
  if (level == "Species") {
    for (row in 1:nrow(table)) {
      genus = table$Genus[row]
      species = table$Species[row]
      tax = paste(genus, species, sep = "_")
      tax_list <- c(tax_list, c(tax))
    }
  }
  return(tax_list)
}

###################################################################################################

dereplicate_taxonomy <- function(tax_vector){
  if (!"collections" %in% installed.packages()) install.packages("collections")
  
  species_list <- c() # list to store final names of OTUs
  species_counts <- collections::dict() # dictionary to help store names of OTUs and their counts
  
  # for each "species" row in taxonomy column
  for(species in tax_vector){
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
  return(species_list)
}


###################################################################################################

AlmostEqual <- function(x, y, tolerance=1e-8) {
  diff <- abs(x - y)
  mag <- pmax(abs(x), abs(y))
  ifelse(mag > tolerance, diff/mag <= tolerance, diff <= tolerance)
}
    
###################################################################################################
    
filter_otus_by_counts_zeros <- function(otu_table, percentage){
    return(otu_table[which(rowMeans(! AlmostEqual(otu_table, 0) ) >= percentage), ])
}

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
    return(otu_table[which(rowSums(otu_table >= min_count) >= col_number), ])
}
