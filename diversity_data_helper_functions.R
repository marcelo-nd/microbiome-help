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
# Function that parses a string on greengenes format and outputs a string in readable format.
# Input: string: "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__acetobutylicum"
# Output: string: "Clostridium acetobutylicum"

greengenes_parser_sp_level <- function(string){
    result_string <- ""
    pieces <- strsplit(string, split  = ";")[[1]]
    len_pieces <- length(pieces)
    
    pieces_counter <- len_pieces
    # while result string is empty.
    while(nchar(result_string) == 0){
        # if we arrived to the end of the categories and not reached taxonomy return the original string, e.g. "Undetermined". To make sure loop stops.
        if(pieces_counter < 2){
            result_string <- string
        }else{
            # lets analyse the last two pieces of the string
            last_piece <- strsplit(pieces[pieces_counter], split  = "__")[[1]]
            ap_piece <- strsplit(pieces[pieces_counter - 1], split  = "__")[[1]]
            # if the last piece has length two (means is has a name on it), "s" means species. Ideal taxzonomy resolution case.
            if(length(last_piece) == 2 && (last_piece[1] == " s" || last_piece[1] == "s")){
                result_string <- paste(ap_piece[2], last_piece[2])
            }else if(length(last_piece) == 2 && last_piece[1] != " s"){ #if last piece has a name but it is not species, add sp
                result_string <- paste(last_piece[2], "sp")
            }else if(length(ap_piece) == 2){ #if antepenultimate piece has a name, add sp
                result_string = paste(ap_piece[2], "sp")
            }
        }
        # we go to the next taxonomy level if we did not find a resolved taxonomy in these levels.
        pieces_counter <- pieces_counter - 1
    }
    return(result_string)
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