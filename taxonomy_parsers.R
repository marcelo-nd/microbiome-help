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

greengenes_parser_strain_level <- function(string){
  result_string <- ""
  pieces <- strsplit(string, split  = ";")[[1]]
  len_pieces <- length(pieces)
  pieces_counter <- len_pieces
  # while result string is empty.
  while(nchar(result_string) == 0){
    # if we arrived to the end of the categories and not reached taxonomy return the original string, e.g. "Undetermined" (to make sure loop stops).
    if(pieces_counter < 2){
      result_string <- string
    }else{
      # lets analyze the last two pieces of the string
      last_piece <- strsplit(pieces[pieces_counter], split  = "__")[[1]]
      ap_piece <- strsplit(pieces[pieces_counter - 1], split  = "__")[[1]]
      if (length(last_piece) == 2 && (last_piece[1] == " n" )) {
        genus_piece <- strsplit(pieces[pieces_counter - 2], split  = "__")[[1]]
        result_string <- paste(genus_piece[2], ap_piece[2], last_piece[2])
      }
      # if the last piece has length two (means is has a name on it), "s" means species. Ideal taxonomy resolution case.
      else if(length(last_piece) == 2 && (last_piece[1] == " s" || last_piece[1] == "s")){
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

#########

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
# Function that parses a string on MIrROR format and outputs a string in readable format.
# Input: string: "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__Clostridium_acetobutylicum"
# Output: string: "Clostridium acetobutylicum"

mirror_parser_sp_level <- function(string){
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
      if(length(last_piece) >= 2 && (last_piece[1] == " s" || last_piece[1] == "s")){
        result_string <- gsub("_", " ", last_piece[2])
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