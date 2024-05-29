# Install and load packages
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

get_palette <- function(nColors = 50){
  colors_vec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442","#0072B2",
    "brown1", "#CC79A7", "olivedrab3", "rosybrown", "darkorange3",
    "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
    "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
    "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
    "lavenderblush4", "lawngreen", "indianred1", "lightblue1", "honeydew4",
    "hotpink", "#e3ae78", "#a23f3f", "#290f76", "#ce7e00",
    "#386857", "#738564", "#e89d56", "#cd541d", "#1a3a46",
    "#9C4A1A", "#ffe599", "#583E26", "#A78B71", "#F7C815",
    "#EC9704", "#4B1E19", "firebrick2", "#C8D2D1", "#14471E",
    "#6279B8", "#DA6A00", "#C0587E", "#FC8B5E", "#FEF4C0",
    "#EA592A", "khaki3", "lavenderblush3", "indianred4", "lightblue",
    "honeydew1", "hotpink4", "ivory3", "#49516F", "#502F4C",
    "#A8C686", "#669BBC", "#29335C", "#E4572E", "#F3A712",
    "#EF5B5B", "#FFBA49", "#20A39E", "#23001E", "#A4A9AD")
  
  #set.seed(1)

  return(colors_vec[sample(1:length(colors_vec), size = nColors)])
  }

barplot_from_feature_table <- function(feature_table){
  # Remove columns (samples) with zero count
  if (ncol(feature_table) > 1) {
    feature_table <- feature_table[, colSums(feature_table != 0) > 0]
  }
  
  feature_table2 <- feature_table
  
  # Generate a column with the names of ASVs/OTUs using rownames.
  feature_table2["species"] <- row.names(feature_table2)
  
  #print(head(feature_table))
  
  # Gather
  feature_table2 <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
  
  
  print(head(feature_table2))
  # Keep order of samples in graph
  feature_table2$sample <- factor(feature_table2$sample, levels = colnames(feature_table))
  
  print(head(feature_table2))
  
  color_palette <- get_palette(nColors = nrow(feature_table))
  
  otu_barplot <- ggplot2::ggplot(feature_table2, ggplot2::aes(x=sample, y=abundance, fill=species)) + 
    ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::scale_fill_manual(values=color_palette) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
          legend.title=ggplot2::element_text(size=14), 
          legend.text=ggplot2::element_text(size=12))
  otu_barplot
  return(otu_barplot)
}


barplot_from_feature_tables <- function(feature_tables, experiments_names, shared_samples = FALSE){
  # Create a barplot from several otu tables. This otutables can correspond to different experiments, runs or treatments.
  # If treatments are used, it can graph samples in a different pane per treatment where each treatment should share the samples.
  # If only different Samples types are to be graphed, samples per treatment can correspond to replicates.
  # feature_tables should be a list (e.g. list(otu_table1, otu_table2, otu_table3))
  # List of feature tables and experiment/treatment names should be in the same order.
  # 1) Clean, join and gather the otu tables.
  for (table in seq(from = 1, to = length(feature_tables), by=1)) {
    #print(head(feature_tables)) # check the list of feature tables
    
    # copy feature table to avoid modifiyng the original table.
    feature_table2 <- feature_tables[[table]]
    
    #print(head(feature_table2)) # check the working feature table
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table2) > 1) {
      feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
    }
    
    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table2["species"] <- row.names(feature_table2)
    
    # Gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table])
    
    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector (order is important)
    
    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data to make a barplot.
    if (table == 1) {
      exp_plot_table <- feature_table_g
    }else{
      exp_plot_table <- rbind(exp_plot_table, feature_table_g)
    }
  }
  
  #print(head(exp_plot_table)) # check gathered table
  
  # Keep order of experiments in graph
  exp_plot_table$experiment <- factor(exp_plot_table$experiment, levels = experiments_names)
  
  # If samples are shared, keep order of samples in graph
  if (shared_samples) {
    exp_plot_table$sample <- factor(exp_plot_table$sample, levels = colnames(feature_table2))
  }
  
  #print(head(exp_plot_table)) # check plot table
  
  #exp_plot_table$experiment <- factor(exp_plot_table$experiment,levels=unique(as.character(exp_plot_table$experiment)))
  #exp_plot_table$species <- factor(exp_plot_table$species,levels = unique(as.character(exp_plot_table$species)))
  
  # Orders the rows of a data frame by the species alphabetical order.
  # to do: select order of species, may chonce between overal abundance, user defined or alphabetical
  exp_plot_table <- exp_plot_table %>%
    dplyr::arrange(species)
  
  # Reorder factor
  #exp_plot_table$species <- forcats::fct_relevel(exp_plot_table$species, after = 0)
  #exp_plot_table$species <- forcats::fct_rev(exp_plot_table$species)
  
  print(head(exp_plot_table)) 
  
  if (shared_samples) {
    otu_barplot <- ggplot(exp_plot_table) +
      geom_bar(aes(x = experiment, y = abundance, fill = species),
               position = position_fill(),
               stat = "identity") + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::scale_fill_manual(values=get_palette()) + # Get color palette
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                     legend.title=ggplot2::element_text(size=14), 
                     legend.text=ggplot2::element_text(size=12)) +
      facet_grid(~sample)
    otu_barplot
    return(otu_barplot)
  } else{
  otu_barplot <- ggplot(exp_plot_table) +
    geom_bar(aes(x = sample, y = abundance, fill = species),
             position = position_fill(),
             stat = "identity") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::scale_fill_manual(values=get_palette()) + # Get color palette
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   legend.title=ggplot2::element_text(size=14), 
                   legend.text=ggplot2::element_text(size=12)) +
    facet_grid(~experiment, scales = "free", space = "free") # this is to remove empty factors
  otu_barplot
  return(otu_barplot)
  }
}
