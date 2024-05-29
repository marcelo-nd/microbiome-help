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

# to do order alphabetically or by overal abundance.
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
  # CreateS a barplot with panels from several otu tables. This otu tables can correspond to different experiments, runs or treatments.
  # Use "shared_samples = FALSE" to treat each otu table as different experiments.Then, each sample in a table is treated as a replicate.
  # In this case, experiment_names are used to name the panels and samples per otu_table are graphed in their corresponding panel.
  # Use "shared_samples = TRUE" to treat each otu table as runs or treatments that share samples. In this case each sample
  # in a table is treated as a sample taken from the runs or treatments. Panels then correspond to the shared samples in each otu_table,
  # and "experiment_names" (treatments/runs) are then graphed inside panels, one sample of the corresponding treatment/run per panel.
  # "feature_tables" should be a list (e.g. list(otu_table1, otu_table2, otu_table3)).
  # Otu tables are DFs in rows = species and cols = samples format. DFs should contain rownames for species and colnames as sample names.
  # experiments_names is a vector of strings containing the names of experiments/treatments/runs.
  # List of feature tables and experiment/treatment names should be in the same order.
  
  #print(head(feature_tables)) # check the list of feature tables
  
  # 1) Clean, join and gather the otu tables.
  for (table in seq(from = 1, to = length(feature_tables), by=1)) {

    # copy feature table to avoid modifying the original table.
    feature_table2 <- feature_tables[[table]]
    
    #print(head(feature_table2)) # check the working feature table
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table2) > 1) {
      feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
    }
    
    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table2["species"] <- row.names(feature_table2)
    
    # Use dplyr gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table]) # check experiment name that corresponds to working feature table.
    
    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector
    
    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data for all tables to make a barplot.
    if (table == 1) {
      exp_plot_table <- feature_table_g
    }else{
      exp_plot_table <- rbind(exp_plot_table, feature_table_g)
    }
  }
  
  #print(head(exp_plot_table)) # check gathered table
  
  # 2) Keep order of treatments and species
  # Keep order of experiments in graph
  exp_plot_table$experiment <- factor(exp_plot_table$experiment, levels = experiments_names)
  
  # If samples are shared, keep order of samples in graph
  if (shared_samples) {
    exp_plot_table$sample <- factor(exp_plot_table$sample, levels = colnames(feature_table2))
  }
  
  #print(head(exp_plot_table)) # check plot table
  
  # 3) Order the rows of a data frame by the species alphabetical order.
  # to do: select order of species, may choose between overall abundance, user defined or alphabetical
  exp_plot_table <- exp_plot_table %>%
    dplyr::arrange(species)
  
  # Reorder factor, will be useful to be able to reorder species in graphs.
  #exp_plot_table$species <- forcats::fct_relevel(exp_plot_table$species, after = 0)
  #exp_plot_table$species <- forcats::fct_rev(exp_plot_table$species)
  
  print(head(exp_plot_table)) 
  
  # 4) Create and return return graph objects
  # if "shared_samples = TRUE" x-axis is "experiment" then, for each experiment a panel is created and all of their samples are graphed within.
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
  } else{ # if "shared_samples = FALSE" x-axis is "sample" then, a panel is created for each sample, and one sample from each of the traetments/runs are graphed within.
  otu_barplot <- ggplot(exp_plot_table) +
    geom_bar(aes(x = sample, y = abundance, fill = species),
             position = position_fill(),
             stat = "identity") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::scale_fill_manual(values=get_palette()) + # Get color palette
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   legend.title=ggplot2::element_text(size=14), 
                   legend.text=ggplot2::element_text(size=12)) +
    facet_grid(~experiment, scales = "free", space = "free") # this is to remove empty factors due to samples being named differently
  otu_barplot
  return(otu_barplot)
  }
}
