# Install and load packages
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

source("C:/Users/marce/Documents/GitHub/microbiome-help/otuTableWrangling.R")

get_palette <- function(nColors = 60){
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

# to do order alphabetically or by overall abundance.
barplot_from_feature_table <- function(feature_table, 
                                       plot_title = "", plot_title_size = 14,
                                       x_axis_text_size = 12, x_axis_title_size = 12,
                                       y_axis_title_size = 12, y_axis_text_size = 12,
                                       legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3,
                                       colour_palette = NULL){
  
  # Remove empty rows (species)
  #feature_table2 <- filter_otus_by_counts_col_counts(feature_table, min_count = 1, col_number = 1)
  feature_table2 <- feature_table
  # Saves species names from row_names
  species <- row.names(feature_table2)
  
  # Remove columns (samples) with zero count
  if (ncol(feature_table2) > 1) {
    feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
  }

  # Generate a column with the names of ASVs/OTUs using rownames.
  feature_table2["species"] <- species
  
  print(head(feature_table))
  
  # Gather
  feature_table2 <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
  
  
  print(head(feature_table2))
  # Keep order of samples in graph
  feature_table2$sample <- factor(feature_table2$sample, levels = colnames(feature_table))
  
  print(head(feature_table2))
  
  if (is.null(colour_palette)) {
    print("Colour pallette generated")
    colour_palette <- get_palette(nColors = nrow(feature_table))
  }
  
  otu_barplot <- ggplot2::ggplot(feature_table2, ggplot2::aes(x=sample, y=abundance, fill=species)) + 
    ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
    ggplot2::scale_fill_manual(values=colour_palette) +
    ggtitle(plot_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                   axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_text_size),
                   axis.title.y = ggplot2::element_text(size=y_axis_title_size),
                   axis.text.y = ggplot2::element_text(size = y_axis_text_size),
                   legend.position=legend_pos,
                   legend.title=ggplot2::element_text(size=legend_title_size),
                   legend.text=ggplot2::element_text(size=legend_text_size),
                   ) +
    guides(fill = guide_legend(ncol = legend_cols))
  otu_barplot
  return(otu_barplot)
}

barplot_from_feature_tables <- function(feature_tables, experiments_names, shared_samples = FALSE,
                                        plot_title = "", plot_title_size = 14,
                                        x_axis_text_size = 12, x_axis_title_size = 12,
                                        y_axis_title_size = 12, y_axis_text_size = 12,
                                        legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3, legend_key_size = 1, 
                                        colour_palette = NULL){
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
  sample_names = c()
  for (table in seq(from = 1, to = length(feature_tables), by=1)) {
    
    # copy feature table to avoid modifying the original table.
    feature_table2 <- feature_tables[[table]]
    
    #print(head(feature_table2)) # check the working feature table
    
    # Remove rows with Zero counts
    feature_table2 <- filter_species_by_col_counts(feature_table2, min_count = 1, col_number = 1)
    
    print(head(feature_table2))
    
    # save names of species
    species_names <- row.names(feature_table2)
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table2) > 1) {
      feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
    }
    
    sample_names <- c(sample_names, colnames(feature_table2))
    
    #print(head(feature_table2))
    
    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table2["species"] <- species_names
    #print(feature_table2$species)
    
    # Use dplyr gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table2, 1:(ncol(feature_table2) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table]) # check experiment name that corresponds to working feature table.
    
    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector
    
    #print(head(feature_table_g))
    
    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data for all tables to make a barplot.
    if (table == 1) {
      exp_plot_table <- feature_table_g
    }else{
      exp_plot_table <- rbind(exp_plot_table, feature_table_g)
    }
  }
  
  #print(sample_names)
  #print(head(exp_plot_table)) # check gathered table
  
  # 2) Keep order of treatments and species
  # Keep order of experiments in graph
  exp_plot_table$experiment <- factor(exp_plot_table$experiment, levels = experiments_names)
  
  #print(head(exp_plot_table))
  
  # If samples are shared, keep order of samples in graph
  if (shared_samples) {
    exp_plot_table$sample <- factor(exp_plot_table$sample, levels = colnames(feature_table2))
  } else{
    #exp_plot_table$sample <- factor(exp_plot_table$sample, levels = sample_names)
  }
  
  
  #print(head(exp_plot_table)) # check plot table
  
  # 3) Order the rows of a data frame by the species alphabetical order.
  # to do: select order of species, may choose between overall abundance, user defined or alphabetical
  exp_plot_table <- exp_plot_table %>%
    dplyr::arrange(species)
  
  # Reorder factor, will be useful to be able to reorder species in graphs.
  #exp_plot_table$species <- forcats::fct_relevel(exp_plot_table$species, after = 0)
  #exp_plot_table$species <- forcats::fct_rev(exp_plot_table$species)
  
  # 4) Create and return return graph objects
  # Check if color palette was passed
  if (is.null(colour_palette)) {
    colour_palette <- get_palette(nColors = length(unique(exp_plot_table$species)))
  }
  
  # if "shared_samples = TRUE" x-axis is "experiment" then, for each experiment a panel is created and all of their samples are graphed within.
  if (shared_samples) {
    otu_barplot <- ggplot(exp_plot_table) +
      geom_bar(aes(x = experiment, y = abundance, fill = species),
               position = position_fill(),
               stat = "identity") + 
      ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                     axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_text_size),
                     axis.title.y = ggplot2::element_text(size=y_axis_title_size),
                     axis.text.y = ggplot2::element_text(size = x_axis_text_size),
                     legend.position=legend_pos,
                     legend.title=ggplot2::element_text(size=legend_title_size),
                     legend.text=ggplot2::element_text(size=legend_text_size)) +
      ggplot2::scale_fill_manual(values=colour_palette) + # Get color palette
      facet_grid(~sample) +
      guides(fill = guide_legend(ncol = legend_cols))
    otu_barplot
    return(otu_barplot)
  } else{ # if "shared_samples = FALSE" x-axis is "sample" then, a panel is created for each sample, and one sample from each of the treatments/runs are graphed within.
    otu_barplot <- ggplot(exp_plot_table) +
      geom_bar(aes(x = sample, y = abundance, fill = species),
               position = position_fill(),
               stat = "identity") +
      ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                     axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_text_size),
                     axis.title.y = ggplot2::element_text(size=y_axis_title_size),
                     axis.text.y = ggplot2::element_text(size = x_axis_text_size),
                     legend.title=ggplot2::element_text(size=legend_title_size),
                     legend.text=ggplot2::element_text(size=legend_text_size),
                     legend.position=legend_pos, legend.key.size = unit(legend_key_size, "cm")) +
      ggplot2::scale_fill_manual(values=colour_palette) + # Get color palette
      facet_grid(~experiment, scales = "free", space = "free") + # this is to remove empty factors due to samples being named differently
      guides(fill = guide_legend(ncol = legend_cols))+
      guides(fill = guide_legend(override.aes = list(size=1)))
    otu_barplot
    return(otu_barplot)
  }
}

barplot_from_feature_table_sorted <- function(feature_table, sort_type = NULL, species_to_sort = NULL,
                                              plot_title = "", plot_title_size = 14,
                                              x_axis_text_size = 12, x_axis_title_size = 12,
                                              y_axis_title_size = 12, y_axis_text_size = 12,
                                              legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3,
                                              colour_palette = NULL){
  # Specify the species of interest
  if (!is.null(species_to_sort) && sort_type == "species_abundance") {
    print("Sort samples by species_abundance")
    
    # Make "Species" column with the rownames 
    df <- feature_table %>% rownames_to_column(var = "Species")
    
    total_abundance <- colSums(df[, -1])
    
    # Filter the row of the species of interest and calculate its proportion with respect to total abundance
    df_proportion <- df %>%
      filter(Species == species_to_sort) %>%
      select(-Species)
    # calculate species of interest proportion
    df_proportion <- df_proportion[1,]/total_abundance
    # Get sample names sorted by the species of interest proportion
    ordered_samples <- df_proportion %>%
      unlist() %>%
      sort(decreasing = TRUE) %>%
      names()
    
    df_long <- df %>%
      pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")
    
    df_long$Sample <- factor(df_long$Sample, levels = ordered_samples)
    
    otu_barplot <- ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species)) + 
      ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
      ggplot2::scale_fill_manual(values=colour_palette) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                     axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_text_size),
                     axis.title.y = ggplot2::element_text(size=y_axis_title_size),
                     axis.text.y = ggplot2::element_text(size = x_axis_text_size),
                     legend.position=legend_pos,
                     legend.title=ggplot2::element_text(size=legend_title_size),
                     legend.text=ggplot2::element_text(size=legend_text_size)) +
      guides(fill = guide_legend(ncol = legend_cols))
    
  }else if (sort_type == "similarity") {
    print("Sort samples by similarity")
    
    ##### Plot with bars ordered by similarity
    
    # transform table
    df2 <- transform_feature_table(feature_table = feature_table, transform_method = "min_max")
    
    # Get the order of samples based on clustering
    ordered_samples_cluster <- order_samples_by_clustering(df2)
    
    df2 <- df2 %>% rownames_to_column(var = "Species")
    
    #print(head(df2))
    
    df_long <- df2 %>%
      pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")
    
    # Update sample factor levels in the long-format data for ggplot
    df_long$Sample <- factor(df_long$Sample, levels = ordered_samples_cluster)

    # Plot ordered by clustering similarity
    otu_barplot <- ggplot2::ggplot(df_long, ggplot2::aes(x=Sample, y=Abundance, fill=Species)) + 
      ggplot2::geom_bar(position="fill", stat="identity", show.legend = TRUE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::scale_fill_manual(values=colour_palette) +
      ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                     axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = x_axis_text_size),
                     axis.title.y = ggplot2::element_text(size=y_axis_title_size),
                     axis.text.y = ggplot2::element_text(size = x_axis_text_size),
                     legend.position="right",
                     legend.title=ggplot2::element_text(size=legend_title_size),
                     legend.text=ggplot2::element_text(size=legend_text_size)) +
      guides(fill = guide_legend(ncol = 1))
    otu_barplot
    
  }else{
    "No valid option chosen"
    #break
  }

  otu_barplot
  return(otu_barplot)
}

barplot_with_replicates <- function(feature_table){
  print(head(feature_table))
  #df_rel <- feature_table
  df_rel <- calculate_relative_abundance(feature_table)
  
  df_rel <- rownames_to_column(df_rel, var = "Species")
  
  ##### Calculate cumulative half addition
  df_hcs <- cumulative_half_addition(df_rel)
  
  ### Now lets put together replicate data
  ### Calculate the means of each sample
  df_means_rel <- df_rel %>%
    # Pivot to long format for easier manipulation
    pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
    # Extract sample IDs without replicate information
    mutate(Sample_ID = sub("_R[0-9]+$", "", Sample)) %>%
    # Group by species and sample to calculate mean abundance
    group_by(Species, Sample_ID) %>%
    summarize(Mean_Abundance = mean(Abundance), .groups = "drop")%>%        # Extract replicate number
    filter(Mean_Abundance != 0)
  
  df_long_dots <- df_hcs %>%
    pivot_longer(cols = -Species, names_to = "Sample", values_to = "Abundance") %>%
    mutate(
      Sample_ID = sub("_R[0-9]+$", "", Sample),  # Extract sample name without replicate info
      Replicate = sub(".*_R", "", Sample)) %>%        # Extract replicate number
    filter(Abundance != 0)
  
  (plot2 <- ggplot(NULL) + 
      ggplot2::geom_bar(data = df_means_rel, aes(x=Sample_ID, y=Mean_Abundance, fill=Species, color = Species), position="fill", stat="identity", show.legend = TRUE) +
      geom_point(data = df_long_dots, aes(x=Sample_ID, y=Abundance, shape = Species), size = 3, color = "white", stroke = 1.3) +
      scale_shape_manual(values = c(0, 1, 17, 3, 4, 5, 6, 20, 8, 9, 10)) +
      ggthemes::theme_tufte()
  )
  
  return(plot2)
}

barplot_w_strain_data <- function(otu_table, strain_data){
  #library(ggpattern)
 # Does a barplot but also showing strain level belonging of identified OTUs
  # Minmax transform
  ot_scree_filtered_norm <- transform_feature_table(otu_table, transform_method = "min_max") # is this necessary?
  # Get clustering order
  ordered_samples_cluster <- order_samples_by_clustering(ot_scree_filtered_norm)
  
  # Generate df long format for Species abundance data
  ot_scree_filtered2 <- ot_scree_filtered %>% rownames_to_column(var = "Species")
  df_otu_long <- ot_scree_filtered2 %>%
    pivot_longer(-Species, names_to = "Sample", values_to = "Abundance")
  
  # Update sample factor levels in the long-format data for ggplot
  df_otu_long$Sample <- factor(df_otu_long$Sample, levels = ordered_samples_cluster)
  
  ##### Now lets work with the strain data
  
  strain_data2 <- strain_data %>%
    mutate(Strain_Number = rep(1:3, times = 10)) # 1, 2, 3 for each species
  # TODO: make a flexible approach for varying strain number.
  
  # Replace `1`s with the strain number for each sample, keeping `0`s unchanged
  for (sample_col in colnames(strain_data2[4:ncol(strain_data2)-1])) {
    strain_data2[[sample_col]] <- ifelse(strain_data2[[sample_col]] == 1, strain_data2$Strain_Number, 0)
  }
  
  # Drop the helper "Strain_Number" column to get the desired output
  df_final <- strain_data2 %>% select(-Strain_Number)
  
  # Collapse the information by species
  df_collapsed <- df_final %>%
    mutate(Species = sapply(strsplit(Species, " "), function(x) paste(x[1:2], collapse = " "))) %>% # Extract species name
    group_by(Species) %>%
    summarise(across(starts_with("SC"), max)) %>% # Take max per sample to represent strain
    ungroup()
  
  # This should not be necessary
  new_row <- c("Staphylococcus aureus", rep(1, ncol(df_collapsed) - 1))
  new_row_df <- as.data.frame(t(new_row), stringsAsFactors = FALSE)
  colnames(new_row_df) <- colnames(df_collapsed)
  #
  
  df_collapsed2 <- rbind(df_collapsed, new_row_df)
  
  ##### Generate table for 
  df_otu_long2 <- df_collapsed2 %>% # this might be df_strain_long
    pivot_longer(-Species, names_to = "Sample", values_to = "Strain")
  
  ### Join the two long dataframes
  
  df_joined <- dplyr::full_join(df_otu_long, df_otu_long2, by = c("Species", "Sample"))
  
  df_joined_filtered <- df_joined %>%
    filter(!is.na(Abundance) & Abundance != 0)
  
  df_joined_filtered <- df_joined_filtered %>%
    filter(!is.na(Strain) & Strain != 0)
  
  df_joined_filtered <- df_joined_filtered %>%
    mutate(Strain = factor(
      case_when(
        Strain == 1 ~ "Strain 1",
        Strain == 2 ~ "Strain 2",
        Strain == 3 ~ "Strain 3"
      )
    ))
  
  # Update sample factor levels in the long-format data for ggplot
  df_joined_filtered$Sample <- factor(df_joined_filtered$Sample, levels = ordered_samples_cluster)
  
  ####### Now plotting
  
  colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                   "lightblue1","brown1", "olivedrab3", "darkorange3", "#23001E","hotpink" )
  
  
  # this is the one, do not touch
  p1 <- ggplot(data = df_joined_filtered, aes(x = Sample, y=Abundance)) +
    ggpattern::geom_bar_pattern(aes(fill = Species, pattern = Strain, pattern_density = Strain),
                     position = "fill",
                     stat="identity",
                     show.legend = TRUE,
                     pattern_color = "white",
                     pattern_fill = "white",
                     pattern_angle = 45,
                     pattern_spacing = 0.025) +
    ggplot2::scale_fill_manual(values=colours_vec) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                   axis.text.y = ggplot2::element_text(size=12),
                   legend.text = element_text(size=12)) +
    guides(pattern = guide_legend(override.aes = list(fill = "black")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
    scale_pattern_density_manual(values = c(0, 0.2, 0.05))
  
  #return(df_collapsed)
  #return(df_collapsed2)
  #return(df_otu_long)
  #return(df_otu_long2)
  return(df_joined_filtered)
  #return(p1)
}

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
  print(df)
}

barplot_w_strain_data2 <- function(feature_table, strain_data){
  #library(ggpattern)
  # Does a barplot but also showing strain level belonging of identified OTUs
  strain_numers_ft <- strain_name2strain_number(feature_table)
  
  # Assuming your dataframe is named df
  df_long <- strain_numers_ft %>%
    rownames_to_column(var = "Species") %>%  # Convert row names to a column
    pivot_longer(
      cols = -Species,  # All columns except "Species" will be pivoted
      names_to = "Sample",
      values_to = "Abundance"
    ) %>%
    mutate(
      Strain = paste0("Strain ", sub(".* ", "", Species)),  # Extract last number as strain
      Species = sub(" \\d+$", "", Species)  # Remove strain number from species name
    )
  
  df_long_filtered <- df_long %>%
    filter(!is.na(Abundance) & Abundance != 0)
  
  df_long_filtered <- df_long_filtered %>%
    filter(!is.na(Strain) & Strain != 0)
  
  if (TRUE) {
    # Minmax transform
    ft_std <- transform_feature_table(feature_table, transform_method = "min_max") # is this necessary?
    # Get clustering order
    ordered_samples_cluster <- order_samples_by_clustering(ft_std)
    # Update sample factor levels in the long-format data for ggplot
    df_long_filtered$Sample <- factor(df_long_filtered$Sample, levels = ordered_samples_cluster)
  }
  
  colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                   "lightblue1","brown1", "olivedrab3", "darkorange3", "#23001E","hotpink" )
  
  
  # this is the one, do not touch
  p1 <- ggplot(data = df_long_filtered, aes(x = Sample, y=Abundance)) +
    ggpattern::geom_bar_pattern(aes(fill = Species, pattern = Strain, pattern_density = Strain),
                                position = "fill",
                                stat="identity",
                                show.legend = TRUE,
                                pattern_color = "white",
                                pattern_fill = "white",
                                pattern_angle = 45,
                                pattern_spacing = 0.025) +
    ggplot2::scale_fill_manual(values=colours_vec) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                   axis.text.y = ggplot2::element_text(size=12),
                   legend.text = element_text(size=12)) +
    guides(pattern = guide_legend(override.aes = list(fill = "black")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
    scale_pattern_density_manual(values = c(0, 0.2, 0.05))
  
  #return(strain_numers_ft)
  #return(df_long)
  #return(df_long_filtered)
  #return(ordered_samples_cluster)
  return(p1)
}
