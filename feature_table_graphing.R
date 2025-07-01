# Install and load packages
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

if (!require("tibble", quietly = TRUE))
  install.packages("tibble")

if (!require("ggpattern", quietly = TRUE))
  install.packages("ggpattern")

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpattern)

source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

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

# todo: strain level
barplot_from_feature_table <- function(feature_table, sort_type = "none", feature_to_sort = NULL, strains = FALSE,
                                       plot_title = "", plot_title_size = 14,
                                       x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                                       y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 90,
                                       legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3,
                                       x_vjust = 0.5, x_hjust = 1, 
                                       colour_palette = NULL){
  ### Step 1. Clean feature table
  # Remove empty rows (features)
  feature_table2 <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1) # why is this not working???
  #feature_table2 <- feature_table
  
  # Remove columns (samples) with zero count
  if (ncol(feature_table2) > 1) {
    feature_table2 <- feature_table2[, colSums(feature_table2 != 0) > 0]
  }
  
  if (isTRUE(strains)) {
    # Convert table with strain names to a strain-number table
    feature_table2 <- strain_name2strain_number(feature_table2)
  }
  
  # Saves species names from row_names
  species <- row.names(feature_table2)
  
  print(head(feature_table2))
  
  ### Step 2. If sorting, determine sample order.
  if (sort_type == "feature_value" && !is.null(feature_to_sort)) {
    print("Sort samples by feature_value")
    # Make "Species" column with the rownames 
    df1 <- feature_table2 %>% rownames_to_column(var = "species")
    
    total_abundance <- colSums(df1[, -1])
    
    # Filter the row of the species of interest and calculate its proportion with respect to total abundance
    df_proportion <- df1 %>%
      filter(species == feature_to_sort) %>%
      select(-species)
    # calculate species of interest proportion
    df_proportion <- df_proportion[1,]/total_abundance
    # Get sample names sorted by the species of interest proportion
    ordered_samples <- df_proportion %>%
      unlist() %>%
      sort(decreasing = TRUE) %>%
      names()
    
  }else if (sort_type == "similarity") {
    print("Sort samples by similarity")
    
    # transform table
    df1 <- transform_feature_table(feature_table = feature_table2, transform_method = "min_max")
    
    # Get the order of samples based on clustering
    ordered_samples <- order_samples_by_clustering(df1)
    
    df1 <- df1 %>% rownames_to_column(var = "species")
    
  }else if (sort_type == "none") {
    print("No sorting chosen")
    df1 <- feature_table2
    ordered_samples <- colnames(feature_table2)
    # Generate a column with the names of ASVs/OTUs using rownames.
    df1["species"] <- species
  }else{
    print("No valid sorting option chosen")
    return()
  }
  
  #print(head(df1))
  #print(species)
  
  ### Step 3. Process features table to ploting table.
  # create the plot table
  plot_df <- df1 %>%
    pivot_longer(-species, names_to = "sample", values_to = "abundance")
  
  # If strain processing has to be done.
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      mutate(
        strain = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)  # Remove strain number from species name
      )
  }
  
  ### Step 4. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    filter(!is.na(abundance) & abundance != 0)
  
  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      filter(!is.na(strain) & strain != 0)
  }
  
  # Factor the "sample" variable so the order of samples is as in "ordered_samples" variable
  plot_df_filtered$sample <- factor(plot_df_filtered$sample, levels = ordered_samples)
  
  print(head(plot_df_filtered))
  
  ### Step 4. Create plot.
  if (is.null(colour_palette)) { # get colour palette
    print("Colour pallette generated")
    nfeatures <- length(unique(plot_df_filtered$species))
    colour_palette <- get_palette(nColors = nfeatures)
    print(colour_palette)
  }
  
  # Create base plot.
  ft_barplot <- ggplot2::ggplot(plot_df_filtered, ggplot2::aes(x=sample, y=abundance, fill=species))
  
  if (isTRUE(strains)) {
    print("strains processing")
    ft_barplot <- ft_barplot + ggpattern::geom_bar_pattern(aes(fill = species2, pattern = strain, pattern_density = strain),
                                                           position = "fill",
                                                           stat="identity",
                                                           show.legend = TRUE,
                                                           pattern_color = "white",
                                                           pattern_fill = "white",
                                                           pattern_angle = 45,
                                                           pattern_spacing = 0.025) +
      ggpattern::scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(pattern = guide_legend(override.aes = list(fill = "black")),
             fill = guide_legend(override.aes = list(pattern = "none")))
  } else{
    print("no strains")
    ft_barplot <- ft_barplot + geom_bar(aes(fill = species),
                                        position = position_fill(),
                                        stat = "identity")
  }
  
  # add theme options
  ft_barplot <- ft_barplot  +
    theme_void() +
    ggplot2::scale_fill_manual(values=colour_palette) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                   axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                   axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, vjust = x_vjust, hjust= x_hjust, size = x_axis_text_size),
                   axis.title.y = ggplot2::element_text(size=y_axis_title_size, angle = y_axis_text_angle, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                   axis.text.y = ggplot2::element_text(size = x_axis_text_size),
                   legend.position=legend_pos,
                   legend.title=ggplot2::element_text(size=legend_title_size),
                   legend.text=ggplot2::element_text(size=legend_text_size)) +
    guides(fill = guide_legend(ncol = legend_cols))
  
  ft_barplot # show plot
  return(ft_barplot) # return plot
}

dendrogram_from_feature_table <- function(df, diss_method = "euclidean", clust_method = "ward.D2",
                                          margin_t = 40, margin_r = 0, margin_b = -25, margin_l = 0){
  # Prepare table
  df_n <- transform_feature_table(df, transform_method = "min_max")
  
  df_n <- df_n %>% rownames_to_column(var = "Species")
  
  df_t <- as.matrix(t(df_n[, -1]))  # Exclude the "Species" column after moving it to row names
  
  #print(head(df_t))
  
  # Perform hierarchical clustering
  d <- dist(df_t, method = diss_method)
  hc <- hclust(d, method = clust_method)
  
  dendo_plot <- ggdendro::ggdendrogram(hc, rotate = 0,
                                       leaf_labels = F) +
    theme(plot.margin = margin(t = margin_t,  # Top margin
                               r = margin_r,  # Right margin
                               b = margin_b,  # Bottom margin
                               l = margin_l)) + # Left margin
    scale_y_continuous(expand = expansion(add = c(1, 1)), labels = NULL) + 
    scale_x_continuous(expand = expansion(add = c(0.35, 0.5)), labels = NULL)
  
  dendo_plot
  return(dendo_plot)
}


barplots_grid <- function(feature_tables, experiments_names, shared_samples = FALSE, strains = FALSE, plot_title = "",
                          plot_title_size = 14, x_axis_text_size = 12, x_axis_title_size = 12, x_axis_text_angle = 0,
                          y_axis_title_size = 12, y_axis_text_size = 12, y_axis_text_angle = 0,
                          legend_pos = "right", legend_title_size = 12, legend_text_size = 12, legend_cols = 3, legend_key_size = 1, 
                          colour_palette = NULL){
  # Creates a grid of Barplots
  
  ### Step 1. Clean, join and gather the otu tables.
  sample_names = c() # to keep track of the sample names
  for (table in seq(from = 1, to = length(feature_tables), by=1)) { # iterate over all the feature tables
    # copy current feature table to avoid modifying the original table.
    feature_table <- feature_tables[[table]]
    
    #print(head(feature_table2)) # check the working feature table
    
    if (isTRUE(strains)) {
      # Convert table with strain names to a strain-number table
      feature_table <- strain_name2strain_number(feature_table)
    }
    
    # Remove rows with Zero counts
    feature_table <- filter_features_by_col_counts(feature_table, min_count = 1, col_number = 1)
    
    #print(head(feature_table2))
    
    # save names of species
    species_names <- row.names(feature_table)
    
    # Remove columns (samples) with zero count
    if (ncol(feature_table) > 1) {
      feature_table <- feature_table[, colSums(feature_table != 0) > 0]
    }
    
    sample_names <- c(sample_names, colnames(feature_table))
    
    #print(head(feature_table2))
    
    # Create a column with the names of ASVs/OTUs using rownames.
    feature_table["species"] <- species_names
    #print(feature_table2$species)
    
    # Use dplyr gather the working feature table.
    feature_table_g <- tidyr::gather(feature_table, 1:(ncol(feature_table) - 1) , key = "sample", value = "abundance")
    
    #print(experiments_names[table]) # check experiment name that corresponds to working feature table.
    
    # Create a column to keep track of from which experiment/treatment the samples come from.
    feature_table_g$experiment <- experiments_names[table] # the experiment name is taken from experiments_names vector
    
    #print(head(feature_table_g))
    
    # rbind the gathered feature tables.
    # Result is exp_plot_table, a table containing in each row species;sample;abundance;experiment data for all tables to make a barplot.
    if (table == 1) {
      plot_df <- feature_table_g
    }else{
      plot_df <- rbind(plot_df, feature_table_g)
    }
  }
  print(sample_names) # check sample_names
  print(head(plot_df)) # check gathered table
  
  ### Step 2. Convert Strain data to a graphing-compatible format.
  # Add strain data column to long dataframe
  if (isTRUE(strains)) {
    plot_df <- plot_df %>%
      mutate(
        strain = paste0("Strain ", sub(".* ", "", species)),  # Extract last number as strain
        species2 = sub(" \\d+$", "", species)  # Remove strain number from species name
      )
  }
  
  print(head(plot_df))
  
  ### Step 3. Clean the long-format table
  plot_df_filtered <- plot_df %>%
    filter(!is.na(abundance) & abundance != 0)
  
  if (isTRUE(strains)) {
    plot_df_filtered <- plot_df_filtered %>%
      filter(!is.na(strain) & strain != 0)
  }
  
  plot_df_filtered$experiment <- factor(plot_df_filtered$experiment, levels = experiments_names)
  
  ### Step 4. Plotting
  # get color palette
  if (is.null(colour_palette)) {
    colour_palette <- get_palette(nColors = length(unique(plot_df$species)))
  }
  
  print(plot_df_filtered) # check final table prevouos to plotting
  
  # Create base plot.
  if (shared_samples) {
    p1 <- ggplot(data = plot_df_filtered, aes(x = experiment, y=abundance)) +
      facet_grid(~sample)
  } else{
    p1 <- ggplot(data = plot_df_filtered, aes(x = sample, y=abundance)) +
      facet_grid(~experiment, scales = "free", space = "free")
  }
  
  # Add elements based on graph type.
  if (isTRUE(strains)) {
    print("strains processing")
    p1 <- p1 + ggpattern::geom_bar_pattern(aes(fill = species2, pattern = strain, pattern_density = strain),
                                           position = "fill",
                                           stat="identity",
                                           show.legend = TRUE,
                                           pattern_color = "white",
                                           pattern_fill = "white",
                                           pattern_angle = 45,
                                           pattern_spacing = 0.025) +
      ggpattern::scale_pattern_manual(values = c("Strain 1" = "none", "Strain 2" = "circle", "Strain 3" = "stripe")) +
      ggpattern::scale_pattern_density_manual(values = c(0, 0.2, 0.1)) +
      guides(pattern = guide_legend(override.aes = list(fill = "black")),
             fill = guide_legend(override.aes = list(pattern = "none")))
  } else{
    print("no strains")
    p1 <- p1 + geom_bar(aes(fill = species),
                        position = position_fill(),
                        stat = "identity")
  }
  
  if (!is.null(colour_palette)) {
    p1 <- p1 + ggplot2::scale_fill_manual(values=colour_palette)
  } else{
    print("Colours vec is null, using standard color palette.")
  }
  
  p1 <- p1 +
    theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = plot_title_size, face = "bold", hjust = 0.5, vjust = 0.5),
                   axis.title.x = ggplot2::element_text(size=x_axis_title_size),
                   axis.text.x = ggplot2::element_text(angle = x_axis_text_angle, vjust = 0.5, hjust=1, size = x_axis_text_size),
                   axis.title.y = ggplot2::element_text(size=y_axis_title_size, angle = 90),
                   axis.text.y = ggplot2::element_text(size = x_axis_text_size, angle = y_axis_text_angle),
                   legend.title=ggplot2::element_text(size=legend_title_size),
                   legend.text=ggplot2::element_text(size=legend_text_size),
                   legend.position=legend_pos, legend.key.size = unit(legend_key_size, "cm")) + 
    guides(fill = guide_legend(ncol = legend_cols))

  # Show plot
  p1
  
  return(p1)
}

feature_table_heatmap <- function(ft1, ft2 = NULL, sig_stars = FALSE, corr_type = "pearson", pval_adjust = TRUE, fdr_correction = "fdr", axis_text_size = 1, stars_size = 0.5, hm_type = "full"){
  ft1_t <- t(ft1)
  ft1_t <- ft1_t[order(row.names(ft1_t)), ] # Ordering by row names
  #print(rownames(ft1_t))
  if (is.null(ft2)) {
    # Get correlation results
    corr_results <- Hmisc::rcorr(x = ft1_t, type = corr_type)
    #corr_mat <- as.data.frame(corr_results$r)
    corr_mat <- corr_results$r
    print(head(corr_mat))
  }else{
    ft2_t <- t(ft2)
    ft2_t <- ft2_t[order(row.names(ft2_t)), ] # Ordering by row names
    corr_results <- Hmisc::rcorr(x = ft1_t, y = ft2_t,type = corr_type)
    #corr_mat <- as.data.frame(corr_results$r)
    corr_mat <- corr_results$r
    
    corr_mat <- corr_mat[1:(nrow(ft1)), (nrow(ft1)+1):ncol(corr_mat)]
  }
  # Calculate significance
  if (sig_stars) {
    pval_mat <- corr_results$P
    
    if (pval_adjust) {
      # Flatten the matrix
      pval_vector <- as.vector(pval_mat)
      
      # Apply FDR correction
      fdr_vector <- p.adjust(pval_vector, method = fdr_correction)
      
      # Reshape back into a matrix
      fdr_matrix <- matrix(fdr_vector, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
      rownames(fdr_matrix) <- rownames(pval_mat)
      colnames(fdr_matrix) <- colnames(pval_mat)
      pval_mat <- fdr_matrix
    }
    
    # keep only the need correlation
    if (!is.null(ft2)) {
      pval_mat <- pval_mat[1:(nrow(ft1)), (nrow(ft1)+1):ncol(pval_mat)]
    }
    print("error1")
    #print(head(corr_mat))
    #return(corr_mat)
    
    #print(pval_mat)
    #return(pval_mat)
    
    print(nrow(corr_mat))
    print(ncol(corr_mat))
    print(nrow(pval_mat))
    print(ncol(pval_mat))
    
    # Remove NAs in p-value matrix
    pval_mat[is.na(pval_mat)] <- 1

    print(head(corr_mat))
    print(head(pval_mat))
    
    corrplot::corrplot(corr = corr_mat, p.mat = pval_mat, method = "circle", tl.cex=axis_text_size, sig.level = 0.05, insig = "label_sig", pch.cex = stars_size, tl.srt = 35, type = hm_type)
  }else{
    corrplot::corrplot(corr = corr_mat, method = "circle", tl.cex = axis_text_size, pch.cex = stars_size, type = hm_type)
  }
}


ft_pca_1 <- function(ft, metadata_table, grouping_col, color_palette = NULL, encircle = FALSE, dist_method = "euclidean",
                     pca_scale = TRUE, pca_center = TRUE){
  # transposing feature table
  ft_t <- t(feature_table)
  if (isTRUE(all.equal(rownames(ft_t), rownames(metadata_table)))) {
    print("Sample names in feature table and metadatable are identical :)")
    
    # 1. Calculate Bray-Curtis dissimilarity
    dist_mat <- vegan::vegdist(ft_t, method = dist_method)
    
    # 2. Convert the distance matrix to a matrix object
    dist_mat <- as.matrix(dist_mat)
    
    print(head(dist_mat))
    
    fia_pca <- PCAtools::pca(dist_mat, scale = pca_scale, center = pca_center, metadata = metadata_table, transposed = TRUE)
    
    if (is.null(color_palette)) {
      color_palette = get_palette(nColors = 60)
    }
    
    p2 <- PCAtools::biplot(fia_pca, showLoadings = TRUE, ntopLoadings = 0, lab = NULL, colby = grouping_col,
                           legendPosition = "right", axisLabSize = 8, legendLabSize = 8, legendIconSize = 2, pointSize = 1.5,
                           colkey = color_palette, encircle = encircle, colLegendTitle = "SynCom")
    
    p2
  }else{
    print("Sample names in feature table and metadatable are not identical")
    #print(all.equal(colnames(ft_t),metadata_table$Sample))
  }
}

ft_pca_2 <- function(feature_table, metadata_table, grouping_col = NULL, p_shape = NULL, dist_method = "euclidean"){
  # Compute distance matrix according to dist_method
  dist_matrix <- vegan::vegdist(feature_table, method = dist_method)
  
  # Perform PCoA
  pcoa_results <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  
  pcoa_df <- as.data.frame(pcoa_results$points)
  colnames(pcoa_df) <- c("PC1", "PC2")  # Rename axes
  pcoa_df$Sample <- rownames(feature_table)  # Add sample names to pcoa_df
  
  metadata_table$Sample <- rownames(metadata_table)  # Add sample names to pcoa_df
  
  #print(pcoa_df$Sample == row.names(metadata_table))
  if (identical(pcoa_df$Sample,row.names(metadata_table))) {
    print("Sample names in feature table and metadatable are identical :)")
    
  } else{
    print("Sample names in feature table and metadatable are not identical")
    return()
  }
  
  pcoa_df <- left_join(pcoa_df, metadata_table, by = c("Sample" = "Sample"))
  colour_palette <- get_palette(nColors = 20)
  #print(colour_palette)
  
  if (is.null(grouping_col) && is.null(p_shape)) {
    ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2")
  }else if(!is.null(grouping_col) && is.null(p_shape)){
    print("Plotting with grouping variable")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }else if(!is.null(grouping_col) && !is.null(p_shape)){
    print("Plotting with two grouping variables")
    ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[grouping_col]], shape = .data[[p_shape]])) +
      geom_point(size = 3) +
      theme_minimal() +
      scale_color_manual(values=colour_palette) +
      labs(title = "PCoA Plot",
           x = "PCoA 1",
           y = "PCoA 2",
           color = "Sample Type") +
      theme(legend.position = "right")  +
      guides(color = guide_legend(ncol = 2))
  }
  
}


# Do PCA plot, prev. fia_pca (DELETE)
ft_pca_1 <- function(feature_table, metadata_table, grouping_col, color_palette = NULL, encircle = FALSE){
  # transposing feature table
  ft_t <- t(feature_table)
  if (isTRUE(all.equal(colnames(ft_t),row.names(metadata_table)))) {
    print("Sample names in feature table and metadatable are identical :)")
    
    fia_pca <- PCAtools::pca(ft_t, scale = TRUE, metadata = metadata_table, transposed = FALSE)
    
    if (is.null(color_palette)) {
      color_palette = get_palette(nColors = 60)
    }
    
    p2 <- PCAtools::biplot(fia_pca, showLoadings = TRUE, ntopLoadings = 0, lab = NULL, colby = grouping_col,
                           legendPosition = "right", axisLabSize = 8, legendLabSize = 8, legendIconSize = 2, pointSize = 1.5,
                           colkey = color_palette, encircle = encircle, colLegendTitle = "SynCom")
    
    p2
  }else{
    print("Sample names in feature table and metadatable are not identical")
    #print(all.equal(colnames(ft_t),metadata_table$Sample))
  }
  
}