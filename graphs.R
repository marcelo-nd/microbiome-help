get_palette <- function(nColors = 50){
  return(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                           "#0072B2","brown1", "#CC79A7", "olivedrab3", "rosybrown",
                           "darkorange3", "blueviolet", "darkolivegreen4", "lightskyblue4", "navajowhite4",
                           "purple4", "springgreen4", "firebrick3", "gold3", "cyan3",
                           "plum", "mediumspringgreen", "blue", "yellow", "#053f73",
                           "#e3ae78", "#a23f3f", "#290f76", "#ce7e00", "#386857",
                           "#738564", "#e89d56", "#cd541d", "#1a3a46", "#ffe599",
                           "#583E26", "#A78B71", "#F7C815", "#EC9704", "#9C4A1A",
                           "firebrick2", "#C8D2D1", "#14471E", "#EE9B01", "#DA6A00",
                           "#4B1E19", "#C0587E", "#FC8B5E", "#EA592A", "#FEF4C0")[1:nColors])
  }

barplot_from_feature_table <- function(feature_table){
  # Remove columns (samples) with zero count
  if (ncol(feature_table) > 1) {
    feature_table <- feature_table[, colSums(feature_table != 0) > 0]
  }
  
  # Generate a column with the names of ASVs/OTUs using rownames.
  feature_table["bacteria"] <- row.names(feature_table)
  
  # Gather
  #soil_table_bp <- gather(soil_table_bp, X13114.control.soil.grass.near.BRF:X13114.shade.23.s010 , key = "sample", value = "abundance")
  feature_table <- gather(feature_table, 1:(ncol(feature_table) - 1) , key = "sample", value = "abundance")
  
  color_palette <- get_palette()
  
  ggplot(feature_table, aes(x=sample, y=abundance, fill=bacteria)) + 
    geom_bar(position="fill", stat="identity", show.legend = TRUE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values=color_palette) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
