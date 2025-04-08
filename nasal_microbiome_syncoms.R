# Load scripts
source("C:/Users/marce/Documents/GitHub/microbiome-help/feature_table_wrangling.R")

### ... ###

otu_table_sctp <- read.csv("E:/SequencingData/SynCom100/TheChampions/emu_results/otu_table.csv", row.names=1)

otu_table_sctp_sorted <- sort_nanopore_table_by_barcodes(df = otu_table_sctp,
                                                         new_names = c("SC4_T1_R1", "SC4_T1_R2", "SC4_T1_R3",
                                                                       "SC4_T2_R1", "SC4_T2_R2", "SC4_T2_R3",
                                                                       "SC4_T3_R1", "SC4_T3_R2", "SC4_T3_R3",
                                                                       "SC4_TF_R1", "SC4_TF_R2", "SC4_TF_R3",
                                                                       "SC7_T1_R1", "SC7_T1_R2", "SC7_T1_R3",
                                                                       "SC7_T2_R1", "SC7_T2_R2", "SC7_T2_R3",
                                                                       "SC7_T3_R1", "SC7_T3_R2", "SC7_T3_R3",
                                                                       "SC7_TF_R1", "SC7_TF_R2", "SC7_TF_R3",
                                                                       "SC9_T1_R1", "SC9_T1_R2", "SC9_T1_R3",
                                                                       "SC9_T2_R1", "SC9_T2_R2", "SC9_T2_R3",
                                                                       "SC9_T3_R1", "SC9_T3_R2", "SC9_T3_R3",
                                                                       "SC9_TF_R1", "SC9_TF_R2", "SC9_TF_R3",
                                                                       "SC10_T1_R1", "SC10_T1_R2", "SC10_T1_R3",
                                                                       "SC10_T2_R1", "SC10_T2_R2", "SC10_T2_R3",
                                                                       "SC10_T3_R1", "SC10_T3_R2", "SC10_T3_R3",
                                                                       "SC10_TF_R1", "SC10_TF_R2", "SC10_TF_R3",
                                                                       "SC11_T1_R1", "SC11_T1_R2", "SC11_T1_R3",
                                                                       "SC11_T2_R1", "SC11_T2_R2", "SC11_T2_R3",
                                                                       "SC11_T3_R1", "SC11_T3_R2", "SC11_T3_R3",
                                                                       "SC11_TF_R1", "SC11_TF_R2", "SC11_TF_R3",
                                                                       "SC12_T1_R1", "SC12_T1_R2", "SC12_T1_R3",
                                                                       "SC12_T2_R1", "SC12_T2_R2", "SC12_T2_R3",
                                                                       "SC12_T3_R1", "SC12_T3_R2", "SC12_T3_R3",
                                                                       "SC12_TF_R1", "SC12_TF_R2", "SC12_TF_R3",
                                                                       "SC13_T1_R1", "SC13_T1_R2", "SC13_T1_R3",
                                                                       "SC13_T2_R1", "SC13_T2_R2", "SC13_T2_R3",
                                                                       "SC13_T3_R1", "SC13_T3_R2", "SC13_T3_R3",
                                                                       "SC13_TF_R1", "SC13_TF_R2", "SC13_TF_R3",
                                                                       "SC14_T1_R1", "SC14_T1_R2", "SC14_T1_R3",
                                                                       "SC14_T2_R1", "SC14_T2_R2", "SC14_T2_R3",
                                                                       "SC14_T3_R1", "SC14_T3_R2", "SC14_T3_R3",
                                                                       "SC14_TF_R1", "SC14_TF_R2", "SC14_TF_R3",
                                                                       "SC20_T1_R1", "SC20_T1_R2", "SC20_T1_R3",
                                                                       "SC20_T2_R1", "SC20_T2_R2", "SC20_T2_R3",
                                                                       "SC20_T3_R1", "SC20_T3_R2", "SC20_T3_R3",
                                                                       "SC20_TF_R1", "SC20_TF_R2", "SC20_TF_R3",
                                                                       "SC23_T1_R1", "SC23_T1_R2", "SC23_T1_R3",
                                                                       "SC23_T2_R1", "SC23_T2_R2", "SC23_T2_R3",
                                                                       "SC23_T3_R1", "SC23_T3_R2", "SC23_T3_R3",
                                                                       "SC23_TF_R1", "SC23_TF_R2", "SC23_TF_R3",
                                                                       "SC24_T1_R1", "SC24_T1_R2", "SC24_T1_R3",
                                                                       "SC24_T2_R1", "SC24_T2_R2", "SC24_T2_R3",
                                                                       "SC24_T3_R1", "SC24_T3_R2", "SC24_T3_R3",
                                                                       "SC24_TF_R1", "SC24_TF_R2", "SC24_TF_R3",
                                                                       "SC25_T1_R1", "SC25_T1_R2", "SC25_T1_R3",
                                                                       "SC25_T2_R1", "SC25_T2_R2", "SC25_T2_R3",
                                                                       "SC25_T3_R1", "SC25_T3_R2", "SC25_T3_R3",
                                                                       "SC25_TF_R1", "SC25_TF_R2", "SC25_TF_R3",
                                                                       "SC26_T1_R1", "SC26_T1_R2", "SC26_T1_R3",
                                                                       "SC26_T2_R1", "SC26_T2_R2", "SC26_T2_R3",
                                                                       "SC26_T3_R1", "SC26_T3_R2", "SC26_T3_R3",
                                                                       "SC26_TF_R1", "SC26_TF_R2", "SC26_TF_R3",
                                                                       "SC28_T1_R1", "SC28_T1_R2", "SC28_T1_R3",
                                                                       "SC28_T2_R1", "SC28_T2_R2", "SC28_T2_R3",
                                                                       "SC28_T3_R1", "SC28_T3_R2", "SC28_T3_R3",
                                                                       "SC28_TF_R1", "SC28_TF_R2", "SC28_TF_R3",
                                                                       "SC32_T1_R1", "SC32_T1_R2", "SC32_T1_R3",
                                                                       "SC32_T2_R1", "SC32_T2_R2", "SC32_T2_R3",
                                                                       "SC32_T3_R1", "SC32_T3_R2", "SC32_T3_R3",
                                                                       "SC32_TF_R1", "SC32_TF_R2", "SC32_TF_R3",
                                                                       "SC36_T1_R1", "SC36_T1_R2", "SC36_T1_R3",
                                                                       "SC36_T2_R1", "SC36_T2_R2", "SC36_T2_R3",
                                                                       "SC36_T3_R1", "SC36_T3_R2", "SC36_T3_R3",
                                                                       "SC36_TF_R1", "SC36_TF_R2", "SC36_TF_R3",
                                                                       "SC42_T1_R1", "SC42_T1_R2", "SC42_T1_R3",
                                                                       "SC42_T2_R1", "SC42_T2_R2", "SC42_T2_R3",
                                                                       "SC42_T3_R1", "SC42_T3_R2", "SC42_T3_R3",
                                                                       "SC42_TF_R1", "SC42_TF_R2", "SC42_TF_R3",
                                                                       "SC43_T1_R1", "SC43_T1_R2", "SC43_T1_R3",
                                                                       "SC43_T2_R1", "SC43_T2_R2", "SC43_T2_R3",
                                                                       "SC43_T3_R1", "SC43_T3_R2", "SC43_T3_R3",
                                                                       "SC43_TF_R1", "SC43_TF_R2", "SC43_TF_R3",
                                                                       "SC47_T1_R1", "SC47_T1_R2", "SC47_T1_R3",
                                                                       "SC47_T2_R1", "SC47_T2_R2", "SC47_T2_R3",
                                                                       "SC47_T3_R1", "SC47_T3_R2", "SC47_T3_R3",
                                                                       "SC47_TF_R1", "SC47_TF_R2", "SC47_TF_R3",
                                                                       "SC53_T1_R1", "SC53_T1_R2", "SC53_T1_R3",
                                                                       "SC53_T2_R1", "SC53_T2_R2", "SC53_T2_R3",
                                                                       "SC53_T3_R1", "SC53_T3_R2", "SC53_T3_R3",
                                                                       "SC53_TF_R1", "SC53_TF_R2", "SC53_TF_R3"))

# Remove species with no counts
otu_table_sctp_filt <- filter_otus_by_counts_col_counts(otu_table_sctp_sorted,
                                                        min_count = 10,
                                                        col_number = 1)

# Remove Anaerococcus octavius, it did not grow on any of the SCs
otu_table_sctp_filt <- otu_table_sctp_filt[!rownames(otu_table_sctp_filt) %in% "Anaerococcus octavius", ]

# Remove Cutibacterium acnes, it did not grow on any of the SCs
otu_table_sctp_filt <- otu_table_sctp_filt[!rownames(otu_table_sctp_filt) %in% "Cutibacterium acnes", ]

# Remove unnassigned reads
otu_table_sctp_filt <- otu_table_sctp_filt[-10,]

# Convert OTU Table to strain level table.

strain_data <- readxl::read_excel(path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/SOPs/1_Nose (HMP) Strains.xlsx", sheet = "SynCom100_2", range = "A1:BA32", col_names = TRUE)

strain_ft <- merge_abundance_by_strain(otu_table_sctp_filt, strain_data)

sc4 <- strain_ft[1:12]
sc7 <- strain_ft[13:24]
sc9 <- strain_ft[25:36]

barplot_w_strain_data2(sc4)
barplot_w_strain_data2(sc7)
barplot_w_strain_data2(sc9)

colours_vec <- c("gold3", "#053f73", "blueviolet", "#CC79A7","#6279B8",
                 "lightblue1", "brown1", "olivedrab3", "darkorange3", "springgreen4")

#colour_palette = colours_vec,

barplot_from_feature_tables2(feature_tables = list(sc4, sc7, sc9),
                            experiments_names = c("SC4", "SC7", "SC9"),
                            x_axis_title_size = 9, x_axis_text_size = 5,
                            y_axis_title_size = 9, y_axis_text_size = 5,
                            legend_pos = "bottom", legend_cols = 2,
                            legend_title_size = 9, legend_text_size = 7,
                            legend_key_size = 0.3)

