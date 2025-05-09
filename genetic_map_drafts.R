library(onemap)

#args <- commandArgs(trailingOnly = TRUE)
#input_file <- args[1]

input_file <- "AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable_non_redundant.raw"

a_Y3088 <- read_onemap(inputfile = input_file)

segreg_test_Y3088 <- test_segregation(a_Y3088)

no_dist_Y3088 <- select_segreg(segreg_test_Y3088, distorted = FALSE, numbers = TRUE)
#dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE)

twopts_Y3088 <- rf_2pts(a_Y3088)

mark_no_dist_Y3088 <- make_seq(twopts_Y3088, c(no_dist_Y3088))

LOD_sug_Y3088 <- suggest_lod(mark_no_dist_Y3088)          # 7.477526

LGs_Y3088 <- group(mark_no_dist_Y3088, LOD=LOD_sug_Y3088)       # 6 groups generated

LGs_35_Y3088 <- group(mark_no_dist_Y3088, LOD=LOD_sug_Y3088, max.rf = 0.35) # based on spruce paper; also 6 groups

LGs_upgma_Y3088 <- group_upgma(mark_no_dist_Y3088, expected.groups = 12, inter = F)     # a different method
#LGs_upgma <- group_upgma(mark_no_dist, expected.groups = 12, inter = T)     # to make the plot interactive
#plot(LGs_upgma)

set_map_fun(type = "kosambi")

##

LG1_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 1)
LG1_upgma_Y3088_rec <- record(LG1_upgma_Y3088, hmm = FALSE)
LG2_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 2)
LG2_upgma_Y3088_rec <- record(LG2_upgma_Y3088, hmm = FALSE)
LG3_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 3)
LG3_upgma_Y3088_rec <- record(LG3_upgma_Y3088, hmm = FALSE)
LG4_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 4)
LG4_upgma_Y3088_rec <- record(LG4_upgma_Y3088, hmm = FALSE)
LG5_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 5)
LG5_upgma_Y3088_rec <- record(LG5_upgma_Y3088, hmm = FALSE)
LG6_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 6)
LG6_upgma_Y3088_rec <- record(LG6_upgma_Y3088, hmm = FALSE)
LG7_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 7)
LG7_upgma_Y3088_rec <- record(LG7_upgma_Y3088, hmm = FALSE)
LG8_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 8)
LG8_upgma_Y3088_rec <- record(LG8_upgma_Y3088, hmm = FALSE)
LG9_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 9)
LG9_upgma_Y3088_rec <- record(LG9_upgma_Y3088, hmm = FALSE)
LG10_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 10)
LG10_upgma_Y3088_rec <- record(LG10_upgma_Y3088, hmm = FALSE)
LG11_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 11)
LG11_upgma_Y3088_rec <- record(LG11_upgma_Y3088, hmm = FALSE)
LG12_upgma_Y3088 <- make_seq(LGs_upgma_Y3088, 12)
LG12_upgma_Y3088_rec <- record(LG12_upgma_Y3088, hmm = FALSE)

rf_graph_table(LG1_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG1_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG1_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG1_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG2_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG2_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG2_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG2_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG3_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG3_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG3_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG3_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG4_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG4_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG4_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG4_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG5_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG5_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG5_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG5_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG6_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG6_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG6_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG6_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG7_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG7_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG7_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG7_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG8_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG8_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG8_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG8_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG9_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG9_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG9_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG9_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG10_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG10_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG10_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG10_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG11_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG11_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG11_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG11_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG12_upgma_Y3088_rec, graph.LOD = TRUE)
ggsave("rf_graph_table_LG12_upgma_Y3088_rec_LOD.pdf", width = 210, height = 210, units = "mm")
rf_graph_table(LG12_upgma_Y3088_rec, graph.LOD = FALSE)
ggsave("rf_graph_table_LG12_upgma_Y3088_rec.pdf", width = 210, height = 210, units = "mm")

### i want to see the chromosome distribution of the markers assigned to 6 groups

vcf_data_Y3088 <- read_tsv("AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable.vcf",
                     comment = "#", col_names = FALSE)

twelve_chromosomes <- paste0("PS_chr", sprintf("%02d", 1:12))

# Create a function to process each LG group
process_LG <- function(LG_num, LGs, seg, vcf_data, twelve_chromosomes) {
  LG_seq <- make_seq(LGs, LG_num)
  markers <- seg$Marker[LG_seq$seq.num]
  
  chrom_df <- vcf_data %>%
    filter(X3 %in% markers) %>%
    count(X1, name = "Counts") %>%
    rename(CHROM = X1) %>%
    mutate(CHROM = if_else(!CHROM %in% twelve_chromosomes, "rest", CHROM)) %>%
    group_by(CHROM) %>%
    summarise(Counts = sum(Counts), .groups = "drop") %>%
    complete(CHROM = c(twelve_chromosomes, "rest"), fill = list(Counts = 0)) %>%
    arrange(factor(CHROM, levels = c(twelve_chromosomes, "rest"))) %>%
    mutate(Group = paste0("LG", LG_num))  # Add LG identifier
  
  return(chrom_df)
}

# Process all LGs and combine into a single dataframe
LG_numbers <- 1:6
chrom_df_all_Y3088_35 <- bind_rows(lapply(LG_numbers, function(LG) process_LG(LG, LGs_35_Y3088, segreg_test_Y3088, vcf_data_Y3088, twelve_chromosomes)))
LG_numbers <- 1:6
chrom_df_all_Y3088 <- bind_rows(lapply(LG_numbers, function(LG) process_LG(LG, LGs_Y3088, segreg_test_Y3088, vcf_data_Y3088, twelve_chromosomes)))
LG_numbers <- 1:12
chrom_df_all_Y3088_upgma <- bind_rows(lapply(LG_numbers, function(LG) process_LG(LG, LGs_upgma_Y3088, segreg_test_Y3088, vcf_data_Y3088, twelve_chromosomes)))
LG_numbers <- 1:13
chrom_df_all_Y3088_upgma_13 <- bind_rows(lapply(LG_numbers, function(LG) process_LG(LG, LGs_upgma_Y3088_13, segreg_test_Y3088, vcf_data_Y3088, twelve_chromosomes)))
LG_numbers <- 1:14
chrom_df_all_Y3088_upgma_14 <- bind_rows(lapply(LG_numbers, function(LG) process_LG(LG, LGs_upgma_Y3088_14, segreg_test_Y3088, vcf_data_Y3088, twelve_chromosomes)))

# Plot stacked bar chart
ggplot(chrom_df_all_Y3088_upgma, aes(x = CHROM, y = Counts, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars by Group
  geom_text(aes(label = Group), 
            position = position_stack(vjust = 0.5),  # Centers text inside each stacked bar
            color = "white", size = 2) +  # Adjust color/size for readability
  theme_minimal() +
  labs(title = "Y3088 LG Groups (UPGMA)", 
       x = "Chromosomes", 
       y = "Counts", 
       fill = "LG Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), legend.position = "none")

ggsave("LGs_Y3088_upgma.pdf")

ggplot(chrom_df_all_Y3088, aes(x = CHROM, y = Counts, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars by Group
  geom_text(aes(label = Group), 
            position = position_stack(vjust = 0.5),  # Centers text inside each stacked bar
            color = "white", size = 2) +  # Adjust color/size for readability
  theme_minimal() +
 labs(title = "Y3088 LG Groups (LOD_sug)", 
      x = "Chromosomes", 
      y = "Counts", 
      fill = "LG Group") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1), 
       plot.title = element_text(hjust = 0.5), legend.position = "none")

ggsave("LGs_Y3088_LOD_sug.pdf")

ggplot(chrom_df_all_Y3088_35, aes(x = CHROM, y = Counts, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars by Group
  geom_text(aes(label = Group), 
            position = position_stack(vjust = 0.5),  # Centers text inside each stacked bar
            color = "white", size = 2) +  # Adjust color/size for readability
  theme_minimal() +
  labs(title = "Y3088 LG Groups (LOD_sug, rf=0.35)", 
       x = "Chromosomes", 
       y = "Counts", 
       fill = "LG Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), legend.position = "none")

ggsave("LGs_Y3088_LOD_sug_rf035.pdf")

## to select which LGs to plot

LG_to_plot <- NULL
chrom_df_all <- chrom_df_all_Y3088_upgma_14


# Apply filtering only if a specific LG is chosen
plot_data <- if (is.null(LG_to_plot)) chrom_df_all else filter(chrom_df_all, Group %in% LG_to_plot)

ggplot(plot_data, aes(x = CHROM, y = Counts, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +  
  geom_text(aes(label = Group), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 2.5) +  
  theme_minimal() +
  labs(title = if (is.null(LG_to_plot)) "Stacked Bar Plot of All LG Groups" else paste("Chromosome Distribution for", paste(LG_to_plot, collapse = ", ")), 
       x = "Chromosomes", 
       y = "Counts", 
       fill = "LG Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # Show legend only for all LGs











