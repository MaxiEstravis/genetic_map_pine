library(onemap)
library(tidyverse)

process_sample_with_bins <- function(sample_id) {
  sample_objects <- list()  # Store all objects for this sample
  sample_objects$prelim <- list()
  sample_objects$LGs_progress <- list()
  
  input_file <- paste0("AxiomGT1.calls_", sample_id, "_HetMarkedMissing_newCoords_newHeader_sorted_variable.raw")
  input_file_noInfo <- paste0("AxiomGT1.calls_", sample_id, "_HetMarkedMissing_newCoords_newHeader_sorted_variable_noInfo.raw")
  
  sample_objects$prelim$a <- read_onemap(inputfile = input_file)
  sample_objects$prelim$a_noInfo <- read_onemap(inputfile = input_file_noInfo)
  
  sample_objects$prelim$bins <- find_bins(sample_objects$prelim$a, exact = FALSE)
  cat("Bins with info found for", sample_id, "\n")
  sample_objects$prelim$bins_noInfo <- find_bins(sample_objects$prelim$a_noInfo, exact = FALSE)
  cat("Bins found for", sample_id, "\n")
  sample_objects$prelim$bins_data <- create_data_bins(sample_objects$prelim$a, sample_objects$prelim$bins)
  cat("Data bins with info created for", sample_id, "\n")
  sample_objects$prelim$bins_data_noInfo <- create_data_bins(sample_objects$prelim$a_noInfo, sample_objects$prelim$bins_noInfo)
  cat("Data bins created for", sample_id, "\n")
  #sample_objects$prelim$segreg_test <- test_segregation(sample_objects$prelim$bins_data)
  sample_objects$prelim$segreg_test <- test_segregation(sample_objects$prelim$bins_data_noInfo)
  cat("Segregation test completed for", sample_id, "\n")
  sample_objects$prelim$no_dist <- select_segreg(sample_objects$prelim$segreg_test, distorted = FALSE, numbers = TRUE)
  #sample_objects$prelim$twopts <- rf_2pts(sample_objects$prelim$bins_data)
  sample_objects$prelim$twopts <- rf_2pts(sample_objects$prelim$bins_data_noInfo)
  cat("Two-point recombination fraction calculated for", sample_id, "\n")
  sample_objects$prelim$mark_no_dist <- make_seq(sample_objects$prelim$twopts, sample_objects$prelim$no_dist)
  sample_objects$prelim$LOD_sug <- suggest_lod(sample_objects$prelim$mark_no_dist)
  sample_objects$LGs_progress$LGs_upgma <- group_upgma(sample_objects$prelim$mark_no_dist, expected.groups = 12, inter = FALSE)
  cat("Markers grouped using UPGMA for", sample_id, "\n")
  
  set_map_fun(type = "kosambi")
  
  for (i in 1:12) {
    LG_name <- paste0("LG", i)
    LG_rec_name <- paste0("LG", i, "_rec")
    LG_rec_map_name <- paste0("LG", i, "_rec_map")
    LG_rec_map_redu_name <- paste0("LG", i, "_rec_map_redu")
    
    sample_objects$LGs_progress[[LG_name]] <- make_seq(sample_objects$LGs_progress$LGs_upgma, i)
    cat("Created sequence object", LG_name, "for", sample_id, "\n")
    sample_objects$LGs_progress[[LG_rec_name]] <- record(sample_objects$LGs_progress[[LG_name]], hmm = FALSE)
    cat("Created ordered object", LG_rec_name, "for", sample_id, "\n")
    sample_objects$LGs_progress[[LG_rec_map_name]] <- onemap::map(sample_objects$LGs_progress[[LG_rec_name]])
    cat("Created mapped object", LG_rec_map_name, "for", sample_id, "\n")
    sample_objects$LGs_progress[[LG_rec_map_redu_name]] <- try({add_redundants(sample_objects$LGs_progress[[LG_rec_map_name]],
                                                                  sample_objects$prelim$a_noInfo,
                                                                  sample_objects$prelim$bins_noInfo)})
    cat("Added redundant markers in", LG_rec_map_redu_name, "for", sample_id, "\n")
    
    rf_graph_table(sample_objects$LGs_progress[[LG_rec_map_name]], mrk.axis = "none")
    ggsave(filename = paste0("rf_graph_plots/rf_graph_table_LG", i, "_", sample_id, "_upgma_rec_map.pdf"), width = 210, height = 210, units = "mm")
    rf_graph_table(sample_objects$LGs_progress[[LG_rec_map_name]], mrk.axis = "none", graph.LOD = TRUE)
    ggsave(filename = paste0("rf_graph_plots/rf_graph_table_LG", i, "_", sample_id, "_upgma_rec_map_LOD.pdf"), width = 210, height = 210, units = "mm")
  }
  
  LG_names <- paste0("LG", 1:12, "_rec_map")
  LG_redu_names <- paste0("LG", 1:12, "_rec_map_redu")
  LG_combined <- as.vector(rbind(LG_names, LG_redu_names))
  group_names <- as.vector(rbind(paste0("LG", 1:12), paste0("LG", 1:12, "redu")))
  
  draw_map2(mget(LG_combined, envir = as.environment(sample_objects$LGs_progress)), 
            group.names = group_names, 
            output = paste0("LGs_", sample_id, "_UPGMA_raw.pdf"))
  
  return(sample_objects)
}

sample_ids <- c("AC1017_strict", "AC1017_loose", "Y3088")
all_samples_with_bins <- list()  # Store all results

for (sample in sample_ids) {
  all_samples_with_bins[[sample]] <- process_sample_with_bins(sample)
}


#### inspect
for (sample in sample_ids) {
  for (i in 1:12) {
    lg_name <- paste0("LG", i)
    cat(paste0(sample, "_LG", i, ":"))
    print(table(all_samples_with_bins[[sample]]$prelim$bins_data$CHROM[all_samples_with_bins[[sample]]$LGs_progress[[lg_name]]$seq.num]))
  }
}

# add marker positions (could be integrated in the big pipeline...)
for (sample in sample_ids) {
  for (i in 1:12) {
    lg_name <- paste0("LG", i, "_rec_map_redu")
    filename <- paste0(sample, "_", lg_name, ".tsv")
    marker_positions <- data.frame(
      id = colnames(all_samples_with_bins[[sample]]$LGs_progress[[lg_name]]$data.name$geno)[all_samples_with_bins[[sample]]$LGs_progress[[lg_name]]$seq.num], 
      pos = c(0, cumsum(kosambi(all_samples_with_bins[[sample]]$LGs_progress[[lg_name]]$seq.rf)))
    )
    gaps <- diff(marker_positions$pos)
    max_gap_index <- which.max(gaps)
    marker1 <- marker_positions$id[max_gap_index]
    marker2 <- marker_positions$id[max_gap_index + 1]
    max_gap <- gaps[max_gap_index]
    cat(sample, "largest gap in", lg_name, ":", max_gap, "between markers", marker1, "and", marker2, "\n")
    
    all_samples_with_bins[[sample]]$LGs_progress[[paste0(lg_name, "_marker_positions")]] <- marker_positions
    write.table(marker_positions, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
    if (file.exists(filename)) {
      cat("File", filename, "written successfully!\n")
    } else {
      cat("There was an error writing the file", filename, "\n")
    }
  }
}


# with the following line we can know all
for (sample in sample_ids) {
  nind <- all_samples_with_bins[[sample]]$prelim$bins$info$n.ind
  nmar <- all_samples_with_bins[[sample]]$prelim$bins$info$n.mar
  nbin <- length(all_samples_with_bins[[sample]]$prelim$bins$bins)
  nlgs <- sum(sapply(1:12, function(i) length(all_samples_with_bins[[sample]]$LGs_progress[[paste0("LG", i, "_rec_map")]]$seq.num)))
  nredu <- sum(sapply(1:12, function(i) length(all_samples_with_bins[[sample]]$LGs_progress[[paste0("LG", i, "_rec_map_redu_marker_positions")]][,1])))
  p <- (nredu * 100) / nmar
  cat(sample, "summary\n")
  cat("Number of individuals:", nind, "\n")
  cat("Number of markers in original dataset:", nmar, "\n")
  cat("Number of bins found:", nbin, "\n")
  cat("Number of non-redundant markers in LGs:", nlgs, "\n")
  cat("Number of total markers in LGs:", nredu, "\n")
  cat("Percentage of original markers used in LGs:", p, "\n\n")
}

#AC1017_strict summary
#Number of individuals: 433 
#Number of markers in original dataset: 11646 
#Number of bins found: 3802 
#Number of non-redundant markers in LGs: 2800 
#Number of total markers in LGs: 9834 
#Percentage of original markers used in LGs: 84.44101 

#AC1017_loose summary
#Number of individuals: 816 
#Number of markers in original dataset: 11220 
#Number of bins found: 4575 
#Number of non-redundant markers in LGs: 2042 
#Number of total markers in LGs: 6514 
#Percentage of original markers used in LGs: 58.05704 

#Y3088 summary
#Number of individuals: 1091 
#Number of markers in original dataset: 12851 
#Number of bins found: 5556 
#Number of non-redundant markers in LGs: 4754 
#Number of total markers in LGs: 11703 
#Percentage of original markers used in LGs: 91.06684


## continues in AC1017_loose_vs_strict_venn.R
