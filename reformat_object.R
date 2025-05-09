### this was useful when i had done the pipeline without the sub-lists. now it's unnecesary but i keep the script for reference

sample_name <- "AC1017_loose"

all_samples_with_bins[[sample_name]]$prelim <- list()

elements_prelim <- c("a", "a_noInfo", "bins", "bins_noInfo", "bins_data", "bins_data_noInfo", "segreg_test", "no_dist", "twopts", "mark_no_dist", "LOD_sug")

for (elem in elements_prelim) {
  all_samples_with_bins[[sample_name]]$prelim[[elem]] <- all_samples_with_bins[[sample_name]][[elem]]
  all_samples_with_bins[[sample_name]][[elem]] <- NULL
  message("Moved ", elem, " to ", sample_name, "$prelim and removed it from the top level.")
}

all_samples_with_bins[[sample_name]]$LGs_progress <- list()

all_samples_with_bins[[sample_name]]$LGs_progress$LGs_upgma <- all_samples_with_bins[[sample_name]]$LGs_upgma
all_samples_with_bins[[sample_name]]$LGs_upgma <- NULL

for (i in 1:12) {
  lg_name <- paste0("LG", i)
  all_samples_with_bins[[sample_name]]$LGs_progress[[lg_name]] <- all_samples_with_bins[[sample_name]][[lg_name]]
  all_samples_with_bins[[sample_name]][[lg_name]] <- NULL
  message("Moved ", lg_name, " to ", sample_name, "$LGs_progress and removed it from the top level.")
  
  lg_rec_name <- paste0("LG", i, "_rec")
  all_samples_with_bins[[sample_name]]$LGs_progress[[lg_rec_name]] <- all_samples_with_bins[[sample_name]][[lg_rec_name]]
  all_samples_with_bins[[sample_name]][[lg_rec_name]] <- NULL
  message("Moved ", lg_rec_name, " to ", sample_name, "$LGs_progress and removed it from the top level.")
  
  lg_rec_map_name <- paste0("LG", i, "_rec_map")
  all_samples_with_bins[[sample_name]]$LGs_progress[[lg_rec_map_name]] <- all_samples_with_bins[[sample_name]][[lg_rec_map_name]]
  all_samples_with_bins[[sample_name]][[lg_rec_map_name]] <- NULL
  message("Moved ", lg_rec_map_name, " to ", sample_name, "$LGs_progress and removed it from the top level.")
  
  lg_rec_map_redu_name <- paste0("LG", i, "_rec_map_redu")
  all_samples_with_bins[[sample_name]]$LGs_progress[[lg_rec_map_redu_name]] <- all_samples_with_bins[[sample_name]][[lg_rec_map_redu_name]]
  all_samples_with_bins[[sample_name]][[lg_rec_map_redu_name]] <- NULL
  message("Moved ", lg_rec_map_redu_name, " to ", sample_name, "$LGs_progress and removed it from the top level.")
  
  lg_rec_map_redu_pos_name <- paste0("LG", i, "_rec_map_redu_marker_positions")
  all_samples_with_bins[[sample_name]]$LGs_progress[[lg_rec_map_redu_pos_name]] <- all_samples_with_bins[[sample_name]][[lg_rec_map_redu_pos_name]]
  all_samples_with_bins[[sample_name]][[lg_rec_map_redu_pos_name]] <- NULL
  message("Moved ", lg_rec_map_redu_pos_name, " to ", sample_name, "$LGs_progress and removed it from the top level.")
}

## continues in genetic_map_outliers.R