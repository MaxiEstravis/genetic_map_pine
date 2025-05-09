library(qs)

qsave(all_samples_with_bins$Y3088, "Y3088_data.qs", preset = "high")
qsave(all_samples_with_bins$AC1017_strict, "AC1017_strict_data.qs", preset = "high")
qsave(all_samples_with_bins$AC1017_loose, "AC1017_loose_data.qs", preset = "high")

qsave(consensus, "consensus_data.qs", preset = "high")

## i want to save smaller things to work on the laptop

qload("AC1017_strict_data.qs")

LGs_rec_map <- list()
LGs_rec_map$LG1_rec_map <- LGs_progress$LG1_rec_map
LGs_rec_map$LG2_rec_map_noSus <- outliers$LG2_rec_map_noSus
LGs_rec_map$LG3_rec_map <- LGs_progress$LG3_rec_map
LGs_rec_map$LG4_rec_map <- LGs_progress$LG4_rec_map
LGs_rec_map$LG5_rec_map <- LGs_progress$LG5_rec_map
LGs_rec_map$LG6_rec_map <- LGs_progress$LG6_rec_map
LGs_rec_map$LG7_rec_map <- LGs_progress$LG7_rec_map
LGs_rec_map$LG8_rec_map <- LGs_progress$LG8_rec_map
LGs_rec_map$LG9_rec_map <- LGs_progress$LG9_rec_map
LGs_rec_map$LG10_rec_map <- LGs_progress$LG10_rec_map
LGs_rec_map$LG11_rec_map <- LGs_progress$LG11_rec_map
LGs_rec_map$LG12_rec_map <- LGs_progress$LG12_rec_map

qsave(LGs_rec_map, "AC1017_strict_LGs_rec_map.qs", preset = "high")
