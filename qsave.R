library(qs)

qsave(all_samples_with_bins$Y3088, "Y3088_data.qs", preset = "high")
qsave(all_samples_with_bins$AC1017_strict, "AC1017_strict_data.qs", preset = "high")
qsave(all_samples_with_bins$AC1017_loose, "AC1017_loose_data.qs", preset = "high")

qsave(consensus, "consensus_data.qs", preset = "high")

