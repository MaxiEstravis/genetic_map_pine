#chr AC1017_strict  Y3088
#PS_chr01  LG1 LG1
#PS_chr02  LG4 LG2
#PS_chr03  LG5 LG3
#PS_chr04  LG3 LG4
#PS_chr05  LG8 LG6
#PS_chr06  LG9 LG8
#PS_chr07  LG10  LG9
#PS_chr08  LG12  LG7
#PS_chr09  LG7 LG10
#PS_chr10  LG6 LG11
#PS_chr11  LG2 LG5
#PS_chr12  LG11  LG12

library(LPmerge)
library(dplyr)

pairings <- list(c("1", "1"), c("4", "2"), c("5", "3"), c("3", "4"), c("8", "6"), c("9", "8"),
                 c("10", "9"), c("12", "7"), c("7", "10"), c("6", "11"), c("2", "5"), c("11", "12"))

map.names <- c("AC1017", "Y3088")
pop.size <- c(433, 1091)

consensus <- list()

consensus$Maps_all <- list()
consensus$weighted_all <- list()
consensus$output_weighted_all <- list()
consensus$rmse_df_all <- list()

for (j in seq_along(pairings)) {
  t <- pairings[[j]]
  
  Maps <- list()
  filenames <- character(2)
  for (i in 1:2) {
    filenames[i] <- paste0(map.names[i], "_LG", t[i], "_def.tsv")
    Maps[[i]] <- read.csv(filenames[i], header = TRUE, sep = "\t")
  }
  names(Maps) <- map.names
  
  consensus$Maps_all[[j]] <- list(maps = Maps, filenames = filenames)
  
  output_weighted <- capture.output({
    weighted <- LPmerge(Maps, max.interval = 1:10, weights = pop.size)
  })
  
  consensus$weighted_all[[j]] <- weighted
  consensus$output_weighted_all[[j]] <- output_weighted
  
  rmse_start_indices <- grep("map\\s+RMSE", output_weighted)
  rmse_blocks <- lapply(rmse_start_indices, function(i) output_weighted[(i + 1):(i + 4)])
  
  rmse_df <- do.call(rbind, lapply(seq_along(rmse_blocks), function(i) {
    lines <- rmse_blocks[[i]]
    split_lines <- strsplit(lines, "\\s+")
    labels <- sapply(split_lines, function(x) x[2])
    values <- as.numeric(sapply(split_lines, function(x) x[3]))
    df_row <- as.data.frame(t(values))
    names(df_row) <- labels
    df_row$max.interval <- i
    df_row[, c("max.interval", "AC1017", "Y3088", "mean", "sd")]
  }))

  consensus$rmse_df_all[[j]] <- rmse_df
  
  message("Finished pair ", j, ": ", filenames[1], " + ", filenames[2])
}

for (i in seq_along(consensus$rmse_df_all)) {
  df <- consensus$rmse_df_all[[i]]
  
  df$mean <- as.numeric(as.character(df$mean))
  df$sd <- as.numeric(as.character(df$sd))
  
  min_mean <- min(df$mean, na.rm = TRUE)
  df_min_mean <- df[df$mean == min_mean, ]
  
  min_sd <- min(df_min_mean$sd, na.rm = TRUE)
  best_rows <- df_min_mean[df_min_mean$sd == min_sd, ]
  
  cat("\n=== Best row(s) for iteration", i, "===\n")
  print(best_rows)
}

#=== Best row(s) for iteration 1 ===
#  max.interval AC1017 Y3088 mean   sd
#9            9   4.38  0.03  2.2 3.07

#=== Best row(s) for iteration 2 ===
#  max.interval AC1017 Y3088 mean   sd
#6            6   3.31  0.49  1.9 1.99

#=== Best row(s) for iteration 3 ===
#  max.interval AC1017 Y3088 mean   sd
#10           10   2.18   1.1 1.64 0.76

#=== Best row(s) for iteration 4 ===
#  max.interval AC1017 Y3088 mean   sd
#8            8   7.45  2.64 5.04 3.41

#=== Best row(s) for iteration 5 ===
#  max.interval AC1017 Y3088 mean   sd
#1            1   5.02  1.24 3.13 2.67

#=== Best row(s) for iteration 6 ===
#  max.interval AC1017 Y3088 mean   sd
#2            2   2.72  0.19 1.45 1.79

#=== Best row(s) for iteration 7 ===
#  max.interval AC1017 Y3088 mean   sd
#3            3   12.8  3.46 8.13 6.61

#=== Best row(s) for iteration 8 ===
#  max.interval AC1017 Y3088 mean   sd
#9            9   1.96  0.05 1.01 1.35

#=== Best row(s) for iteration 9 ===
#  max.interval AC1017 Y3088 mean   sd
#9            9   3.03  0.98 2.01 1.45

#=== Best row(s) for iteration 10 ===
#  max.interval AC1017 Y3088 mean   sd
#6            6   3.64  0.02 1.83 2.55
#7            7   3.64  0.03 1.83 2.55
#8            8   3.63  0.03 1.83 2.55
#9            9   3.63  0.03 1.83 2.55

#=== Best row(s) for iteration 11 ===
#  max.interval AC1017 Y3088 mean   sd
#5            5   4.76  0.71 2.74 2.86

#=== Best row(s) for iteration 12 ===
#  max.interval AC1017 Y3088 mean   sd
#4            4   3.29  1.25 2.27 1.45

consensus$consensus_map <- list()
for (i in 1:12) {
  LG_name <- paste0("LG", i)
  consensus$consensus_map[[LG_name]] <- list()
}

consensus$consensus_map$LG1$map <- consensus$weighted_all[[1]][9]
consensus$consensus_map$LG1$filenames <- consensus$Maps_all[[1]]$filenames

consensus$consensus_map$LG2$map <- consensus$weighted_all[[2]][6]
consensus$consensus_map$LG2$filenames <- consensus$Maps_all[[2]]$filenames

consensus$consensus_map$LG3$map <- consensus$weighted_all[[3]][10]
consensus$consensus_map$LG3$filenames <- consensus$Maps_all[[3]]$filenames

consensus$consensus_map$LG4$map <- consensus$weighted_all[[4]][8]
consensus$consensus_map$LG4$filenames <- consensus$Maps_all[[4]]$filenames

consensus$consensus_map$LG5$map <- consensus$weighted_all[[5]][1]
consensus$consensus_map$LG5$filenames <- consensus$Maps_all[[5]]$filenames

consensus$consensus_map$LG6$map <- consensus$weighted_all[[6]][2]
consensus$consensus_map$LG6$filenames <- consensus$Maps_all[[6]]$filenames

consensus$consensus_map$LG7$map <- consensus$weighted_all[[7]][3]
consensus$consensus_map$LG7$filenames <- consensus$Maps_all[[7]]$filenames

consensus$consensus_map$LG8$map <- consensus$weighted_all[[8]][9]
consensus$consensus_map$LG8$filenames <- consensus$Maps_all[[8]]$filenames

consensus$consensus_map$LG9$map <- consensus$weighted_all[[9]][9]
consensus$consensus_map$LG9$filenames <- consensus$Maps_all[[9]]$filenames

consensus$consensus_map$LG10$map <- consensus$weighted_all[[10]][6]
consensus$consensus_map$LG10$filenames <- consensus$Maps_all[[10]]$filenames

consensus$consensus_map$LG11$map <- consensus$weighted_all[[11]][5]
consensus$consensus_map$LG11$filenames <- consensus$Maps_all[[11]]$filenames

consensus$consensus_map$LG12$map <- consensus$weighted_all[[12]][4]
consensus$consensus_map$LG12$filenames <- consensus$Maps_all[[12]]$filenames

# now to save...

for (i in 1:12) {
  lg <- paste0("LG", i)
  filename <- paste0("consensus_LG", i, ".tsv")
  write.table(consensus$consensus_map[[lg]]$map, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  if (file.exists(filename)) {
    cat("File", filename, "written successfully!\n")
  } else {
    cat("There was an error writing the file", filename, "\n")
  }
}


