library(qs)

# Load your data
## same then for AC1017
qload("Y3088_data.qs")

for (lg in 1:12) {
  lg_name <- paste0("LG", lg)
  LG_names <- colnames(LGs_progress[[lg_name]]$data.name$geno)[LGs_progress[[lg_name]]$seq.num]
  
  for (i in 1:100) {
    sampled_names <- sample(LG_names, 100)
    file_name <- sprintf("Y3088_%s_r%d_names.txt", lg_name, i)
    write.table(sampled_names, file = file_name, sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("Generated", file_name, "\n")
  }
}

## bash
# RAW_FILE="AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable_noInfo.raw"
# for lg in {1..12}; do for i in {1..100}; do names_file="Y3088_LG${lg}_r${i}_names.txt"; raw_file="Y3088_LG${lg}_r${i}_noInfo.raw"; echo 'data type outcross' > "$raw_file"; echo '1091 100 0 0 0' >> "$raw_file"; head -3 "$RAW_FILE" | tail -1 >> "$raw_file"; grep -Fwf "$names_file" "$RAW_FILE" >> "$raw_file"; echo "Generated $raw_file"; done; done

# except for LG3, 4 and 9! they will give huge numbers if they sample the suspect markers

LG3_noSus_names <- colnames(outliers$LG3_rec_map_noSus$data.name$geno)[outliers$LG3_rec_map_noSus$seq.num]
LG4_noSus_names <- colnames(outliers$LG4_rec_map_noSus$data.name$geno)[outliers$LG4_rec_map_noSus$seq.num]
LG9_noSus_names <- colnames(outliers$LG9_rec_map_noSus$data.name$geno)[outliers$LG9_rec_map_noSus$seq.num]

for (i in 1:100) {
  sampled_names <- sample(LG3_noSus_names, 100)
  file_name <- sprintf("Y3088_LG3_noSus_r%d_names.txt", i)
  write.table(sampled_names, file = file_name, sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Generated", file_name, "\n")
}

for (i in 1:100) {
  sampled_names <- sample(LG4_noSus_names, 100)
  file_name <- sprintf("Y3088_LG4_noSus_r%d_names.txt", i)
  write.table(sampled_names, file = file_name, sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Generated", file_name, "\n")
}

for (i in 1:100) {
  sampled_names <- sample(LG9_noSus_names, 100)
  file_name <- sprintf("Y3088_LG9_noSus_r%d_names.txt", i)
  write.table(sampled_names, file = file_name, sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Generated", file_name, "\n")
}

## bash
# RAW_FILE="AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable_noInfo.raw"
# for i in {1..100}; do names_file="Y3088_LG3_noSus_r${i}_names.txt"; raw_file="Y3088_LG3_noSus_r${i}_noInfo.raw"; echo 'data type outcross' > "$raw_file"; echo '1091 100 0 0 0' >> "$raw_file"; head -3 "$RAW_FILE" | tail -1 >> "$raw_file"; grep -Fwf "$names_file" "$RAW_FILE" >> "$raw_file"; echo "Generated $raw_file"; done
# for i in {1..100}; do names_file="Y3088_LG4_noSus_r${i}_names.txt"; raw_file="Y3088_LG4_noSus_r${i}_noInfo.raw"; echo 'data type outcross' > "$raw_file"; echo '1091 100 0 0 0' >> "$raw_file"; head -3 "$RAW_FILE" | tail -1 >> "$raw_file"; grep -Fwf "$names_file" "$RAW_FILE" >> "$raw_file"; echo "Generated $raw_file"; done
# for i in {1..100}; do names_file="Y3088_LG9_noSus_r${i}_names.txt"; raw_file="Y3088_LG9_noSus_r${i}_noInfo.raw"; echo 'data type outcross' > "$raw_file"; echo '1091 100 0 0 0' >> "$raw_file"; head -3 "$RAW_FILE" | tail -1 >> "$raw_file"; grep -Fwf "$names_file" "$RAW_FILE" >> "$raw_file"; echo "Generated $raw_file"; done

library(onemap)
library(dplyr)

for (lg_num in 1:12) {
  lg_name <- paste0("LG", lg_num)
  pattern <- paste0("Y3088_", lg_name, "_r[0-9]+_noInfo\\.raw$")
  
  raw_files <- list.files(pattern = pattern)
  map_lengths <- numeric(length(raw_files))
  
  cat(sprintf("\n===== Processing %s (%d files) =====\n", lg_name, length(raw_files)))
  
  for (i in seq_along(raw_files)) {
    input_file <- raw_files[i]
    cat(sprintf("Processing %s (%d of %d)...\n", input_file, i, length(raw_files)))
    
    a_noInfo <- read_onemap(inputfile = input_file, verbose = FALSE)
    bins_noInfo <- find_bins(a_noInfo, exact = FALSE)
    bins_data_noInfo <- create_data_bins(a_noInfo, bins_noInfo)
    segreg_test <- test_segregation(bins_data_noInfo)
    no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE)
    twopts <- rf_2pts(bins_data_noInfo, verbose = FALSE)
    mark_no_dist <- make_seq(twopts, no_dist)
    LOD_sug <- suggest_lod(mark_no_dist)
    LGs_upgma <- group_upgma(mark_no_dist, expected.groups = 1, inter = FALSE)
    
    set_map_fun(type = "kosambi")
    LG <- make_seq(LGs_upgma, 1)
    LG_rec <- record(LG, hmm = FALSE)
    LG_rec_map <- onemap::map(LG_rec)
    LG_rec_map_len <- dplyr::last(cumsum(kosambi(LG_rec_map$seq.rf)))
    
    map_lengths[i] <- LG_rec_map_len
    cat(sprintf("  Map length: %.2f cM\n", LG_rec_map_len))
  }
  
  out_file <- paste0("Y3088_", lg_name, "_map_lengths.txt")
  write.table(data.frame(File = raw_files, Length_cM = map_lengths),
              file = out_file, row.names = FALSE, sep = "\t", quote = FALSE)
}

## do separately for noSus files changing "pattern" object.
### same for AC1017