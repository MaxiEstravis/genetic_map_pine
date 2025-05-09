library(onemap)

input_file_Y3088 <- "AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable.raw"
input_file_AC1017_strict <- "AxiomGT1.calls_AC1017_strict_HetMarkedMissing_newCoords_newHeader_sorted_variable.raw"
input_file_AC1017_loose <- "AxiomGT1.calls_AC1017_loose_HetMarkedMissing_newCoords_newHeader_sorted_variable.raw"

a_Y3088 <- read_onemap(inputfile = input_file_Y3088)
a_AC1017_strict <- read_onemap(inputfile = input_file_AC1017_strict)
a_AC1017_loose <- read_onemap(inputfile = input_file_AC1017_loose)

bins_a_Y3088 <- find_bins(a_Y3088, exact = FALSE)
bins_a_AC1017_strict <- find_bins(a_AC1017_strict, exact = FALSE)
bins_a_AC1017_loose <- find_bins(a_AC1017_loose, exact = FALSE)

bins_data_a_AC1017_strict <- create_data_bins(a_AC1017_strict, bins_a_AC1017_strict)

write_onemap_raw(bins_data_a_AC1017_strict, file.name = "AxiomGT1.calls_AC1017_strict_HetMarkedMissing_newCoords_newHeader_sorted_variable_non_redundant.raw")


output_file <- "bins_table_AC1017_loose.txt"
file_conn <- file(output_file, open = "wt")

# Loop through the elements of bins_a$bins
for (i in seq_along(bins_a_AC1017_loose$bins)) {
  if (i == 1) {
    # Write the first element with the header
    write.table(bins_a_AC1017_loose$bins[[i]], file = file_conn, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  } else {
    # Append the rest without the header
    write.table(bins_a_AC1017_loose$bins[[i]], file = file_conn, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
  }
}

# Close the connection
close(file_conn)
