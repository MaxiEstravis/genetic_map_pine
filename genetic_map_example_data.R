onemap_example_out <- read_onemap(inputfile = "onemap_example_out.raw")
vcf_example_out <- onemap_read_vcfR(vcf = system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"), parent1 = "P1", parent2 = "P2", cross = "outcross")
vcf_filtered <- filter_missing(vcf_example_out, threshold = 0.25)
comb_example <- combine_onemap(onemap_example_out, vcf_example_out)
bins <- find_bins(comb_example, exact = FALSE)
bins_example <- create_data_bins(comb_example, bins)
segreg_test <- test_segregation(bins_example)
dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE)
no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE)
twopts <- rf_2pts(bins_example)
mark_no_dist <- make_seq(twopts, c(no_dist))
LOD_sug <- suggest_lod(mark_no_dist)
LGs <- group(mark_no_dist, LOD=LOD_sug)

LG2 <- make_seq(LGs, 2)
LG2_rcd <- rcd(LG2, hmm = F)
LG2_map <- onemap::map(LG2_rcd)

add_redundants(LG2_map, comb_example, bins)   # works!


LG1 <- make_seq(LGs, 1)
LG1_rcd <- rcd(LG1, hmm = F)
LG1_map <- onemap::map(LG1_rcd)

add_redundants(LG1_map, comb_example, bins)   # works!


###########

vvcf_example_out <- onemap_read_vcfR(vcf = "vcf_example_out.vcf.gz", parent1 = "P1", parent2 = "P2", cross = "outcross")
vvcf_filtered <- filter_missing(vvcf_example_out, threshold = 0.25)
#oonemap_example_out <- read_onemap(inputfile = "onemap_example_out.raw")
#bbins <- find_bins(oonemap_example_out, exact = FALSE)
#bbins_example <- create_data_bins(oonemap_example_out, bbins)
bbins <- find_bins(vvcf_filtered, exact = FALSE)
bbins_example <- create_data_bins(vvcf_filtered, bbins)
ssegreg_test <- test_segregation(bbins_example)
ddist <- select_segreg(ssegreg_test, distorted = TRUE, numbers = TRUE)
nno_dist <- select_segreg(ssegreg_test, distorted = FALSE, numbers = TRUE)
ttwopts <- rf_2pts(bbins_example)
mmark_no_dist <- make_seq(ttwopts, c(nno_dist))
LLOD_sug <- suggest_lod(mmark_no_dist)
LLGs <- group(mmark_no_dist, LOD=LLOD_sug)

LLG1 <- make_seq(LLGs, 1)
LLG1_rcd <- rcd(LLG1, hmm = F)
LLG1_map <- onemap::map(LLG1_rcd)

LLG2 <- make_seq(LLGs, 2)
LLG2_rcd <- rcd(LLG2, hmm = F)
LLG2_map <- onemap::map(LLG2_rcd)

LLG3 <- make_seq(LLGs, 3)
LLG3_rcd <- rcd(LLG3, hmm = F)
LLG3_map <- onemap::map(LLG3_rcd)

LLG4 <- make_seq(LLGs, 4)
LLG4_rcd <- rcd(LLG4, hmm = F)
LLG4_map <- onemap::map(LLG4_rcd)

LLG1_map_redu <- add_redundants(LLG1_map, vvcf_filtered, bbins)  # works!
LLG2_map_redu <- add_redundants(LLG2_map, vvcf_filtered, bbins)  # Error in do.call(c, mks) : second argument must be a list
LLG3_map_redu <- add_redundants(LLG3_map, vvcf_filtered, bbins)  # Error in do.call(c, mks) : second argument must be a list
LLG4_map_redu <- add_redundants(LLG4_map, vvcf_filtered, bbins)  # works!




