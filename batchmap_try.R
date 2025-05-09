## with example data:

input_file <- system.file("example/sim7.5k.txt.gz",package = "BatchMap")
outcross <- read.outcross2(input_file)
bins <- find.bins(outcross, exact = FALSE)
outcross_clean <- create.data.bins(outcross, bins)
twopt_table <- rf.2pts(outcross_clean)
linkage_groups <- BatchMap::group(make.seq(input.obj = twopt_table, "all"), LOD = 12)

library(BatchMap)

input_file_BM <- "AxiomGT1.calls_Y3088_HetMarkedMissing_newCoords_newHeader_sorted_variable_non_redundant_BM.raw"
# the input file for BatchMap has a slightly different format

a_BM <- read.outcross2(input_file_BM)
bins_BM <- find.bins(a_BM, exact = FALSE)
outcross_clean <- create.data.bins(a_BM, bins_BM)

####
> outcross_clean
This is an object of class 'outcross'
No. individuals:    1091 
No. markers:        5471 
Segregation types:
  D1.11:	5471
No. traits:         0 
####  # barely better than the input, which was already binned by OneMap; i continue with this

suggested_LOD <- suggest.lod(a)     # i think this is a bug, if i write suggest.lod(outcross_clean) or suggest.lod(a_BM) it tells me it's not a onemap object...i know! it's "outcross" object but seems to be the same!
# 7.60937412249278
twopt_table <- rf.2pts(outcross_clean, LOD = suggested_LOD)
linkage_groups <- BatchMap::group(make.seq(input.obj = twopt_table, "all"))
linkage_groups_2 <- BatchMap::group(make.seq(input.obj = twopt_table, "all"), LOD = suggested_LOD, max.rf = 0.35)

