## after deciding to keep AC1017_strict, i inspect the LGs for suspicious markers

# for Y3088, LG4

## the first 5 markers are clearly misplaced
#############
#> all_samples_with_bins$Y3088$LGs_progress$LG4_rec_map
#4373 AX-601663619              0.00           a |  | b       o |  | o 
#5464 AX-601707602             75.71           a |  | b       o |  | o 
#4364 AX-613169263            104.03           b |  | a       o |  | o 
#2036 AX-600608066            105.86           b |  | a       o |  | o 
#1454 AX-602744081            125.42           a |  | b       o |  | o
#[...]
#############

all_samples_with_bins$Y3088$outliers <- list()

all_samples_with_bins$Y3088$outliers$suspects_LG4 <- drop_marker(all_samples_with_bins$Y3088$LGs_progress$LG4_rec_map,
                                                                       c(4373, 5464, 4364, 2036, 1454))

all_samples_with_bins$Y3088$outliers$LG4_rec_map_noSus <- onemap::map(all_samples_with_bins$Y3088$outliers$suspects_LG4)

all_samples_with_bins$Y3088$outliers$LG4_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG4_rec_map_noSus, 4373)
all_samples_with_bins$Y3088$outliers$LG4_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG4_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG4_test_seq$LOD)))
all_samples_with_bins$Y3088$outliers$LG4_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG4_extended, 5464)
all_samples_with_bins$Y3088$outliers$LG4_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG4_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG4_test_seq$LOD)))
all_samples_with_bins$Y3088$outliers$LG4_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG4_extended, 4364)
all_samples_with_bins$Y3088$outliers$LG4_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG4_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG4_test_seq$LOD)))
all_samples_with_bins$Y3088$outliers$LG4_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG4_extended, 2036)
all_samples_with_bins$Y3088$outliers$LG4_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG4_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG4_test_seq$LOD)))
all_samples_with_bins$Y3088$outliers$LG4_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG4_extended, 1454)
all_samples_with_bins$Y3088$outliers$LG4_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG4_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG4_test_seq$LOD)))

draw_map2(all_samples_with_bins$Y3088$LGs_progress$LG4_rec_map, 
          all_samples_with_bins$Y3088$LGs_progress$LG4_rec_map_redu,
          all_samples_with_bins$Y3088$outliers$LG4_rec_map_noSus, 
          all_samples_with_bins$Y3088$outliers$LG4_extended,
          main = "Y3088", group.names = c("LG4", "LG4_redu", "LG4_noSus", "LG4_extended"),
          tag = c("AX-601663619", "AX-601707602", "AX-613169263", "AX-600608066", "AX-602744081"), pos = F, cex.label = 0.3,
          output = "Y3088_LG4s.pdf")

############ they were re-mapped in the same way. i will remove them


# for Y3088, LG3

#############
#> all_samples_with_bins$Y3088$LGs_progress$LG3_rec_map
#[...]
#1118 AX-614120084            187.43           b |  | a       o |  | o 
#############

all_samples_with_bins$Y3088$outliers$suspects_LG3 <- drop_marker(all_samples_with_bins$Y3088$LGs_progress$LG3_rec_map, 
                                                                 c(1118))

all_samples_with_bins$Y3088$outliers$LG3_rec_map_noSus <- onemap::map(all_samples_with_bins$Y3088$outliers$suspects_LG3)

all_samples_with_bins$Y3088$outliers$LG3_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG3_rec_map_noSus, 1118)
all_samples_with_bins$Y3088$outliers$LG3_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG3_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG3_test_seq$LOD)))

draw_map2(all_samples_with_bins$Y3088$LGs_progress$LG3_rec_map, 
          all_samples_with_bins$Y3088$LGs_progress$LG3_rec_map_redu, 
          all_samples_with_bins$Y3088$outliers$LG3_rec_map_noSus, 
          all_samples_with_bins$Y3088$outliers$LG3_extended,
          main = "Y3088", group.names = c("LG3", "LG3_redu", "LG3_noSus", "LG3_extended"),
          tag = c("AX-614120084"), pos = F, cex.label = 0.3,
          output = "Y3088_LG3s.pdf")


############ they were re-mapped in the same way. i will remove them



## for Y3088, LG9
### a not-so-noticeable outlier in the end...let's see

#############
#> all_samples_with_bins$Y3088$LGs_progress$LG9_rec_map
#[...]
#3212 AX-601718530            111.92           b |  | a       o |  | o 
#############

all_samples_with_bins$Y3088$outliers$suspects_LG9 <- drop_marker(all_samples_with_bins$Y3088$LGs_progress$LG9_rec_map, 
                                                                 c(3212))

all_samples_with_bins$Y3088$outliers$LG9_rec_map_noSus <- onemap::map(all_samples_with_bins$Y3088$outliers$suspects_LG9)

all_samples_with_bins$Y3088$outliers$LG9_test_seq <- try_seq(all_samples_with_bins$Y3088$outliers$LG9_rec_map_noSus, 3212)
all_samples_with_bins$Y3088$outliers$LG9_extended <- make_seq(all_samples_with_bins$Y3088$outliers$LG9_test_seq, 
                                                              which.min(abs(all_samples_with_bins$Y3088$outliers$LG9_test_seq$LOD)))

draw_map2(all_samples_with_bins$Y3088$LGs_progress$LG9_rec_map, 
          all_samples_with_bins$Y3088$LGs_progress$LG9_rec_map_redu, 
          all_samples_with_bins$Y3088$outliers$LG9_rec_map_noSus, 
          all_samples_with_bins$Y3088$outliers$LG9_extended,
          main = "Y3088", group.names = c("LG9", "LG9_redu", "LG9_noSus", "LG9_extended"),
          tag = c("AX-601718530"), pos = F, cex.label = 0.3,
          output = "Y3088_LG9s.pdf")


# for AC1017_strict, LG2

#############
#> all_samples_with_bins$AC1017_strict$LGs_progress$LG2_rec_map
#[...]
#3317 AX-599941399            162.50           b |  | a       o |  | o 
#1372 AX-600301767            169.44           b |  | a       o |  | o 
#70 AX-388338409            169.44           b |  | a       o |  | o 
#129 AX-613129619            170.71           a |  | b       o |  | o 
#2075 AX-602614412            181.32           b |  | a       o |  | o 
#767 AX-602055318            189.47           b |  | a       o |  | o 
#479 AX-601666419            216.86           b |  | a       o |  | o 
#146 AX-601101069            265.41           b |  | a       o |  | o 
#3685 AX-601609311            321.17           b |  | a       o |  | o 
#2952 AX-613550901            365.70           b |  | a       o |  | o  
#############

all_samples_with_bins$AC1017_strict$outliers <- list()

all_samples_with_bins$AC1017_strict$outliers$suspects_LG2 <- drop_marker(all_samples_with_bins$AC1017_strict$LGs_progress$LG2_rec_map, 
                                                                         c(3317, 1372, 70, 129, 2075, 767, 479, 146, 3685, 2952))

# i map again without the sus (i will add the definite to the list at the end)
all_samples_with_bins$AC1017_strict$outliers$LG2_rec_map_noSus <- onemap::map(all_samples_with_bins$AC1017_strict$outliers$suspects_LG2)

all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_rec_map_noSus, 3317)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 1372)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 70)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 129)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 2075)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq,
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 767)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 479)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 146)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 3685)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))
all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq <- try_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_extended, 2952)
all_samples_with_bins$AC1017_strict$outliers$LG2_extended <- make_seq(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq, 
                                       which.min(abs(all_samples_with_bins$AC1017_strict$outliers$LG2_test_seq$LOD)))

draw_map2(all_samples_with_bins$AC1017_strict$LGs_progress$LG2_rec_map,
          all_samples_with_bins$AC1017_strict$LGs_progress$LG2_rec_map_redu,
          all_samples_with_bins$AC1017_strict$outliers$LG2_rec_map_noSus,
          all_samples_with_bins$AC1017_strict$outliers$LG2_extended,
          main = "AC1017_strict", group.names = c("LG2", "LG2_redu", "LG2_noSus", "LG2_extended"),
          tag = c("AX-599941399", "AX-613129619", "AX-600301767", "AX-388338409", "AX-602055318",
                  "AX-602614412", "AX-613550901", "AX-601609311", "AX-601666419", "AX-601101069"),
          pos = F, cex.label = 0.3,
          output = "AC1017_strict_LG2s.pdf")

############ they were re-mapped in the same way. i will remove them





#### making the definite LGs without outliers:

all_samples_with_bins$Y3088$LGs_def <- list()
all_samples_with_bins$AC1017_strict$LGs_def <- list()


all_samples_with_bins$Y3088$LGs_def$LG4 <- add_redundants(all_samples_with_bins$Y3088$outliers$LG4_rec_map_noSus,
                                                          all_samples_with_bins$Y3088$prelim$a_noInfo,
                                                          all_samples_with_bins$Y3088$prelim$bins_noInfo)
all_samples_with_bins$Y3088$LGs_def$LG3 <- add_redundants(all_samples_with_bins$Y3088$outliers$LG3_rec_map_noSus,
                                               all_samples_with_bins$Y3088$prelim$a_noInfo,
                                               all_samples_with_bins$Y3088$prelim$bins_noInfo)
all_samples_with_bins$Y3088$LGs_def$LG9 <- add_redundants(all_samples_with_bins$Y3088$outliers$LG9_rec_map_noSus,
                                               all_samples_with_bins$Y3088$prelim$a_noInfo,
                                               all_samples_with_bins$Y3088$prelim$bins_noInfo)

all_samples_with_bins$AC1017_strict$LGs_def$LG2 <- add_redundants(all_samples_with_bins$AC1017_strict$outliers$LG2_rec_map_noSus,
                                                                  all_samples_with_bins$AC1017_strict$prelim$a_noInfo,
                                                                  all_samples_with_bins$AC1017_strict$prelim$bins_noInfo)


sample_name <- "Y3088"
except_numbers <- c(3, 4, 9)

for (i in 1:12) {
  if (i %in% except_numbers) next
  lg_rec_map_redu_name <- paste0("LG", i, "_rec_map_redu")
  new_name <- paste0("LG", i)
  all_samples_with_bins[[sample_name]]$LGs_def[[new_name]] <- all_samples_with_bins[[sample_name]]$LGs_progress[[lg_rec_map_redu_name]]
  message("Copied ", lg_rec_map_redu_name, " to ", sample_name, "$LGs_def")
}


##

sample_ids <- c("AC1017_strict", "Y3088")

for (sample in sample_ids) {
  for (i in 1:12) {
    lg_name <- paste0("LG", i)
    filename <- paste0(sample, "_", lg_name, "_def.tsv")
    marker_positions <- data.frame(
      id = colnames(all_samples_with_bins[[sample]]$LGs_def[[lg_name]]$data.name$geno)[all_samples_with_bins[[sample]]$LGs_def[[lg_name]]$seq.num], 
      pos = c(0, cumsum(kosambi(all_samples_with_bins[[sample]]$LGs_def[[lg_name]]$seq.rf)))
    )
    gaps <- diff(marker_positions$pos)
    max_gap_index <- which.max(gaps)
    marker1 <- marker_positions$id[max_gap_index]
    marker2 <- marker_positions$id[max_gap_index + 1]
    max_gap <- gaps[max_gap_index]
    cat(sample, "largest gap in", lg_name, ":", max_gap, "between markers", marker1, "and", marker2, "\n")
    
    all_samples_with_bins[[sample]]$LGs_def[[paste0(lg_name, "_marker_positions")]] <- marker_positions
    write.table(marker_positions, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
    if (file.exists(filename)) {
      cat("File", filename, "written successfully!\n")
    } else {
      cat("There was an error writing the file", filename, "\n")
    }
  }
}

#AC1017_strict largest gap in LG1 : 4.978856 between markers AX-601646580 and AX-611396146 
#AC1017_strict largest gap in LG2 : 10.04337 between markers AX-612486918 and AX-612670125 
#AC1017_strict largest gap in LG3 : 4.865226 between markers AX-600089946 and AX-602688305 
#AC1017_strict largest gap in LG4 : 3.469809 between markers AX-388291743 and AX-601276436 
#AC1017_strict largest gap in LG5 : 2.542605 between markers AX-600629136 and AX-601630018 
#AC1017_strict largest gap in LG6 : 3.006583 between markers AX-388338706 and AX-600575500 
#AC1017_strict largest gap in LG7 : 4.166718 between markers AX-613278935 and AX-600726253 
#AC1017_strict largest gap in LG8 : 2.542605 between markers AX-601855776 and AX-599797360 
#AC1017_strict largest gap in LG9 : 2.371375 between markers AX-600687531 and AX-601580805 
#AC1017_strict largest gap in LG10 : 4.413774 between markers AX-606723910 and AX-613966730 
#AC1017_strict largest gap in LG11 : 4.475772 between markers AX-600428750 and AX-388046379 
#AC1017_strict largest gap in LG12 : 3.136985 between markers AX-612482189 and AX-389881949 
#Y3088 largest gap in LG1 : 3.847066 between markers AX-601610175 and AX-600414952 
#Y3088 largest gap in LG2 : 6.265555 between markers AX-612418925 and AX-387860822 
#Y3088 largest gap in LG3 : 2.654802 between markers AX-606034332 and AX-606072772 
#Y3088 largest gap in LG4 : 4.779286 between markers AX-600444235 and AX-601881132 
#Y3088 largest gap in LG5 : 3.304556 between markers AX-601710147 and AX-601615194 
#Y3088 largest gap in LG6 : 2.109757 between markers AX-601351462 and AX-600973078 
#Y3088 largest gap in LG7 : 7.29233 between markers AX-600430751 and AX-606249668 
#Y3088 largest gap in LG8 : 7.24821 between markers AX-612838135 and AX-388286590 
#Y3088 largest gap in LG9 : 2.093904 between markers AX-389541845 and AX-613762226 
#Y3088 largest gap in LG10 : 4.597611 between markers AX-613604270 and AX-601610282 
#Y3088 largest gap in LG11 : 2.752578 between markers AX-600591070 and AX-601605713 
#Y3088 largest gap in LG12 : 4.595981 between markers AX-606316295 and AX-601294220 

#### at this stage i change the names of the .tsv files removing the "strict" and just leaving "AC1017"

for (sample in sample_ids) {
  nind <- all_samples_with_bins[[sample]]$prelim$bins$info$n.ind
  nmar <- all_samples_with_bins[[sample]]$prelim$bins$info$n.mar
  nbin <- length(all_samples_with_bins[[sample]]$prelim$bins$bins)
  nlgs <- sum(sapply(1:12, function(i) length(all_samples_with_bins[[sample]]$LGs_progress[[paste0("LG", i, "_rec_map")]]$seq.num)))
  nredu <- sum(sapply(1:12, function(i) length(all_samples_with_bins[[sample]]$LGs_def[[paste0("LG", i, "_marker_positions")]][,1])))
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
#Number of total markers in LGs: 9789 
#Percentage of original markers used in LGs: 84.05461 

#Y3088 summary
#Number of individuals: 1091 
#Number of markers in original dataset: 12851 
#Number of bins found: 5556 
#Number of non-redundant markers in LGs: 4754 
#Number of total markers in LGs: 11696 
#Percentage of original markers used in LGs: 91.01237 


for (sample in sample_ids) {
  for (i in 1:12) {
  nlgs <- length(all_samples_with_bins[[sample]]$LGs_progress[[paste0("LG", i, "_rec_map")]]$seq.num)
  nredu <- length(all_samples_with_bins[[sample]]$LGs_def[[paste0("LG", i, "_marker_positions")]][,1])
  cat(paste0(sample, " LG", i, " summary\n"))
  cat("Number of bins:", nlgs, "\n")
  cat("Number of markers:", nredu, "\n")
  }
}

##AC1017_strict LG1 summary
#Number of bins: 239 
#Number of markers: 796 
##AC1017_strict LG2 summary
#Number of bins: 160 
#Number of markers: 487 
##AC1017_strict LG3 summary
#Number of bins: 223 
#Number of markers: 737 
##AC1017_strict LG4 summary
#Number of bins: 213 
#Number of markers: 744 
##AC1017_strict LG5 summary
#Number of bins: 278 
#Number of markers: 959 
##AC1017_strict LG6 summary
#Number of bins: 213 
#Number of markers: 717 
##AC1017_strict LG7 summary
#Number of bins: 228 
#Number of markers: 812 
##AC1017_strict LG8 summary
#Number of bins: 265 
#Number of markers: 855 
##AC1017_strict LG9 summary
#Number of bins: 231 
#Number of markers: 869 
##AC1017_strict LG10 summary
#Number of bins: 261 
#Number of markers: 1021 
##AC1017_strict LG11 summary
#Number of bins: 235 
#Number of markers: 852 
##AC1017_strict LG12 summary
#Number of bins: 254 
#Number of markers: 940 
##Y3088 LG1 summary
#Number of bins: 332 
#Number of markers: 665 
##Y3088 LG2 summary
#Number of bins: 376 
#Number of markers: 909 
##Y3088 LG3 summary
#Number of bins: 426 
#Number of markers: 1086 
##Y3088 LG4 summary
#Number of bins: 409 
#Number of markers: 904 
##Y3088 LG5 summary
#Number of bins: 389 
#Number of markers: 992 
##Y3088 LG6 summary
#Number of bins: 409 
#Number of markers: 988 
##Y3088 LG7 summary
#Number of bins: 441 
#Number of markers: 1185 
##Y3088 LG8 summary
#Number of bins: 367 
#Number of markers: 912 
##Y3088 LG9 summary
#Number of bins: 428 
#Number of markers: 1152 
##Y3088 LG10 summary
#Number of bins: 458 
#Number of markers: 1063 
##Y3088 LG11 summary
#Number of bins: 346 
#Number of markers: 850 
##Y3088 LG12 summary
#Number of bins: 373 
#Number of markers: 990 

## continues in consensus.R
