library(tidyverse)

# Load SNP marker data
snp_df <- read_tsv("pine_SNParray_to_genome.tsv", col_names = c("SNP_ID", "SNP_Chromosome", "SNP_Position"))

# Load primer data
primer_df <- read_tsv("primers_2003_vs_Pinsy01_chromosomes.psl", col_names = paste0("V", 1:20)) %>%
  transmute(
    Primer = V10,
    Chromosome = V14,
    Start = as.numeric(V16),
    End = as.numeric(V17),
    Midpoint = floor((Start + End) / 2)
  )

# Compute closest SNP for each primer
closest_matches <- primer_df %>%
  rowwise() %>%
  mutate(
    Closest_SNP = list({
      snp_subset <- snp_df %>% filter(SNP_Chromosome == Chromosome)
      if (nrow(snp_subset) == 0) {
        tibble(SNP_ID = NA, SNP_Position = NA, Distance = NA)
      } else {
        snp_subset %>%
          mutate(Distance = abs(SNP_Position - Midpoint)) %>%
          arrange(Distance) %>%
          slice(1)
      }
    })
  ) %>%
  unnest(Closest_SNP, names_sep = "_")

# Save the results
closest_matches %>%
  select(Primer, Chromosome, Midpoint, SNP_ID = Closest_SNP_SNP_ID, SNP_Position = Closest_SNP_SNP_Position, Distance = Closest_SNP_Distance) %>%
  write_tsv("closest_snps.tsv")
