dat = readRDS("data/Created/processed_tiling_newC0.rds")
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")


########################################################################
# Tiling
########################################################################

## Train data:

distances = 1:10

nucleotides_combinations_df = expand.grid(mono1=nucleotides, mono2=nucleotides)
nucleotides_combinations = split(nucleotides_combinations_df, seq(nrow(nucleotides_combinations_df)))

n_gap_di = map2(rep(nucleotides_combinations, length(distances)), 
                rep(distances, each=16), 
                ~stri_count(dat$x50mer, 
                            regex=paste0("(?=(",.x[1,1],".{",.y,"}",.x[1,2],"))"))) %>%
  bind_cols()

colnames(n_gap_di) = paste0(rep(paste(nucleotides_combinations_df$mono1, nucleotides_combinations_df$mono2, sep="_"), length(distances)), 
                            "_gap", rep(distances, each=16))


## Test data:

n_gap_di_test = map2(rep(nucleotides_combinations, length(distances)), 
                     rep(distances, each=16), 
                     ~stri_count(dat_test$x50mer, 
                                 regex=paste0("(?=(",.x[1,1],".{",.y,"}",.x[1,2],"))"))) %>%
  bind_cols()

colnames(n_gap_di_test) = paste0(rep(paste(nucleotides_combinations_df$mono1, nucleotides_combinations_df$mono2, sep="_"), length(distances)), 
                                 "_gap", rep(distances, each=16))


########################################################################
# Random
########################################################################

n_gap_di_random_all = map2(rep(nucleotides_combinations, length(distances)), 
                     rep(distances, each=16), 
                     ~stri_count(dat_random_all$x50mer, 
                                 regex=paste0("(?=(",.x[1,1],".{",.y,"}",.x[1,2],"))"))) %>%
  bind_cols()

colnames(n_gap_di_random_all) = paste0(rep(paste(nucleotides_combinations_df$mono1, nucleotides_combinations_df$mono2, sep="_"), length(distances)), 
                                 "_gap", rep(distances, each=16))
