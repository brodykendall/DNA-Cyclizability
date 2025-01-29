library(tidyverse)
library(TmCalculator)

cycle5 <- read_csv("cycle5.txt")

dat_105 <- cycle5

free_energy <- read_csv("data/Created/kim_min_free_energy.csv")
free_energy <- free_energy[1:82368,]

dat_105 <- dat_105 %>%
  bind_cols(free_energy)

source("scripts/functions/sequence-functions.R")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
tetranucleotides <- gtools::permutations(n = 4, r = 4, v = nucleotides,
                                         repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
pentanucleotides <- gtools::permutations(n = 4, r = 5, v = nucleotides,
                                         repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

dat_105 %>% map_df(~ sum(is.na(.x)))

dat_105 <- dat_105 %>%
  dplyr::select(x50mer = Sequence,
                C0,
                grna_energy, grna_scaffold_energy) %>%
  mutate(x41mer = substring(x50mer, 5, 45))

dat_105_melting <- dat_105 %>%
  group_by(1:n()) %>%
  transmute(tm1 = Tm_Wallace(x41mer)$Tm,
            tm2 = Tm_Wallace(substring(x41mer, 1, 4))$Tm,
            tm3 = Tm_Wallace(substring(x41mer, 5, 12))$Tm,
            tm4 = Tm_Wallace(substring(x41mer, 36, 40))$Tm) %>%
  ungroup() %>%
  select(-`1:n()`)

sequence_1_df_105 <- dat_105 %>%
  pull(x50mer) %>%
  sequence_df()

sequence_2_df_105 <- dat_105 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_3_df_105 <- dat_105 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 3, 2)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_4_df_105 <- dat_105 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 4, 3)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_5_df_105 <- dat_105 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 5, 4)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

mono_int_indices_105 = c(9,10,11)
for(i in mono_int_indices_105) {
  assign(paste0("interaction_", i, "_df_105"), dat_105 %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap2(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

di_int_indices_105 = c(8,9,10)
# di_int_indices_105 = c(8)
for(i in di_int_indices_105) {
  assign(paste0("di_interaction_", i, "_df_105"), dat_105 %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap3(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_105), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

sequence_1_factor_105 <- sequence_1_df_105 %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))

sequence_2_factor_105 <- sequence_2_df_105 %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_3_factor_105 <- sequence_3_df_105 %>%
  map_df(~ factor(.x, levels = trinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tri")) %>%
  select(-X49tri, -X50tri)

sequence_4_factor_105 <- sequence_4_df_105 %>%
  map_df(~ factor(.x, levels = tetranucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tetra")) %>%
  select(-X48tetra, -X49tetra, -X50tetra)

sequence_5_factor_105 <- sequence_5_df_105 %>%
  map_df(~ factor(.x, levels = pentanucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "penta")) %>%
  select(-X47penta, -X48penta, -X49penta, -X50penta)

factor_seq_105 <- sequence_1_factor_105 %>%
  bind_cols(sequence_2_factor_105, sequence_3_factor_105, sequence_4_factor_105, sequence_5_factor_105)

for(i in mono_int_indices_105) {
  assign(paste0("interaction_", i, "_factor_105"), get(paste0("interaction_", i, "_df_105")) %>%
           map_df(~ factor(.x, levels = dinucleotides)) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "int", i)))
}
# interaction_49_factor <- interaction_49_factor %>%
#   rename(X1int49 = .int49)

for(i in di_int_indices_105) {
  assign(paste0("di_interaction_", i, "_factor_105"), get(paste0("di_interaction_", i, "_df_105")) %>%
           map_df(~ factor(.x, levels = tetranucleotides)) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "di_int", i)))
}

factor_int_105 <- interaction_9_factor_105
for(i in mono_int_indices_105[-1]) {
  factor_int_105 <- factor_int_105 %>%
    bind_cols(get(paste0("interaction_", i, "_factor_105")))
}

factor_di_int_105 <- di_interaction_8_factor_105
for(i in di_int_indices_105[-1]) {
  factor_di_int_105 <- factor_di_int_105 %>%
    bind_cols(get(paste0("di_interaction_", i, "_factor_105")))
}

library(stringi)

nc_1_105 <- nucleotides %>%
  map(~ str_count(dat_105$x50mer, .x)) %>%
  bind_cols()

nc_2_105 <- dinucleotides %>%
  map(~ str_count(dat_105$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

nc_3_105 <- trinucleotides %>%
  map(~ str_count(dat_105$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

colnames(nc_1_105) <- nucleotides
colnames(nc_2_105) <- dinucleotides
colnames(nc_3_105) <- trinucleotides

for(i in mono_int_indices_105) {
  assign(paste0("nc_int", i, "_105"), dinucleotides %>%
           map(~ stri_count(dat_105$x50mer, regex = paste0(substr(.x,1,1),
                                                           paste0('.{', i-1, '}'),
                                                           substr(.x,2,2)))) %>%
           bind_cols())
}

colnames(nc_int9_105) <- paste0(dinucleotides, "int9")
colnames(nc_int10_105) <- paste0(dinucleotides, "int10")
colnames(nc_int11_105) <- paste0(dinucleotides, "int11")

#GC count
gc_count_105 <- nc_1_105 %>%
  transmute(gc_count_105 = G + C)

# dat_seq <- factor_seq %>%
#   bind_cols(nc_1,
#             nc_2,
#             nc_3,
#             gc_count,
#             dat)

dat_105 <- dat_105 %>%
  bind_cols(factor_seq_105,
            factor_int_105,
            factor_di_int_105,
            nc_1_105,
            nc_2_105,
            nc_3_105, 
            gc_count_105,
            dat_105_melting)

for(i in mono_int_indices_105) {
  dat_105 <- dat_105 %>%
    bind_cols(get(paste0("nc_int", i, "_105")))
}

set.seed(50)

train_indices_105 = sample(1:nrow(dat_105), nrow(dat_105)*.9, replace=FALSE)

dat_train_105 = dat_105[train_indices_105,]
dat_test_105 = dat_105[-train_indices_105,]
dat_105 = dat_train_105

saveRDS(dat_105, "data/Created/processed_105_int.rds")
saveRDS(dat_test_105, "data/Created/processed_105_int_test.rds")
