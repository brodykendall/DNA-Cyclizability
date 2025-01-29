library(tidyverse)
library(TmCalculator)

cycle5 <- read_csv("cycle5.txt")

dat_AT <- cycle5

free_energy <- read_csv("data/Created/kim_min_free_energy.csv")
free_energy <- free_energy[1:82368,]

dat_AT <- dat_AT %>%
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

dat_AT %>% map_df(~ sum(is.na(.x)))

dat_AT <- dat_AT %>%
  dplyr::select(x50mer = Sequence,
                C0,
                grna_energy, grna_scaffold_energy) %>%
  mutate(x41mer = substring(x50mer, 5, 45))

dat_AT_melting <- dat_AT %>%
  group_by(1:n()) %>%
  transmute(tm1 = Tm_Wallace(x41mer)$Tm,
            tm2 = Tm_Wallace(substring(x41mer, 1, 4))$Tm,
            tm3 = Tm_Wallace(substring(x41mer, 5, 12))$Tm,
            tm4 = Tm_Wallace(substring(x41mer, 36, 40))$Tm) %>%
  ungroup() %>%
  select(-`1:n()`)

sequence_1_df_AT <- dat_AT %>%
  pull(x50mer) %>%
  sequence_df_AT()

sequence_2_df_AT <- dat_AT %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  map(~contain_AT(.x)) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_3_df_AT <- dat_AT %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 3, 2)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  map(~contain_AT(.x)) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_4_df_AT <- dat_AT %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 4, 3)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  map(~contain_AT(.x)) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_5_df_AT <- dat_AT %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 5, 4)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  map(~contain_AT(.x)) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

mono_int_indices_AT = c(9,10,11)
for(i in mono_int_indices_AT) {
  assign(paste0("interaction_", i, "_df_AT"), dat_AT %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap2(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           map(~contain_AT(.x)) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

di_int_indices_AT = c(8,9,10)
for(i in di_int_indices_AT) {
  assign(paste0("di_interaction_", i, "_df_AT"), dat_AT %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap3(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           map(~contain_AT(.x)) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_AT), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

sequence_1_factor_AT <- sequence_1_df_AT %>%
  map_df(~ factor(.x, levels = c(0,1))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))

sequence_2_factor_AT <- sequence_2_df_AT %>%
  map_df(~ factor(.x, levels = c(0,1))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_3_factor_AT <- sequence_3_df_AT %>%
  map_df(~ factor(.x, levels = c(0,1))) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tri")) %>%
  select(-X49tri, -X50tri)

sequence_4_factor_AT <- sequence_4_df_AT %>%
  map_df(~ factor(.x, levels = c(0,1))) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tetra")) %>%
  select(-X48tetra, -X49tetra, -X50tetra)

sequence_5_factor_AT <- sequence_5_df_AT %>%
  map_df(~ factor(.x, levels = c(0,1))) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "penta")) %>%
  select(-X47penta, -X48penta, -X49penta, -X50penta)

factor_seq_AT <- sequence_1_factor_AT %>%
  bind_cols(sequence_2_factor_AT, sequence_3_factor_AT, sequence_4_factor_AT, sequence_5_factor_AT)

for(i in mono_int_indices_AT) {
  assign(paste0("interaction_", i, "_factor_AT"), get(paste0("interaction_", i, "_df_AT")) %>%
           map_df(~ factor(.x, levels = c(0,1))) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "int", i)))
}
# interaction_49_factor <- interaction_49_factor %>%
#   rename(X1int49 = .int49)

for(i in di_int_indices_AT) {
  assign(paste0("di_interaction_", i, "_factor_AT"), get(paste0("di_interaction_", i, "_df_AT")) %>%
           map_df(~ factor(.x, levels = c(0,1))) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "di_int", i)))
}

factor_int_AT <- interaction_9_factor_AT
for(i in mono_int_indices_AT[-1]) {
  factor_int_AT <- factor_int_AT %>%
    bind_cols(get(paste0("interaction_", i, "_factor_AT")))
}

factor_di_int_AT <- di_interaction_8_factor_AT
for(i in di_int_indices_AT[-1]) {
  factor_di_int_AT <- factor_di_int_AT %>%
    bind_cols(get(paste0("di_interaction_", i, "_factor_AT")))
}

library(stringi)

nc_1 <- nucleotides %>%
  map(~ str_count(dat_AT$x50mer, .x)) %>%
  bind_cols()

nc_2 <- dinucleotides %>%
  map(~ str_count(dat_AT$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

nc_3 <- trinucleotides %>%
  map(~ str_count(dat_AT$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

colnames(nc_1) <- nucleotides
colnames(nc_2) <- dinucleotides
colnames(nc_3) <- trinucleotides

## FIXME

nc_1_AT = nc_1 %>%
  mutate(nc1_AorT = rowSums(select(., (contains("A")) | (contains("T")))),
         nc1_no_AorT = rowSums(select(., -((contains("A")) | (contains("T")))))) %>%
  select(-all_of(nucleotides))

nc_2_AT = nc_2 %>%
  mutate(nc2_AorT = rowSums(select(., (contains("A")) | (contains("T")))),
         nc2_no_AorT = rowSums(select(., -((contains("A")) | (contains("T")))))) %>%
  select(-all_of(dinucleotides))

nc_3_AT = nc_3 %>%
  mutate(nc3_AorT = rowSums(select(., (contains("A")) | (contains("T")))),
         nc3_no_AorT = rowSums(select(., -((contains("A")) | (contains("T")))))) %>%
  select(-all_of(trinucleotides))

#FIXME

for(i in mono_int_indices_AT) {
  assign(paste0("nc_int", i), dinucleotides %>%
           map(~ stri_count(dat_AT$x50mer, regex = paste0(substr(.x,1,1),
                                                           paste0('.{', i-1, '}'),
                                                           substr(.x,2,2)))) %>%
           bind_cols())
}

colnames(nc_int9) <- paste0(dinucleotides, "int9")
colnames(nc_int10) <- paste0(dinucleotides, "int10")
colnames(nc_int11) <- paste0(dinucleotides, "int11")

nc_int9_AT = nc_int9 %>%
  mutate(nc_1_int9_AT = rowSums(select(., (contains("A")) | (contains("T", ignore.case = FALSE)))),
         nc_1_int9_no_AorT = rowSums(select(., -((contains("A")) | (contains("T", ignore.case = FALSE)))))) %>%
  select(-all_of(paste0(dinucleotides, "int9")))

nc_int10_AT = nc_int10 %>%
  mutate(nc_1_int10_AT = rowSums(select(., (contains("A")) | (contains("T", ignore.case = FALSE)))),
         nc_1_int10_no_AorT = rowSums(select(., -((contains("A")) | (contains("T", ignore.case = FALSE)))))) %>%
  select(-all_of(paste0(dinucleotides, "int10")))

nc_int11_AT = nc_int11 %>%
  mutate(nc_1_int11_AT = rowSums(select(., (contains("A")) | (contains("T", ignore.case = FALSE)))),
         nc_1_int11_no_AorT = rowSums(select(., -((contains("A")) | (contains("T", ignore.case = FALSE)))))) %>%
  select(-all_of(paste0(dinucleotides, "int11")))


# dat_seq <- factor_seq %>%
#   bind_cols(nc_1,
#             nc_2,
#             nc_3,
#             gc_count,
#             dat)

dat_AT <- dat_AT %>%
  bind_cols(factor_seq_AT,
            factor_int_AT,
            factor_di_int_AT,
            nc_1_AT,
            nc_2_AT,
            nc_3_AT,
            dat_AT_melting)

for(i in mono_int_indices_AT) {
  dat_AT <- dat_AT %>%
    bind_cols(get(paste0("nc_int", i, "_AT")))
}

set.seed(50)

train_indices_AT = sample(1:nrow(dat_AT), nrow(dat_AT)*.9, replace=FALSE)

dat_train_AT = dat_AT[train_indices_AT,]
dat_test_AT = dat_AT[-train_indices_AT,]
dat_AT = dat_train_AT

saveRDS(dat_AT, "data/Created/processed_AT_int.rds")
saveRDS(dat_test_AT, "data/Created/processed_AT_int_test.rds")
