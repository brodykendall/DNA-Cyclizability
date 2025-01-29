library(tidyverse)
library(TmCalculator)

cycle3 <- read_csv("cycle3.txt")

dat_random_library <- cycle3

free_energy_random_library <- read_csv("data/Created/kim_min_free_energy_random_library.csv")
# free_energy_random_library <- free_energy_random_library[1:12472,]

dat_random_library <- dat_random_library %>%
  bind_cols(free_energy_random_library)

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

dat_random_library %>% map_df(~ sum(is.na(.x)))

dat_random_library <- dat_random_library %>%
  dplyr::select(x50mer = Sequence,
                C0,
                grna_energy, grna_scaffold_energy) %>%
  mutate(x41mer = substring(x50mer, 5, 45))

dat_melting_random_library <- dat_random_library %>%
  group_by(1:n()) %>%
  transmute(tm1 = Tm_Wallace(x41mer)$Tm,
            tm2 = Tm_Wallace(substring(x41mer, 1, 4))$Tm,
            tm3 = Tm_Wallace(substring(x41mer, 5, 12))$Tm,
            tm4 = Tm_Wallace(substring(x41mer, 36, 40))$Tm) %>%
  ungroup() %>%
  select(-`1:n()`)

sequence_1_df_random_library <- dat_random_library %>%
  pull(x50mer) %>%
  sequence_df()

sequence_2_df_random_library <- dat_random_library %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_3_df_random_library <- dat_random_library %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 3, 2)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_4_df_random_library <- dat_random_library %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 4, 3)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_5_df_random_library <- dat_random_library %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 5, 4)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

for(i in 3:49) {
  assign(paste0("interaction_", i, "_df_random_library"), dat_random_library %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap2(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

for(i in 3:47) {
  assign(paste0("di_interaction_", i, "_df_random_library"), dat_random_library %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap3(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat_random_library), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

sequence_1_factor_random_library <- sequence_1_df_random_library %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))

sequence_2_factor_random_library <- sequence_2_df_random_library %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_3_factor_random_library <- sequence_3_df_random_library %>%
  map_df(~ factor(.x, levels = trinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tri")) %>%
  select(-X49tri, -X50tri)

sequence_4_factor_random_library <- sequence_4_df_random_library %>%
  map_df(~ factor(.x, levels = tetranucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tetra")) %>%
  select(-X48tetra, -X49tetra, -X50tetra)

sequence_5_factor_random_library <- sequence_5_df_random_library %>%
  map_df(~ factor(.x, levels = pentanucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "penta")) %>%
  select(-X47penta, -X48penta, -X49penta, -X50penta)

factor_seq_random_library <- sequence_1_factor_random_library %>%
  bind_cols(sequence_2_factor_random_library, sequence_3_factor_random_library,
            sequence_4_factor_random_library, sequence_5_factor_random_library)

for(i in 3:49) {
  assign(paste0("interaction_", i, "_factor_random_library"), get(paste0("interaction_", i, "_df_random_library")) %>%
           map_df(~ factor(.x, levels = dinucleotides)) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "int", i)))
}
interaction_49_factor_random_library <- interaction_49_factor_random_library %>%
  rename(X1int49 = .int49)

factor_int_random_library <- interaction_3_factor_random_library
for(i in 4:49) {
  factor_int_random_library <- factor_int_random_library %>%
    bind_cols(get(paste0("interaction_", i, "_factor_random_library")))
}


for(i in 3:47) {
  assign(paste0("di_interaction_", i, "_factor_random_library"), 
         get(paste0("di_interaction_", i, "_df_random_library")) %>%
           map_df(~ factor(.x, levels = tetranucleotides)) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "di_int", i)))
}
di_interaction_47_factor_random_library <- di_interaction_47_factor_random_library %>%
  rename(X1di_int47 = .di_int47)

factor_di_int_random_library <- di_interaction_3_factor_random_library
for(i in 4:47) {
  factor_di_int_random_library <- factor_di_int_random_library %>%
    bind_cols(get(paste0("di_interaction_", i, "_factor_random_library")))
}

library(stringi)

nc_1_random_library <- nucleotides %>%
  map(~ str_count(dat_random_library$x50mer, .x)) %>%
  bind_cols()

nc_2_random_library<- dinucleotides %>%
  map(~ str_count(dat_random_library$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

nc_3_random_library <- trinucleotides %>%
  map(~ str_count(dat_random_library$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

colnames(nc_1_random_library) <- nucleotides
colnames(nc_2_random_library) <- dinucleotides
colnames(nc_3_random_library) <- trinucleotides

# for(i in 3:49) {
#   assign(paste0("nc_int", i), dinucleotides %>%
#            map(~ stri_count(dat$x50mer, regex = paste0(substr(.x,1,1), 
#                                                        paste0('.{', i-1, '}'), 
#                                                        substr(.x,2,2)))) %>%
#            bind_cols())
# }

# colnames(nc_int3) <- paste0(dinucleotides, "int3")
# colnames(nc_int4) <- paste0(dinucleotides, "int4")
# colnames(nc_int5) <- paste0(dinucleotides, "int5")
# colnames(nc_int6) <- paste0(dinucleotides, "int6")
# colnames(nc_int7) <- paste0(dinucleotides, "int7")
# colnames(nc_int3) <- paste0(dinucleotides, "int8")
# colnames(nc_int9) <- paste0(dinucleotides, "int9")
# colnames(nc_int10) <- paste0(dinucleotides, "int10")
# colnames(nc_int11) <- paste0(dinucleotides, "int11")
# colnames(nc_int12) <- paste0(dinucleotides, "int12")
# colnames(nc_int13) <- paste0(dinucleotides, "int13")
# colnames(nc_int14) <- paste0(dinucleotides, "int14")
# colnames(nc_int15) <- paste0(dinucleotides, "int15")
# colnames(nc_int16) <- paste0(dinucleotides, "int16")
# colnames(nc_int17) <- paste0(dinucleotides, "int17")
# colnames(nc_int18) <- paste0(dinucleotides, "int18")
# colnames(nc_int19) <- paste0(dinucleotides, "int19")
# colnames(nc_int10) <- paste0(dinucleotides, "int20")
# colnames(nc_int21) <- paste0(dinucleotides, "int21")
# colnames(nc_int22) <- paste0(dinucleotides, "int22")
# colnames(nc_int23) <- paste0(dinucleotides, "int23")
# colnames(nc_int24) <- paste0(dinucleotides, "int24")
# colnames(nc_int25) <- paste0(dinucleotides, "int25")
# colnames(nc_int26) <- paste0(dinucleotides, "int26")
# colnames(nc_int27) <- paste0(dinucleotides, "int27")
# colnames(nc_int28) <- paste0(dinucleotides, "int28")
# colnames(nc_int29) <- paste0(dinucleotides, "int29")
# colnames(nc_int30) <- paste0(dinucleotides, "int30")
# colnames(nc_int31) <- paste0(dinucleotides, "int31")
# colnames(nc_int32) <- paste0(dinucleotides, "int32")
# colnames(nc_int33) <- paste0(dinucleotides, "int33")
# colnames(nc_int34) <- paste0(dinucleotides, "int34")
# colnames(nc_int35) <- paste0(dinucleotides, "int35")
# colnames(nc_int36) <- paste0(dinucleotides, "int36")
# colnames(nc_int37) <- paste0(dinucleotides, "int37")
# colnames(nc_int38) <- paste0(dinucleotides, "int38")
# colnames(nc_int39) <- paste0(dinucleotides, "int39")
# colnames(nc_int40) <- paste0(dinucleotides, "int40")
# colnames(nc_int41) <- paste0(dinucleotides, "int41")
# colnames(nc_int42) <- paste0(dinucleotides, "int42")
# colnames(nc_int43) <- paste0(dinucleotides, "int43")
# colnames(nc_int44) <- paste0(dinucleotides, "int44")
# colnames(nc_int45) <- paste0(dinucleotides, "int45")
# colnames(nc_int46) <- paste0(dinucleotides, "int46")
# colnames(nc_int47) <- paste0(dinucleotides, "int47")
# colnames(nc_int48) <- paste0(dinucleotides, "int48")
# colnames(nc_int49) <- paste0(dinucleotides, "int49")

# for(i in 3:47) {
#   assign(paste0("nc_di_int", i), tetranucleotides %>%
#            map(~ stri_count(dat$x50mer, regex = paste0(substr(.x,1,2), 
#                                                        paste0('.{', i-1, '}'), 
#                                                        substr(.x,3,4)))) %>%
#            bind_cols())
# }

# colnames(nc_di_int3) <- paste0(tetranucleotides, "int3")
# colnames(nc_di_int4) <- paste0(tetranucleotides, "int4")
# colnames(nc_di_int5) <- paste0(tetranucleotides, "int5")
# colnames(nc_di_int6) <- paste0(tetranucleotides, "int6")
# colnames(nc_di_int7) <- paste0(tetranucleotides, "int7")
# colnames(nc_di_int8) <- paste0(tetranucleotides, "int8")
# colnames(nc_di_int9) <- paste0(tetranucleotides, "int9")
# colnames(nc_di_int10) <- paste0(tetranucleotides, "int10")
# colnames(nc_di_int11) <- paste0(tetranucleotides, "int11")
# colnames(nc_di_int12) <- paste0(tetranucleotides, "int12")
# colnames(nc_di_int13) <- paste0(tetranucleotides, "int13")
# colnames(nc_di_int14) <- paste0(tetranucleotides, "int14")
# colnames(nc_di_int15) <- paste0(tetranucleotides, "int15")
# colnames(nc_di_int16) <- paste0(tetranucleotides, "int16")
# colnames(nc_di_int17) <- paste0(tetranucleotides, "int17")
# colnames(nc_di_int18) <- paste0(tetranucleotides, "int18")
# colnames(nc_di_int19) <- paste0(tetranucleotides, "int19")
# colnames(nc_di_int20) <- paste0(tetranucleotides, "int20")
# colnames(nc_di_int21) <- paste0(tetranucleotides, "int21")
# colnames(nc_di_int22) <- paste0(tetranucleotides, "int22")
# colnames(nc_di_int23) <- paste0(tetranucleotides, "int23")
# colnames(nc_di_int24) <- paste0(tetranucleotides, "int24")
# colnames(nc_di_int25) <- paste0(tetranucleotides, "int25")
# colnames(nc_di_int26) <- paste0(tetranucleotides, "int26")
# colnames(nc_di_int27) <- paste0(tetranucleotides, "int27")
# colnames(nc_di_int28) <- paste0(tetranucleotides, "int28")
# colnames(nc_di_int29) <- paste0(tetranucleotides, "int29")
# colnames(nc_di_int30) <- paste0(tetranucleotides, "int30")
# colnames(nc_di_int31) <- paste0(tetranucleotides, "int31")
# colnames(nc_di_int32) <- paste0(tetranucleotides, "int32")
# colnames(nc_di_int33) <- paste0(tetranucleotides, "int33")
# colnames(nc_di_int34) <- paste0(tetranucleotides, "int34")
# colnames(nc_di_int35) <- paste0(tetranucleotides, "int35")
# colnames(nc_di_int36) <- paste0(tetranucleotides, "int36")
# colnames(nc_di_int37) <- paste0(tetranucleotides, "int37")
# colnames(nc_di_int38) <- paste0(tetranucleotides, "int38")
# colnames(nc_di_int39) <- paste0(tetranucleotides, "int39")
# colnames(nc_di_int40) <- paste0(tetranucleotides, "int40")
# colnames(nc_di_int41) <- paste0(tetranucleotides, "int41")
# colnames(nc_di_int42) <- paste0(tetranucleotides, "int42")
# colnames(nc_di_int43) <- paste0(tetranucleotides, "int43")
# colnames(nc_di_int44) <- paste0(tetranucleotides, "int44")
# colnames(nc_di_int45) <- paste0(tetranucleotides, "int45")
# colnames(nc_di_int46) <- paste0(tetranucleotides, "int46")
# colnames(nc_di_int47) <- paste0(tetranucleotides, "int47")

#GC count
gc_count_random_library <- nc_1_random_library %>%
  transmute(gc_count = G + C)

dat_seq_random_library <- factor_seq_random_library %>%
  bind_cols(nc_1_random_library,
            nc_2_random_library,
            nc_3_random_library,
            gc_count_random_library,
            dat_random_library)

dat_random_library <- dat_random_library %>%
  bind_cols(factor_seq_random_library,
            factor_int_random_library,
            factor_di_int_random_library,
            nc_1_random_library,
            nc_2_random_library,
            nc_3_random_library, 
            gc_count_random_library,
            dat_melting_random_library)

# for(i in 3:49) {
#   dat <- dat %>%
#     bind_cols(get(paste0("nc_int", i)))
# }

# for(i in 3:47) {
#   dat <- dat %>%
#     bind_cols(get(paste0("nc_di_int", i)))
# }

set.seed(50)

train_indices_random_library = sample(1:nrow(dat_random_library), 
                                      nrow(dat_random_library)*.9, replace=FALSE)

dat_train_random_library = dat_random_library[train_indices_random_library,]
dat_test_random_library = dat_random_library[-train_indices_random_library,]
dat_seq_random_library = dat_seq_random_library[train_indices_random_library,]
dat_random_library = dat_train_random_library

saveRDS(dat_random_library, "data/Created/processed_random_library.rds")
saveRDS(dat_test_random_library, "data/Created/processed_test_random_library.rds")
saveRDS(dat_seq_random_library, "data/Created/processed_seq_random_library.rds")
