library(tidyverse)
library(TmCalculator)

cycle5 <- read_csv("cycle5.txt")

dat <- cycle5

free_energy <- read_csv("data/Created/kim_min_free_energy.csv")
free_energy <- free_energy[1:82368,]

dat <- dat %>%
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

dat %>% map_df(~ sum(is.na(.x)))

dat <- dat %>%
  dplyr::select(x50mer = Sequence,
                C0,
                grna_energy, grna_scaffold_energy) %>%
  mutate(x41mer = substring(x50mer, 5, 45))

dat_melting <- dat %>%
  group_by(1:n()) %>%
  transmute(tm1 = Tm_Wallace(x41mer)$Tm,
            tm2 = Tm_Wallace(substring(x41mer, 1, 4))$Tm,
            tm3 = Tm_Wallace(substring(x41mer, 5, 12))$Tm,
            tm4 = Tm_Wallace(substring(x41mer, 36, 40))$Tm) %>%
  ungroup() %>%
  select(-`1:n()`)

sequence_1_df <- dat %>%
  pull(x50mer) %>%
  sequence_df()

sequence_2_df <- dat %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_3_df <- dat %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 3, 2)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_4_df <- dat %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 4, 3)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_5_df <- dat %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 5, 4)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

for(i in 3:49) {
  assign(paste0("interaction_", i, "_df"), dat %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap2(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

for(i in 3:47) {
  assign(paste0("di_interaction_", i, "_df"), dat %>%
           pull(x50mer) %>%
           map(~unlist(strsplit(.x, ""))) %>%
           map(~splitWithOverlap3(.x, i+1, i)) %>%
           map_depth(2, paste0, collapse = "") %>%
           map(unlist) %>%
           unlist() %>%
           {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
           data.frame(stringsAsFactors = FALSE))
}

sequence_1_factor <- sequence_1_df %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))

sequence_2_factor <- sequence_2_df %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_3_factor <- sequence_3_df %>%
  map_df(~ factor(.x, levels = trinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tri")) %>%
  select(-X49tri, -X50tri)

sequence_4_factor <- sequence_4_df %>%
  map_df(~ factor(.x, levels = tetranucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "tetra")) %>%
  select(-X48tetra, -X49tetra, -X50tetra)

sequence_5_factor <- sequence_5_df %>%
  map_df(~ factor(.x, levels = pentanucleotides)) %>%
  as.data.frame() %>%
  rename_all(~ paste0(.x, "penta")) %>%
  select(-X47penta, -X48penta, -X49penta, -X50penta)

factor_seq <- sequence_1_factor %>%
  bind_cols(sequence_2_factor, sequence_3_factor, sequence_4_factor, sequence_5_factor)

for(i in 3:49) {
  assign(paste0("interaction_", i, "_factor"), get(paste0("interaction_", i, "_df")) %>%
           map_df(~ factor(.x, levels = dinucleotides)) %>%
           as.data.frame() %>%
           rename_all(~paste0(.x, "int", i)))
}
interaction_49_factor <- interaction_49_factor %>%
  rename(X1int49 = .int49)

factor_int <- interaction_3_factor
for(i in 4:49) {
  factor_int <- factor_int %>%
    bind_cols(get(paste0("interaction_", i, "_factor")))
}

library(stringi)

nc_1 <- nucleotides %>%
  map(~ str_count(dat$x50mer, .x)) %>%
  bind_cols()

nc_2 <- dinucleotides %>%
  map(~ str_count(dat$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

nc_3 <- trinucleotides %>%
  map(~ str_count(dat$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

colnames(nc_1) <- nucleotides
colnames(nc_2) <- dinucleotides
colnames(nc_3) <- trinucleotides

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

#GC count
gc_count <- nc_1 %>%
  transmute(gc_count = G + C)

dat_seq <- factor_seq %>%
  bind_cols(nc_1,
            nc_2,
            nc_3,
            gc_count,
            dat)

dat <- dat %>%
  bind_cols(factor_seq,
            factor_int,
            nc_1,
            nc_2,
            nc_3, 
            gc_count,
            dat_melting)

# for(i in 3:49) {
#   dat <- dat %>%
#     bind_cols(get(paste0("nc_int", i)))
# }

set.seed(50)

train_indices = sample(1:nrow(dat), nrow(dat)*.9, replace=FALSE)

dat_train = dat[train_indices,]
dat_test = dat[-train_indices,]
dat = dat_train

saveRDS(dat, "data/Created/processed_single_int_only.rds")
saveRDS(dat_test, "data/Created/processed_single_int_only_test.rds")