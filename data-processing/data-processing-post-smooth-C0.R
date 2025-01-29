library(tidyverse)

# dat_chrv_2 <- read.csv("data/Created/yeast_chrV_ir_lstm_tiling_post_smoothed_matched.csv")
dat_chrv_2 <- read.csv("data/Created/cycle6_ir_lstm_cn_tiling_post_smoothed_matched.csv")
# Remove irrelevant columns:
# dat_chrv_2 = dat_chrv_2[,c(1:17, 25:31)]
# # Remove first sequence:
# dat_chrv_2 = dat_chrv_2[-1,]
# # Remove final sequence:
# dat_chrv_2 = dat_chrv_2[-nrow(dat_chrv_2),]

# dat_chrv_2 = dat_chrv_2 %>%
#   rename(x50mer = Sequence,
#          C26 = n.26,
#          C29 = n.29,
#          C31 = n.31,
#          pred_tiling_post_smooth = model_ir_lstm_post_smoothed_C0_chrv_pred_chrV,
#          model_ir_lstm_tiling_pred_chrV = model_ir_lstm_tiling_pred_chrV,
#          ir_lstm_post_smooth_1bp = model_ir_lstm_tiling_pred_chrV_post_smooth_1bp,
#          ir_lstm_post_smooth_2bp = model_ir_lstm_tiling_pred_chrV_post_smooth_2bp,
#          ir_lstm_post_smooth_3bp = model_ir_lstm_tiling_pred_chrV_post_smooth_3bp,
#          ir_lstm_post_smooth_4bp = model_ir_lstm_tiling_pred_chrV_post_smooth_4bp,
#          ir_lstm_post_smooth_4bp_weighted = model_ir_lstm_tiling_pred_chrV_post_smooth_4bp_weighted,
#          ir_lstm_post_smooth_5bp = model_ir_lstm_tiling_pred_chrV_post_smooth_5bp,
#          DNAcycP_post_smooth_1bp = DNAcycP_pred_chrV_post_smooth_1bp,
#          DNAcycP_post_smooth_2bp = DNAcycP_pred_chrV_post_smooth_2bp,
#          DNAcycP_post_smooth_3bp = DNAcycP_pred_chrV_post_smooth_3bp,
#          DNAcycP_post_smooth_4bp = DNAcycP_pred_chrV_post_smooth_4bp,
#          DNAcycP_post_smooth_4bp_weighted = DNAcycP_pred_chrV_post_smooth_4bp_weighted,
#          DNAcycP_post_smooth_5bp = DNAcycP_pred_chrV_post_smooth_5bp,
#   )

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

dat_chrv_2 %>% map_df(~ sum(is.na(.x)))


sequence_1_df <- dat_chrv_2 %>%
  pull(x50mer) %>%
  sequence_df()

sequence_2_df <- dat_chrv_2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv_2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_3_df <- dat_chrv_2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 3, 2)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv_2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_4_df <- dat_chrv_2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 4, 3)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv_2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

sequence_5_df <- dat_chrv_2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 5, 4)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv_2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)


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


library(stringi)

nc_1 <- nucleotides %>%
  map(~ str_count(dat_chrv_2$x50mer, .x)) %>%
  bind_cols()

nc_2 <- dinucleotides %>%
  map(~ str_count(dat_chrv_2$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

nc_3 <- trinucleotides %>%
  map(~ str_count(dat_chrv_2$x50mer, paste0("(?=",.x,")"))) %>%
  bind_cols()

colnames(nc_1) <- nucleotides
colnames(nc_2) <- dinucleotides
colnames(nc_3) <- trinucleotides

gc_count <- nc_1 %>%
  transmute(gc_count = G + C)

dat_chrv_2_seq <- factor_seq %>%
  bind_cols(nc_1,
            nc_2,
            nc_3,
            gc_count,
            dat_chrv_2)

dat_chrv_2_final <- dat_chrv_2 %>%
  bind_cols(factor_seq,
            nc_1,
            nc_2,
            nc_3, 
            gc_count)


set.seed(50)

train_indices = sample(1:nrow(dat_chrv_2_final), nrow(dat_chrv_2_final)*.9, replace=FALSE)

dat_chrv_2_train = dat_chrv_2_final[train_indices,]
dat_chrv_2_test = dat_chrv_2_final[-train_indices,]
dat_chrv_2_seq = dat_chrv_2_seq[train_indices,]
# dat_final = dat_train

saveRDS(dat_chrv_2_train, "data/Created/processed_chrV_post_smooth_C0.rds")
saveRDS(dat_chrv_2_test, "data/Created/processed_chrV_post_smooth_C0_test.rds")
saveRDS(dat_chrv_2_seq, "data/Created/processed_chrV_post_smooth_C0_seq.rds")
