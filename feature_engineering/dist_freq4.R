library(tidyverse)
library(lightgbm)
library(stringi)
source("scripts/functions/lgbm-helper-functions.R")
library(Matrix)
library(glmnet)

dat = readRDS("data/Created/processed_tiling_newC0.rds")
y = dat$C0_new
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
y_test = dat_test$C0_new
dat_random = readRDS("data/Created/processed_random_newC0.rds")
y_random = dat_random$C0_new
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")
y_random_test = dat_random_test$C0_new
dat_random_all = rbind(dat_random, dat_random_test)
y_random_all = c(y_random, y_random_test)


nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

distances = 2:48

dist_freq_AorTdi_CorGdi_bid2 = map(distances, 
                                   ~stri_count(dat$x50mer, 
                                               regex=paste0("(?=([AT]{2}.{",.x-2,"}[CG]{2})|([CG]{2}.{",.x-2,"}[AT]{2}))"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_CorGdi_bid2) = paste0("AorT_CorG_dist", distances)

dist_freq_AorTdi_AorTdi_bid2 = map(distances, 
                                   ~stri_count(dat$x50mer, 
                                               regex=paste0("(?=[AT]{2}.{",.x-2,"}[AT]{2})"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_AorTdi_bid2) = paste0("AorT_AorT_dist", distances)

dist_freq_CorGdi_CorGdi_bid2 = map(distances, 
                                   ~stri_count(dat$x50mer, 
                                               regex=paste0("(?=[CG]{2}.{",.x-2,"}[CG]{2})"))) %>%
  bind_cols()

colnames(dist_freq_CorGdi_CorGdi_bid2) = paste0("CorG_CorG_dist", distances)

saveRDS(dist_freq_AorTdi_AorTdi_bid2, "data/Created/tiling_dist_freq_AorTdi_AorTdi_bid2.rds")
saveRDS(dist_freq_AorTdi_CorGdi_bid2, "data/Created/tiling_dist_freq_AorTdi_CorGdi_bid2.rds")
saveRDS(dist_freq_CorGdi_CorGdi_bid2, "data/Created/tiling_dist_freq_CorGdi_CorGdi_bid2.rds")

# Test data:
dist_freq_AorTdi_CorGdi_bid2_test = map(distances, 
                                        ~stri_count(dat_test$x50mer, 
                                                    regex=paste0("(?=([AT]{2}.{",.x-2,"}[CG]{2})|([CG]{2}.{",.x-2,"}[AT]{2}))"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_CorGdi_bid2_test) = paste0("AorT_CorG_dist", distances)

dist_freq_AorTdi_AorTdi_bid2_test = map(distances, 
                                        ~stri_count(dat_test$x50mer, 
                                                    regex=paste0("(?=[AT]{2}.{",.x-2,"}[AT]{2})"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_AorTdi_bid2_test) = paste0("AorT_AorT_dist", distances)

dist_freq_CorGdi_CorGdi_bid2_test = map(distances, 
                                        ~stri_count(dat_test$x50mer, 
                                                    regex=paste0("(?=[CG]{2}.{",.x-2,"}[CG]{2})"))) %>%
  bind_cols()

colnames(dist_freq_CorGdi_CorGdi_bid2_test) = paste0("CorG_CorG_dist", distances)

saveRDS(dist_freq_AorTdi_AorTdi_bid2_test, "data/Created/tiling_dist_freq_AorTdi_AorTdi_bid2_test.rds")
saveRDS(dist_freq_AorTdi_CorGdi_bid2_test, "data/Created/tiling_dist_freq_AorTdi_CorGdi_bid2_test.rds")
saveRDS(dist_freq_CorGdi_CorGdi_bid2_test, "data/Created/tiling_dist_freq_CorGdi_CorGdi_bid2_test.rds")


# Random Library:
dist_freq_AorTdi_CorGdi_bid2_random_all = map(distances, 
                                              ~stri_count(dat_random_all$x50mer, 
                                                          regex=paste0("(?=([AT]{2}.{",.x-2,"}[CG]{2})|([CG]{2}.{",.x-2,"}[AT]{2}))"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_CorGdi_bid2_random_all) = paste0("AorT_CorG_dist", distances)

dist_freq_AorTdi_AorTdi_bid2_random_all = map(distances, 
                                              ~stri_count(dat_random_all$x50mer, 
                                                          regex=paste0("(?=[AT]{2}.{",.x-2,"}[AT]{2})"))) %>%
  bind_cols()

colnames(dist_freq_AorTdi_AorTdi_bid2_random_all) = paste0("AorT_AorT_dist", distances)

dist_freq_CorGdi_CorGdi_bid2_random_all = map(distances, 
                                              ~stri_count(dat_random_all$x50mer, 
                                                          regex=paste0("(?=[CG]{2}.{",.x-2,"}[CG]{2})"))) %>%
  bind_cols()

colnames(dist_freq_CorGdi_CorGdi_bid2_random_all) = paste0("CorG_CorG_dist", distances)

saveRDS(dist_freq_AorTdi_AorTdi_bid2_random_all, "data/Created/random_all_dist_freq_AorTdi_AorTdi_bid2.rds")
saveRDS(dist_freq_AorTdi_CorGdi_bid2_random_all, "data/Created/random_all_dist_freq_AorTdi_CorGdi_bid2.rds")
saveRDS(dist_freq_CorGdi_CorGdi_bid2_random_all, "data/Created/random_all_dist_freq_CorGdi_CorGdi_bid2.rds")



