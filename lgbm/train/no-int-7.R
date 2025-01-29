library(tidyverse)
library(lightgbm)

kfold_obj <- readRDS("data/Created/train-10-fold-ratio.rds")

set.seed(50)

source("scripts/functions/lgbm-helper-functions.R")

nucleotides <- c("A", "C", "G", "T")
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")
ps4 <- paste0("X", 1:47, "tetra")
ps5 <- paste0("X", 1:46, "penta")

dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

thermo_feat <-c("tm1", "tm2", "tm3", "tm4", "grna_energy", "grna_scaffold_energy")

parameter_list <- list(learning_rate = c(0.1),
                       objective = c("regression"),
                       min_data_in_leaf = 100,
                       max_depth = c(12),
                       num_leaves = c(1000),
                       lambda_l2 = c(1),
                       boosting = c("gbdt")) %>%
  cross()


model_cols <- c(ps1, ps2, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count", "ratio_score_2nd")

subset_eval = parameter_list %>%
  map(~ lgb_cv_eval(kfold_obj = kfold_obj,
                    feature_cols = model_cols,
                    parameters=.x,
                    early_stop = 50,
                    n_rounds = 500,
                    feature_importance_threshold = 0))

subset_eval[[1]]$results

saveRDS(subset_eval[[1]]$results, "lgbm/results/train-no-int-7-results.rds")

saveRDS(subset_eval, "lgbm/results/train-no-int-7-subset-eval.rds")
