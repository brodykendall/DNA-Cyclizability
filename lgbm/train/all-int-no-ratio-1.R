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

ps_int1 = character()
mono_some_int_indices = 3:49
for(i in mono_some_int_indices) {
  ps_int1 <- c(ps_int1, paste0("X", 1:(50-i), "int", i))
}

ps_int2 = character()
di_some_int_indices = 3:47
for(i in di_some_int_indices) {
  ps_int2 <- c(ps_int2, paste0("X", 1:(48-i), "di_int", i))
}

parameter_list <- list(learning_rate = c(0.1),
                       objective = c("regression"),
                       min_data_in_leaf = 100,
                       max_bin = 120,
                       max_depth = c(4),
                       num_leaves = c(20),
                       feature_fraction = 0.5,
                       bagging_fraction = 0.4,
                       lambda_l2 = c(1),
                       boosting = c("gbdt")) %>%
  cross()


model_cols <- c(ps1, ps2, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, 
                thermo_feat, "gc_count", ps_int1, ps_int2)

subset_eval = parameter_list %>%
  map(~ lgb_cv_eval(kfold_obj = kfold_obj,
                    feature_cols = model_cols,
                    parameters=.x,
                    early_stop = 50,
                    n_rounds = 6000,
                    feature_importance_threshold = 0.001))

subset_eval[[1]]$results

saveRDS(subset_eval[[1]]$results, "lgbm/results/train-all-int-no-ratio-1-results.rds")

saveRDS(subset_eval, "lgbm/results/train-all-int-no-ratio-1-subset-eval.rds")
