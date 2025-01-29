#Hyperparameter search for LGBM model
#Kim dataset
#For use with array job
library(tidyverse)
library(modelr)
library(lightgbm)

#Lgbm helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

nucleotides <- c("A", "C", "G", "T")

#Features
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")

mono_int_indices = c(6,9,10,12,15,21)
for(i in mono_int_indices) {
  assign(paste0("int", i), paste0("X", 1:(50-i), "int", i))
}

di_int_indices = c(3,4,5,8,9,10,11,14,20)
for(i in di_int_indices) {
  assign(paste0("di_int", i), paste0("X", 1:(48-i), "di_int", i))
}

dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

thermo_feat <-c("tm1", "tm2", "tm3", "tm4", "grna_energy", "grna_scaffold_energy")

model_cols <- c(ps1, ps2, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count")

for(i in mono_int_indices) {
  model_cols <- c(model_cols, get(paste0("int", i)))
}

for(i in di_int_indices) {
  model_cols <- c(model_cols, get(paste0("di_int", i)))
}

kfold_obj <- readRDS("data/Created/train-10-fold-network-int.rds")

#Grid for random hyperparameter search
#2,772,000 combinations
# parameter_list <- list(learning_rate = c(0.01, 0.1),
#                        objective = c("regression"),
#                        min_data = c(0, 1, 5, 10, 20, 50, 100, 200, 300),
#                        max_bin = c(60, 120, 255, 300, 500),
#                        max_depth = c(2, 3, 4, 5, 6, 12, -1),
#                        num_leaves = c(2, 5, 10, 20, 30, 60, 120, 250, 500, 1000, 2000),
#                        feature_fraction = seq(0.1, 1, by = 0.1),
#                        bagging_fraction = seq(0.1, 1, by = 0.1),
#                        lambda_l2 = c(0, 0.1, 1, 10),
#                        boosting = c("gbdt"),
#                        num_threads = 14) %>%
# cross()

param_init <- list(learning_rate = c(0.1),
                   objective = c("regression"),
                   min_data = c(0),
                   max_bin = c(500),
                   max_depth = c(12),
                   num_leaves = c(500),
                   feature_fraction = c(0.5),
                   bagging_fraction = c(0.5),
                   lambda_l2 = c(1),
                   boosting = c("gbdt"),
                   num_threads = 14)

#Random subset for evaluation
# parameter_subset <- parameter_list[sample.int(n = length(parameter_list), size = 10000, replace = FALSE)]

# rm(parameter_list)

# #20 jobs, 300 iterations each
# #Subset for current job
# parameter_subset <- parameter_subset[(arg_index*300-299):(arg_index*300)]

# #20 jobs, 10 iterations each
# #Subset for current job
# parameter_subset <- parameter_subset[(arg_index*10-9):(arg_index*10)]

#20 jobs, 100 iterations each
#Subset for current job
# parameter_subset <- parameter_subset[(arg_index*100-99):(arg_index*100)]

# parameter_subset_1 <- parameter_subset[1]

# subset_eval <- param_init %>%
#   map(~ lgb_cv_eval(kfold_obj = kfold_obj,
#                     feature_cols = model_cols,
#                     parameters = .x,
#                     n_rounds = 6000,
#                     early_stop = 50,
#                     feature_importance_threshold = 0))

init_eval <- lgb_cv_eval(kfold_obj = kfold_obj,
                    feature_cols = model_cols,
                    parameters = param_init,
                    n_rounds = 6000,
                    early_stop = 50,
                    feature_importance_threshold = 0)

write_csv(init_eval$results, "tuning/hyperparameter-results-feat-sel-0-no-ratio.csv")
