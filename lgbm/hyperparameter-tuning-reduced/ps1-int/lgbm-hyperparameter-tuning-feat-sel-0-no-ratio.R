#Hyperparameter search for LGBM model
#Kim dataset
#For use with array job
library(tidyverse)
library(modelr)
library(lightgbm)

#Argument indicates which subset to evaluate in job
args <- commandArgs(TRUE)
arg_index <- as.integer(args[1])

#Lgbm helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

nucleotides <- c("A", "C", "G", "T")

#Features
ps1 <- paste0("X", 1:50, "mono")

ps_int1 = character()
mono_some_int_indices = 3:49
for(i in mono_some_int_indices) {
  ps_int1 <- c(ps_int1, paste0("X", 1:(50-i), "int", i))
}

model_cols <- c(ps1, ps_int1)

kfold_obj <- readRDS("data/Created/train-10-fold-ratio.rds")

#Grid for random hyperparameter search
#2,772,000 combinations
parameter_list <- list(learning_rate = c(0.01, 0.1),
                       objective = c("regression"),
                       min_data = c(0, 1, 5, 10, 20, 50, 100, 200, 300),
                       max_bin = c(60, 120, 255, 300, 500),
                       max_depth = c(2, 3, 4, 5, 6, 12, -1),
                       num_leaves = c(2, 5, 10, 20, 30, 60, 120, 250, 500, 1000, 2000),
                       feature_fraction = seq(0.1, 1, by = 0.1),
                       bagging_fraction = seq(0.1, 1, by = 0.1),
                       lambda_l2 = c(0, 0.1, 1, 10),
                       boosting = c("gbdt"),
                       num_threads = 14) %>%
  cross()

#Random subset for evaluation
parameter_subset <- parameter_list[sample.int(n = length(parameter_list), size = 10000, replace = FALSE)]

rm(parameter_list)

# #20 jobs, 300 iterations each
# #Subset for current job
# parameter_subset <- parameter_subset[(arg_index*300-299):(arg_index*300)]

#20 jobs, 10 iterations each
#Subset for current job
parameter_subset <- parameter_subset[(arg_index*10-9):(arg_index*10)]

#20 jobs, 100 iterations each
#Subset for current job
# parameter_subset <- parameter_subset[(arg_index*100-99):(arg_index*100)]

subset_eval <- parameter_subset %>%
  map(~ lgb_cv_eval(kfold_obj = kfold_obj,
                    feature_cols = model_cols,
                    parameters = .x,
                    n_rounds = 6000,
                    early_stop = 50,
                    feature_importance_threshold = 0))

subset_results <- subset_eval %>%
  map(~ .x$results) %>%
  bind_rows()

write_csv(subset_results,
          paste0("tuning-reduced/ps1-int/hyperparameter-results-feat-sel-0-no-ratio-", arg_index, ".csv"))
