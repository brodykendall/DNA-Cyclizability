library(tidyverse)
library(lightgbm)

kfold_obj <- readRDS("data/Created/train-10-fold-ratio.rds")

set.seed(50)

source("scripts/functions/lgbm-helper-functions.R")

ps1 <- paste0("X", 1:50, "mono")

ps_int1 = character()
mono_some_int_indices = 3:49
for(i in mono_some_int_indices) {
  ps_int1 <- c(ps_int1, paste0("X", 1:(50-i), "int", i))
}

parameter_list <- list(learning_rate = c(0.1),
                       objective = c("regression"),
                       max_bin = 120,
                       min_data_in_leaf = 100,
                       max_depth = c(4),
                       num_leaves = c(20),
                       feature_fraction = 0.5,
                       bagging_fraction = 0.4,
                       lambda_l2 = c(1),
                       boosting = c("gbdt")) %>%
  cross()


model_cols <- c(ps1, "ratio_score_2nd", ps_int1)

subset_eval = parameter_list %>%
  map(~ lgb_cv_eval(kfold_obj = kfold_obj,
                    feature_cols = model_cols,
                    parameters=.x,
                    early_stop = 50,
                    n_rounds = 6000,
                    feature_importance_threshold = 0.001))

subset_eval[[1]]$results

saveRDS(subset_eval, "lgbm/results/train-ps1-int-only-1-subset-eval.rds")
