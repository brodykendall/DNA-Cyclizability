# Model Validation

#Sequence functions
source("scripts/functions/sequence-functions.R")

#Markov feature functions
source("scripts/functions/markov-functions.R")

#LGBM helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

tuning_files_all <- list.files("./tuning-reduced/ps1-int", full.names = TRUE)

tuning_files_feat_sel_0 <- grep("hyperparameter-results-feat-sel-0-[0-9]+.csv", tuning_files_all, value = TRUE)
tuning_files_feat_sel_0_no_ratio <- grep("hyperparameter-results-feat-sel-0-no-ratio-[0-9]+.csv", tuning_files_all, value = TRUE)
tuning_files_feat_sel_001 <- grep("hyperparameter-results-feat-sel-001-[0-9]+.csv", tuning_files_all, value = TRUE)
tuning_files_feat_sel_001_no_ratio <- grep("hyperparameter-results-feat-sel-001-no-ratio-[0-9]+.csv", tuning_files_all, value = TRUE)
tuning_files_no_feat_sel <- grep("hyperparameter-results-no-feat-sel-[0-9]+.csv", tuning_files_all, value = TRUE)
tuning_files_no_feat_sel_no_ratio <- grep("hyperparameter-results-no-feat-sel-no-ratio-[0-9]+.csv", tuning_files_all, value = TRUE)

tuning_results_feat_sel_0 <- tuning_files_feat_sel_0 %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "yes",
         feature_selection_level = "0")

tuning_results_feat_sel_0_no_ratio <- tuning_files_feat_sel_0_no_ratio %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "no",
         feature_selection_level = "0")

tuning_results_feat_sel_001 <- tuning_files_feat_sel_001 %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "yes",
         feature_selection_level = "001")

tuning_results_feat_sel_001_no_ratio <- tuning_files_feat_sel_001_no_ratio %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "no",
         feature_selection_level = "001")

tuning_results_no_feat_sel <- tuning_files_no_feat_sel %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "yes",
         feature_selection_level = "none")

tuning_results_no_feat_sel_no_ratio <- tuning_files_no_feat_sel_no_ratio %>%
  map(read_csv, col_types = cols(
    best_iterations = col_character(),
    feat_freq = col_character())) %>%
  bind_rows() %>%
  mutate(inc_ratio = "no",
         feature_selection_level = "none")

tuning_results <- tuning_results_feat_sel_0 %>%
  bind_rows(tuning_results_feat_sel_0_no_ratio,
            tuning_results_feat_sel_001,
            tuning_results_feat_sel_001_no_ratio,
            tuning_results_no_feat_sel,
            tuning_results_no_feat_sel_no_ratio) %>%
  arrange(desc(avg_pearson))

top_parameters <- tuning_results %>%
  filter(avg_pearson == max(tuning_results$avg_pearson))

top_parameters
#ratio included, feature selection level none

optim_parameters <- list(learning_rate = 0.01,
                         objective = "regression",
                         min_data = 50,
                         max_bin = 500,
                         max_depth = 6,
                         num_leaves = 500,
                         feature_fraction = 0.1,
                         bagging_fraction = 0.1,
                         lambda_l2 = 1,
                         boosting = "gbdt",
                         num_threads = 14)

iterations <- ceiling(top_parameters$avg_iterations)

ps1 <- paste0("X", 1:50, "mono")

ps_int1 = character()
mono_some_int_indices = 3:49
for(i in mono_some_int_indices) {
  ps_int1 <- c(ps_int1, paste0("X", 1:(50-i), "int", i))
}

model_cols <- c(ps1, ps_int1, "ratio_score_2nd")

write_csv(tibble(features = model_cols), "model/reduced/ps1-int/all-features.csv")
write_csv(tibble(iterations = iterations), "model/reduced/ps1-int/iterations.csv")
saveRDS(optim_parameters, "model/reduced/ps1-int/optim_parameters.rds")

kfold_obj <- readRDS("data/Created/train-10-fold-ratio.rds")

kfold_obj_with_predictions <- lgb_cv_eval_2(kfold_obj = kfold_obj,
                                            feature_cols = model_cols,
                                            parameters = optim_parameters,
                                            n_rounds = 6000,
                                            early_stop = 50,
                                            feature_importance_threshold = NA)

boostmec_confidence_intervals <- kfold_obj_with_predictions %>%
  mutate(test_set_n = map(test, ~nrow(.x))) %>%
  select(pearson_cor, test_set_n) %>%
  mutate_all(unlist) %>%
  mutate(lower_bound = tanh(atanh(pearson_cor) - qnorm(0.975)/sqrt(test_set_n - 3)),
         upper_bound = tanh(atanh(pearson_cor) + qnorm(0.975)/sqrt(test_set_n - 3)))

boostmec_confidence_intervals

mean(boostmec_confidence_intervals$pearson_cor)
