library(tidyverse)
library(lightgbm)

dat <- readRDS("data/Created/processed_ratio.rds")
dat_test <- readRDS("data/Created/processed_test_ratio.rds")

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
                       min_data_in_leaf = 3000,
                       max_depth = c(7),
                       num_leaves = c(50),
                       lambda_l2 = c(1),
                       boosting = c("gbdt")) %>%
  cross()


model_cols <- c(ps1, ps2, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count", "ratio_score_2nd")

initial_train_set <- dat[,model_cols]
initial_train_rules <- lgb.convert_with_rules(initial_train_set)
initial_train_data <- initial_train_rules$data
initial_train_lgb_obj <- lgb.Dataset(data = as.matrix(initial_train_data),
                                     label = dat$C0,
                                     categorical_feature = handle_categorical(model_cols))
initial_mod = lgb.train(params = parameter_list[[1]],
                        data = initial_train_lgb_obj,
                        nrounds = 500)
feature_importance <- as.data.frame(lgb.importance(initial_mod))
features <- feature_importance[feature_importance$Gain > 0.1, "Feature"]
train_set <- dat[, features]
train_rules <- lgb.convert_with_rules(train_set)
train_lgb_obj <- lgb.Dataset(data = as.matrix(train_rules$data),
                             label = dat$C0,
                             categorical_feature = handle_categorical(colnames(train_rules$data)))
lgbm_model <- lgb.train(params = parameter_list[[1]],
                        data = train_lgb_obj,
                        nrounds = 24)

test_set <- dat_test[, features]
test_rules <- lgb.convert_with_rules(test_set)
test_lgb_obj <- lgb.Dataset(data = as.matrix(test_rules$data),
                            label = dat_test$C0,
                            categorical_feature = handle_categorical(colnames(test_rules$data)))

test.pred = predict(lgbm_model, 
                    as.matrix(test_set))

result = cor(test.pred, dat_test$C0)
saveRDS(result, "lgbm/results/test-no-int-1-cor.rds")