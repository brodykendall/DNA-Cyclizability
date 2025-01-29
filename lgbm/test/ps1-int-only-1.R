library(tidyverse)
library(lightgbm)

dat <- readRDS("data/Created/processed_ratio.rds")
dat_test <- readRDS("data/Created/processed_test_ratio.rds")

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
features <- feature_importance[feature_importance$Gain > 0.001, "Feature"]
train_set <- dat[, features]
train_rules <- lgb.convert_with_rules(train_set)
train_lgb_obj <- lgb.Dataset(data = as.matrix(train_rules$data),
                             label = dat$C0,
                             categorical_feature = handle_categorical(colnames(train_rules$data)))
lgbm_model <- lgb.train(params = parameter_list[[1]],
                        data = train_lgb_obj,
                        nrounds = 1589)

test_set <- dat_test[, features]
test_rules <- lgb.convert_with_rules(test_set)
test_lgb_obj <- lgb.Dataset(data = as.matrix(test_rules$data),
                            label = dat_test$C0,
                            categorical_feature = handle_categorical(colnames(test_rules$data)))

test.pred = predict(lgbm_model, 
                    as.matrix(test_set))

result = cor(test.pred, dat_test$C0)
saveRDS(result, "lgbm/results/test-ps1-int-only-1-cor.rds")
