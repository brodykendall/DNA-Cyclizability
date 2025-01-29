# Load data:
dat = readRDS("data/Created/processed_tiling_newC0.rds")
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
dat_random = readRDS("data/Created/processed_random_newC0.rds")
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")

library(mgcv)
library(tidyverse)





########################################################################
# IR_LSTM_NEWC0
########################################################################

dat_first_conv = readRDS("data/Created/tiling_tiling_first_conv_output_train.rds")
dat_first_conv_test = readRDS("data/Created/tiling_tiling_first_conv_output_test.rds")
dat_first_conv_random_test = readRDS("data/Created/tiling_random_first_conv_output_test.rds")

dat_first_conv$Y = dat$C0_new
dat_first_conv_test$Y = dat_test$C0_new

lm.model_first_conv = lm(Y~., data = dat_first_conv)

# Prediction accuracy with train data:
cor(lm.model_first_conv$fitted.values, dat$C0_new, method="pearson")
# 0.4238779
mean((lm.model_first_conv$fitted.values - dat$C0_new)^2)
# 0.1988036

# Predict values for test dataset:
lm.model_first_conv.pred = predict(lm.model_first_conv, newdata = dat_first_conv_test)

# Prediction accuracy with test data:
cor(lm.model_first_conv.pred, dat_test$C0_new, method="pearson")
# 0.377503
mean((lm.model_first_conv.pred - dat_test$C0_new)^2)
# 0.2074552

# Predict values for random test dataset:
lm.model_first_conv.pred_random = predict(lm.model_first_conv, newdata = dat_first_conv_random_test)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.pred_random, dat_random_test$C0_new, method="pearson")
# 0.3590254
mean((lm.model_first_conv.pred_random - dat_random_test$C0_new)^2)
# 0.138732




# First in parallel line 1
# Load data:
dat_second_conv = readRDS("data/Created/tiling_tiling_second_conv_output_train.rds")
dat_second_conv_test = readRDS("data/Created/tiling_tiling_second_conv_output_test.rds")

dat_second_conv_random_test = readRDS("data/Created/tiling_random_second_conv_output_test.rds")


dat_second_conv$Y = dat$C0_new
dat_second_conv_test$Y = dat_test$C0_new

lm.model_second_conv = lm(Y~., data = dat_second_conv)

# Prediction accuracy with train data:
cor(lm.model_second_conv$fitted.values, dat$C0_new, method="pearson")
# 0.4897099
mean((lm.model_second_conv$fitted.values - dat$C0_new)^2)
# 0.184228

# Predict values for test dataset:
lm.model_second_conv.pred = predict(lm.model_second_conv, newdata = dat_second_conv_test)

# Prediction accuracy with test data:
cor(lm.model_second_conv.pred, dat_test$C0_new, method="pearson")
# 0.4589663
mean((lm.model_second_conv.pred - dat_test$C0_new)^2)
# 0.1904921

# Predict values for random test dataset:
lm.model_second_conv.pred_random = predict(lm.model_second_conv, newdata = dat_second_conv_random_test)

# Prediction accuracy with random test data:
cor(lm.model_second_conv.pred_random, dat_random_test$C0_new, method="pearson")
# 0.4505374
mean((lm.model_second_conv.pred_random - dat_random_test$C0_new)^2)
# 0.1233665






# Second in parallel line 1
# Load data:
dat_third_conv = readRDS("data/Created/tiling_tiling_third_conv_output_train.rds")
dat_third_conv_test = readRDS("data/Created/tiling_tiling_third_conv_output_test.rds")

dat_third_conv_random_test = readRDS("data/Created/tiling_random_third_conv_output_test.rds")

dat_third_conv$Y = dat$C0_new
dat_third_conv_test$Y = dat_test$C0_new

lm.model_third_conv = lm(Y~., data = dat_third_conv)

# Prediction accuracy with train data:
cor(lm.model_third_conv$fitted.values, dat$C0_new, method="pearson")
# 0.5073504
mean((lm.model_third_conv$fitted.values - dat$C0_new)^2)
# 0.1799655

# Predict values for test dataset:
lm.model_third_conv.pred = predict(lm.model_third_conv, newdata = dat_third_conv_test)

# Prediction accuracy with test data:
cor(lm.model_third_conv.pred, dat_test$C0_new, method="pearson")
# 0.4278708
mean((lm.model_third_conv.pred - dat_test$C0_new)^2)
# 0.2021132

# Predict values for random test dataset:
lm.model_third_conv.pred_random = predict(lm.model_third_conv, newdata = dat_third_conv_random_test)

# Prediction accuracy with random test data:
cor(lm.model_third_conv.pred_random, dat_random_test$C0_new, method="pearson")
# 0.477897
mean((lm.model_third_conv.pred_random - dat_random_test$C0_new)^2)
# 0.1190561




# First in parallel line 2
# Load data:
dat_fourth_conv = readRDS("data/Created/tiling_tiling_fourth_conv_output_train.rds")
dat_fourth_conv_test = readRDS("data/Created/tiling_tiling_fourth_conv_output_test.rds")

dat_fourth_conv_random_test = readRDS("data/Created/tiling_random_fourth_conv_output_test.rds")

dat_fourth_conv$Y = dat$C0_new
dat_fourth_conv_test$Y = dat_test$C0_new

lm.model_fourth_conv = lm(Y~., data = dat_fourth_conv)

# Prediction accuracy with train data:
cor(lm.model_fourth_conv$fitted.values, dat$C0_new, method="pearson")
# 0.6864899
mean((lm.model_fourth_conv$fitted.values - dat$C0_new)^2)
# 0.1281363

# Predict values for test dataset:
lm.model_fourth_conv.pred = predict(lm.model_fourth_conv, newdata = dat_fourth_conv_test)

# Prediction accuracy with test data:
cor(lm.model_fourth_conv.pred, dat_test$C0_new, method="pearson")
# 0.6871318
mean((lm.model_fourth_conv.pred - dat_test$C0_new)^2)
# 0.1272171

# Predict values for random test dataset:
lm.model_fourth_conv.pred_random = predict(lm.model_fourth_conv, newdata = dat_fourth_conv_random_test)

# Prediction accuracy with random test data:
cor(lm.model_fourth_conv.pred_random, dat_random_test$C0_new, method="pearson")
# 0.6767955
mean((lm.model_fourth_conv.pred_random - dat_random_test$C0_new)^2)
# 0.08465991







# Second in parallel line 2
# Load data:
dat_fifth_conv = readRDS("data/Created/tiling_tiling_fifth_conv_output_train.rds")
dat_fifth_conv_test = readRDS("data/Created/tiling_tiling_fifth_conv_output_test.rds")

dat_fifth_conv_random_test = readRDS("data/Created/tiling_random_fifth_conv_output_test.rds")

dat_fifth_conv$Y = dat$C0_new
dat_fifth_conv_test$Y = dat_test$C0_new

lm.model_fifth_conv = lm(Y~., data = dat_fifth_conv)

# Prediction accuracy with train data:
cor(lm.model_fifth_conv$fitted.values, dat$C0_new, method="pearson")
# 0.8901451
mean((lm.model_fifth_conv$fitted.values - dat$C0_new)^2)
# 0.05032127

# Predict values for test dataset:
lm.model_fifth_conv.pred = predict(lm.model_fifth_conv, newdata = dat_fifth_conv_test)

# Prediction accuracy with test data:
cor(lm.model_fifth_conv.pred, dat_test$C0_new, method="pearson")
# 0.8924833 (!!!!)
mean((lm.model_fifth_conv.pred - dat_test$C0_new)^2)
# 0.04905821

# Predict values for random test dataset:
lm.model_fifth_conv.pred_random = predict(lm.model_fifth_conv, newdata = dat_fifth_conv_random_test)

# Prediction accuracy with random test data:
cor(lm.model_fifth_conv.pred_random, dat_random_test$C0_new, method="pearson")
# 0.9025493
mean((lm.model_fifth_conv.pred_random - dat_random_test$C0_new)^2)
# 0.02930911






# Combine outputs from first layer with 2 parallel lines
# Load data:
dat_add = readRDS("data/Created/tiling_tiling_add_output_train.rds")
dat_add_test = readRDS("data/Created/tiling_tiling_add_output_test.rds")

dat_add_random_test = readRDS("data/Created/tiling_random_add_output_test.rds")

dat_add$Y = dat$C0_new
dat_add_test$Y = dat_test$C0_new

lm.model_add = lm(Y~., data = dat_add)

# Prediction accuracy with train data:
cor(lm.model_add$fitted.values, dat$C0_new, method="pearson")
# 0.9096397
mean((lm.model_add$fitted.values - dat$C0_new)^2)
# 0.04181829

# Predict values for test dataset:
lm.model_add.pred = predict(lm.model_add, newdata = dat_add_test)

# Prediction accuracy with test data:
cor(lm.model_add.pred, dat_test$C0_new, method="pearson")
# 0.907292
mean((lm.model_add.pred - dat_test$C0_new)^2)
# 0.04261492

# Predict values for random test dataset:
lm.model_add.pred_random = predict(lm.model_add, newdata = dat_add_random_test)

# Prediction accuracy with random test data:
cor(lm.model_add.pred_random, dat_random_test$C0_new, method="pearson")
# 0.9195114
mean((lm.model_add.pred_random - dat_random_test$C0_new)^2)
# 0.02438631










# LSTM output
# Load data:
dat_lstm = readRDS("data/Created/tiling_tiling_lstm_output_train.rds")
dat_lstm_test = readRDS("data/Created/tiling_tiling_lstm_output_test.rds")

dat_lstm_random_test = readRDS("data/Created/tiling_random_lstm_output_test.rds")

dat_lstm$Y = dat$C0_new
dat_lstm_test$Y = dat_test$C0_new

lm.model_lstm = lm(Y~., data = dat_lstm)

# Prediction accuracy with train data:
cor(lm.model_lstm$fitted.values, dat$C0_new, method="pearson")
# 0.7405922
mean((lm.model_lstm$fitted.values - dat$C0_new)^2)
# 0.1094251

# Predict values for test dataset:
lm.model_lstm.pred = predict(lm.model_lstm, newdata = dat_lstm_test)

# Prediction accuracy with test data:
cor(lm.model_lstm.pred, dat_test$C0_new, method="pearson")
# 0.7421592
mean((lm.model_lstm.pred - dat_test$C0_new)^2)
# 0.1082678

# Predict values for random test dataset:
lm.model_lstm.pred_random = predict(lm.model_lstm, newdata = dat_lstm_random_test)

# Prediction accuracy with random test data:
cor(lm.model_lstm.pred_random, dat_random_test$C0_new, method="pearson")
# 0.7727764
mean((lm.model_lstm.pred_random - dat_random_test$C0_new)^2)
# 0.06249908












########################################################################
# IR_LSTM_NEWC0 + Fourier
########################################################################

# Set up models with fourier features included
dat_first_conv_fourier_sq.1 = bind_cols(dat_first_conv, Xone_AT_fourier_cos^2, Xone_AT_fourier_sin^2)

lm.model_first_conv_fourier_sq.1 = lm(Y~., data = dat_first_conv_fourier_sq.1)

# Prediction accuracy with train data:
cor(lm.model_first_conv_fourier_sq.1$fitted.values, dat$C0_new, method="pearson")
# 0.5850556
mean((lm.model_first_conv_fourier_sq.1$fitted.values - dat$C0_new)^2)
# 0.1593938

# Set up test data
dat_first_conv_test.fourier_sq.1 = bind_cols(dat_first_conv_test, Xone_AT_fourier_cos_test^2, Xone_AT_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv_test.fourier_sq.1.pred = predict(lm.model_first_conv_fourier_sq.1, newdata = dat_first_conv_test.fourier_sq.1)

# Prediction accuracy with test data:
cor(lm.model_first_conv_test.fourier_sq.1.pred, dat_test$C0_new, method="pearson")
# 0.5582504
mean((lm.model_first_conv_test.fourier_sq.1.pred - dat_test$C0_new)^2)
# 0.1661409

# Set up random test data
dat_first_conv_random_test.fourier_sq.1 = bind_cols(dat_first_conv_random_test, Xone_AT_random_fourier_cos_test^2, Xone_AT_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.fourier_sq.1.pred_random = predict(lm.model_first_conv_fourier_sq.1, newdata = dat_first_conv_random_test.fourier_sq.1)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.fourier_sq.1.pred_random, dat_random_test$C0_new, method="pearson")
# 0.5257918
mean((lm.model_first_conv.fourier_sq.1.pred_random - dat_random_test$C0_new)^2)
# 0.1182379



#FIXME
# First in parallel line 2
# Set up models with fourier features included
dat_fourth_conv_fourier_sq.1 = bind_cols(dat_fourth_conv, Xone_AT_fourier_cos^2, Xone_AT_fourier_sin^2)

lm.model_fourth_conv_fourier_sq.1 = lm(Y~., data = dat_fourth_conv_fourier_sq.1)

# Prediction accuracy with train data:
cor(lm.model_fourth_conv_fourier_sq.1$fitted.values, dat$C0_new, method="pearson")
# 0.7195142
mean((lm.model_fourth_conv_fourier_sq.1$fitted.values - dat$C0_new)^2)
# 0.1168836

# Set up test data
dat_fourth_conv_test.fourier_sq.1 = bind_cols(dat_fourth_conv_test, Xone_AT_fourier_cos_test^2, Xone_AT_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_fourth_conv_test.fourier_sq.1.pred = predict(lm.model_fourth_conv_fourier_sq.1, newdata = dat_fourth_conv_test.fourier_sq.1)

# Prediction accuracy with test data:
cor(lm.model_fourth_conv_test.fourier_sq.1.pred, dat_test$C0_new, method="pearson")
# 0.718429
mean((lm.model_fourth_conv_test.fourier_sq.1.pred - dat_test$C0_new)^2)
# 0.1166149

# Set up random test data
dat_fourth_conv_random_test.fourier_sq.1 = bind_cols(dat_fourth_conv_random_test, Xone_AT_random_fourier_cos_test^2, Xone_AT_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_fourth_conv.fourier_sq.1.pred_random = predict(lm.model_fourth_conv_fourier_sq.1, newdata = dat_fourth_conv_random_test.fourier_sq.1)

# Prediction accuracy with random test data:
cor(lm.model_fourth_conv.fourier_sq.1.pred_random, dat_random_test$C0_new, method="pearson")
# 0.7047076
mean((lm.model_fourth_conv.fourier_sq.1.pred_random - dat_random_test$C0_new)^2)
# 0.07981403













########################################################################
# CONV_ONLY
########################################################################

dat_first_conv.conv_only = readRDS("data/Created/tiling_tiling_first_conv_output_train_conv_only.rds")
dat_first_conv_test.conv_only = readRDS("data/Created/tiling_tiling_first_conv_output_test_conv_only.rds")
dat_first_conv_random_test.conv_only = readRDS("data/Created/tiling_random_first_conv_output_test_conv_only.rds")

dat_first_conv.conv_only$Y = dat$C0_new
dat_first_conv_test.conv_only$Y = dat_test$C0_new

lm.model_first_conv.conv_only = lm(Y~., data = dat_first_conv.conv_only)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only$fitted.values, dat$C0_new, method="pearson")
# 0.3871514
mean((lm.model_first_conv.conv_only$fitted.values - dat$C0_new)^2)
# 0.2060222

# Predict values for test dataset:
lm.model_first_conv.conv_only.pred = predict(lm.model_first_conv.conv_only, newdata = dat_first_conv_test.conv_only)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only.pred, dat_test$C0_new, method="pearson")
# 0.3682102
mean((lm.model_first_conv.conv_only.pred - dat_test$C0_new)^2)
# 0.2085247

# Predict values for random test dataset:
lm.model_first_conv.conv_only.pred_random = predict(lm.model_first_conv.conv_only, newdata = dat_first_conv_random_test.conv_only)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only.pred_random, dat_random_test$C0_new, method="pearson")
# 0.3528178
mean((lm.model_first_conv.conv_only.pred_random - dat_random_test$C0_new)^2)
# 0.1381328










########################################################################
# CONV_ONLY + Fourier
########################################################################

# Set up models with fourier features included
dat_first_conv.conv_only_fourier_sq.1 = bind_cols(dat_first_conv.conv_only, Xone_AT_fourier_cos^2, Xone_AT_fourier_sin^2)

lm.model_first_conv.conv_only_fourier_sq.1 = lm(Y~., data = dat_first_conv.conv_only_fourier_sq.1)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_fourier_sq.1$fitted.values, dat$C0_new, method="pearson")
# 0.560766
mean((lm.model_first_conv.conv_only_fourier_sq.1$fitted.values - dat$C0_new)^2)
# 0.1661386

# Set up test data
dat_first_conv_test.conv_only_fourier_sq.1 = bind_cols(dat_first_conv_test.conv_only, Xone_AT_fourier_cos_test^2, Xone_AT_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_fourier_sq.1.pred = predict(lm.model_first_conv.conv_only_fourier_sq.1, newdata = dat_first_conv_test.conv_only_fourier_sq.1)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_fourier_sq.1.pred, dat_test$C0_new, method="pearson")
# 0.5499149
mean((lm.model_first_conv.conv_only_fourier_sq.1.pred - dat_test$C0_new)^2)
# 0.1681648

# Set up random test data
dat_first_conv_random_test.conv_only_fourier_sq.1 = bind_cols(dat_first_conv_random_test.conv_only, Xone_AT_random_fourier_cos_test^2, Xone_AT_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_fourier_sq.1.pred_random = predict(lm.model_first_conv.conv_only_fourier_sq.1, newdata = dat_first_conv_random_test.conv_only_fourier_sq.1)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_fourier_sq.1.pred_random, dat_random_test$C0_new, method="pearson")
# 0.5137014
mean((lm.model_first_conv.conv_only_fourier_sq.1.pred_random - dat_random_test$C0_new)^2)
# 0.1190525









########################################################################
# CONV_ONLY_2
########################################################################

dat_first_conv.conv_only_2 = readRDS("data/Created/tiling_tiling_first_conv_output_train_conv_only_2.rds")
dat_first_conv_test.conv_only_2 = readRDS("data/Created/tiling_tiling_first_conv_output_test_conv_only_2.rds")
dat_first_conv_random_test.conv_only_2 = readRDS("data/Created/tiling_random_first_conv_output_test_conv_only_2.rds")

dat_first_conv.conv_only_2$Y = dat$C0_new
dat_first_conv_test.conv_only_2$Y = dat_test$C0_new

lm.model_first_conv.conv_only_2 = lm(Y~., data = dat_first_conv.conv_only_2)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 0.3885772
mean((lm.model_first_conv.conv_only_2$fitted.values - dat$C0_new)^2)
# 0.2057541

# Predict values for test dataset:
lm.model_first_conv.conv_only_2.pred = predict(lm.model_first_conv.conv_only_2, newdata = dat_first_conv_test.conv_only_2)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 0.3783012
mean((lm.model_first_conv.conv_only_2.pred - dat_test$C0_new)^2)
# 0.2066198

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2.pred_random = predict(lm.model_first_conv.conv_only_2, newdata = dat_first_conv_random_test.conv_only_2)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2.pred_random, dat_random_test$C0_new, method="pearson")
# 0.3758269
mean((lm.model_first_conv.conv_only_2.pred_random - dat_random_test$C0_new)^2)
# 0.1341676


dat_first_conv.conv_only_2.lgb_train_set = dat_first_conv.conv_only_2 %>%
  select(-Y)
dat_first_conv.conv_only_2.lgb_train_rules = lgb.convert_with_rules(data = dat_first_conv.conv_only_2.lgb_train_set)
dat_first_conv.conv_only_2.lgb_train_data = lgb.Dataset(data = as.matrix(dat_first_conv.conv_only_2.lgb_train_rules$data),
                                                        label = dat_first_conv.conv_only_2$Y)
# dat_first_conv.conv_only_2.lgb_obj = lightgbm(data = dat_first_conv.conv_only_2.lgb_train_data,
#                                               params = list(learning_rate = c(0.1),
#                                                             objective = c("regression"),
#                                                             min_data_in_leaf = 500,
#                                                             max_depth = c(7),
#                                                             num_leaves = c(1000),
#                                                             lambda_l2 = c(1),
#                                                             boosting = c("gbdt")),
#                                               nrounds = 500)

dat_first_conv.conv_only_2.lgb_test_set = dat_first_conv_test.conv_only_2 %>%
  select(-Y)
dat_first_conv.conv_only_2.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_first_conv.conv_only_2.lgb_test_set, rules = dat_first_conv.conv_only_2.lgb_train_rules$rules)$data)
dat_first_conv.conv_only_2.lgb_test_data = lgb.Dataset(data = dat_first_conv.conv_only_2.lgb_test_matrix,
                                                       label = dat_test_first_conv.conv_only_2$C0_new,
                                                       categorical_feature = handle_categorical(colnames(dat_first_conv.conv_only_2.lgb_test_matrix)))

dat_first_conv.conv_only_2.lgb_obj = lightgbm(data = dat_first_conv.conv_only_2.lgb_train_data,
                                              params = list(learning_rate = c(0.1),
                                                            objective = c("regression"),
                                                            min_data_in_leaf = 100,
                                                            max_depth = c(5),
                                                            num_leaves = c(1000),
                                                            lambda_l2 = c(1),
                                                            boosting = c("gbdt")),
                                              valids = list(valid = dat_first_conv.conv_only_2.lgb_test_data),
                                              nrounds = 4000,
                                              early_stopping_rounds = 50)


dat_first_conv.conv_only_2.lgb_train_matrix = as.matrix(lgb.convert_with_rules(data = dat_first_conv.conv_only_2.lgb_train_set, rules = dat_first_conv.conv_only_2.lgb_train_rules$rules)$data)
dat_first_conv.conv_only_2.lgb_train_pred = predict(dat_first_conv.conv_only_2.lgb_obj, data = dat_first_conv.conv_only_2.lgb_train_matrix)

cor(dat_first_conv.conv_only_2.lgb_train_pred, dat$C0_new, method="pearson")
# 0.8171923 (nround = 4000)

dat_first_conv.conv_only_2.lgb_test_pred = predict(dat_first_conv.conv_only_2.lgb_obj, data = dat_first_conv.conv_only_2.lgb_test_matrix)

cor(dat_first_conv.conv_only_2.lgb_test_pred, dat_test$C0_new, method="pearson")
# 0.6947794 (nrounds = 4000)


poly2.model_fist_conv.conv_only_2.features = colnames(dat_first_conv.conv_only_2)[which(apply(dat_first_conv.conv_only_2, 2, function(col){
  length(unique(col))
}) > 2)]
poly2.model_fist_conv.conv_only_2.formula = as.formula(paste0("Y~",paste(paste0("poly(", poly2.model_fist_conv.conv_only_2.features[-length(poly2.model_fist_conv.conv_only_2.features)], ", 2)"), collapse = "+")))

poly2.model_first_conv.conv_only_2 = lm(poly2.model_fist_conv.conv_only_2.formula, data=dat_first_conv.conv_only_2)

# Prediction accuracy with train data:
cor(poly2.model_first_conv.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 
mean((poly2.model_first_conv.conv_only_2$fitted.values - dat$C0_new)^2)
# 

# Predict values for test dataset:
poly2.model_first_conv.conv_only_2.pred = predict(poly2.model_first_conv.conv_only_2, newdata = dat_first_conv_test.conv_only_2)

# Prediction accuracy with test data:
cor(poly2.model_first_conv.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 
mean((poly2.model_first_conv.conv_only_2.pred - dat_test$C0_new)^2)
# 

# Predict values for random test dataset:
poly2.model_first_conv.conv_only_2.pred_random = predict(poly2.model_first_conv.conv_only_2, newdata = dat_first_conv_random_test.conv_only_2)

# Prediction accuracy with random test data:
cor(poly2.model_first_conv.conv_only_2.pred_random, dat_random_test$C0_new, method="pearson")
# 
mean((poly2.model_first_conv.conv_only_2.pred_random - dat_random_test$C0_new)^2)
# 







#FIXME
gam.model_fist_conv.conv_only_2.features = colnames(dat_first_conv.conv_only_2)[which(apply(dat_first_conv.conv_only_2, 2, function(col){
  length(unique(col))
}) > 2)]
gam.model_fist_conv.conv_only_2.formula = as.formula(paste0("Y~",paste(paste0("s(", gam.model_fist_conv.conv_only_2.features, ")"), collapse = "+")))

gam.model_first_conv.conv_only_2 = bam(gam.model_fist_conv.conv_only_2.formula, data=dat_first_conv.conv_only_2)

# Prediction accuracy with train data:
cor(gam.model_first_conv.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 
mean((gam.model_first_conv.conv_only_2$fitted.values - dat$C0_new)^2)
# 

# Predict values for test dataset:
gam.model_first_conv.conv_only_2.pred = predict(gam.model_first_conv.conv_only_2, newdata = dat_first_conv_test.conv_only_2)

# Prediction accuracy with test data:
cor(gam.model_first_conv.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 
mean((gam.model_first_conv.conv_only_2.pred - dat_test$C0_new)^2)
# 

# Predict values for random test dataset:
gam.model_first_conv.conv_only_2.pred_random = predict(gam.model_first_conv.conv_only_2, newdata = dat_first_conv_random_test.conv_only_2)

# Prediction accuracy with random test data:
cor(gam.model_first_conv.conv_only_2.pred_random, dat_random_test$C0_new, method="pearson")
# 
mean((gam.model_first_conv.conv_only_2.pred_random - dat_random_test$C0_new)^2)
# 







dat_max_pool.conv_only_2 = readRDS("data/Created/tiling_tiling_max_pool_output_train_conv_only_2.rds")
dat_max_pool_test.conv_only_2 = readRDS("data/Created/tiling_tiling_max_pool_output_test_conv_only_2.rds")

dat_max_pool.conv_only_2$Y = dat$C0_new
dat_max_pool_test.conv_only_2$Y = dat_test$C0_new

lm.model_max_pool.conv_only_2 = lm(Y~., data = dat_max_pool.conv_only_2)

# Prediction accuracy with train data:
cor(lm.model_max_pool.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 0.3556862
mean((lm.model_max_pool.conv_only_2$fitted.values - dat$C0_new)^2)
# 0.2116867

# Predict values for test dataset:
lm.model_max_pool.conv_only_2.pred = predict(lm.model_max_pool.conv_only_2, newdata = dat_max_pool_test.conv_only_2)

# Prediction accuracy with test data:
cor(lm.model_max_pool.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 0.3534044
mean((lm.model_max_pool.conv_only_2.pred - dat_test$C0_new)^2)
# 0.2109547



dat_batch_norm.conv_only_2 = readRDS("data/Created/tiling_tiling_batch_norm_output_train_conv_only_2.rds")
dat_batch_norm_test.conv_only_2 = readRDS("data/Created/tiling_tiling_batch_norm_output_test_conv_only_2.rds")

dat_batch_norm.conv_only_2$Y = dat$C0_new
dat_batch_norm_test.conv_only_2$Y = dat_test$C0_new

lm.model_batch_norm.conv_only_2 = lm(Y~., data = dat_batch_norm.conv_only_2)

# Prediction accuracy with train data:
cor(lm.model_batch_norm.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 0.3556862
mean((lm.model_batch_norm.conv_only_2$fitted.values - dat$C0_new)^2)
# 0.2116867

# Predict values for test dataset:
lm.model_batch_norm.conv_only_2.pred = predict(lm.model_batch_norm.conv_only_2, newdata = dat_batch_norm_test.conv_only_2)

# Prediction accuracy with test data:
cor(lm.model_batch_norm.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 0.3534044
mean((lm.model_batch_norm.conv_only_2.pred - dat_test$C0_new)^2)
# 0.2109547



dat_dense.conv_only_2 = readRDS("data/Created/tiling_tiling_dense_output_train_conv_only_2.rds")
dat_dense_test.conv_only_2 = readRDS("data/Created/tiling_tiling_dense_output_test_conv_only_2.rds")

dat_dense.conv_only_2 = data.frame(dat_dense.conv_only_2)
dat_dense_test.conv_only_2 = data.frame(dat_dense_test.conv_only_2)

dat_dense.conv_only_2$Y = dat$C0_new
dat_dense_test.conv_only_2$Y = dat_test$C0_new

lm.model_dense.conv_only_2 = lm(Y~., data = dat_dense.conv_only_2)

# Prediction accuracy with train data:
cor(lm.model_dense.conv_only_2$fitted.values, dat$C0_new, method="pearson")
# 0.01231283
mean((lm.model_dense.conv_only_2$fitted.values - dat$C0_new)^2)
# 0.2423099

# Predict values for test dataset:
lm.model_dense.conv_only_2.pred = predict(lm.model_dense.conv_only_2, newdata = dat_dense_test.conv_only_2)

# Prediction accuracy with test data:
cor(lm.model_dense.conv_only_2.pred, dat_test$C0_new, method="pearson")
# 
mean((lm.model_dense.conv_only_2.pred - dat_test$C0_new)^2)
# 











########################################################################
# CONV_ONLY_2 + Fourier
########################################################################

# Set up models with fourier features included
dat_first_conv.conv_only_2_fourier_sq.1 = bind_cols(dat_first_conv.conv_only_2, Xone_AT_fourier_cos^2, Xone_AT_fourier_sin^2)

lm.model_first_conv.conv_only_2_fourier_sq.1 = lm(Y~., data = dat_first_conv.conv_only_2_fourier_sq.1)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.1$fitted.values, dat$C0_new, method="pearson")
# 0.5589223
mean((lm.model_first_conv.conv_only_2_fourier_sq.1$fitted.values - dat$C0_new)^2)
# 0.166639

# Set up test data
dat_first_conv_test.conv_only_2_fourier_sq.1 = bind_cols(dat_first_conv_test.conv_only_2, Xone_AT_fourier_cos_test^2, Xone_AT_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.1.pred = predict(lm.model_first_conv.conv_only_2_fourier_sq.1, newdata = dat_first_conv_test.conv_only_2_fourier_sq.1)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.1.pred, dat_test$C0_new, method="pearson")
# 0.5500898
mean((lm.model_first_conv.conv_only_2_fourier_sq.1.pred - dat_test$C0_new)^2)
# 0.1681023

# Set up random test data
dat_first_conv_random_test.conv_only_2_fourier_sq.1 = bind_cols(dat_first_conv_random_test.conv_only_2, Xone_AT_random_fourier_cos_test^2, Xone_AT_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.1.pred_random = predict(lm.model_first_conv.conv_only_2_fourier_sq.1, newdata = dat_first_conv_random_test.conv_only_2_fourier_sq.1)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.1.pred_random, dat_random_test$C0_new, method="pearson")
# 0.5148304
mean((lm.model_first_conv.conv_only_2_fourier_sq.1.pred_random - dat_random_test$C0_new)^2)
# 0.119102

dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_set = dat_first_conv.conv_only_2_fourier_sq.1 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_rules = lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_set)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_data = lgb.Dataset(data = as.matrix(dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_rules$data),
                                                                     label = dat_first_conv.conv_only_2_fourier_sq.1$Y)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_obj = lightgbm(data = dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_data,
                                                           params = list(learning_rate = c(0.1),
                                                                         objective = c("regression"),
                                                                         min_data_in_leaf = 100,
                                                                         max_depth = c(7),
                                                                         num_leaves = c(1000),
                                                                         lambda_l2 = c(1),
                                                                         boosting = c("gbdt")),
                                                           nrounds = 2000)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_set = dat_first_conv_test.conv_only_2_fourier_sq.1 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_data = as.matrix(lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_set, rules = dat_first_conv.conv_only_2_fourier_sq.1.lgb_train_rules$rules)$data)
dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_pred = predict(dat_first_conv.conv_only_2_fourier_sq.1.lgb_obj, data = dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_data)

cor(dat_first_conv.conv_only_2_fourier_sq.1.lgb_test_pred, dat_test$C0_new, method="pearson")
# 0.5453537 (nrounds = 50)
# 0.6178642 (nrounds = 200)
# 0.6555547 (nrounds = 500)
# 0.6749467 (nround = 1000)
# 0.6850614 (nround = 2000)



# Set up models with fourier features included
dat_first_conv.conv_only_2_fourier_sq.2 = bind_cols(dat_first_conv.conv_only_2, Xtwo_AT_fourier_cos^2, Xtwo_AT_fourier_sin^2)

lm.model_first_conv.conv_only_2_fourier_sq.2 = lm(Y~., data = dat_first_conv.conv_only_2_fourier_sq.2)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.2$fitted.values, dat$C0_new, method="pearson")
# 0.5334231
mean((lm.model_first_conv.conv_only_2_fourier_sq.2$fitted.values - dat$C0_new)^2)
# 0.1733892

# Set up test data
dat_first_conv_test.conv_only_2_fourier_sq.2 = bind_cols(dat_first_conv_test.conv_only_2, Xtwo_AT_fourier_cos_test^2, Xtwo_AT_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.2.pred = predict(lm.model_first_conv.conv_only_2_fourier_sq.2, newdata = dat_first_conv_test.conv_only_2_fourier_sq.2)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.2.pred, dat_test$C0_new, method="pearson")
# 0.5317166
mean((lm.model_first_conv.conv_only_2_fourier_sq.2.pred - dat_test$C0_new)^2)
# 0.172868

# Set up random test data
dat_first_conv_random_test.conv_only_2_fourier_sq.2 = bind_cols(dat_first_conv_random_test.conv_only_2, Xtwo_AT_random_fourier_cos_test^2, Xtwo_AT_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.2.pred_random = predict(lm.model_first_conv.conv_only_2_fourier_sq.2, newdata = dat_first_conv_random_test.conv_only_2_fourier_sq.2)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.2.pred_random, dat_random_test$C0_new, method="pearson")
# 0.4779607
mean((lm.model_first_conv.conv_only_2_fourier_sq.2.pred_random - dat_random_test$C0_new)^2)
# 0.1204049

dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_set = dat_first_conv.conv_only_2_fourier_sq.2 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_rules = lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_set)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_data = lgb.Dataset(data = as.matrix(dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_rules$data),
                                                                     label = dat_first_conv.conv_only_2_fourier_sq.2$Y)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_obj = lightgbm(data = dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_data,
                                                           params = list(learning_rate = c(0.1),
                                                                         objective = c("regression"),
                                                                         min_data_in_leaf = 100,
                                                                         max_depth = c(7),
                                                                         num_leaves = c(1000),
                                                                         lambda_l2 = c(1),
                                                                         boosting = c("gbdt")),
                                                           nrounds = 2000)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_set = dat_first_conv_test.conv_only_2_fourier_sq.2 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_data = as.matrix(lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_set, rules = dat_first_conv.conv_only_2_fourier_sq.2.lgb_train_rules$rules)$data)
dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_pred = predict(dat_first_conv.conv_only_2_fourier_sq.2.lgb_obj, data = dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_data)

cor(dat_first_conv.conv_only_2_fourier_sq.2.lgb_test_pred, dat_test$C0_new, method="pearson")
# 0.6849322 (nround 2000)


# Set up models with fourier features included
dat_first_conv.conv_only_2_fourier_sq.3 = bind_cols(dat_first_conv.conv_only_2, Xone_CG_fourier_cos^2, Xone_CG_fourier_sin^2)

lm.model_first_conv.conv_only_2_fourier_sq.3 = lm(Y~., data = dat_first_conv.conv_only_2_fourier_sq.3)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.3$fitted.values, dat$C0_new, method="pearson")
# 0.5589203
mean((lm.model_first_conv.conv_only_2_fourier_sq.3$fitted.values - dat$C0_new)^2)
# 0.1666395

# Set up test data
dat_first_conv_test.conv_only_2_fourier_sq.3 = bind_cols(dat_first_conv_test.conv_only_2, Xone_CG_fourier_cos_test^2, Xone_CG_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.3.pred = predict(lm.model_first_conv.conv_only_2_fourier_sq.3, newdata = dat_first_conv_test.conv_only_2_fourier_sq.3)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.3.pred, dat_test$C0_new, method="pearson")
# 0.5501346
mean((lm.model_first_conv.conv_only_2_fourier_sq.3.pred - dat_test$C0_new)^2)
# 0.1680903

# Set up random test data
dat_first_conv_random_test.conv_only_2_fourier_sq.3 = bind_cols(dat_first_conv_random_test.conv_only_2, Xone_CG_random_fourier_cos_test^2, Xone_CG_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.3.pred_random = predict(lm.model_first_conv.conv_only_2_fourier_sq.3, newdata = dat_first_conv_random_test.conv_only_2_fourier_sq.3)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.3.pred_random, dat_random_test$C0_new, method="pearson")
# 0.5147154
mean((lm.model_first_conv.conv_only_2_fourier_sq.3.pred_random - dat_random_test$C0_new)^2)
# 0.1191346



# Set up models with fourier features included
dat_first_conv.conv_only_2_fourier_sq.4 = bind_cols(dat_first_conv.conv_only_2, Xtwo_CG_fourier_cos^2, Xtwo_CG_fourier_sin^2)

lm.model_first_conv.conv_only_2_fourier_sq.4 = lm(Y~., data = dat_first_conv.conv_only_2_fourier_sq.4)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.4$fitted.values, dat$C0_new, method="pearson")
# 0.4207978
mean((lm.model_first_conv.conv_only_2_fourier_sq.4$fitted.values - dat$C0_new)^2)
# 0.1994341

# Set up test data
dat_first_conv_test.conv_only_2_fourier_sq.4 = bind_cols(dat_first_conv_test.conv_only_2, Xtwo_CG_fourier_cos_test^2, Xtwo_CG_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.4.pred = predict(lm.model_first_conv.conv_only_2_fourier_sq.4, newdata = dat_first_conv_test.conv_only_2_fourier_sq.4)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.4.pred, dat_test$C0_new, method="pearson")
# 0.4061322
mean((lm.model_first_conv.conv_only_2_fourier_sq.4.pred - dat_test$C0_new)^2)
# 0.2013736

# Set up random test data
dat_first_conv_random_test.conv_only_2_fourier_sq.4 = bind_cols(dat_first_conv_random_test.conv_only_2, Xtwo_CG_random_fourier_cos_test^2, Xtwo_CG_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.4.pred_random = predict(lm.model_first_conv.conv_only_2_fourier_sq.4, newdata = dat_first_conv_random_test.conv_only_2_fourier_sq.4)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.4.pred_random, dat_random_test$C0_new, method="pearson")
# 0.4210632
mean((lm.model_first_conv.conv_only_2_fourier_sq.4.pred_random - dat_random_test$C0_new)^2)
# 0.1297992



# Set up models with fourier features included
dat_first_conv.conv_only_2_fourier_sq.5 = bind_cols(dat_first_conv.conv_only_2, Xtwo_AT_fourier_cos^2, Xtwo_AT_fourier_sin^2, 
                                                    Xtwo_CG_fourier_cos^2, Xtwo_CG_fourier_sin^2)

lm.model_first_conv.conv_only_2_fourier_sq.5 = lm(Y~., data = dat_first_conv.conv_only_2_fourier_sq.5)

# Prediction accuracy with train data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.5$fitted.values, dat$C0_new, method="pearson")
# 0.5382848
mean((lm.model_first_conv.conv_only_2_fourier_sq.5$fitted.values - dat$C0_new)^2)
# 0.1721265

# Set up test data
dat_first_conv_test.conv_only_2_fourier_sq.5 = bind_cols(dat_first_conv_test.conv_only_2, Xtwo_AT_fourier_cos_test^2, Xtwo_AT_fourier_sin_test^2, 
                                                         Xtwo_CG_fourier_cos_test^2, Xtwo_CG_fourier_sin_test^2)

# Predict values for test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.5.pred = predict(lm.model_first_conv.conv_only_2_fourier_sq.5, newdata = dat_first_conv_test.conv_only_2_fourier_sq.5)

# Prediction accuracy with test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.5.pred, dat_test$C0_new, method="pearson")
# 0.5348621
mean((lm.model_first_conv.conv_only_2_fourier_sq.5.pred - dat_test$C0_new)^2)
# 0.1720627

# Set up random test data
dat_first_conv_random_test.conv_only_2_fourier_sq.5 = bind_cols(dat_first_conv_random_test.conv_only_2, Xtwo_AT_random_fourier_cos_test^2, Xtwo_AT_random_fourier_sin_test^2,
                                                                Xtwo_CG_random_fourier_cos_test^2, Xtwo_CG_random_fourier_sin_test^2)

# Predict values for random test dataset:
lm.model_first_conv.conv_only_2_fourier_sq.5.pred_random = predict(lm.model_first_conv.conv_only_2_fourier_sq.5, newdata = dat_first_conv_random_test.conv_only_2_fourier_sq.5)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.conv_only_2_fourier_sq.5.pred_random, dat_random_test$C0_new, method="pearson")
# 0.4956422
mean((lm.model_first_conv.conv_only_2_fourier_sq.5.pred_random - dat_random_test$C0_new)^2)
# 0.1180861

dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_set = dat_first_conv.conv_only_2_fourier_sq.5 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_rules = lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_set)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_data = lgb.Dataset(data = as.matrix(dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_rules$data),
                                                                     label = dat_first_conv.conv_only_2_fourier_sq.5$Y)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_obj = lightgbm(data = dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_data,
                                                           params = list(learning_rate = c(0.1),
                                                                         objective = c("regression"),
                                                                         min_data_in_leaf = 100,
                                                                         max_depth = c(7),
                                                                         num_leaves = c(1000),
                                                                         lambda_l2 = c(1),
                                                                         boosting = c("gbdt")),
                                                           nrounds = 2000)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_set = dat_first_conv_test.conv_only_2_fourier_sq.5 %>%
  select(-Y)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_data = as.matrix(lgb.convert_with_rules(data = dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_set, rules = dat_first_conv.conv_only_2_fourier_sq.5.lgb_train_rules$rules)$data)
dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_pred = predict(dat_first_conv.conv_only_2_fourier_sq.5.lgb_obj, data = dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_data)

cor(dat_first_conv.conv_only_2_fourier_sq.5.lgb_test_pred, dat_test$C0_new, method="pearson")
# 0.6735351 (nround = 2000)










########################################################################
# 2CONV_ONLY
########################################################################

dat_first_conv.2conv_only = readRDS("data/Created/tiling_tiling_first_conv_output_train_2conv_only.rds")
dat_first_conv_test.2conv_only = readRDS("data/Created/tiling_tiling_first_conv_output_test_2conv_only.rds")
dat_first_conv_random_test.2conv_only = readRDS("data/Created/tiling_random_first_conv_output_test_2conv_only.rds")

dat_first_conv.2conv_only$Y = dat$C0_new
dat_first_conv_test.2conv_only$Y = dat_test$C0_new

lm.model_first_conv.2conv_only = lm(Y~., data = dat_first_conv.2conv_only)

# gam.model_first_conv.2conv_only = gam(as.formula(paste0("Y ~ s(", setdiff(colnames(dat_first_conv.2conv_only), "Y") %>%
#                                                           paste0(collapse = ") + s("), ")")), data=dat_first_conv.2conv_only)


# Prediction accuracy with train data:
cor(lm.model_first_conv.2conv_only$fitted.values, dat$C0_new, method="pearson")
# 0.4290705
mean((lm.model_first_conv.2conv_only$fitted.values - dat$C0_new)^2)
# 0.1977302

# Predict values for test dataset:
lm.model_first_conv.2conv_only.pred = predict(lm.model_first_conv.2conv_only, newdata = dat_first_conv_test.2conv_only)

# Prediction accuracy with test data:
cor(lm.model_first_conv.2conv_only.pred, dat_test$C0_new, method="pearson")
# 0.3758825
mean((lm.model_first_conv.2conv_only.pred - dat_test$C0_new)^2)
# 0.2079272

# Predict values for random test dataset:
lm.model_first_conv.2conv_only.pred_random = predict(lm.model_first_conv.2conv_only, newdata = dat_first_conv_random_test.2conv_only)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.2conv_only.pred_random, dat_random_test$C0_new, method="pearson")
# 0.3464437
mean((lm.model_first_conv.2conv_only.pred_random - dat_random_test$C0_new)^2)
# 0.1415889





dat_second_conv.2conv_only = readRDS("data/Created/tiling_tiling_second_conv_output_train_2conv_only.rds")
dat_second_conv_test.2conv_only = readRDS("data/Created/tiling_tiling_second_conv_output_test_2conv_only.rds")
dat_second_conv_random_test.2conv_only = readRDS("data/Created/tiling_random_second_conv_output_test_2conv_only.rds")

dat_second_conv.2conv_only$Y = dat$C0_new
dat_second_conv_test.2conv_only$Y = dat_test$C0_new

lm.model_second_conv.2conv_only = lm(Y~., data = dat_second_conv.2conv_only)

# Prediction accuracy with train data:
cor(lm.model_second_conv.2conv_only$fitted.values, dat$C0_new, method="pearson")
# 0.77974
mean((lm.model_second_conv.2conv_only$fitted.values - dat$C0_new)^2)
# 0.09500122

# Predict values for test dataset:
lm.model_second_conv.2conv_only.pred = predict(lm.model_second_conv.2conv_only, newdata = dat_second_conv_test.2conv_only)

# Prediction accuracy with test data:
cor(lm.model_second_conv.2conv_only.pred, dat_test$C0_new, method="pearson")
# 0.7795399
mean((lm.model_second_conv.2conv_only.pred - dat_test$C0_new)^2)
# 0.09455835

# Predict values for random test dataset:
lm.model_second_conv.2conv_only.pred_random = predict(lm.model_second_conv.2conv_only, newdata = dat_second_conv_random_test.2conv_only)

# Prediction accuracy with random test data:
cor(lm.model_second_conv.2conv_only.pred_random, dat_random_test$C0_new, method="pearson")
# 0.7847345
mean((lm.model_second_conv.2conv_only.pred_random - dat_random_test$C0_new)^2)
# 0.0590556









########################################################################
# 1CONV_LSTM
########################################################################

dat_first_conv.1conv_lstm = readRDS("data/Created/tiling_tiling_first_conv_output_train_1conv_lstm.rds")
dat_first_conv_test.1conv_lstm = readRDS("data/Created/tiling_tiling_first_conv_output_test_1conv_lstm.rds")
dat_first_conv_random_test.1conv_lstm = readRDS("data/Created/tiling_random_first_conv_output_test_1conv_lstm.rds")

dat_first_conv.1conv_lstm$Y = dat$C0_new
dat_first_conv_test.1conv_lstm$Y = dat_test$C0_new

lm.model_first_conv.1conv_lstm = lm(Y~., data = dat_first_conv.1conv_lstm)

# Prediction accuracy with train data:
cor(lm.model_first_conv.1conv_lstm$fitted.values, dat$C0_new, method="pearson")
# 0.3527634
mean((lm.model_first_conv.1conv_lstm$fitted.values - dat$C0_new)^2)
# 0.2121885

# Predict values for test dataset:
lm.model_first_conv.1conv_lstm.pred = predict(lm.model_first_conv.1conv_lstm, newdata = dat_first_conv_test.1conv_lstm)

# Prediction accuracy with test data:
cor(lm.model_first_conv.1conv_lstm.pred, dat_test$C0_new, method="pearson")
# 0.3384393
mean((lm.model_first_conv.1conv_lstm.pred - dat_test$C0_new)^2)
# 0.2135239

# Predict values for random test dataset:
lm.model_first_conv.1conv_lstm.pred_random = predict(lm.model_first_conv.1conv_lstm, newdata = dat_first_conv_random_test.1conv_lstm)

# Prediction accuracy with random test data:
cor(lm.model_first_conv.1conv_lstm.pred_random, dat_random_test$C0_new, method="pearson")
# 0.3020655
mean((lm.model_first_conv.1conv_lstm.pred_random - dat_random_test$C0_new)^2)
# 0.1426499






dat_lstm.1conv_lstm = readRDS("data/Created/tiling_tiling_lstm_output_train_1conv_lstm.rds")
dat_lstm_test.1conv_lstm = readRDS("data/Created/tiling_tiling_lstm_output_test_1conv_lstm.rds")
dat_lstm_random_test.1conv_lstm = readRDS("data/Created/tiling_random_lstm_output_test_1conv_lstm.rds")

dat_lstm.1conv_lstm$Y = dat$C0_new
dat_lstm_test.1conv_lstm$Y = dat_test$C0_new

lm.model_lstm.1conv_lstm = lm(Y~., data = dat_lstm.1conv_lstm)

# Prediction accuracy with train data:
cor(lm.model_lstm.1conv_lstm$fitted.values, dat$C0_new, method="pearson")
# 0.8206614
mean((lm.model_lstm.1conv_lstm$fitted.values - dat$C0_new)^2)
# 0.07912978

# Predict values for test dataset:
lm.model_lstm.1conv_lstm.pred = predict(lm.model_lstm.1conv_lstm, newdata = dat_lstm_test.1conv_lstm)

# Prediction accuracy with test data:
cor(lm.model_lstm.1conv_lstm.pred, dat_test$C0_new, method="pearson")
# 0.8268062
mean((lm.model_lstm.1conv_lstm.pred - dat_test$C0_new)^2)
# 0.07625315

# Predict values for random test dataset:
lm.model_lstm.1conv_lstm.pred_random = predict(lm.model_lstm.1conv_lstm, newdata = dat_lstm_random_test.1conv_lstm)

# Prediction accuracy with random test data:
cor(lm.model_lstm.1conv_lstm.pred_random, dat_random_test$C0_new, method="pearson")
# 0.8337224
mean((lm.model_lstm.1conv_lstm.pred_random - dat_random_test$C0_new)^2)
# 0.04879125

