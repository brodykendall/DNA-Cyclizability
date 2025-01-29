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
dat_tiling_all = rbind(dat, dat_test)
y_tiling_all = c(y, y_test)
dat_random = readRDS("data/Created/processed_random_newC0.rds")
y_random = dat_random$C0_new
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")
y_random_test = dat_random_test$C0_new
dat_random_all = rbind(dat_random, dat_random_test)
y_random_all = c(y_random, y_random_test)
dat_chrV = readRDS("data/Created/processed_chrV_newC0.rds")
y_chrV = dat_chrV$C0_new
dat_chrV_test = readRDS("data/Created/processed_chrV_test_newC0.rds")
y_chrV_test = dat_chrV_test$C0_new
dat_chrV_all = rbind(dat_chrV, dat_chrV_test)
y_chrV_all = c(y_chrV, y_chrV_test)
dat_yeast = readRDS("data/Created/processed_yeast_newC0.rds")
y_yeast = dat_yeast$C0_new
dat_yeast_test = readRDS("data/Created/processed_yeast_test_newC0.rds")
y_yeast_test = dat_yeast_test$C0_new
dat_yeast_all = rbind(dat_yeast, dat_yeast_test)
y_yeast_all = c(y_yeast, y_yeast_test)

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

# distances = c(4,5,6, 9,10,11, 14,15,16, 19,20,21, 24,25,26, 29,30,31, 34,35,36, 39,40,41, 44,45,46)
distances = 2:48

# Tiling Library:
dist_freq_di_AorTdi_bid2 = map2(rep(dinucleotides, length(distances)), 
                                rep(distances, each=16), 
                                ~stri_count(dat$x50mer, 
                                            regex=paste0("(?=(",.x,".{",.y-2,"}[AT]{2})|([AT]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_AorTdi_bid2) = paste0(rep(dinucleotides, length(distances)), 
                                            "_AorT_dist", rep(distances, each=16))
saveRDS(dist_freq_di_AorTdi_bid2, "data/Created/tiling_dist_freq_di_AorTdi_bid2.rds")

dist_freq_di_CorGdi_bid2 = map2(rep(dinucleotides, length(distances)), 
                                rep(distances, each=16), 
                                ~stri_count(dat$x50mer, 
                                            regex=paste0("(?=(",.x,".{",.y-2,"}[CG]{2})|([CG]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_CorGdi_bid2) = paste0(rep(dinucleotides, length(distances)), 
                                            "_CorG_dist", rep(distances, each=16))
saveRDS(dist_freq_di_CorGdi_bid2, "data/Created/tiling_dist_freq_di_CorGdi_bid2.rds")


# Test data:
dist_freq_di_AorTdi_bid2_test = map2(rep(dinucleotides, length(distances)), 
                                     rep(distances, each=16), 
                                     ~stri_count(dat_test$x50mer, 
                                                 regex=paste0("(?=(",.x,".{",.y-2,"}[AT]{2})|([AT]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_AorTdi_bid2_test) = paste0(rep(dinucleotides, length(distances)), 
                                                 "_AorT_dist", rep(distances, each=16))
saveRDS(dist_freq_di_AorTdi_bid2_test, "data/Created/tiling_dist_freq_di_AorTdi_bid2_test.rds")

dist_freq_di_CorGdi_bid2_test = map2(rep(dinucleotides, length(distances)), 
                                     rep(distances, each=16), 
                                     ~stri_count(dat_test$x50mer, 
                                                 regex=paste0("(?=(",.x,".{",.y-2,"}[CG]{2})|([CG]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_CorGdi_bid2_test) = paste0(rep(dinucleotides, length(distances)), 
                                                 "_CorG_dist", rep(distances, each=16))
saveRDS(dist_freq_di_CorGdi_bid2_test, "data/Created/tiling_dist_freq_di_CorGdi_bid2_test.rds")

dist_freq_di_AorTdi_bid2_tiling_all = rbind(dist_freq_di_AorTdi_bid2, 
                                            dist_freq_di_AorTdi_bid2_test)
dist_freq_di_CorGdi_bid2_tiling_all = rbind(dist_freq_di_CorGdi_bid2, 
                                            dist_freq_di_CorGdi_bid2_test)

# Random Library:
dist_freq_di_AorTdi_bid2_random_all = map2(rep(dinucleotides, length(distances)), 
                                           rep(distances, each=16), 
                                           ~stri_count(dat_random_all$x50mer, 
                                                       regex=paste0("(?=(",.x,".{",.y-2,"}[AT]{2})|([AT]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_AorTdi_bid2_random_all) = paste0(rep(dinucleotides, length(distances)), 
                                                       "_AorT_dist", rep(distances, each=16))
saveRDS(dist_freq_di_AorTdi_bid2_random_all, "data/Created/random_all_dist_freq_di_AorTdi_bid2.rds")


dist_freq_di_CorGdi_bid2_random_all = map2(rep(dinucleotides, length(distances)), 
                                           rep(distances, each=16), 
                                           ~stri_count(dat_random_all$x50mer, 
                                                       regex=paste0("(?=(",.x,".{",.y-2,"}[CG]{2})|([CG]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_CorGdi_bid2_random_all) = paste0(rep(dinucleotides, length(distances)), 
                                                       "_CorG_dist", rep(distances, each=16))
saveRDS(dist_freq_di_CorGdi_bid2_random_all, "data/Created/random_all_dist_freq_di_CorGdi_bid2.rds")



# ChrV Library:
dist_freq_di_AorTdi_bid2_chrV_all = map2(rep(dinucleotides, length(distances)), 
                                         rep(distances, each=16), 
                                         ~stri_count(dat_chrV_all$x50mer, 
                                                     regex=paste0("(?=(",.x,".{",.y-2,"}[AT]{2})|([AT]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_AorTdi_bid2_chrV_all) = paste0(rep(dinucleotides, length(distances)), 
                                                     "_AorT_dist", rep(distances, each=16))
saveRDS(dist_freq_di_AorTdi_bid2_chrV_all, "data/Created/chrV_all_dist_freq_di_AorTdi_bid2.rds")


dist_freq_di_CorGdi_bid2_chrV_all = map2(rep(dinucleotides, length(distances)), 
                                         rep(distances, each=16), 
                                         ~stri_count(dat_chrV_all$x50mer, 
                                                     regex=paste0("(?=(",.x,".{",.y-2,"}[CG]{2})|([CG]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_CorGdi_bid2_chrV_all) = paste0(rep(dinucleotides, length(distances)), 
                                                     "_CorG_dist", rep(distances, each=16))
saveRDS(dist_freq_di_CorGdi_bid2_chrV_all, "data/Created/chrV_all_dist_freq_di_CorGdi_bid2.rds")



# Yeast Library:
dist_freq_di_AorTdi_bid2_yeast_all = map2(rep(dinucleotides, length(distances)), 
                                          rep(distances, each=16), 
                                          ~stri_count(dat_yeast_all$x50mer, 
                                                      regex=paste0("(?=(",.x,".{",.y-2,"}[AT]{2})|([AT]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_AorTdi_bid2_yeast_all) = paste0(rep(dinucleotides, length(distances)), 
                                                      "_AorT_dist", rep(distances, each=16))
saveRDS(dist_freq_di_AorTdi_bid2_yeast_all, "data/Created/yeast_all_dist_freq_di_AorTdi_bid2.rds")


dist_freq_di_CorGdi_bid2_yeast_all = map2(rep(dinucleotides, length(distances)), 
                                          rep(distances, each=16), 
                                          ~stri_count(dat_yeast_all$x50mer, 
                                                      regex=paste0("(?=(",.x,".{",.y-2,"}[CG]{2})|([CG]{2}.{",.y-2,"}",.x,"))"))) %>%
  bind_cols()
colnames(dist_freq_di_CorGdi_bid2_yeast_all) = paste0(rep(dinucleotides, length(distances)), 
                                                      "_CorG_dist", rep(distances, each=16))
saveRDS(dist_freq_di_CorGdi_bid2_yeast_all, "data/Created/yeast_all_dist_freq_di_CorGdi_bid2.rds")








# Simple Linear Model Training on Tiling train set:
dat_temp_dist = cbind(dist_freq_di_AorTdi_bid2, dist_freq_di_CorGdi_bid2, dat%>%select(C0_new))
temp_dist_lm = lm(C0_new~., data=dat_temp_dist)

cor(temp_dist_lm$fitted.values, y)
# distances = 5 10 15 20 25 30 35 40 45
# 0.5847679
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.6600189
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.6892242
# distances = 2-48
# 0.7324237

# Prediction on tiling test set:
dat_test_temp_dist = cbind(dist_freq_di_AorTdi_bid2_test, dist_freq_di_CorGdi_bid2_test,
                           dat_test%>%select(C0_new))
temp_dist_pred = predict(temp_dist_lm, dat_test_temp_dist)

cor(temp_dist_pred, y_test)
# distances = 5 10 15 20 25 30 35 40 45
# 0.5821323
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.6600768
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.6890315
# distances = 2-48
# 0.7304727




# Prediction on entire random library:
dat_random_all_temp_dist = cbind(dist_freq_di_AorTdi_bid2_random_all, 
                                 dist_freq_di_CorGdi_bid2_random_all,
                                 dat_random_all%>%select(C0_new))
random_all_temp_dist_pred = predict(temp_dist_lm, dat_random_all_temp_dist)

cor(random_all_temp_dist_pred, y_random_all)
# distances = 2-48, modeled on train set
# 0.6490919

# Prediction on entire chrV library:
dat_chrV_all_temp_dist = cbind(dist_freq_di_AorTdi_bid2_chrV_all, 
                               dist_freq_di_CorGdi_bid2_chrV_all,
                               dat_chrV_all%>%select(C0_new))
chrV_all_temp_dist_pred = predict(temp_dist_lm, dat_chrV_all_temp_dist)

cor(chrV_all_temp_dist_pred, y_chrV_all)
# distances = 2-48
# 0.5893484

# Prediction on entire yeast library:
dat_yeast_all_temp_dist = cbind(dist_freq_di_AorTdi_bid2_yeast_all, 
                                dist_freq_di_CorGdi_bid2_yeast_all,
                                dat_yeast_all%>%select(C0_new))
yeast_all_temp_dist_pred = predict(temp_dist_lm, dat_yeast_all_temp_dist)

cor(yeast_all_temp_dist_pred, y_yeast_all)
# distances = 2-48
# 0.665529





# Simple Linear Model Training on entire Tiling library:
dat_tiling_all_temp_dist = cbind(dist_freq_di_AorTdi_bid2_tiling_all, 
                                 dist_freq_di_CorGdi_bid2_tiling_all, 
                                 dat_tiling_all%>%select(C0_new))
# Correlation on Training data:
temp_dist_tiling_all_lm = lm(C0_new~., data=dat_tiling_all_temp_dist)

cor(temp_dist_tiling_all_lm$fitted.values, y_tiling_all)
# distances = 2-48
# 0.7328266

# Correlation on Test data (other libraries):
random_all_temp_dist_pred_tiling_all = predict(temp_dist_tiling_all_lm, dat_random_all_temp_dist)

cor(random_all_temp_dist_pred_tiling_all, y_random_all)
# distances = 2-48, modeled on entire tiling library
# 0.6502421

chrV_all_temp_dist_pred_tiling_all = predict(temp_dist_tiling_all_lm, dat_chrV_all_temp_dist)

cor(chrV_all_temp_dist_pred_tiling_all, y_chrV_all)
# distances = 2-48, modeled on entire tiling library
# 0.589719

yeast_all_temp_dist_pred_tiling_all = predict(temp_dist_tiling_all_lm, dat_yeast_all_temp_dist)

cor(yeast_all_temp_dist_pred_tiling_all, y_yeast_all)
# distances = 2-48, modeled on entire tiling library
# 0.6658656






# Simple Linear Model Training on entire Random library:
temp_dist_random_all_lm = lm(C0_new~., data=dat_random_all_temp_dist)

# Correlation on Training data:
cor(temp_dist_random_all_lm$fitted.values, y_random_all)
# distances = 2-48
# 0.7379272

# Correlation on Test data (other libraries):
tiling_all_temp_dist_pred_random_all = predict(temp_dist_random_all_lm, dat_tiling_all_temp_dist)

cor(tiling_all_temp_dist_pred_random_all, y_tiling_all)
# distances = 2-48
# 0.6523197

chrV_all_temp_dist_pred_random_all = predict(temp_dist_random_all_lm, dat_chrV_all_temp_dist)

cor(chrV_all_temp_dist_pred_random_all, y_chrV_all)
# distances = 2-48
# 0.5405849

yeast_all_temp_dist_pred_random_all = predict(temp_dist_random_all_lm, dat_yeast_all_temp_dist)

cor(yeast_all_temp_dist_pred_random_all, y_yeast_all)
# distances = 2-48
# 0.6206823




# Simple Linear Model Training on entire chrV library:
temp_dist_chrV_all_lm = lm(C0_new~., data=dat_chrV_all_temp_dist)

# Correlation on Training data:
cor(temp_dist_chrV_all_lm$fitted.values, y_chrV_all)
# distances = 2-48
# 0.6149188

# Correlation on Test data (other libraries):
tiling_all_temp_dist_pred_chrV_all = predict(temp_dist_chrV_all_lm, dat_tiling_all_temp_dist)

cor(tiling_all_temp_dist_pred_chrV_all, y_tiling_all)
# distances = 2-48
# 0.7052136

random_all_temp_dist_pred_chrV_all = predict(temp_dist_chrV_all_lm, dat_random_all_temp_dist)

cor(random_all_temp_dist_pred_chrV_all, y_random_all)
# distances = 2-48
# 0.6429841

yeast_all_temp_dist_pred_chrV_all = predict(temp_dist_chrV_all_lm, dat_yeast_all_temp_dist)

cor(yeast_all_temp_dist_pred_chrV_all, y_yeast_all)
# distances = 2-48
# 0.6642024




# Simple Linear Model Training on entire yeast library:
temp_dist_yeast_all_lm = lm(C0_new~., data=dat_yeast_all_temp_dist)

# Correlation on Training data:
cor(temp_dist_yeast_all_lm$fitted.values, y_yeast_all)
# distances = 2-48
# 0.710464

# Correlation on Test data (other libraries):
tiling_all_temp_dist_pred_yeast_all = predict(temp_dist_yeast_all_lm, dat_tiling_all_temp_dist)

cor(tiling_all_temp_dist_pred_yeast_all, y_tiling_all)
# distances = 2-48
# 0.6869731

random_all_temp_dist_pred_yeast_all = predict(temp_dist_yeast_all_lm, dat_random_all_temp_dist)

cor(random_all_temp_dist_pred_yeast_all, y_random_all)
# distances = 2-48
# 0.6238325

chrV_all_temp_dist_pred_yeast_all = predict(temp_dist_yeast_all_lm, dat_chrV_all_temp_dist)

cor(chrV_all_temp_dist_pred_yeast_all, y_chrV_all)
# distances = 2-48
# 0.5721629







# LightGBM
set.seed(50)
dat_temp_dist.lgb_train_set = dat_temp_dist %>%
  select(-C0_new)
dat_temp_dist.lgb_train_rules = lgb.convert_with_rules(data = dat_temp_dist.lgb_train_set)
dat_temp_dist.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp_dist.lgb_train_rules$data),
                                           label = dat_temp_dist$C0_new,
                                           categorical_feature = handle_categorical(colnames(dat_temp_dist.lgb_train_rules$data)))
dat_temp_dist.lgb_test_set = dat_test_temp_dist %>%
  select(-C0_new)
dat_temp_dist.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp_dist.lgb_test_set, rules = dat_temp_dist.lgb_train_rules$rules)$data)
dat_temp_dist.lgb_test_data = lgb.Dataset(data = dat_temp_dist.lgb_test_matrix,
                                          label = dat_test_temp_dist$C0_new,
                                          categorical_feature = handle_categorical(colnames(dat_temp_dist.lgb_test_matrix)))
dat_temp_dist.lgb_obj = lightgbm(data = dat_temp_dist.lgb_train_data,
                                 params = list(learning_rate = c(0.1),
                                               objective = c("regression"),
                                               min_data_in_leaf = 100,
                                               max_depth = c(5),
                                               num_leaves = c(1000),
                                               lambda_l2 = c(1),
                                               boosting = c("gbdt")),
                                 valids = list(valid = dat_temp_dist.lgb_test_data),
                                 nrounds = 8000,
                                 early_stopping_rounds = 50)

dat_temp_dist.lgb_train_pred = predict(dat_temp_dist.lgb_obj, data = as.matrix(dat_temp_dist.lgb_train_rules$data))
cor(dat_temp_dist.lgb_train_pred, y, method="pearson")
# distances = 2-48
# 0.9863918 (nrounds = 5436)

# Prediciton on tiling test set:
dat_temp_dist.lgb_test_pred = predict(dat_temp_dist.lgb_obj, data = dat_temp_dist.lgb_test_matrix)

cor(dat_temp_dist.lgb_test_pred, y_test, method="pearson")
# distances = 2-48
# 0.7693199

mean((dat_temp_dist.lgb_test_pred - y_test)^2)
# distances = 2-48
# 0.09856109

# Prediction on entire random library:
dat_random_all_temp_dist.lgb_test_set = dat_random_all_temp_dist %>%
  select(-C0_new)
dat_random_all_temp_dist.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_random_all_temp_dist.lgb_test_set, 
                                                                            rules = dat_temp_dist.lgb_train_rules$rules)$data)

dat_random_all_temp_dist.lgb_test_pred = predict(dat_temp_dist.lgb_obj, 
                                                 data = dat_random_all_temp_dist.lgb_test_matrix)

cor(dat_random_all_temp_dist.lgb_test_pred, y_random_all)
# distances = 2-48, modeled on train set
# 0.6059586


# Lasso:
X_temp_dist = dat_temp_dist %>% select(-C0_new)
x_temp_dist = sparse.model.matrix(~., data=X_temp_dist)

set.seed(50)
la.model_temp_dist = cv.glmnet(x_temp_dist, y,
                               alpha=1, family="gaussian")

saveRDS(la.model_temp_dist, "model/lasso_dist_freq_di_AorTdi_bid2.rds")

la.temp_dist_best_lambda = la.model_temp_dist$lambda.min
la.temp_dist_best_lambda

la.model_temp_dist.train_pred = predict(la.model_temp_dist, 
                                        newx = x_temp_dist, 
                                        s=la.temp_dist_best_lambda)

cor(la.model_temp_dist.train_pred, y, method="pearson")
# 0.7318081
mean((la.model_temp_dist.train_pred - y)^2)
# 0.1125691

# Prepare test data:
X_temp_dist_test = dat_test_temp_dist %>% select(-C0_new)
x_temp_dist_test = sparse.model.matrix(~., data=X_temp_dist_test)

# Predict values for test dataset: (lambda.min)
la.model_temp_dist.pred = predict(la.model_temp_dist, 
                                  newx = x_temp_dist_test,
                                  s=la.temp_dist_best_lambda)

cor(la.model_temp_dist.pred, y_test, method="pearson")
# 0.7302937
mean((la.model_temp_dist.pred - y_test)^2)
# 0.1124889

# Predict values for test dataset: (lambda.1se)
la.model_temp_dist.pred = predict(la.model_temp_dist, 
                                  newx = x_temp_dist_test,
                                  s=la.model_all_dist_freq_di_AorTdi_bid$lambda.1se)

cor(la.model_temp_dist.pred, y_test, method="pearson")
# 0.7256155
mean((la.model_temp_dist.pred - y_test)^2)
# 0.1142668

# Prediction on entire random library: (lambda.min)
X_random_all_temp_dist = dat_random_all_temp_dist %>% select(-C0_new)
x_random_all_temp_dist = sparse.model.matrix(~., data=X_random_all_temp_dist)

la.model_temp_dist.random_all_pred = predict(la.model_temp_dist,
                                             newx = x_random_all_temp_dist,
                                             s=la.temp_dist_best_lambda)

cor(la.model_temp_dist.random_all_pred, y_random_all)
# distances = 2-48, modeled on train set
# 0.648069
mean((la.model_temp_dist.random_all_pred - y_random_all)^2)
# 0.08598576

# Prediction on entire random library: (lambda.1se)
la.model_temp_dist.random_all_pred = predict(la.model_temp_dist,
                                             newx = x_random_all_temp_dist,
                                             s=la.model_all_dist_freq_di_AorTdi_bid$lambda.1se)

cor(la.model_temp_dist.random_all_pred, y_random_all)
# distances = 2-48, modeled on train set
# 0.6440565
mean((la.model_temp_dist.random_all_pred - y_random_all)^2)
# 0.08595627



# Add counts of trinucleotides:
dat_temp_dist2 = cbind(dist_freq_di_AorTdi_bid2, dist_freq_di_CorGdi_bid2, 
                       dat%>%select(all_of(trinucleotides), C0_new))
temp_dist2_lm = lm(C0_new~., data=dat_temp_dist2)

cor(temp_dist2_lm$fitted.values, y)
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.6745137
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.70292
# distances = 2-48
# 0.7416589

dat_test_temp_dist2 = cbind(dist_freq_di_AorTdi_bid2_test, dist_freq_di_CorGdi_bid2_test,
                            dat_test%>%select(all_of(trinucleotides), C0_new))
temp_dist2_pred = predict(temp_dist2_lm, dat_test_temp_dist2)

cor(temp_dist2_pred, y_test)
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.673619
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.7017362
# distances = 2-48
# 0.7395515




set.seed(50)
dat_temp_dist2.lgb_train_set = dat_temp_dist2 %>%
  select(-C0_new)
dat_temp_dist2.lgb_train_rules = lgb.convert_with_rules(data = dat_temp_dist2.lgb_train_set)
dat_temp_dist2.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp_dist2.lgb_train_rules$data),
                                            label = dat_temp_dist2$C0_new,
                                            categorical_feature = handle_categorical(colnames(dat_temp_dist2.lgb_train_rules$data)))
dat_temp_dist2.lgb_test_set = dat_test_temp_dist2 %>%
  select(-C0_new)
dat_temp_dist2.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp_dist2.lgb_test_set, rules = dat_temp_dist2.lgb_train_rules$rules)$data)
dat_temp_dist2.lgb_test_data = lgb.Dataset(data = dat_temp_dist2.lgb_test_matrix,
                                           label = dat_test_temp_dist2$C0_new,
                                           categorical_feature = handle_categorical(colnames(dat_temp_dist2.lgb_test_matrix)))
dat_temp_dist2.lgb_obj = lightgbm(data = dat_temp_dist2.lgb_train_data,
                                  params = list(learning_rate = c(0.1),
                                                objective = c("regression"),
                                                min_data_in_leaf = 100,
                                                max_depth = c(5),
                                                num_leaves = c(1000),
                                                lambda_l2 = c(1),
                                                boosting = c("gbdt")),
                                  valids = list(valid = dat_temp_dist2.lgb_test_data),
                                  nrounds = 8000,
                                  early_stopping_rounds = 50)

dat_temp_dist2.lgb_train_pred = predict(dat_temp_dist2.lgb_obj, data = as.matrix(dat_temp_dist2.lgb_train_rules$data))
cor(dat_temp_dist2.lgb_train_pred, y, method="pearson")
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.9583576 (nrounds = 3198)
# distances = 2-48
# 0.9931754 (nrounds = 7395)

dat_temp_dist2.lgb_test_pred = predict(dat_temp_dist2.lgb_obj, data = dat_temp_dist2.lgb_test_matrix)

cor(dat_temp_dist2.lgb_test_pred, y_test, method="pearson")
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.7452596 (nrounds = 3198)
# distances = 2-48
# 0.7733564 (nrounds = 7395)

mean((dat_temp_dist2.lgb_test_pred - y_test)^2)
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.1072083 (nrounds = 3198)
# distances = 2-48
# 0.09697333 (nrounds = 7395)






X_temp_dist2 = dat_temp_dist2 %>% select(-C0_new)
x_temp_dist2 = sparse.model.matrix(~., data=X_temp_dist2)

set.seed(50)
la.model_temp_dist2 = cv.glmnet(x_temp_dist2, y,
                                alpha=1, family="gaussian")

saveRDS(la.model_temp_dist2, "model/lasso_dist_freq_di_AorTdi_bid2_2.rds")

la.temp_dist2_best_lambda = la.model_temp_dist2$lambda.min
la.temp_dist2_best_lambda

la.model_temp_dist2.train_pred = predict(la.model_temp_dist2, 
                                         newx = x_temp_dist2, 
                                         s=la.temp_dist2_best_lambda)

cor(la.model_temp_dist2.train_pred, y, method="pearson")
# 0.741043
mean((la.model_temp_dist2.train_pred - y)^2)
# 0.1092695

# Prepare test data:
X_temp_dist2_test = dat_test_temp_dist2 %>% select(-C0_new)
x_temp_dist2_test = sparse.model.matrix(~., data=X_temp_dist2_test)

# Predict values for test dataset:
la.model_temp_dist2.pred = predict(la.model_temp_dist2, 
                                   newx = x_temp_dist2_test,
                                   s=la.temp_dist2_best_lambda)

cor(la.model_temp_dist2.pred, y_test, method="pearson")
# 0.7391268
mean((la.model_temp_dist2.pred - y_test)^2)
# 0.1093595







stri_count(dat$x50mer[6], regex = "(?=GG.{3}[AT]{2})") + 
  stri_count(dat$x50mer[6], regex = "(?=[AT]{2}.{3}GG)")


stri_count(dat$x50mer[1], regex = "(?=(GG.{3}[AT]{2})|([AT]{2}.{3}GG))")
