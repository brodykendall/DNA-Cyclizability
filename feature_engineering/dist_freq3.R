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

distances = 2:48

dinucleotides_combinations_df = expand.grid(di1=dinucleotides, di2=dinucleotides) %>%
  filter(as.character(di1) <= as.character(di2))
dinucleotides_combinations = split(dinucleotides_combinations_df, seq(nrow(dinucleotides_combinations_df)))

dist_freq_di_di_bid2 = map2(rep(dinucleotides_combinations, length(distances)), 
                            rep(distances, each=136), 
                            ~stri_count(dat$x50mer, 
                                        regex=paste0("(?=(",.x[1,1],".{",.y-2,"}",.x[1,2],")|(",.x[1,2],".{",.y-2,"}",.x[1,1],"))"))) %>%
  bind_cols()

colnames(dist_freq_di_di_bid2) = paste0(rep(paste(dinucleotides_combinations_df$di1, dinucleotides_combinations_df$di2, sep="_"), length(distances)), 
                                        "_dist", rep(distances, each=136))

saveRDS(dist_freq_di_di_bid2, "data/Created/tiling_dist_freq_di_di_bid2.rds")



# Test data:
dist_freq_di_di_bid2_test = map2(rep(dinucleotides_combinations, length(distances)), 
                                 rep(distances, each=136), 
                                 ~stri_count(dat_test$x50mer, 
                                             regex=paste0("(?=(",.x[1,1],".{",.y-2,"}",.x[1,2],")|(",.x[1,2],".{",.y-2,"}",.x[1,1],"))"))) %>%
  bind_cols()

colnames(dist_freq_di_di_bid2_test) = paste0(rep(paste(dinucleotides_combinations_df$di1, dinucleotides_combinations_df$di2, sep="_"), length(distances)), 
                                             "_dist", rep(distances, each=136))

saveRDS(dist_freq_di_di_bid2_test, "data/Created/tiling_dist_freq_di_di_bid2_test.rds")



# Random Library:
dist_freq_di_di_bid2_random_all = map2(rep(dinucleotides_combinations, length(distances)), 
                                       rep(distances, each=136), 
                                       ~stri_count(dat_random_all$x50mer, 
                                                   regex=paste0("(?=(",.x[1,1],".{",.y-2,"}",.x[1,2],")|(",.x[1,2],".{",.y-2,"}",.x[1,1],"))"))) %>%
  bind_cols()

colnames(dist_freq_di_di_bid2_random_all) = paste0(rep(paste(dinucleotides_combinations_df$di1, dinucleotides_combinations_df$di2, sep="_"), length(distances)), 
                                                   "_dist", rep(distances, each=136))

saveRDS(dist_freq_di_di_bid2_random_all, "data/Created/random_all_dist_freq_di_di_bid2.rds")



# ChrV Library:
dist_freq_di_di_bid2_chrV_all = map2(rep(dinucleotides_combinations, length(distances)), 
                                       rep(distances, each=136), 
                                       ~stri_count(dat_chrV_all$x50mer, 
                                                   regex=paste0("(?=(",.x[1,1],".{",.y-2,"}",.x[1,2],")|(",.x[1,2],".{",.y-2,"}",.x[1,1],"))"))) %>%
  bind_cols()

colnames(dist_freq_di_di_bid2_chrV_all) = paste0(rep(paste(dinucleotides_combinations_df$di1, dinucleotides_combinations_df$di2, sep="_"), length(distances)), 
                                                   "_dist", rep(distances, each=136))

saveRDS(dist_freq_di_di_bid2_chrV_all, "data/Created/chrV_all_dist_freq_di_di_bid2.rds")



# Yeast Library:
dist_freq_di_di_bid2_yeast_all = map2(rep(dinucleotides_combinations, length(distances)), 
                                       rep(distances, each=136), 
                                       ~stri_count(dat_yeast_all$x50mer, 
                                                   regex=paste0("(?=(",.x[1,1],".{",.y-2,"}",.x[1,2],")|(",.x[1,2],".{",.y-2,"}",.x[1,1],"))"))) %>%
  bind_cols()

colnames(dist_freq_di_di_bid2_yeast_all) = paste0(rep(paste(dinucleotides_combinations_df$di1, dinucleotides_combinations_df$di2, sep="_"), length(distances)), 
                                                   "_dist", rep(distances, each=136))

saveRDS(dist_freq_di_di_bid2_yeast_all, "data/Created/yeast_all_dist_freq_di_di_bid2.rds")





dat_temp_dist_di_di = cbind(dist_freq_di_di_bid2, dat%>%select(C0_new))
temp_dist_di_di_lm = lm(C0_new~., data=dat_temp_dist_di_di)

cor(temp_dist_di_di_lm$fitted.values, y)
# distances = 5 10 15 20 25 30 35 40 45
# 0.6240238
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.712149
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.7466421
# distances = 2-48
# 0.8011437

dat_test_temp_dist_di_di = cbind(dist_freq_di_di_bid2_test, 
                                 dat_test%>%select(C0_new))
temp_dist_di_di_pred = predict(temp_dist_di_di_lm, dat_test_temp_dist_di_di)

cor(temp_dist_di_di_pred, y_test)
# distances = 5 10 15 20 25 30 35 40 45
# 0.6163852
# distances = 5  6 10 11 15 16 20 21 25 26 30 31 35 36 40 41 45 46
# 0.6979923
# distances = 4  5  6  9 10 11 14 15 16 19 20 21 24 25 26 29 30 31 34 35 36 39 40 41 44 45 46
# 0.7250886
# distances = 2-48
# 0.7666998





set.seed(50)
dat_temp_dist_di_di.lgb_train_set = dat_temp_dist_di_di %>%
  select(-C0_new)
dat_temp_dist_di_di.lgb_train_rules = lgb.convert_with_rules(data = dat_temp_dist_di_di.lgb_train_set)
dat_temp_dist_di_di.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp_dist_di_di.lgb_train_rules$data),
                                                 label = dat_temp_dist_di_di$C0_new,
                                                 categorical_feature = handle_categorical(colnames(dat_temp_dist_di_di.lgb_train_rules$data)))
dat_temp_dist_di_di.lgb_test_set = dat_test_temp_dist_di_di %>%
  select(-C0_new)
dat_temp_dist_di_di.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp_dist_di_di.lgb_test_set, rules = dat_temp_dist_di_di.lgb_train_rules$rules)$data)
dat_temp_dist_di_di.lgb_test_data = lgb.Dataset(data = dat_temp_dist_di_di.lgb_test_matrix,
                                                label = dat_test_temp_dist_di_di$C0_new,
                                                categorical_feature = handle_categorical(colnames(dat_temp_dist_di_di.lgb_test_matrix)))
dat_temp_dist_di_di.lgb_obj = lightgbm(data = dat_temp_dist_di_di.lgb_train_data,
                                       params = list(learning_rate = c(0.1),
                                                     objective = c("regression"),
                                                     min_data_in_leaf = 100,
                                                     max_depth = c(5),
                                                     num_leaves = c(1000),
                                                     lambda_l2 = c(1),
                                                     boosting = c("gbdt")),
                                       valids = list(valid = dat_temp_dist_di_di.lgb_test_data),
                                       nrounds = 8000,
                                       early_stopping_rounds = 50)

dat_temp_dist_di_di.lgb_train_pred = predict(dat_temp_dist_di_di.lgb_obj, data = as.matrix(dat_temp_dist_di_di.lgb_train_rules$data))
cor(dat_temp_dist_di_di.lgb_train_pred, y, method="pearson")
# distances = 2-48
# 0.9839357 (nrounds = 5506)

dat_temp_dist_di_di.lgb_test_pred = predict(dat_temp_dist_di_di.lgb_obj, data = dat_temp_dist_di_di.lgb_test_matrix)

cor(dat_temp_dist_di_di.lgb_test_pred, y_test, method="pearson")
# distances = 2-48
# 0.7977623

mean((dat_temp_dist_di_di.lgb_test_pred - y_test)^2)
# distances = 2-48
# 0.08789541








# Add counts of trinucleotides:
dat_temp_dist_di_di2 = cbind(dist_freq_di_di_bid2, 
                             dat%>%select(all_of(trinucleotides), C0_new))
temp_dist_di_di2_lm = lm(C0_new~., data=dat_temp_dist_di_di2)

cor(temp_dist_di_di2_lm$fitted.values, y)
# distances = 2-48
# 0.8025348

dat_test_temp_dist_di_di2 = cbind(dist_freq_di_di_bid2_test, 
                                  dat_test%>%select(all_of(trinucleotides), C0_new))
temp_dist_di_di2_pred = predict(temp_dist_di_di2_lm, dat_test_temp_dist_di_di2)

cor(temp_dist_di_di2_pred, y_test)
# distances = 2-48
# 0.768055
