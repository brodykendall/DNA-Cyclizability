library(tidyverse)
library(lightgbm)
library(markovchain)
library(modelr)
library(readxl)
# library(janitor)
# library(ggpubr)

#Sequence functions
source("scripts/functions/sequence-functions.R")

#Markov feature functions
source("scripts/functions/markov-functions.R")

#LGBM helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

optim_parameters <- list(learning_rate = 0.1,
                         objective = "regression",
                         min_data = 200,
                         max_bin = 120,
                         max_depth = 4,
                         num_leaves = 500,
                         feature_fraction = 1,
                         bagging_fraction = 0.3,
                         lambda_l2 = 1,
                         boosting = "gbdt",
                         num_threads = 14)

iterations <- 903

nucleotides <- c("A", "C", "G", "T")
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

write_csv(tibble(features = model_cols), "model/all-features-network-int.csv")
write_csv(tibble(iterations = iterations), "model/iterations-network-int.csv")
saveRDS(optim_parameters, "model/optim_parameters-network-int.rds")

kfold_obj <- readRDS("data/Created/train-10-fold-network-int.rds")

kfold_obj_with_predictions <- lgb_cv_eval_2(kfold_obj = kfold_obj,
                                            feature_cols = model_cols,
                                            parameters = optim_parameters,
                                            n_rounds = 6000,
                                            early_stop = 50)


png(file="interpret_1.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`1`[[1]])
dev.off()
png(file="interpret_2.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`2`[[1]])
dev.off()
png(file="interpret_3.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`3`[[1]])
dev.off()
png(file="interpret_4.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`4`[[1]])
dev.off()
png(file="interpret_5.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`5`[[1]])
dev.off()
png(file="interpret_6.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`6`[[1]])
dev.off()
png(file="interpret_7.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`7`[[1]])
dev.off()
png(file="interpret_8.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`8`[[1]])
dev.off()
png(file="interpret_9.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`9`[[1]])
dev.off()
png(file="interpret_10.png")
lgb.plot.interpretation(kfold_obj_with_predictions$lgb_interpret$`10`[[1]])
dev.off()

# saveRDS(kfold_obj_with_predictions, "data/Created/kfold-obj-with-predictions-network-int.rds")

boostmec_confidence_intervals <- kfold_obj_with_predictions %>%
  mutate(test_set_n = map(test, ~nrow(.x))) %>%
  select(pearson_cor, test_set_n) %>%
  mutate_all(unlist) %>%
  mutate(lower_bound = tanh(atanh(pearson_cor) - qnorm(0.975)/sqrt(test_set_n - 3)),
         upper_bound = tanh(atanh(pearson_cor) + qnorm(0.975)/sqrt(test_set_n - 3)))

write_csv(boostmec_confidence_intervals, "results/validation_confidence_intervals.csv")

lgb_cv_shap <- kfold_obj_with_predictions %>%
  select(lgbm_mod, test_mat)

saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`1`, "data/Created/lgb-cv-lgbm-mod-1")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`2`, "data/Created/lgb-cv-lgbm-mod-2")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`3`, "data/Created/lgb-cv-lgbm-mod-3")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`4`, "data/Created/lgb-cv-lgbm-mod-4")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`5`, "data/Created/lgb-cv-lgbm-mod-5")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`6`, "data/Created/lgb-cv-lgbm-mod-6")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`7`, "data/Created/lgb-cv-lgbm-mod-7")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`8`, "data/Created/lgb-cv-lgbm-mod-8")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`9`, "data/Created/lgb-cv-lgbm-mod-9")
saveRDS.lgb.Booster(lgb_cv_shap$lgbm_mod$`10`, "data/Created/lgb-cv-lgbm-mod-10")


saveRDS(lgb_cv_shap, "data/Created/lgb-cv-shap-network-int.rds")

library(SHAPforxgboost)

# lgb_cv_shap = kfold_obj_with_predictions %>%
#   mutate(test_set = map2(test, features, ~ .x[, .y]),
#          test_rules = map(test_set, ~ lgb.convert_with_rules(data = .x)),
#          test_data = map(test_rules, ~.x$data),
#          shap = map2(lgbm_mod, test_data, ~ shap.prep(xgb_model = .x, X_train = as.matrix(.y))))
# 
# shap.plot.summary.wrap1(model = lgb_cv_shap$lgbm_mod$`1`, as.matrix(lgb_cv_shap$test_data$`1`), top_n=20)

# png(file="shap_1.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`1`, kfold_obj_with_predictions$test_mat$`1`, top_n=20)
# dev.off()
# png(file="shap_2.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`2`, kfold_obj_with_predictions$test_mat$`2`, top_n=20)
# dev.off()
# png(file="shap_3.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`3`, kfold_obj_with_predictions$test_mat$`3`, top_n=20)
# dev.off()
# png(file="shap_4.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`4`, kfold_obj_with_predictions$test_mat$`4`, top_n=20)
# dev.off()
# png(file="shap_5.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`5`, kfold_obj_with_predictions$test_mat$`5`, top_n=20)
# dev.off()
# png(file="shap_6.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`6`, kfold_obj_with_predictions$test_mat$`6`, top_n=20)
# dev.off()
# png(file="shap_7.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`7`, kfold_obj_with_predictions$test_mat$`7`, top_n=20)
# dev.off()
# png(file="shap_8.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`8`, kfold_obj_with_predictions$test_mat$`8`, top_n=20)
# dev.off()
# png(file="shap_9.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`9`, kfold_obj_with_predictions$test_mat$`9`, top_n=20)
# dev.off()
# png(file="shap_10.png")
# shap.plot.summary.wrap1(model = kfold_obj_with_predictions$lgbm_mod$`10`, kfold_obj_with_predictions$test_mat$`10`, top_n=20)
# dev.off()

lgb_cv_lgbm_mod_1 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-1")
lgb_cv_lgbm_mod_2 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-2")
lgb_cv_lgbm_mod_3 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-3")
lgb_cv_lgbm_mod_4 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-4")
lgb_cv_lgbm_mod_5 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-5")
lgb_cv_lgbm_mod_6 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-6")
lgb_cv_lgbm_mod_7 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-7")
lgb_cv_lgbm_mod_8 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-8")
lgb_cv_lgbm_mod_9 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-9")
lgb_cv_lgbm_mod_10 <- readRDS.lgb.Booster("data/Created/lgb-cv-lgbm-mod-10")


png(file="shap_1.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_1, lgb_cv_shap$test_mat$`1`, top_n=20)
dev.off()
png(file="shap_2.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_2, lgb_cv_shap$test_mat$`2`, top_n=20)
dev.off()
png(file="shap_3.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_3, lgb_cv_shap$test_mat$`3`, top_n=20)
dev.off()
png(file="shap_4.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_4, lgb_cv_shap$test_mat$`4`, top_n=20)
dev.off()
png(file="shap_5.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_5, lgb_cv_shap$test_mat$`5`, top_n=20)
dev.off()
png(file="shap_6.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_6, lgb_cv_shap$test_mat$`6`, top_n=20)
dev.off()
png(file="shap_7.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_7, lgb_cv_shap$test_mat$`7`, top_n=20)
dev.off()
png(file="shap_8.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_8, lgb_cv_shap$test_mat$`8`, top_n=20)
dev.off()
png(file="shap_9.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_9, lgb_cv_shap$test_mat$`9`, top_n=20)
dev.off()
png(file="shap_10.png")
shap.plot.summary.wrap1(model = lgb_cv_lgbm_mod_10, lgb_cv_shap$test_mat$`10`, top_n=20)
dev.off()



