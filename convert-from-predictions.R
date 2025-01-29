library(tidyverse)
library(lightgbm)
library(markovchain)
library(modelr)
library(readxl)

#Sequence functions
source("scripts/functions/sequence-functions.R")

#Markov feature functions
source("scripts/functions/markov-functions.R")

#LGBM helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

kfold_obj_with_predictions <- readRDS("data/Created/kfold-obj-with-predictions-network-int.rds")

boostmec_confidence_intervals <- kfold_obj_with_predictions %>%
  mutate(test_set_n = map(test, ~nrow(.x))) %>%
  select(pearson_cor, test_set_n) %>%
  mutate_all(unlist) %>%
  mutate(lower_bound = tanh(atanh(pearson_cor) - qnorm(0.975)/sqrt(test_set_n - 3)),
         upper_bound = tanh(atanh(pearson_cor) + qnorm(0.975)/sqrt(test_set_n - 3)))

write_csv(boostmec_confidence_intervals, "results/validation_confidence_intervals.csv")

lgb_cv_shap <- kfold_obj_with_predictions %>%
  select(lgbm_mod, test_mat)

saveRDS(lgb_cv_shap, "data/Created/lgb-cv-shap-network-int.rds")

