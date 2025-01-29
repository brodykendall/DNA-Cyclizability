library(tidyverse)
library(lightgbm)
library(markovchain)
library(modelr)

#Sequence functions
source("scripts/functions/sequence-functions.R")

#Markov feature functions
source("scripts/functions/markov-functions.R")

#LGBM helper functions
source("scripts/functions/lgbm-helper-functions.R")

#Set seed for reproducibility
set.seed(50)

dat <- readRDS("data/Created/processed_network_int.rds")

#Features

model_cols <- read_csv("model/all-features-network-int.csv") %>% pull(features)

#Parameters

optim_iter <- read_csv("model/iterations-network-int.csv") %>% pull(iterations)

optim_parameters <- readRDS("model/optim_parameters-network-int.rds")

cutoffs <- quantile(dat$C0, c(0.2, 0.8))


dat_efficient <- dat[dat$C0 > cutoffs[2],]
dat_inefficient <- dat[dat$C0 < cutoffs[1],]

eff_mat_2nd <- position_specific_score_matrices(dat_efficient, order = 2)
ineff_mat_2nd <- position_specific_score_matrices(dat_inefficient, order = 2)
ratio_mat_2nd <- map2(eff_mat_2nd, ineff_mat_2nd, `-`)

dat$ratio_score_2nd <- pssm_sum_score(dat, ratio_mat_2nd)

saveRDS(ratio_mat_2nd, "model/ratio-matrices-order-2.rds")

initial_train_set <- dat[,model_cols]
initial_train_rules <- lgb.convert_with_rules(initial_train_set)
initial_train_data <- initial_train_rules$data
initial_train_lgb_obj <- lgb.Dataset(data = as.matrix(initial_train_data),
                                     label = dat$C0,
                                     categorical_feature = handle_categorical(model_cols))

initial_mod = lgb.train(params = optim_parameters,
                        data = initial_train_lgb_obj,
                        nrounds = 500)

feature_importance <- as.data.frame(lgb.importance(initial_mod))
features <- feature_importance[feature_importance$Gain > 0, "Feature"]
train_set <- dat[, features]
train_rules <- lgb.convert_with_rules(train_set)
train_lgb_obj <- lgb.Dataset(data = as.matrix(train_rules$data),
                             label = dat$C0,
                             categorical_feature = handle_categorical(colnames(train_rules$data)))
lgbm_model <- lgb.train(params = optim_parameters,
                        data = train_lgb_obj,
                        nrounds = optim_iter)

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps_original <- c(ps1, ps2)

top_features_plot <- feature_importance %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()

ggsave("figures/top-features-network-int.png", top_features_plot, width = 6, height = 4)

#Save selected features
write_csv(tibble(features = features), "model/selected_features-network-int.csv")

#Save lgbm data rules
saveRDS(train_rules, "model/training_rules-network-int.rds")

#Save model
lgb.save(lgbm_model, "model/model-network-int.txt")