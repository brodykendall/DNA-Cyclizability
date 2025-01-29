# Model 2:

#Set seed for reproducibility
set.seed(50)

dat <- readRDS("data/Created/processed_ratio.rds")

#Features:
model_cols <- read_csv("model/reduced/ps1-int/all-features.csv") %>% pull(features)

#Parameters:
optim_iter <- read_csv("model/reduced/ps1-int/iterations.csv") %>% pull(iterations)
optim_parameters <- readRDS("model/reduced/ps1-int/optim_parameters.rds")

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
# features <- feature_importance[feature_importance$Gain > 0, "Feature"]
features <- feature_importance[, "Feature"]
train_set <- dat[, features]
train_rules <- lgb.convert_with_rules(train_set)
train_lgb_obj <- lgb.Dataset(data = as.matrix(train_rules$data),
                             label = dat$C0,
                             categorical_feature = handle_categorical(colnames(train_rules$data)))
lgbm_model <- lgb.train(params = optim_parameters,
                        data = train_lgb_obj,
                        nrounds = iterations)

top_features_plot <- feature_importance %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:20) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()

ggsave("figures/reduced/ps1-int/top-features.png", top_features_plot, width = 6, height = 4)

#Save lgbm data rules
saveRDS(train_rules, "model/reduced/ps1-int/training_rules.rds")

#Save model
lgb.save(lgbm_model, "model/reduced/ps1-int/model.txt")
