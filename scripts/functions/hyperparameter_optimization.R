# Hyperparameter Optimization

library(rBayesianOptimization)

parameter_bounds <- list(learning_rate = c(0.01, 0.1),
                         min_data = as.integer(c(0, 300)),
                         max_depth = as.integer(c(-1, 12)),
                         max_bin = as.integer(c(60, 500)),
                         num_leaves = as.integer(c(2, 2000)),
                         feature_fraction = c(0.1, 1),
                         bagging_fraction = c(0.1, 1),
                         lambda_l2 = c(0, 10))

bayes_lgb_optimize = function(train_X, train_y, test_X, test_y, init_points, n_iter) {
  bayes.lgb_train_rules = lgb.convert_with_rules(data = train_X)
  bayes.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = test_X, rules = bayes.lgb_train_rules$rules)$data)

  bayes_func = function(learning_rate, min_data, max_depth, num_leaves, max_bin,
                        feature_fraction, bagging_fraction, lambda_l2) {
    bayes.lgb_train_data = lgb.Dataset(data = as.matrix(bayes.lgb_train_rules$data),
                                       label = train_y,
                                       feature_pre_filter = FALSE,
                                       free_raw_data = FALSE,
                                       categorical_feature = handle_categorical(colnames(bayes.lgb_train_rules$data)))
    bayes.lgb_test_data = lgb.Dataset(data = bayes.lgb_test_matrix,
                                      label = test_y,
                                      feature_pre_filter = FALSE,
                                      categorical_feature = handle_categorical(colnames(bayes.lgb_test_matrix)))
    
    bayes.lgb_obj = lightgbm(data = bayes.lgb_train_data,
                             params = list(learning_rate = learning_rate,
                                           objective = c("regression"),
                                           min_data = min_data,
                                           max_depth = max_depth,
                                           num_leaves = num_leaves,
                                           feature_fraction = feature_fraction,
                                           bagging_fraction = bagging_fraction,
                                           lambda_l2 = lambda_l2,
                                           boosting = c("gbdt")),
                             valids = list(valid = bayes.lgb_test_data),
                             nrounds = 8000,
                             early_stopping_rounds = 50)
    
    bayes.lgb_test_pred = predict(bayes.lgb_obj, data = bayes.lgb_test_matrix)
    
    list(Score = cor(bayes.lgb_test_pred, test_y, method="pearson"),
         Pred = bayes.lgb_test_pred)
  }
  
  return(BayesianOptimization(bayes_func,
                       bounds = parameter_bounds,
                       init_points = init_points,
                       n_iter = n_iter))
  
}
