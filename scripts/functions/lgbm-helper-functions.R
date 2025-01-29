#Lightgbm helper functions
#lightgbm, tidyverse, modelr required

# #User defined evaluation function
# #Spearman correlation
# eval_spearman <- function(preds, dtrain){
#   
#   labels <- getinfo(dtrain, "label")
#   spearman_cor <- cor(preds, labels, method = "spearman")
#   return(list(name = "Spearman", value = spearman_cor, higher_better = TRUE))
#   
# }

#User defined evaluation function
#Pearson correlation
eval_pearson <- function(preds, dtrain){
  
  labels <- getinfo(dtrain, "label")
  pearson_cor <- cor(preds, labels, method = "pearson")
  return(list(name = "Pearson", value = pearson_cor, higher_better = TRUE))
  
}

#Helper function to denote index of categorical columns in data
# handle_categorical <- function(feat_vec){
#   
#   cat_matches <- grep("mono|di|tri|tetra|penta", feat_vec)
#   
#   if(length(cat_matches) == 0) {
#     return(NULL)
#   } else {
#     return(cat_matches)
#   }
# }

handle_categorical <- function(feat_vec){
  
  cat_matches <- which(startsWith(feat_vec, "X"))
  
  if(length(cat_matches) == 0) {
    return(NULL)
  } else {
    return(cat_matches)
  }
}

#Function to perform cross-validation of a LightGBM model
#given columns and parameters. If importance_threshold is provided, feature selection is applied beforehand with fixed number of iterations
#Full model is then fitted using early stopping and max n_rounds iterations
# lgb_cv_eval <- function(kfold_obj, feature_cols, parameters, n_rounds, early_stop, feature_importance_threshold = 0.001){
#   
#   #For each fold
#   #Set up initial train pools
#   #Train initial model for feature selection
#   #Subset features and create updated pools
#   #Train model
#   #Make predictions and evaluate
#   if(!is.na(feature_importance_threshold)){
#     
#     kfold_obj <- kfold_obj %>%
#       mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
#              initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
#              initial_train_data = map(initial_train_rules, ~.x$data),
#              initial_train_lgb_obj = map2(initial_train_data, train,
#                                           ~ lgb.Dataset(data = as.matrix(.x),
#                                                         label = .y$C0,
#                                                         categorical_feature = handle_categorical(feature_cols))),
#              initial_mod = map(initial_train_lgb_obj, ~ lgb.train(params = parameters,
#                                                                   data = .x,
#                                                                   nrounds = 500)),
#              feature_importance = map(initial_mod, ~ as.data.frame(lgb.importance(.x))),
#              features = map(feature_importance, ~ .x[.x$Gain > feature_importance_threshold, "Feature"]),
#              train_set = map2(train, features, ~ .x[, .y]),
#              test_set = map2(test, features, ~ .x[, .y]),
#              train_rules = map(train_set, ~ lgb.convert_with_rules(data = .x)),
#              train_lgb_obj = map2(train_rules, train,
#                                   ~ lgb.Dataset(data = as.matrix(.x$data),
#                                                 label = .y$C0,
#                                                 categorical_feature = handle_categorical(colnames(.x$data)))),
#              test_mat = map2(train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y, rules = .x$rules)$data)),
#              test_lgb_obj = pmap(list(train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
#              lgbm_mod = map2(train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
#                                                                       data = .x,
#                                                                       valids = list(valid = .y),
#                                                                       nrounds = n_rounds,
#                                                                       early_stopping_rounds = early_stop)),
#              test_pred = map2(lgbm_mod, test_mat, predict),
#              spearman_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "spearman")),
#              best_iter = map(lgbm_mod, ~ .x$best_iter))
#   } else {
#     
#     kfold_obj <- kfold_obj %>%
#       mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
#              initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
#              initial_train_data = map(initial_train_rules, ~.x$data),
#              initial_train_lgb_obj = map2(initial_train_data, train,
#                                           ~ lgb.Dataset(data = as.matrix(.x),
#                                                         label = .y$C0,
#                                                         categorical_feature = handle_categorical(feature_cols))),
#              test_set = map(test, ~.x[, feature_cols]),
#              test_mat = map2(initial_train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y,
#                                                                                                rules = .x$rules)$data)),
#              test_lgb_obj = pmap(list(initial_train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
#              lgbm_mod = map2(initial_train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
#                                                                               data = .x,
#                                                                               valids = list(valid = .y),
#                                                                               nrounds = n_rounds,
#                                                                               early_stopping_rounds = early_stop)),
#              feature_importance = map(lgbm_mod, ~ "Not calculated"),
#              features = map(feature_importance, ~ feature_cols),
#              test_pred = map2(lgbm_mod, test_mat, predict),
#              spearman_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "spearman")),
#              best_iter = map(lgbm_mod, ~ .x$best_iter))
#   }
#   
#   avg_spearman <- kfold_obj %>%
#     pull(spearman_cor) %>%
#     unlist() %>%
#     mean()
#   
#   spearman_sd <- kfold_obj %>%
#     pull(spearman_cor) %>%
#     unlist() %>%
#     sd()
#   
#   best_iterations <- kfold_obj %>%
#     pull(best_iter) %>%
#     unlist() %>%
#     paste(collapse = ",")
#   
#   avg_iterations <- kfold_obj %>%
#     pull(best_iter) %>%
#     unlist() %>%
#     mean()
#   
#   selected_features <- kfold_obj %>%
#     pull(features) %>%
#     unlist() %>%
#     table() %>%
#     as.data.frame(stringsAsFactors = FALSE) %>%
#     {.$`.`} %>%
#     paste(collapse = ",")
#   
#   feature_frequency <- kfold_obj %>%
#     pull(features) %>%
#     unlist() %>%
#     table() %>%
#     as.data.frame(stringsAsFactors = FALSE) %>%
#     {.$Freq} %>%
#     paste(collapse = ",")
#   
#   #Return results as df
#   results_df <- data.frame(
#     avg_spearman = avg_spearman,
#     sd = spearman_sd,
#     avg_iterations = avg_iterations,
#     best_iterations = best_iterations,
#     parameters = paste(names(parameters), parameters, sep = ": ", collapse = ", "),
#     select_feat = selected_features,
#     feat_freq = feature_frequency,
#     stringsAsFactors = FALSE)
#   
#   #Return feature importance folds
#   feature_importance <- kfold_obj %>%
#     select(feature_importance)
#   
#   #Return together as list
#   return(list(results = results_df, feature_importance = feature_importance))
#   
# }

#Function to perform cross-validation of a LightGBM model
#given columns and parameters. If importance_threshold is provided, feature selection is applied beforehand with fixed number of iterations
#Full model is then fitted using early stopping and max n_rounds iterations
#Updated with Pearson correlation
lgb_cv_eval <- function(kfold_obj, feature_cols, parameters, n_rounds, early_stop, feature_importance_threshold = 0.001){
  
  #For each fold
  #Set up initial train pools
  #Train initial model for feature selection
  #Subset features and create updated pools
  #Train model
  #Make predictions and evaluate
  if(!is.na(feature_importance_threshold)){
    
    kfold_obj <- kfold_obj %>%
      mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
             initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
             initial_train_data = map(initial_train_rules, ~.x$data),
             initial_train_lgb_obj = map2(initial_train_data, train,
                                          ~ lgb.Dataset(data = as.matrix(.x),
                                                        label = .y$C0,
                                                        categorical_feature = handle_categorical(feature_cols))),
             initial_mod = map(initial_train_lgb_obj, ~ lgb.train(params = parameters,
                                                                  data = .x,
                                                                  nrounds = 500)),
             feature_importance = map(initial_mod, ~ as.data.frame(lgb.importance(.x))),
             features = map(feature_importance, ~ .x[.x$Gain > feature_importance_threshold, "Feature"]),
             train_set = map2(train, features, ~ .x[, .y]),
             test_set = map2(test, features, ~ .x[, .y]),
             train_rules = map(train_set, ~ lgb.convert_with_rules(data = .x)),
             train_lgb_obj = map2(train_rules, train,
                                  ~ lgb.Dataset(data = as.matrix(.x$data),
                                                label = .y$C0,
                                                categorical_feature = handle_categorical(colnames(.x$data)))),
             test_mat = map2(train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y, rules = .x$rules)$data)),
             test_lgb_obj = pmap(list(train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
             lgbm_mod = map2(train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
                                                                      data = .x,
                                                                      valids = list(valid = .y),
                                                                      nrounds = n_rounds,
                                                                      early_stopping_rounds = early_stop)),
             test_pred = map2(lgbm_mod, test_mat, predict),
             pearson_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "pearson")),
             best_iter = map(lgbm_mod, ~ .x$best_iter))
  } else {
    
    kfold_obj <- kfold_obj %>%
      mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
             initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
             initial_train_data = map(initial_train_rules, ~.x$data),
             initial_train_lgb_obj = map2(initial_train_data, train,
                                          ~ lgb.Dataset(data = as.matrix(.x),
                                                        label = .y$C0,
                                                        categorical_feature = handle_categorical(feature_cols))),
             test_set = map(test, ~.x[, feature_cols]),
             test_mat = map2(initial_train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y,
                                                                                               rules = .x$rules)$data)),
             test_lgb_obj = pmap(list(initial_train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
             lgbm_mod = map2(initial_train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
                                                                              data = .x,
                                                                              valids = list(valid = .y),
                                                                              nrounds = n_rounds,
                                                                              early_stopping_rounds = early_stop)),
             feature_importance = map(lgbm_mod, ~ "Not calculated"),
             features = map(feature_importance, ~ feature_cols),
             test_pred = map2(lgbm_mod, test_mat, predict),
             pearson_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "pearson")),
             best_iter = map(lgbm_mod, ~ .x$best_iter))
  }
  
  avg_pearson <- kfold_obj %>%
    pull(pearson_cor) %>%
    unlist() %>%
    mean()
  
  pearson_sd <- kfold_obj %>%
    pull(pearson_cor) %>%
    unlist() %>%
    sd()
  
  best_iterations <- kfold_obj %>%
    pull(best_iter) %>%
    unlist() %>%
    paste(collapse = ",")
  
  avg_iterations <- kfold_obj %>%
    pull(best_iter) %>%
    unlist() %>%
    mean()
  
  selected_features <- kfold_obj %>%
    pull(features) %>%
    unlist() %>%
    table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    {.$`.`} %>%
    paste(collapse = ",")
  
  feature_frequency <- kfold_obj %>%
    pull(features) %>%
    unlist() %>%
    table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    {.$Freq} %>%
    paste(collapse = ",")
  
  #Return results as df
  results_df <- data.frame(
    avg_pearson = avg_pearson,
    sd = pearson_sd,
    avg_iterations = avg_iterations,
    best_iterations = best_iterations,
    parameters = paste(names(parameters), parameters, sep = ": ", collapse = ", "),
    select_feat = selected_features,
    feat_freq = feature_frequency,
    stringsAsFactors = FALSE)
  
  #Return feature importance folds
  feature_importance <- kfold_obj %>%
    select(feature_importance)
  
  lgbm_mod <- kfold_obj %>%
    select(lgbm_mod) %>%
    unlist()
  
  features <- kfold_obj %>%
    pull(features)
  
  
  #Return together as list
  return(list(results = results_df, feature_importance = feature_importance, lgbm_mod = lgbm_mod,
              features = features))
  
  
}


# #Abridged version of lgb_cv_eval(), returns kfold_obj after adding predictions and evaluation
# lgb_cv_eval_2 <- function(kfold_obj, feature_cols, parameters, n_rounds, early_stop, feature_importance_threshold = 0.001){
#   
#   #For each fold
#   #Set up initial train pools
#   #Train initial model for feature selection
#   #Subset features and create updated pools
#   #Train model
#   #Make predictions and evaluate
#   if(!is.na(feature_importance_threshold)){
#     
#     kfold_obj <- kfold_obj %>%
#       mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
#              initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
#              initial_train_data = map(initial_train_rules, ~.x$data),
#              initial_train_lgb_obj = map2(initial_train_data, train,
#                                           ~ lgb.Dataset(data = as.matrix(.x),
#                                                         label = .y$C0,
#                                                         categorical_feature = handle_categorical(feature_cols))),
#              initial_mod = map(initial_train_lgb_obj, ~ lgb.train(params = parameters,
#                                                                   data = .x,
#                                                                   nrounds = 500)),
#              feature_importance = map(initial_mod, ~ as.data.frame(lgb.importance(.x))),
#              features = map(feature_importance, ~ .x[.x$Gain > feature_importance_threshold, "Feature"]),
#              train_set = map2(train, features, ~ .x[, .y]),
#              test_set = map2(test, features, ~ .x[, .y]),
#              train_rules = map(train_set, ~ lgb.convert_with_rules(data = .x)),
#              train_lgb_obj = map2(train_rules, train,
#                                   ~ lgb.Dataset(data = as.matrix(.x$data),
#                                                 label = .y$C0,
#                                                 categorical_feature = handle_categorical(colnames(.x$data)))),
#              test_mat = map2(train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y, rules = .x$rules)$data)),
#              test_lgb_obj = pmap(list(train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
#              lgbm_mod = map2(train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
#                                                                       data = .x,
#                                                                       valids = list(valid = .y),
#                                                                       nrounds = n_rounds,
#                                                                       early_stopping_rounds = early_stop)),
#              test_pred = map2(lgbm_mod, test_mat, predict),
#              spearman_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "spearman")),
#              best_iter = map(lgbm_mod, ~ .x$best_iter))
#   } else {
#     
#     kfold_obj <- kfold_obj %>%
#       mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
#              initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
#              initial_train_data = map(initial_train_rules, ~.x$data),
#              initial_train_lgb_obj = map2(initial_train_data, train,
#                                           ~ lgb.Dataset(data = as.matrix(.x),
#                                                         label = .y$C0,
#                                                         categorical_feature = handle_categorical(feature_cols))),
#              test_set = map(test, ~.x[, feature_cols]),
#              test_mat = map2(initial_train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y,
#                                                                                                rules = .x$rules)$data)),
#              test_lgb_obj = pmap(list(initial_train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
#              lgbm_mod = map2(initial_train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
#                                                                               data = .x,
#                                                                               valids = list(valid = .y),
#                                                                               nrounds = n_rounds,
#                                                                               early_stopping_rounds = early_stop)),
#              feature_importance = map(lgbm_mod, ~ "Not calculated"),
#              features = map(feature_importance, ~ feature_cols),
#              test_pred = map2(lgbm_mod, test_mat, predict),
#              spearman_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "spearman")),
#              best_iter = map(lgbm_mod, ~ .x$best_iter))
#   }
#   
#   return(kfold_obj)
# }

#Abridged version of lgb_cv_eval(), returns kfold_obj after adding predictions and evaluation
#Updated with Pearson correlation
lgb_cv_eval_2 <- function(kfold_obj, feature_cols, parameters, n_rounds, early_stop, feature_importance_threshold = 0.001){
  
  #For each fold
  #Set up initial train pools
  #Train initial model for feature selection
  #Subset features and create updated pools
  #Train model
  #Make predictions and evaluate
  if(!is.na(feature_importance_threshold)){
    
    kfold_obj <- kfold_obj %>%
      mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
             initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
             initial_train_data = map(initial_train_rules, ~.x$data),
             initial_train_lgb_obj = map2(initial_train_data, train,
                                          ~ lgb.Dataset(data = as.matrix(.x),
                                                        label = .y$C0,
                                                        categorical_feature = handle_categorical(feature_cols))),
             initial_mod = map(initial_train_lgb_obj, ~ lgb.train(params = parameters,
                                                                  data = .x,
                                                                  nrounds = 500)),
             feature_importance = map(initial_mod, ~ as.data.frame(lgb.importance(.x))),
             features = map(feature_importance, ~ .x[.x$Gain > feature_importance_threshold, "Feature"]),
             train_set = map2(train, features, ~ .x[, .y]),
             test_set = map2(test, features, ~ .x[, .y]),
             train_rules = map(train_set, ~ lgb.convert_with_rules(data = .x)),
             train_lgb_obj = map2(train_rules, train,
                                  ~ lgb.Dataset(data = as.matrix(.x$data),
                                                label = .y$C0,
                                                categorical_feature = handle_categorical(colnames(.x$data)))),
             test_mat = map2(train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y, rules = .x$rules)$data)),
             test_lgb_obj = pmap(list(train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
             lgbm_mod = map2(train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
                                                                      data = .x,
                                                                      valids = list(valid = .y),
                                                                      nrounds = n_rounds,
                                                                      early_stopping_rounds = early_stop)),
             test_pred = map2(lgbm_mod, test_mat, predict),
             pearson_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "pearson")),
             best_iter = map(lgbm_mod, ~ .x$best_iter),
             lgb_interpret = map2(lgbm_mod, test_mat, ~lgb.interprete(.x, .y, idxset = nrow(.y))))
  } else {
    
    kfold_obj <- kfold_obj %>%
      mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
             initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
             initial_train_data = map(initial_train_rules, ~.x$data),
             initial_train_lgb_obj = map2(initial_train_data, train,
                                          ~ lgb.Dataset(data = as.matrix(.x),
                                                        label = .y$C0,
                                                        categorical_feature = handle_categorical(feature_cols))),
             test_set = map(test, ~.x[, feature_cols]),
             test_mat = map2(initial_train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y,
                                                                                               rules = .x$rules)$data)),
             test_lgb_obj = pmap(list(initial_train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
             lgbm_mod = map2(initial_train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
                                                                              data = .x,
                                                                              valids = list(valid = .y),
                                                                              nrounds = n_rounds,
                                                                              early_stopping_rounds = early_stop)),
             feature_importance = map(lgbm_mod, ~ "Not calculated"),
             features = map(feature_importance, ~ feature_cols),
             test_pred = map2(lgbm_mod, test_mat, predict),
             pearson_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "pearson")),
             best_iter = map(lgbm_mod, ~ .x$best_iter),
             lgb_interpret = map2(lgbm_mod, test_mat, ~lgb.interprete(.x, .y, idxset = nrow(.y))))
  }
  
  return(kfold_obj)
}

##FIXME
# lgb_cv_interprete <- function(kfold_obj, feature_cols){
#   
#   #For each fold
#   #Set up initial train pools
#   #Train initial model for feature selection
#   #Subset features and create updated pools
#   #Train model
#   #Make predictions and evaluate
#   
#   kfold_obj <- kfold_obj %>%
#     mutate(initial_train_set = map(train, ~ .x[, feature_cols]),
#            initial_train_rules = map(initial_train_set, ~ lgb.convert_with_rules(data = .x)),
#            initial_train_data = map(initial_train_rules, ~.x$data),
#            initial_train_lgb_obj = map2(initial_train_data, train,
#                                         ~ lgb.Dataset(data = as.matrix(.x),
#                                                       label = .y$C0,
#                                                       categorical_feature = handle_categorical(feature_cols))),
#            test_set = map(test, ~.x[, feature_cols]),
#            test_mat = map2(initial_train_rules, test_set, ~ as.matrix(lgb.convert_with_rules(data = .y,
#                                                                                              rules = .x$rules)$data)),
#            test_lgb_obj = pmap(list(initial_train_lgb_obj, test_mat, test), ~ lgb.Dataset.create.valid(..1, data = ..2, label = ..3$C0)),
#            lgbm_mod = map2(initial_train_lgb_obj, test_lgb_obj, ~ lgb.train(params = parameters,
#                                                                             data = .x,
#                                                                             valids = list(valid = .y),
#                                                                             nrounds = n_rounds,
#                                                                             early_stopping_rounds = early_stop)),
#            feature_importance = map(lgbm_mod, ~ "Not calculated"),
#            features = map(feature_importance, ~ feature_cols),
#            test_pred = map2(lgbm_mod, test_mat, predict),
#            pearson_cor = map2(test_pred, test, ~cor(.x, .y$C0, method = "pearson")),
#            best_iter = map(lgbm_mod, ~ .x$best_iter))
#   
#   avg_pearson <- kfold_obj %>%
#     pull(pearson_cor) %>%
#     unlist() %>%
#     mean()
#   
#   pearson_sd <- kfold_obj %>%
#     pull(pearson_cor) %>%
#     unlist() %>%
#     sd()
#   
#   best_iterations <- kfold_obj %>%
#     pull(best_iter) %>%
#     unlist() %>%
#     paste(collapse = ",")
#   
#   avg_iterations <- kfold_obj %>%
#     pull(best_iter) %>%
#     unlist() %>%
#     mean()
#   
#   selected_features <- kfold_obj %>%
#     pull(features) %>%
#     unlist() %>%
#     table() %>%
#     as.data.frame(stringsAsFactors = FALSE) %>%
#     {.$`.`} %>%
#     paste(collapse = ",")
#   
#   feature_frequency <- kfold_obj %>%
#     pull(features) %>%
#     unlist() %>%
#     table() %>%
#     as.data.frame(stringsAsFactors = FALSE) %>%
#     {.$Freq} %>%
#     paste(collapse = ",")
#   
#   #Return results as df
#   results_df <- data.frame(
#     avg_pearson = avg_pearson,
#     sd = pearson_sd,
#     avg_iterations = avg_iterations,
#     best_iterations = best_iterations,
#     parameters = paste(names(parameters), parameters, sep = ": ", collapse = ", "),
#     select_feat = selected_features,
#     feat_freq = feature_frequency,
#     stringsAsFactors = FALSE)
#   
#   #Return feature importance folds
#   feature_importance <- kfold_obj %>%
#     select(feature_importance)
#   
#   lgbm_model <- kfold_obj %>%
#     select(lgbm_model) %>%
#     unlist()
#   
#   
#   #Return together as list
#   return(list(results = results_df, feature_importance = feature_importance, lgbm_model = lgbm_model))
#   
# }
