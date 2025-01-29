train_basic_1_subset_eval = readRDS("lgbm/results/train-basic-1-subset-eval.rds")
train_basic_2_subset_eval = readRDS("lgbm/results/train-basic-2-subset-eval.rds")
train_basic_3_subset_eval = readRDS("lgbm/results/train-basic-3-subset-eval.rds")
train_basic_4_subset_eval = readRDS("lgbm/results/train-basic-4-subset-eval.rds")

train_no_int_1_subset_eval = readRDS("lgbm/results/train-no-int-1-subset-eval.rds")
train_no_int_2_subset_eval = readRDS("lgbm/results/train-no-int-2-subset-eval.rds")
train_no_int_3_subset_eval = readRDS("lgbm/results/train-no-int-3-subset-eval.rds")
train_no_int_4_subset_eval = readRDS("lgbm/results/train-no-int-4-subset-eval.rds")
train_no_int_5_subset_eval = readRDS("lgbm/results/train-no-int-5-subset-eval.rds")
train_no_int_6_subset_eval = readRDS("lgbm/results/train-no-int-6-subset-eval.rds")
train_no_int_7_subset_eval = readRDS("lgbm/results/train-no-int-7-subset-eval.rds")
train_no_int_8_subset_eval = readRDS("lgbm/results/train-no-int-8-subset-eval.rds")

train_ps1_int_1_subset_eval = readRDS("lgbm/results/train-ps1-int-1-subset-eval.rds")
train_ps1_int_2_subset_eval = readRDS("lgbm/results/train-ps1-int-2-subset-eval.rds")
train_ps1_int_3_subset_eval = readRDS("lgbm/results/train-ps1-int-3-subset-eval.rds")
train_ps1_int_4_subset_eval = readRDS("lgbm/results/train-ps1-int-4-subset-eval.rds")
train_ps1_int_5_subset_eval = readRDS("lgbm/results/train-ps1-int-5-subset-eval.rds")
train_ps1_int_6_subset_eval = readRDS("lgbm/results/train-ps1-int-6-subset-eval.rds")
train_ps1_int_7_subset_eval = readRDS("lgbm/results/train-ps1-int-7-subset-eval.rds")
train_ps1_int_8_subset_eval = readRDS("lgbm/results/train-ps1-int-8-subset-eval.rds")

train_ps1_int_only_1_subset_eval = readRDS("lgbm/results/train-ps1-int-only-1-subset-eval.rds")

train_ps1_int_no_ratio_2_subset_eval = readRDS("lgbm/results/train-ps1-int-no-ratio-2-subset-eval.rds")
train_ps1_int_no_ratio_3_subset_eval = readRDS("lgbm/results/train-ps1-int-no-ratio-3-subset-eval.rds")

train_ps2_int_3_subset_eval = readRDS("lgbm/results/train-ps2-int-3-subset-eval.rds")

train_ps2_int_no_ratio_1_subset_eval = readRDS("lgbm/results/train-ps2-int-no-ratio-1-subset-eval.rds")

train_all_int_1_subset_eval = readRDS("lgbm/results/train-all-int-1-subset-eval.rds")
train_all_int_2_results = readRDS("lgbm/results/train-all-int-2-results.rds")

train_all_int_no_ratio_1_subset_eval = readRDS("lgbm/results/train-all-int-no-ratio-1-subset-eval.rds")





test_basic_1_cor = readRDS("lgbm/results/test-basic-1-cor.rds")
test_basic_2_cor = readRDS("lgbm/results/test-basic-2-cor.rds")
test_basic_3_cor = readRDS("lgbm/results/test-basic-3-cor.rds")
test_basic_4_cor = readRDS("lgbm/results/test-basic-4-cor.rds")

test_no_int_1_cor = readRDS("lgbm/results/test-no-int-1-cor.rds")
test_no_int_2_cor = readRDS("lgbm/results/test-no-int-2-cor.rds")
test_no_int_3_cor = readRDS("lgbm/results/test-no-int-3-cor.rds")
test_no_int_4_cor = readRDS("lgbm/results/test-no-int-4-cor.rds")
test_no_int_5_cor = readRDS("lgbm/results/test-no-int-5-cor.rds")
test_no_int_6_cor = readRDS("lgbm/results/test-no-int-6-cor.rds")

test_ps1_int_1_cor = readRDS("lgbm/results/test-ps1-int-1-cor.rds")
test_ps1_int_2_cor = readRDS("lgbm/results/test-ps1-int-2-cor.rds")
test_ps1_int_3_cor = readRDS("lgbm/results/test-ps1-int-3-cor.rds")
test_ps1_int_4_cor = readRDS("lgbm/results/test-ps1-int-4-cor.rds")
test_ps1_int_5_cor = readRDS("lgbm/results/test-ps1-int-5-cor.rds")
test_ps1_int_6_cor = readRDS("lgbm/results/test-ps1-int-6-cor.rds")

test_ps1_int_no_ratio_1_cor = readRDS("lgbm/results/test-ps1-int-no-ratio-1-cor.rds")

test_ps2_int_1_cor = readRDS("lgbm/results/test-ps2-int-1-cor.rds")

test_all_int_1_cor = readRDS("lgbm/results/test-all-int-1-cor.rds")
test_all_int_2_cor = readRDS("lgbm/results/test-all-int-2-cor.rds")
test_all_int_3_cor = readRDS("lgbm/results/test-all-int-3-cor.rds")



# Parameters:
#   learning_rate = 0.1
#   min_data = 300
#   max_bin = 60
#   max_depth = 2
#   num_leaves = 30
#   feature_fraction = 0.5
#   bagging_fraction = 0.5
#   n_rounds = 50
#   feature_importance_threshold = 0.1
train_basic_1_subset_eval[[1]]$results
# Avg pearson: 0.3475196, avg_iterations: 28.9
test_basic_1_cor
# 0.5062941

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 5000
#   max_depth = 7
#   num_leaves = 50
#   n_rounds = 100
#   feature_importance_threshold = 0.1
train_basic_2_subset_eval[[1]]$results
# Avg pearson: 0.3445856, avg_iterations: 23.8
test_basic_2_cor
# 0.4974211

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 3000
#   max_depth = 7
#   num_leaves = 50
#   n_rounds = 100
#   feature_importance_threshold = 0.1
train_basic_3_subset_eval[[1]]$results
# Avg pearson: 0.3460465, avg_iterations: 23.7
test_basic_3_cor
# 0.5026839

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 1000
#   max_depth = 7
#   num_leaves = 100
#   n_rounds = 100
#   feature_importance_threshold = 0.1
train_basic_4_subset_eval[[1]]$results
# Avg pearson: 0.3474027, avg_iterations: 25.7
test_basic_4_cor
# 0.5071405




# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 3000
#   max_depth = 7
#   num_leaves = 50
#   n_rounds = 100
#   feature_importance_threshold = 0.1
train_no_int_1_subset_eval[[1]]$results
# Avg pearson: 0.3460465, avg_iterations: 23.7
test_no_int_1_cor
# 0.5026839

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 100
#   feature_importance_threshold = 0.1
train_no_int_2_subset_eval[[1]]$results
# Avg pearson: 0.3473869, avg_iterations: 21.5
test_no_int_2_cor
# 0.5110459

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 100
#   feature_importance_threshold = 0.001
train_no_int_3_subset_eval[[1]]$results
# Avg pearson: 0.3852064, avg_iterations: 70.9
test_no_int_3_cor
# 0.5049713

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_no_int_4_subset_eval[[1]]$results
# Avg pearson: 0.3853302, avg_iterations: 73.6
test_no_int_4_cor
# 0.5049713

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 12
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_no_int_5_subset_eval[[1]]$results
# Avg pearson: 0.3843692, avg_iterations: 54
test_no_int_5_cor
# 0.5050136

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_no_int_6_subset_eval[[1]]$results
# Avg pearson: 0.3868883, avg_iterations: 77.7
test_no_int_6_cor
# 0.5046858

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 12
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0
train_no_int_7_subset_eval[[1]]$results
# Avg pearson: 0.3832211, avg_iterations: 52.9
test_no_int_7_cor
# 

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 2000
#   feature_importance_threshold = 0
train_no_int_8_subset_eval[[1]]$results
# Avg pearson: 0.3991063, avg_iterations: 962.2
test_no_int_8_cor
# 





# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_ps1_int_1_subset_eval[[1]]$results
# Avg pearson: 0.3849654, avg_iterations: 69.7
test_ps1_int_1_cor
# 0.503721

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0
train_ps1_int_2_subset_eval[[1]]$results
# Avg pearson: 0.3857793, avg_iterations: 77
test_ps1_int_2_cor
# 0.5053998

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0
train_ps1_int_3_subset_eval[[1]]$results
# Avg pearson: 0.3948072, avg_iterations: 499.6
test_ps1_int_3_cor
# 0.5026737

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 2000
#   feature_importance_threshold = 0
train_ps1_int_4_subset_eval[[1]]$results
# Avg pearson: 0.4128242, avg_iterations: 1754.9
test_ps1_int_4_cor
# 0.5043449

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 5000
#   feature_importance_threshold = 0
train_ps1_int_5_subset_eval[[1]]$results
# Avg pearson: 0.4211332, avg_iterations: 2934.8
test_ps1_int_5_cor
# 0.5049293

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 12
#   num_leaves = 1000
#   n_rounds = 5000
#   feature_importance_threshold = 0
train_ps1_int_6_subset_eval[[1]]$results
# Avg pearson: 0.4250708, avg_iterations: 2874.2
test_ps1_int_6_cor
# 0.5065852

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 12
#   num_leaves = 3000
#   n_rounds = 5000
#   feature_importance_threshold = 0
train_ps1_int_7_subset_eval[[1]]$results
# Avg pearson: 0.4250708, avg_iterations: 2874.2 (same as ps1_int_6, not interesting)
test_ps1_int_7_cor
# 

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 1000
#   feature_importance_threshold = 0.001
train_ps1_int_8_subset_eval[[1]]$results
# Avg pearson: 0.3968968, avg_iterations: 824
test_ps1_int_8_cor
# 




# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_bin = 120
#   max_depth = 4
#   num_leaves = 20
#   feature_fraction = 0.5
#   bagging_fraction = 0.4
#   lambda_l2 = 1
#   n_rounds = 6000
#   feature_importance_threshold = 0.001
train_ps1_int_only_1_subset_eval[[1]]$results
# Avg pearson: 0.5445895, avg_iterations: 1588.8
test_ps1_int_only_1_cor
# 0.5083728





# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 1000
#   feature_importance_threshold = 0.001
train_ps1_int_no_ratio_1_subset_eval[[1]]$results
# Avg pearson: 0.4147751, avg_iterations: 997.8
test_ps1_int_no_ratio_1_cor
# 0.2742038

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 3000
#   feature_importance_threshold = 0.001
train_ps1_int_no_ratio_2_subset_eval[[1]]$results
# Avg pearson: 0.4132008, avg_iterations: 1334.6
test_ps1_int_no_ratio_2_cor
# 

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_bin = 120
#   max_depth = 4
#   num_leaves = 20
#   feature_fraction = 0.5
#   bagging_fraction = 0.4
#   lambda_l2 = 1
#   n_rounds = 6000
#   feature_importance_threshold = 0.001
train_ps1_int_no_ratio_3_subset_eval[[1]]$results
# Avg pearson: 0.385635, avg_iterations: 190.2
test_ps1_int_no_ratio_3_cor




# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 2000
#   feature_importance_threshold = 0
train_ps2_int_1_subset_eval[[1]]$results
# Avg pearson: 0.5425753, avg_iterations: 1999.9
test_ps2_int_1_cor
# 0.4990151

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 5000
#   feature_importance_threshold = 0
train_ps2_int_2_subset_eval[[1]]$results
# Avg pearson: , avg_iterations: 
test_ps2_int_2_cor
# 

# Parameters:
#   learning_rate = 0.01
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 2000
#   feature_importance_threshold = 0
train_ps2_int_3_subset_eval[[1]]$results
# Avg pearson: 0.5427445, avg_iterations: 2000
test_ps2_int_3_cor





# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_bin = 120
#   max_depth = 4
#   num_leaves = 20
#   feature_fraction = 0.5
#   bagging_fraction = 0.4
#   lambda_l2 = 1
#   n_rounds = 6000
#   feature_importance_threshold = 0.001
train_ps2_int_no_ratio_1_subset_eval[[1]]$results
# Avg pearson: 0.5351102, avg_iterations: 929.1
test_ps2_int_no_ratio_1_cor





# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_all_int_1_subset_eval[[1]]$results
# Avg pearson: 0.4939855, avg_iterations: 418.5
test_all_int_1_cor
# 0.5041774

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 7
#   num_leaves = 80
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_all_int_2_results
# Avg pearson: 0.4925252, avg_iterations: 429.9
test_all_int_2_cor
# 0.5030119

# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_depth = 12
#   num_leaves = 1000
#   n_rounds = 500
#   feature_importance_threshold = 0.001
train_all_int_3_subset_eval[[1]]$results
# Avg pearson: 0.4554265, avg_iterations: 260.1
test_all_int_3_cor
# 0.5046731





# Parameters:
#   learning_rate = 0.1
#   min_data_in_leaf = 100
#   max_bin = 120
#   max_depth = 4
#   num_leaves = 20
#   feature_fraction = 0.5
#   bagging_fraction = 0.4
#   lambda_l2 = 1
#   n_rounds = 6000
#   feature_importance_threshold = 0.001
train_all_int_no_ratio_1_subset_eval[[1]]$results
# Avg pearson: 0.5336693, avg_iterations: 845.2
test_all_int_no_ratio_1_cor
# 



