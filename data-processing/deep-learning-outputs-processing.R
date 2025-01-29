tiling_nrow = 82368
set.seed(50)
tiling_train_indices = sample(1:tiling_nrow, tiling_nrow*.9, replace=FALSE)

random_nrow = 12472
set.seed(50)
random_train_indices = sample(1:random_nrow, random_nrow*.9, replace=FALSE)

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")





########################################################################
# IR_LSTM_NEWC0
########################################################################

tiling_tiling_first_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train = tiling_tiling_first_conv_output[tiling_train_indices,]
tiling_tiling_first_conv_output_test = tiling_tiling_first_conv_output[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train, "data/Created/tiling_tiling_first_conv_output_train.rds")
saveRDS(tiling_tiling_first_conv_output_test, "data/Created/tiling_tiling_first_conv_output_test.rds")

tiling_random_first_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train = tiling_random_first_conv_output[random_train_indices,]
tiling_random_first_conv_output_test = tiling_random_first_conv_output[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train, "data/Created/tiling_random_first_conv_output_train.rds")
saveRDS(tiling_random_first_conv_output_test, "data/Created/tiling_random_first_conv_output_test.rds")



# First in parallel line 1
tiling_tiling_second_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_second_conv_output", row.names = 1)
tiling_tiling_second_conv_output_train = tiling_tiling_second_conv_output[tiling_train_indices,]
tiling_tiling_second_conv_output_test = tiling_tiling_second_conv_output[-tiling_train_indices,]

saveRDS(tiling_tiling_second_conv_output_train, "data/Created/tiling_tiling_second_conv_output_train.rds")
saveRDS(tiling_tiling_second_conv_output_test, "data/Created/tiling_tiling_second_conv_output_test.rds")

tiling_random_second_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_second_conv_output", row.names = 1)
tiling_random_second_conv_output_train = tiling_random_second_conv_output[random_train_indices,]
tiling_random_second_conv_output_test = tiling_random_second_conv_output[-random_train_indices,]

saveRDS(tiling_random_second_conv_output_train, "data/Created/tiling_random_second_conv_output_train.rds")
saveRDS(tiling_random_second_conv_output_test, "data/Created/tiling_random_second_conv_output_test.rds")




# Second in parallel line 1
tiling_tiling_third_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_third_conv_output", row.names = 1)
tiling_tiling_third_conv_output_train = tiling_tiling_third_conv_output[tiling_train_indices,]
tiling_tiling_third_conv_output_test = tiling_tiling_third_conv_output[-tiling_train_indices,]

saveRDS(tiling_tiling_third_conv_output_train, "data/Created/tiling_tiling_third_conv_output_train.rds")
saveRDS(tiling_tiling_third_conv_output_test, "data/Created/tiling_tiling_third_conv_output_test.rds")

tiling_random_third_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_third_conv_output", row.names = 1)
tiling_random_third_conv_output_train = tiling_random_third_conv_output[random_train_indices,]
tiling_random_third_conv_output_test = tiling_random_third_conv_output[-random_train_indices,]

saveRDS(tiling_random_third_conv_output_train, "data/Created/tiling_random_third_conv_output_train.rds")
saveRDS(tiling_random_third_conv_output_test, "data/Created/tiling_random_third_conv_output_test.rds")





# First in parallel line 2
tiling_tiling_fourth_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_fourth_conv_output", row.names = 1)
tiling_tiling_fourth_conv_output_train = tiling_tiling_fourth_conv_output[tiling_train_indices,]
tiling_tiling_fourth_conv_output_test = tiling_tiling_fourth_conv_output[-tiling_train_indices,]

saveRDS(tiling_tiling_fourth_conv_output_train, "data/Created/tiling_tiling_fourth_conv_output_train.rds")
saveRDS(tiling_tiling_fourth_conv_output_test, "data/Created/tiling_tiling_fourth_conv_output_test.rds")

tiling_random_fourth_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_fourth_conv_output", row.names = 1)
tiling_random_fourth_conv_output_train = tiling_random_fourth_conv_output[random_train_indices,]
tiling_random_fourth_conv_output_test = tiling_random_fourth_conv_output[-random_train_indices,]

saveRDS(tiling_random_fourth_conv_output_train, "data/Created/tiling_random_fourth_conv_output_train.rds")
saveRDS(tiling_random_fourth_conv_output_test, "data/Created/tiling_random_fourth_conv_output_test.rds")





# Second in parallel line 2
tiling_tiling_fifth_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_fifth_conv_output", row.names = 1)
tiling_tiling_fifth_conv_output_train = tiling_tiling_fifth_conv_output[tiling_train_indices,]
tiling_tiling_fifth_conv_output_test = tiling_tiling_fifth_conv_output[-tiling_train_indices,]

saveRDS(tiling_tiling_fifth_conv_output_train, "data/Created/tiling_tiling_fifth_conv_output_train.rds")
saveRDS(tiling_tiling_fifth_conv_output_test, "data/Created/tiling_tiling_fifth_conv_output_test.rds")

tiling_random_fifth_conv_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_fifth_conv_output", row.names = 1)
tiling_random_fifth_conv_output_train = tiling_random_fifth_conv_output[random_train_indices,]
tiling_random_fifth_conv_output_test = tiling_random_fifth_conv_output[-random_train_indices,]

saveRDS(tiling_random_fifth_conv_output_train, "data/Created/tiling_random_fifth_conv_output_train.rds")
saveRDS(tiling_random_fifth_conv_output_test, "data/Created/tiling_random_fifth_conv_output_test.rds")





# Combine outputs from first layer and 2 parallel lines
tiling_tiling_add_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_add_output", row.names = 1)
tiling_tiling_add_output_train = tiling_tiling_add_output[tiling_train_indices,]
tiling_tiling_add_output_test = tiling_tiling_add_output[-tiling_train_indices,]

saveRDS(tiling_tiling_add_output_train, "data/Created/tiling_tiling_add_output_train.rds")
saveRDS(tiling_tiling_add_output_test, "data/Created/tiling_tiling_add_output_test.rds")

tiling_random_add_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_add_output", row.names = 1)
tiling_random_add_output_train = tiling_random_add_output[random_train_indices,]
tiling_random_add_output_test = tiling_random_add_output[-random_train_indices,]

saveRDS(tiling_random_add_output_train, "data/Created/tiling_random_add_output_train.rds")
saveRDS(tiling_random_add_output_test, "data/Created/tiling_random_add_output_test.rds")






# LSTM layer output
tiling_tiling_lstm_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_tiling_lstm_output", row.names = 1)
tiling_tiling_lstm_output_train = tiling_tiling_lstm_output[tiling_train_indices,]
tiling_tiling_lstm_output_test = tiling_tiling_lstm_output[-tiling_train_indices,]

saveRDS(tiling_tiling_lstm_output_train, "data/Created/tiling_tiling_lstm_output_train.rds")
saveRDS(tiling_tiling_lstm_output_test, "data/Created/tiling_tiling_lstm_output_test.rds")

tiling_random_lstm_output = read.csv("benchmarks/deep-learning/ir_lstm_newC0_tiling_random_lstm_output", row.names = 1)
tiling_random_lstm_output_train = tiling_random_lstm_output[random_train_indices,]
tiling_random_lstm_output_test = tiling_random_lstm_output[-random_train_indices,]

saveRDS(tiling_random_lstm_output_train, "data/Created/tiling_random_lstm_output_train.rds")
saveRDS(tiling_random_lstm_output_test, "data/Created/tiling_random_lstm_output_test.rds")















########################################################################
# CONV_ONLY
########################################################################

tiling_tiling_first_conv_output.conv_only = read.csv("benchmarks/deep-learning/conv_only_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only = tiling_tiling_first_conv_output.conv_only[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only = tiling_tiling_first_conv_output.conv_only[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only, "data/Created/tiling_tiling_first_conv_output_train_conv_only.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only, "data/Created/tiling_tiling_first_conv_output_test_conv_only.rds")

tiling_random_first_conv_output.conv_only = read.csv("benchmarks/deep-learning/conv_only_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only = tiling_random_first_conv_output.conv_only[random_train_indices,]
tiling_random_first_conv_output_test.conv_only = tiling_random_first_conv_output.conv_only[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only, "data/Created/tiling_random_first_conv_output_train_conv_only.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only, "data/Created/tiling_random_first_conv_output_test_conv_only.rds")




tiling_first_conv_kernels.conv_only = read.csv("benchmarks/deep-learning/conv_only_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only) = as.vector(outer(nucleotides, 1:3, paste0))
saveRDS(tiling_first_conv_kernels.conv_only, "data/Created/tiling_first_conv_kernels_conv_only.rds")

tiling_first_conv_biases.conv_only = read.csv("benchmarks/deep-learning/conv_only_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only, "data/Created/tiling_first_conv_biases_conv_only.rds")








########################################################################
# CONV_ONLY_2
########################################################################

tiling_tiling_first_conv_output.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_2) = as.vector(outer(paste0("feature", 1:16), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_2 = tiling_tiling_first_conv_output.conv_only_2[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_2 = tiling_tiling_first_conv_output.conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_2, "data/Created/tiling_tiling_first_conv_output_train_conv_only_2.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_2, "data/Created/tiling_tiling_first_conv_output_test_conv_only_2.rds")

tiling_tiling_max_pool_output.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_tiling_max_pool_output", row.names = 1)
# colnames(tiling_tiling_max_pool_output.conv_only_2) = as.vector(outer(paste0("feature", 1:16), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_max_pool_output_train.conv_only_2 = tiling_tiling_max_pool_output.conv_only_2[tiling_train_indices,]
tiling_tiling_max_pool_output_test.conv_only_2 = tiling_tiling_max_pool_output.conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_max_pool_output_train.conv_only_2, "data/Created/tiling_tiling_max_pool_output_train_conv_only_2.rds")
saveRDS(tiling_tiling_max_pool_output_test.conv_only_2, "data/Created/tiling_tiling_max_pool_output_test_conv_only_2.rds")


tiling_tiling_batch_norm_output.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_tiling_batch_norm_output", row.names = 1)
# colnames(tiling_tiling_batch_norm_output.conv_only_2) = as.vector(outer(paste0("feature", 1:16), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_batch_norm_output_train.conv_only_2 = tiling_tiling_batch_norm_output.conv_only_2[tiling_train_indices,]
tiling_tiling_batch_norm_output_test.conv_only_2 = tiling_tiling_batch_norm_output.conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_batch_norm_output_train.conv_only_2, "data/Created/tiling_tiling_batch_norm_output_train_conv_only_2.rds")
saveRDS(tiling_tiling_batch_norm_output_test.conv_only_2, "data/Created/tiling_tiling_batch_norm_output_test_conv_only_2.rds")


tiling_tiling_dense_output.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_tiling_dense_output", row.names = 1)
# colnames(tiling_tiling_dense_output.conv_only_2) = as.vector(outer(paste0("feature", 1:16), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_dense_output_train.conv_only_2 = tiling_tiling_dense_output.conv_only_2[tiling_train_indices,]
tiling_tiling_dense_output_test.conv_only_2 = tiling_tiling_dense_output.conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_dense_output_train.conv_only_2, "data/Created/tiling_tiling_dense_output_train_conv_only_2.rds")
saveRDS(tiling_tiling_dense_output_test.conv_only_2, "data/Created/tiling_tiling_dense_output_test_conv_only_2.rds")



tiling_random_first_conv_output.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_2) = as.vector(outer(paste0("feature", 1:16), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_2 = tiling_random_first_conv_output.conv_only_2[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_2 = tiling_random_first_conv_output.conv_only_2[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_2, "data/Created/tiling_random_first_conv_output_train_conv_only_2.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_2, "data/Created/tiling_random_first_conv_output_test_conv_only_2.rds")




tiling_first_conv_kernels.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_2) = as.vector(outer(nucleotides, 1:3, paste0))
rownames(tiling_first_conv_kernels.conv_only_2) = paste0("feature", 1:16)
saveRDS(tiling_first_conv_kernels.conv_only_2, "data/Created/tiling_first_conv_kernels_conv_only_2.rds")

tiling_first_conv_biases.conv_only_2 = read.csv("benchmarks/deep-learning/conv_only_2_tiling_first_conv_biases", row.names = 1)
rownames(tiling_first_conv_biases.conv_only_2) = paste0("feature", 1:16)
saveRDS(tiling_first_conv_biases.conv_only_2, "data/Created/tiling_first_conv_biases_conv_only_2.rds")










########################################################################
# CONV_ONLY_3
########################################################################

tiling_tiling_first_conv_output.conv_only_3 = read.csv("benchmarks/deep-learning/conv_only_3_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_3) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:46, "_", 5:50), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_3 = tiling_tiling_first_conv_output.conv_only_3[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_3 = tiling_tiling_first_conv_output.conv_only_3[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_3, "data/Created/tiling_tiling_first_conv_output_train_conv_only_3.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_3, "data/Created/tiling_tiling_first_conv_output_test_conv_only_3.rds")

tiling_random_first_conv_output.conv_only_3 = read.csv("benchmarks/deep-learning/conv_only_3_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_3) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:46, "_", 5:50), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_3 = tiling_random_first_conv_output.conv_only_3[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_3 = tiling_random_first_conv_output.conv_only_3[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_3, "data/Created/tiling_random_first_conv_output_train_conv_only_3.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_3, "data/Created/tiling_random_first_conv_output_test_conv_only_3.rds")




tiling_first_conv_kernels.conv_only_3 = read.csv("benchmarks/deep-learning/conv_only_3_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_3) = as.vector(outer(nucleotides, 1:5, paste0))
rownames(tiling_first_conv_kernels.conv_only_3) = paste0("feature", 1:48)
saveRDS(tiling_first_conv_kernels.conv_only_3, "data/Created/tiling_first_conv_kernels_conv_only_3.rds")

tiling_first_conv_biases.conv_only_3 = read.csv("benchmarks/deep-learning/conv_only_3_tiling_first_conv_biases", row.names = 1)
rownames(tiling_first_conv_biases.conv_only_3) = paste0("feature", 1:48)
saveRDS(tiling_first_conv_biases.conv_only_3, "data/Created/tiling_first_conv_biases_conv_only_3.rds")










########################################################################
# CONV_ONLY_4
########################################################################

tiling_tiling_first_conv_output.conv_only_4 = read.csv("benchmarks/deep-learning/conv_only_4_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_4) = as.vector(outer(paste0("feature", 1:4), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_4 = tiling_tiling_first_conv_output.conv_only_4[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_4 = tiling_tiling_first_conv_output.conv_only_4[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_4, "data/Created/tiling_tiling_first_conv_output_train_conv_only_4.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_4, "data/Created/tiling_tiling_first_conv_output_test_conv_only_4.rds")

tiling_random_first_conv_output.conv_only_4 = read.csv("benchmarks/deep-learning/conv_only_4_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_4) = as.vector(outer(paste0("feature", 1:4), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_4 = tiling_random_first_conv_output.conv_only_4[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_4 = tiling_random_first_conv_output.conv_only_4[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_4, "data/Created/tiling_random_first_conv_output_train_conv_only_4.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_4, "data/Created/tiling_random_first_conv_output_test_conv_only_4.rds")




tiling_first_conv_kernels.conv_only_4 = read.csv("benchmarks/deep-learning/conv_only_4_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_4) = as.vector(outer(nucleotides, 1:3, paste0))
rownames(tiling_first_conv_kernels.conv_only_4) = paste0("feature", 1:4)
saveRDS(tiling_first_conv_kernels.conv_only_4, "data/Created/tiling_first_conv_kernels_conv_only_4.rds")

tiling_first_conv_biases.conv_only_4 = read.csv("benchmarks/deep-learning/conv_only_4_tiling_first_conv_biases", row.names = 1)
rownames(tiling_first_conv_biases.conv_only_4) = paste0("feature", 1:4)
saveRDS(tiling_first_conv_biases.conv_only_4, "data/Created/tiling_first_conv_biases_conv_only_4.rds")









########################################################################
# CONV_ONLY_5
########################################################################

tiling_tiling_first_conv_output.conv_only_5 = read.csv("benchmarks/deep-learning/conv_only_5_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_5 = tiling_tiling_first_conv_output.conv_only_5[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_5 = tiling_tiling_first_conv_output.conv_only_5[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_5, "data/Created/tiling_tiling_first_conv_output_train_conv_only_5.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_5, "data/Created/tiling_tiling_first_conv_output_test_conv_only_5.rds")

tiling_random_first_conv_output.conv_only_5 = read.csv("benchmarks/deep-learning/conv_only_5_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:48, "_", 3:50), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_5 = tiling_random_first_conv_output.conv_only_5[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_5 = tiling_random_first_conv_output.conv_only_5[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_5, "data/Created/tiling_random_first_conv_output_train_conv_only_5.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_5, "data/Created/tiling_random_first_conv_output_test_conv_only_5.rds")




tiling_first_conv_kernels.conv_only_5 = read.csv("benchmarks/deep-learning/conv_only_5_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_5) = as.vector(outer(nucleotides, 1:3, paste0))
rownames(tiling_first_conv_kernels.conv_only_5) = paste0("feature", 1:8)
saveRDS(tiling_first_conv_kernels.conv_only_5, "data/Created/tiling_first_conv_kernels_conv_only_5.rds")

tiling_first_conv_biases.conv_only_5 = read.csv("benchmarks/deep-learning/conv_only_5_tiling_first_conv_biases", row.names = 1)
rownames(tiling_first_conv_biases.conv_only_5) = paste0("feature", 1:8)
saveRDS(tiling_first_conv_biases.conv_only_5, "data/Created/tiling_first_conv_biases_conv_only_5.rds")








########################################################################
# 2CONV_ONLY
########################################################################

tiling_tiling_first_conv_output.2conv_only = read.csv("benchmarks/deep-learning/2conv_only_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train.2conv_only = tiling_tiling_first_conv_output.2conv_only[tiling_train_indices,]
tiling_tiling_first_conv_output_test.2conv_only = tiling_tiling_first_conv_output.2conv_only[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.2conv_only, "data/Created/tiling_tiling_first_conv_output_train_2conv_only.rds")
saveRDS(tiling_tiling_first_conv_output_test.2conv_only, "data/Created/tiling_tiling_first_conv_output_test_2conv_only.rds")

tiling_random_first_conv_output.2conv_only = read.csv("benchmarks/deep-learning/2conv_only_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train.2conv_only = tiling_random_first_conv_output.2conv_only[random_train_indices,]
tiling_random_first_conv_output_test.2conv_only = tiling_random_first_conv_output.2conv_only[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.2conv_only, "data/Created/tiling_random_first_conv_output_train_2conv_only.rds")
saveRDS(tiling_random_first_conv_output_test.2conv_only, "data/Created/tiling_random_first_conv_output_test_2conv_only.rds")



tiling_tiling_second_conv_output.2conv_only = read.csv("benchmarks/deep-learning/2conv_only_tiling_tiling_second_conv_output", row.names = 1)
tiling_tiling_second_conv_output_train.2conv_only = tiling_tiling_second_conv_output.2conv_only[tiling_train_indices,]
tiling_tiling_second_conv_output_test.2conv_only = tiling_tiling_second_conv_output.2conv_only[-tiling_train_indices,]

saveRDS(tiling_tiling_second_conv_output_train.2conv_only, "data/Created/tiling_tiling_second_conv_output_train_2conv_only.rds")
saveRDS(tiling_tiling_second_conv_output_test.2conv_only, "data/Created/tiling_tiling_second_conv_output_test_2conv_only.rds")

tiling_random_second_conv_output.2conv_only = read.csv("benchmarks/deep-learning/2conv_only_tiling_random_second_conv_output", row.names = 1)
tiling_random_second_conv_output_train.2conv_only = tiling_random_second_conv_output.2conv_only[random_train_indices,]
tiling_random_second_conv_output_test.2conv_only = tiling_random_second_conv_output.2conv_only[-random_train_indices,]

saveRDS(tiling_random_second_conv_output_train.2conv_only, "data/Created/tiling_random_second_conv_output_train_2conv_only.rds")
saveRDS(tiling_random_second_conv_output_test.2conv_only, "data/Created/tiling_random_second_conv_output_test_2conv_only.rds")






########################################################################
# 2CONV_ONLY_2
########################################################################

tiling_tiling_first_conv_output.2conv_only_2 = read.csv("benchmarks/deep-learning/2conv_only_2_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train.2conv_only_2 = tiling_tiling_first_conv_output.2conv_only_2[tiling_train_indices,]
tiling_tiling_first_conv_output_test.2conv_only_2 = tiling_tiling_first_conv_output.2conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.2conv_only_2, "data/Created/tiling_tiling_first_conv_output_train_2conv_only_2.rds")
saveRDS(tiling_tiling_first_conv_output_test.2conv_only_2, "data/Created/tiling_tiling_first_conv_output_test_2conv_only_2.rds")

tiling_random_first_conv_output.2conv_only_2 = read.csv("benchmarks/deep-learning/2conv_only_2_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train.2conv_only_2 = tiling_random_first_conv_output.2conv_only_2[random_train_indices,]
tiling_random_first_conv_output_test.2conv_only_2 = tiling_random_first_conv_output.2conv_only_2[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.2conv_only_2, "data/Created/tiling_random_first_conv_output_train_2conv_only_2.rds")
saveRDS(tiling_random_first_conv_output_test.2conv_only_2, "data/Created/tiling_random_first_conv_output_test_2conv_only_2.rds")



tiling_tiling_second_conv_output.2conv_only_2 = read.csv("benchmarks/deep-learning/2conv_only_2_tiling_tiling_second_conv_output", row.names = 1)
tiling_tiling_second_conv_output_train.2conv_only_2 = tiling_tiling_second_conv_output.2conv_only_2[tiling_train_indices,]
tiling_tiling_second_conv_output_test.2conv_only_2 = tiling_tiling_second_conv_output.2conv_only_2[-tiling_train_indices,]

saveRDS(tiling_tiling_second_conv_output_train.2conv_only_2, "data/Created/tiling_tiling_second_conv_output_train_2conv_only_2.rds")
saveRDS(tiling_tiling_second_conv_output_test.2conv_only_2, "data/Created/tiling_tiling_second_conv_output_test_2conv_only_2.rds")

tiling_random_second_conv_output.2conv_only_2 = read.csv("benchmarks/deep-learning/2conv_only_2_tiling_random_second_conv_output", row.names = 1)
tiling_random_second_conv_output_train.2conv_only_2 = tiling_random_second_conv_output.2conv_only_2[random_train_indices,]
tiling_random_second_conv_output_test.2conv_only_2 = tiling_random_second_conv_output.2conv_only_2[-random_train_indices,]

saveRDS(tiling_random_second_conv_output_train.2conv_only_2, "data/Created/tiling_random_second_conv_output_train_2conv_only_2.rds")
saveRDS(tiling_random_second_conv_output_test.2conv_only_2, "data/Created/tiling_random_second_conv_output_test_2conv_only_2.rds")













########################################################################
# 2CONV_ONLY_3
########################################################################

tiling_tiling_first_conv_output.2conv_only_3 = read.csv("benchmarks/deep-learning/2conv_only_3_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train.2conv_only_3 = tiling_tiling_first_conv_output.2conv_only_3[tiling_train_indices,]
tiling_tiling_first_conv_output_test.2conv_only_3 = tiling_tiling_first_conv_output.2conv_only_3[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.2conv_only_3, "data/Created/tiling_tiling_first_conv_output_train_2conv_only_3.rds")
saveRDS(tiling_tiling_first_conv_output_test.2conv_only_3, "data/Created/tiling_tiling_first_conv_output_test_2conv_only_3.rds")

tiling_random_first_conv_output.2conv_only_3 = read.csv("benchmarks/deep-learning/2conv_only_3_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train.2conv_only_3 = tiling_random_first_conv_output.2conv_only_3[random_train_indices,]
tiling_random_first_conv_output_test.2conv_only_3 = tiling_random_first_conv_output.2conv_only_3[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.2conv_only_3, "data/Created/tiling_random_first_conv_output_train_2conv_only_3.rds")
saveRDS(tiling_random_first_conv_output_test.2conv_only_3, "data/Created/tiling_random_first_conv_output_test_2conv_only_3.rds")



tiling_tiling_second_conv_output.2conv_only_3 = read.csv("benchmarks/deep-learning/2conv_only_3_tiling_tiling_second_conv_output", row.names = 1)
tiling_tiling_second_conv_output_train.2conv_only_3 = tiling_tiling_second_conv_output.2conv_only_3[tiling_train_indices,]
tiling_tiling_second_conv_output_test.2conv_only_3 = tiling_tiling_second_conv_output.2conv_only_3[-tiling_train_indices,]

saveRDS(tiling_tiling_second_conv_output_train.2conv_only_3, "data/Created/tiling_tiling_second_conv_output_train_2conv_only_3.rds")
saveRDS(tiling_tiling_second_conv_output_test.2conv_only_3, "data/Created/tiling_tiling_second_conv_output_test_2conv_only_3.rds")

tiling_random_second_conv_output.2conv_only_3 = read.csv("benchmarks/deep-learning/2conv_only_3_tiling_random_second_conv_output", row.names = 1)
tiling_random_second_conv_output_train.2conv_only_3 = tiling_random_second_conv_output.2conv_only_3[random_train_indices,]
tiling_random_second_conv_output_test.2conv_only_3 = tiling_random_second_conv_output.2conv_only_3[-random_train_indices,]

saveRDS(tiling_random_second_conv_output_train.2conv_only_3, "data/Created/tiling_random_second_conv_output_train_2conv_only_3.rds")
saveRDS(tiling_random_second_conv_output_test.2conv_only_3, "data/Created/tiling_random_second_conv_output_test_2conv_only_3.rds")













########################################################################
# 2CONV_ONLY_4
########################################################################

tiling_tiling_first_conv_output.2conv_only_4 = read.csv("benchmarks/deep-learning/2conv_only_4_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train.2conv_only_4 = tiling_tiling_first_conv_output.2conv_only_4[tiling_train_indices,]
tiling_tiling_first_conv_output_test.2conv_only_4 = tiling_tiling_first_conv_output.2conv_only_4[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.2conv_only_4, "data/Created/tiling_tiling_first_conv_output_train_2conv_only_4.rds")
saveRDS(tiling_tiling_first_conv_output_test.2conv_only_4, "data/Created/tiling_tiling_first_conv_output_test_2conv_only_4.rds")

tiling_random_first_conv_output.2conv_only_4 = read.csv("benchmarks/deep-learning/2conv_only_4_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train.2conv_only_4 = tiling_random_first_conv_output.2conv_only_4[random_train_indices,]
tiling_random_first_conv_output_test.2conv_only_4 = tiling_random_first_conv_output.2conv_only_4[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.2conv_only_4, "data/Created/tiling_random_first_conv_output_train_2conv_only_4.rds")
saveRDS(tiling_random_first_conv_output_test.2conv_only_4, "data/Created/tiling_random_first_conv_output_test_2conv_only_4.rds")



tiling_tiling_second_conv_output.2conv_only_4 = read.csv("benchmarks/deep-learning/2conv_only_4_tiling_tiling_second_conv_output", row.names = 1)
tiling_tiling_second_conv_output_train.2conv_only_4 = tiling_tiling_second_conv_output.2conv_only_4[tiling_train_indices,]
tiling_tiling_second_conv_output_test.2conv_only_4 = tiling_tiling_second_conv_output.2conv_only_4[-tiling_train_indices,]

saveRDS(tiling_tiling_second_conv_output_train.2conv_only_4, "data/Created/tiling_tiling_second_conv_output_train_2conv_only_4.rds")
saveRDS(tiling_tiling_second_conv_output_test.2conv_only_4, "data/Created/tiling_tiling_second_conv_output_test_2conv_only_4.rds")

tiling_random_second_conv_output.2conv_only_4 = read.csv("benchmarks/deep-learning/2conv_only_4_tiling_random_second_conv_output", row.names = 1)
tiling_random_second_conv_output_train.2conv_only_4 = tiling_random_second_conv_output.2conv_only_4[random_train_indices,]
tiling_random_second_conv_output_test.2conv_only_4 = tiling_random_second_conv_output.2conv_only_4[-random_train_indices,]

saveRDS(tiling_random_second_conv_output_train.2conv_only_4, "data/Created/tiling_random_second_conv_output_train_2conv_only_4.rds")
saveRDS(tiling_random_second_conv_output_test.2conv_only_4, "data/Created/tiling_random_second_conv_output_test_2conv_only_4.rds")













########################################################################
# 1CONV_LSTM
########################################################################

tiling_tiling_first_conv_output.1conv_lstm = read.csv("benchmarks/deep-learning/1conv_lstm_tiling_tiling_first_conv_output", row.names = 1)
tiling_tiling_first_conv_output_train.1conv_lstm = tiling_tiling_first_conv_output.1conv_lstm[tiling_train_indices,]
tiling_tiling_first_conv_output_test.1conv_lstm = tiling_tiling_first_conv_output.1conv_lstm[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.1conv_lstm, "data/Created/tiling_tiling_first_conv_output_train_1conv_lstm.rds")
saveRDS(tiling_tiling_first_conv_output_test.1conv_lstm, "data/Created/tiling_tiling_first_conv_output_test_1conv_lstm.rds")

tiling_random_first_conv_output.1conv_lstm = read.csv("benchmarks/deep-learning/1conv_lstm_tiling_random_first_conv_output", row.names = 1)
tiling_random_first_conv_output_train.1conv_lstm = tiling_random_first_conv_output.1conv_lstm[random_train_indices,]
tiling_random_first_conv_output_test.1conv_lstm = tiling_random_first_conv_output.1conv_lstm[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.1conv_lstm, "data/Created/tiling_random_first_conv_output_train_1conv_lstm.rds")
saveRDS(tiling_random_first_conv_output_test.1conv_lstm, "data/Created/tiling_random_first_conv_output_test_1conv_lstm.rds")





tiling_tiling_lstm_output.1conv_lstm = read.csv("benchmarks/deep-learning/1conv_lstm_tiling_tiling_lstm_output", row.names = 1)
tiling_tiling_lstm_output_train.1conv_lstm = tiling_tiling_lstm_output.1conv_lstm[tiling_train_indices,]
tiling_tiling_lstm_output_test.1conv_lstm = tiling_tiling_lstm_output.1conv_lstm[-tiling_train_indices,]

saveRDS(tiling_tiling_lstm_output_train.1conv_lstm, "data/Created/tiling_tiling_lstm_output_train_1conv_lstm.rds")
saveRDS(tiling_tiling_lstm_output_test.1conv_lstm, "data/Created/tiling_tiling_lstm_output_test_1conv_lstm.rds")

tiling_random_lstm_output.1conv_lstm = read.csv("benchmarks/deep-learning/1conv_lstm_tiling_random_lstm_output", row.names = 1)
tiling_random_lstm_output_train.1conv_lstm = tiling_random_lstm_output.1conv_lstm[random_train_indices,]
tiling_random_lstm_output_test.1conv_lstm = tiling_random_lstm_output.1conv_lstm[-random_train_indices,]

saveRDS(tiling_random_lstm_output_train.1conv_lstm, "data/Created/tiling_random_lstm_output_train_1conv_lstm.rds")
saveRDS(tiling_random_lstm_output_test.1conv_lstm, "data/Created/tiling_random_lstm_output_test_1conv_lstm.rds")







########################################################################
# CONV_ONLY_DI
########################################################################

tiling_tiling_first_conv_output.conv_only_di = read.csv("benchmarks/deep-learning/conv_only_di_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_di) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:29, "_", 21:49), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_di = tiling_tiling_first_conv_output.conv_only_di[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_di = tiling_tiling_first_conv_output.conv_only_di[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_di, "data/Created/tiling_tiling_first_conv_output_train_conv_only_di.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_di, "data/Created/tiling_tiling_first_conv_output_test_conv_only_di.rds")

tiling_random_first_conv_output.conv_only_di = read.csv("benchmarks/deep-learning/conv_only_di_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_di) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:29, "_", 21:49), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_di = tiling_random_first_conv_output.conv_only_di[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_di = tiling_random_first_conv_output.conv_only_di[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_di, "data/Created/tiling_random_first_conv_output_train_conv_only_di.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_di, "data/Created/tiling_random_first_conv_output_test_conv_only_di.rds")




tiling_first_conv_kernels.conv_only_di = read.csv("benchmarks/deep-learning/conv_only_di_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_di) = as.vector(outer(dinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_di, "data/Created/tiling_first_conv_kernels_conv_only_di.rds")

tiling_first_conv_biases.conv_only_di = read.csv("benchmarks/deep-learning/conv_only_di_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_di, "data/Created/tiling_first_conv_biases_conv_only_di.rds")









########################################################################
# CONV_ONLY_DI_5
########################################################################

tiling_tiling_first_conv_output.conv_only_di_5 = read.csv("benchmarks/deep-learning/conv_only_di_5_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_di_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:29, "_", 21:49), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_di_5 = tiling_tiling_first_conv_output.conv_only_di_5[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_di_5 = tiling_tiling_first_conv_output.conv_only_di_5[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_di_5, "data/Created/tiling_tiling_first_conv_output_train_conv_only_di_5.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_di_5, "data/Created/tiling_tiling_first_conv_output_test_conv_only_di_5.rds")

tiling_random_first_conv_output.conv_only_di_5 = read.csv("benchmarks/deep-learning/conv_only_di_5_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_di_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:29, "_", 21:49), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_di_5 = tiling_random_first_conv_output.conv_only_di_5[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_di_5 = tiling_random_first_conv_output.conv_only_di_5[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_di_5, "data/Created/tiling_random_first_conv_output_train_conv_only_di_5.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_di_5, "data/Created/tiling_random_first_conv_output_test_conv_only_di_5.rds")




tiling_first_conv_kernels.conv_only_di_5 = read.csv("benchmarks/deep-learning/conv_only_di_5_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_di_5) = as.vector(outer(dinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_di_5, "data/Created/tiling_first_conv_kernels_conv_only_di_5.rds")

tiling_first_conv_biases.conv_only_di_5 = read.csv("benchmarks/deep-learning/conv_only_di_5_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_di_5, "data/Created/tiling_first_conv_biases_conv_only_di_5.rds")







########################################################################
# CONV_ONLY_TRI
########################################################################

tiling_tiling_first_conv_output.conv_only_tri = read.csv("benchmarks/deep-learning/conv_only_tri_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri = tiling_tiling_first_conv_output.conv_only_tri[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri = tiling_tiling_first_conv_output.conv_only_tri[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri.rds")

tiling_random_first_conv_output.conv_only_tri = read.csv("benchmarks/deep-learning/conv_only_tri_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri = tiling_random_first_conv_output.conv_only_tri[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri = tiling_random_first_conv_output.conv_only_tri[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri, "data/Created/tiling_random_first_conv_output_train_conv_only_tri.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri, "data/Created/tiling_random_first_conv_output_test_conv_only_tri.rds")




tiling_first_conv_kernels.conv_only_tri = read.csv("benchmarks/deep-learning/conv_only_tri_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri) = as.vector(outer(trinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri, "data/Created/tiling_first_conv_kernels_conv_only_tri.rds")

tiling_first_conv_biases.conv_only_tri = read.csv("benchmarks/deep-learning/conv_only_tri_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri, "data/Created/tiling_first_conv_biases_conv_only_tri.rds")







########################################################################
# CONV_ONLY_TRI_2
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_2 = read.csv("benchmarks/deep-learning/conv_only_tri_2_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_2) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:18, "_", 31:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_2 = tiling_tiling_first_conv_output.conv_only_tri_2[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_2 = tiling_tiling_first_conv_output.conv_only_tri_2[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_2, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_2.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_2, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_2.rds")

tiling_random_first_conv_output.conv_only_tri_2 = read.csv("benchmarks/deep-learning/conv_only_tri_2_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_2) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:18, "_", 31:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_2 = tiling_random_first_conv_output.conv_only_tri_2[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_2 = tiling_random_first_conv_output.conv_only_tri_2[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_2, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_2.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_2, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_2.rds")




tiling_first_conv_kernels.conv_only_tri_2 = read.csv("benchmarks/deep-learning/conv_only_tri_2_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_2) = as.vector(outer(trinucleotides, 1:31, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_2, "data/Created/tiling_first_conv_kernels_conv_only_tri_2.rds")

tiling_first_conv_biases.conv_only_tri_2 = read.csv("benchmarks/deep-learning/conv_only_tri_2_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_2, "data/Created/tiling_first_conv_biases_conv_only_tri_2.rds")







########################################################################
# CONV_ONLY_TRI_3
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_3 = read.csv("benchmarks/deep-learning/conv_only_tri_3_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_3) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:8, "_", 41:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_3 = tiling_tiling_first_conv_output.conv_only_tri_3[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_3 = tiling_tiling_first_conv_output.conv_only_tri_3[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_3, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_3.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_3, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_3.rds")

tiling_random_first_conv_output.conv_only_tri_3 = read.csv("benchmarks/deep-learning/conv_only_tri_3_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_3) = as.vector(outer(paste0("feature", 1:48), paste0("position", 1:8, "_", 41:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_3 = tiling_random_first_conv_output.conv_only_tri_3[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_3 = tiling_random_first_conv_output.conv_only_tri_3[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_3, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_3.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_3, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_3.rds")




tiling_first_conv_kernels.conv_only_tri_3 = read.csv("benchmarks/deep-learning/conv_only_tri_3_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_3) = as.vector(outer(trinucleotides, 1:41, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_3, "data/Created/tiling_first_conv_kernels_conv_only_tri_3.rds")

tiling_first_conv_biases.conv_only_tri_3 = read.csv("benchmarks/deep-learning/conv_only_tri_3_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_3, "data/Created/tiling_first_conv_biases_conv_only_tri_3.rds")








########################################################################
# CONV_ONLY_TRI_4
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_4 = read.csv("benchmarks/deep-learning/conv_only_tri_4_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_4) = as.vector(outer(paste0("feature", 1:4), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_4 = tiling_tiling_first_conv_output.conv_only_tri_4[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_4 = tiling_tiling_first_conv_output.conv_only_tri_4[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_4, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_4.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_4, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_4.rds")

tiling_random_first_conv_output.conv_only_tri_4 = read.csv("benchmarks/deep-learning/conv_only_tri_4_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_4) = as.vector(outer(paste0("feature", 1:4), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_4 = tiling_random_first_conv_output.conv_only_tri_4[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_4 = tiling_random_first_conv_output.conv_only_tri_4[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_4, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_4.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_4, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_4.rds")




tiling_first_conv_kernels.conv_only_tri_4 = read.csv("benchmarks/deep-learning/conv_only_tri_4_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_4) = as.vector(outer(trinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_4, "data/Created/tiling_first_conv_kernels_conv_only_tri_4.rds")

tiling_first_conv_biases.conv_only_tri_4 = read.csv("benchmarks/deep-learning/conv_only_tri_4_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_4, "data/Created/tiling_first_conv_biases_conv_only_tri_4.rds")







########################################################################
# CONV_ONLY_TRI_5
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_5 = read.csv("benchmarks/deep-learning/conv_only_tri_5_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_5 = tiling_tiling_first_conv_output.conv_only_tri_5[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_5 = tiling_tiling_first_conv_output.conv_only_tri_5[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_5, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_5.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_5, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_5.rds")

tiling_random_first_conv_output.conv_only_tri_5 = read.csv("benchmarks/deep-learning/conv_only_tri_5_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_5) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_5 = tiling_random_first_conv_output.conv_only_tri_5[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_5 = tiling_random_first_conv_output.conv_only_tri_5[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_5, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_5.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_5, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_5.rds")




tiling_first_conv_kernels.conv_only_tri_5 = read.csv("benchmarks/deep-learning/conv_only_tri_5_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_5) = as.vector(outer(trinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_5, "data/Created/tiling_first_conv_kernels_conv_only_tri_5.rds")

tiling_first_conv_biases.conv_only_tri_5 = read.csv("benchmarks/deep-learning/conv_only_tri_5_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_5, "data/Created/tiling_first_conv_biases_conv_only_tri_5.rds")







########################################################################
# CONV_ONLY_TRI_7
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_7 = read.csv("benchmarks/deep-learning/conv_only_tri_7_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_7) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:18, "_", 31:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_7 = tiling_tiling_first_conv_output.conv_only_tri_7[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_7 = tiling_tiling_first_conv_output.conv_only_tri_7[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_7, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_7.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_7, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_7.rds")

tiling_random_first_conv_output.conv_only_tri_7 = read.csv("benchmarks/deep-learning/conv_only_tri_7_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_7) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:18, "_", 31:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_7 = tiling_random_first_conv_output.conv_only_tri_7[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_7 = tiling_random_first_conv_output.conv_only_tri_7[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_7, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_7.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_7, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_7.rds")




tiling_first_conv_kernels.conv_only_tri_7 = read.csv("benchmarks/deep-learning/conv_only_tri_7_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_7) = as.vector(outer(trinucleotides, 1:31, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_7, "data/Created/tiling_first_conv_kernels_conv_only_tri_7.rds")

tiling_first_conv_biases.conv_only_tri_7 = read.csv("benchmarks/deep-learning/conv_only_tri_7_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_7, "data/Created/tiling_first_conv_biases_conv_only_tri_7.rds")










########################################################################
# CONV_ONLY_TRI_5_FOURIER_RESID
########################################################################

tiling_tiling_first_conv_output.conv_only_tri_5_fourier_resid = read.csv("benchmarks/deep-learning/conv_only_tri_5_fourier_resid_tiling_tiling_first_conv_output", row.names = 1)
colnames(tiling_tiling_first_conv_output.conv_only_tri_5_fourier_resid) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_tiling_first_conv_output_train.conv_only_tri_5_fourier_resid = tiling_tiling_first_conv_output.conv_only_tri_5_fourier_resid[tiling_train_indices,]
tiling_tiling_first_conv_output_test.conv_only_tri_5_fourier_resid = tiling_tiling_first_conv_output.conv_only_tri_5_fourier_resid[-tiling_train_indices,]

saveRDS(tiling_tiling_first_conv_output_train.conv_only_tri_5_fourier_resid, "data/Created/tiling_tiling_first_conv_output_train_conv_only_tri_5_fourier_resid.rds")
saveRDS(tiling_tiling_first_conv_output_test.conv_only_tri_5_fourier_resid, "data/Created/tiling_tiling_first_conv_output_test_conv_only_tri_5_fourier_resid.rds")

tiling_random_first_conv_output.conv_only_tri_5_fourier_resid = read.csv("benchmarks/deep-learning/conv_only_tri_5_fourier_resid_tiling_random_first_conv_output", row.names = 1)
colnames(tiling_random_first_conv_output.conv_only_tri_5_fourier_resid) = as.vector(outer(paste0("feature", 1:8), paste0("position", 1:28, "_", 21:48), paste, sep="_"))
tiling_random_first_conv_output_train.conv_only_tri_5_fourier_resid = tiling_random_first_conv_output.conv_only_tri_5_fourier_resid[random_train_indices,]
tiling_random_first_conv_output_test.conv_only_tri_5_fourier_resid = tiling_random_first_conv_output.conv_only_tri_5_fourier_resid[-random_train_indices,]

saveRDS(tiling_random_first_conv_output_train.conv_only_tri_5_fourier_resid, "data/Created/tiling_random_first_conv_output_train_conv_only_tri_5_fourier_resid.rds")
saveRDS(tiling_random_first_conv_output_test.conv_only_tri_5_fourier_resid, "data/Created/tiling_random_first_conv_output_test_conv_only_tri_5_fourier_resid.rds")




tiling_first_conv_kernels.conv_only_tri_5_fourier_resid = read.csv("benchmarks/deep-learning/conv_only_tri_5_fourier_resid_tiling_first_conv_kernels", row.names = 1)
colnames(tiling_first_conv_kernels.conv_only_tri_5_fourier_resid) = as.vector(outer(trinucleotides, 1:21, paste0))
saveRDS(tiling_first_conv_kernels.conv_only_tri_5_fourier_resid, "data/Created/tiling_first_conv_kernels_conv_only_tri_5_fourier_resid.rds")

tiling_first_conv_biases.conv_only_tri_5_fourier_resid = read.csv("benchmarks/deep-learning/conv_only_tri_5_fourier_resid_tiling_first_conv_biases", row.names = 1)
saveRDS(tiling_first_conv_biases.conv_only_tri_5_fourier_resid, "data/Created/tiling_first_conv_biases_conv_only_tri_5_fourier_resid.rds")





