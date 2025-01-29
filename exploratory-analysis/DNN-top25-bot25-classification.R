library(keras)
# library(mlbench)
# library(dplyr)
# library(magrittr)
# library(neuralnet)

dat <- readRDS("data/Created/processed_ratio.rds")
dat_test <- readRDS("data/Created/processed_test_ratio.rds")

dat_q1_class = dat_q1 %>%
  select(-C0) %>%
  select( all_of(ps1)) %>%
  mutate(high_cyc = F)
dat_q4_class = dat_q4 %>%
  select(-C0) %>%
  select( all_of(ps1)) %>%
  mutate(high_cyc = T)
dat_q1_q4 = bind_rows(dat_q1_class, dat_q4_class)

# One-hot encoding
x_q1_q4 = sparse.model.matrix(high_cyc ~ ., data = dat_q1_q4)

model <- keras_model_sequential()
model %>%
  layer_embedding(input_dim = 151, output_dim = 50) %>%
  layer_simple_rnn(units = 50) %>% 
  layer_dense(units = 1, activation = "relu")

