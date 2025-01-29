library(tidyverse)
library(glmnet)

dat_105 = readRDS("data/Created/processed_105_int.rds")

########## 10.5bp interaction

nucleotides <- c("A", "C", "G", "T")
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")
ps4 <- paste0("X", 1:47, "tetra")
ps5 <- paste0("X", 1:46, "penta")

mono_int_indices_105 = c(9,10,11)
for(i in mono_int_indices_105) {
  assign(paste0("int", i), paste0("X", 1:(50-i), "int", i))
}

di_int_indices_105 = c(8,9,10)
for(i in di_int_indices_105) {
  assign(paste0("di_int", i), paste0("X", 1:(48-i), "di_int", i))
}

dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

thermo_feat <-c("tm1", "tm2", "tm3", "tm4", "grna_energy", "grna_scaffold_energy")

model_cols_105 <- c(ps1, ps2, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count_105")

for(i in mono_int_indices_105) {
  model_cols_105 <- c(model_cols_105, get(paste0("int", i)))
}

for(i in di_int_indices_105) {
  model_cols_105 <- c(model_cols_105, get(paste0("di_int", i)))
}

for(i in mono_int_indices_105) {
  model_cols_105 <- c(model_cols_105, paste0(dinucleotides, "int", i))
}

X_105 = dat_105 %>%
  select(all_of(model_cols_105))

x_105 = sparse.model.matrix(~., data=X_105)
y_105 = dat_105$C0

set.seed(50)
la.model_105 = cv.glmnet(x_105, y_105, alpha=1, family="gaussian")
ri.model_105 = cv.glmnet(x_105, y_105, alpha=0, family="gaussian")

plot(la.model_105)
plot(ri.model_105)

la.105_best_lambda = la.model_105$lambda.min
la.105_best_lambda

ri.105_best_lambda = ri.model_105$lambda.min
ri.105_best_lambda

la.model_105.train_pred = predict(la.model_105, newx = x_105, s=la.105_best_lambda)
ri.model_105.train_pred = predict(ri.model_105, newx = x_105, s=ri.105_best_lambda)

cor(la.model_105.train_pred, dat_105$C0, method="pearson")
# 0.5870294

cor(ri.model_105.train_pred, dat_105$C0, method="pearson")
# 0.7163692

dat_test_105 = readRDS("data/Created/processed_105_int_test.rds")

X_105_test = dat_test_105 %>%
  select(all_of(model_cols_105))

x_105_test = sparse.model.matrix(~., data=X_105_test)

la.model_105.pred = predict(la.model_105, newx = x_105_test, s=la.105_best_lambda)
ri.model_105.pred = predict(ri.model_105, newx = x_105_test, s=ri.105_best_lambda)

cor(la.model_105.pred, dat_test_105$C0, method="pearson")
# 0.435631 with up to ps3
# 0.4506208 with up to ps4
# 0.4608161 with up to ps5

cor(ri.model_105.pred, dat_test_105$C0, method="pearson")
# 0.4081583 with up to ps3
# 0.4278197 with up to ps4
# 0.4314025 with up to ps5

saveRDS(la.model_105, "model/lasso-105-int")
saveRDS(ri.model_105, "model/ridge-105-int")