library(tidyverse)
library(glmnet)

dat = readRDS("data/Created/processed.rds")
dat_105 = readRDS("data/Created/processed_105_int.rds")

########## 10.5bp interaction

nucleotides <- c("A", "C", "G", "T")
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")

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

model_cols_105 <- c(ps1, ps2, ps3, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count_105")

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

x_105 = model.matrix(~., data=X_105)
y_105 = dat_105$C0

la.model_cv_105 = cv.glmnet(x_105, y_105, alpha=1, family="gaussian")
plot(lamodel_cv_105)

la.model_105 = glmnet(x_105, y_105, alpha=1, model="gaussian")

saveRDS(la.model_105, "model/lasso-105-int")


la.model_cv_105 = cv.glmnet(x_105, y_105, alpha=1, family="gaussian")

saveRDS(la.model_cv_105, "model/lasso-cv-105-int")
# plot(lamodel_cv_105)

x_105.sparse = sparse.model.matrix(~., data=X_105)


#find optimal lambda value that minimizes test MSE
# best_lambda = cv_model_105$lambda.min
# best_lambda
# 
# #produce plot of test MSE by lambda value
# plot(cv_model_105) 