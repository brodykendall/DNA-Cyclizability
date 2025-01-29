library(tidyverse)
library(glmnet)

dat_no_int = readRDS("data/Created/processed_105_int.rds")

dat_no_int = dat_no_int %>%
  rename(gc_count_no_int = gc_count_105)

nucleotides <- c("A", "C", "G", "T")
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")
ps4 <- paste0("X", 1:47, "tetra")
ps5 <- paste0("X", 1:46, "penta")

dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

thermo_feat <-c("tm1", "tm2", "tm3", "tm4", "grna_energy", "grna_scaffold_energy")

model_cols_no_int <- c(ps1, ps2, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count_no_int")

X_no_int = dat_no_int %>%
  select(all_of(model_cols_no_int))

x_no_int = sparse.model.matrix(~., data=X_no_int)
y_no_int = dat_no_int$C0

set.seed(50)

elnet.model_no_int.cv = cva.glmnet(x_no_int, y_no_int, family="gaussian")
elnet.model_no_int.cv

plot(elnet.model_no_int.cv)

elnet.no_int_best_alpha = elnet.model_no_int.cv$alpha[[7]]
elnet.model_no_int = elnet.model_no_int.cv$modlist[[7]]

plot(elnet.model_no_int)

elnet.no_int_best_lambda = elnet.model_no_int$lambda.min
elnet.no_int_best_lambda

elnet.model_no_int.train_pred = predict(elnet.model_no_int, newx = x_no_int, s=elnet.no_int_best_lambda)

cor(elnet.model_no_int.train_pred, dat_no_int$C0, method="pearson")
# 0.4987124

dat_test_no_int = readRDS("data/Created/processed_105_int_test.rds")

dat_test_no_int = dat_test_no_int %>%
  rename(gc_count_no_int = gc_count_105)

X_no_int_test = dat_test_no_int %>%
  select(all_of(model_cols_no_int))

x_no_int_test = sparse.model.matrix(~., data=X_no_int_test)

elnet.model_no_int.pred = predict(elnet.model_no_int, newx = x_no_int_test, s=elnet.no_int_best_lambda)

cor(elnet.model_no_int.pred, dat_test_no_int$C0, method="pearson")
# 0.3973511

saveRDS(elnet.model_no_int, "model/elnet-no-int")
saveRDS(elnet.model_no_int.cv, "model/elnet-no-int-cv")





model_cols_all_int.1 <- c(ps1)

X_all_int.1 = dat_no_int %>%
  select(all_of(model_cols_all_int.1))

x_all_int.1 = sparse.model.matrix(~.^2, data=X_all_int.1)
y_all_int = dat_no_int$C0

set.seed(50)
elnet.model_all_int.cv.1 = cva.glmnet(x_all_int.1, y_all_int, family="gaussian")
elnet.model_all_int.cv.1

plot(elnet.model_all_int.cv.1)

elnet.all_int_best_alpha.1 = elnet.model_all_int.cv.1$alpha[[8]]
elnet.model_all_int.1 = elnet.model_all_int.cv.1$modlist[[8]]

plot(elnet.model_all_int.1)

elnet.all_int_best_lambda.1 = elnet.model_all_int.1$lambda.min
elnet.all_int_best_lambda.1

elnet.model_all_int.train_pred.1 = predict(elnet.model_all_int.1, newx = x_all_int.1, s=elnet.all_int_best_lambda.1)

cor(elnet.model_all_int.train_pred.1, dat_no_int$C0, method="pearson")
# 0.698277

dat_test_all_int = readRDS("data/Created/processed_105_int_test.rds")

dat_test_all_int = dat_test_all_int %>%
  rename(gc_count_no_int = gc_count_105)

X_all_int_test.1 = dat_test_all_int %>%
  select(all_of(model_cols_all_int.1))

x_all_int_test.1 = sparse.model.matrix(~.^2, data=X_all_int_test.1)

elnet.model_all_int.pred.1 = predict(elnet.model_all_int.1, newx = x_all_int_test.1, s=elnet.all_int_best_lambda.1)

cor(elnet.model_all_int.pred.1, dat_test_all_int$C0, method="pearson")
# 0.5913607

saveRDS(elnet.model_all_int.1, "model/elnet-all-int-1")
saveRDS(elnet.model_all_int.cv.1, "model/elnet-all-int-cv-1")





