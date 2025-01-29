library(tidyverse)
library(glmnet)

# Load data:
dat = readRDS("data/Created/processed_ratio.rds")
dat_AT = readRDS("data/Created/processed_AT_int.rds")
dat_AT = dat_AT %>% setNames(paste0(names(.), "_AT"))
dat_comb = bind_cols(dat, dat_AT)

dat_test = readRDS("data/Created/processed_test_ratio.rds")
dat_test_AT = readRDS("data/Created/processed_AT_int_test.rds")
dat_test_AT = dat_test_AT %>% setNames(paste0(names(.), "_AT"))
dat_test_comb = bind_cols(dat_test, dat_test_AT)

# Position specific variable names 
# (e.g. X1di represents the dinucleotide in the 1st position of the sequence):
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")
ps4 <- paste0("X", 1:47, "tetra")
ps5 <- paste0("X", 1:46, "penta")

ps1_AT <- paste0("X", 1:50, "mono_AT")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")
trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")


# Set up a few basic models:
model_cols_all_int.1 <- c(ps1)
X_all_int.1 = dat %>%
  select(all_of(model_cols_all_int.1))
x_all_int.1 = sparse.model.matrix(~.^2, data=X_all_int.1)
y_all_int = dat$C0

model_cols_all_int.2 <- c(ps2)
X_all_int.2 = dat %>%
  select(all_of(model_cols_all_int.2))
x_all_int.2 = sparse.model.matrix(~.^2, data=X_all_int.2)

x_all_int.12.dimnames<-list(union(rownames(x_all_int.1),rownames(x_all_int.2)),c(colnames(x_all_int.1),colnames(x_all_int.2)))
x_all_int.12<- Matrix(0,nrow=length(x_all_int.12.dimnames[[1]]),ncol=length(x_all_int.12.dimnames[[2]]),dimnames=x_all_int.12.dimnames)
x_all_int.12[rownames(x_all_int.1),colnames(x_all_int.1)]<-x_all_int.1
x_all_int.12[rownames(x_all_int.2),colnames(x_all_int.2)]<-x_all_int.2

# Train some basic models:
set.seed(50)
la.model_all_int.1 = cv.glmnet(x_all_int.1, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.1 = cv.glmnet(x_all_int.1, y_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_all_int.2 = cv.glmnet(x_all_int.2, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.2 = cv.glmnet(x_all_int.2, y_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_all_int.12 = cv.glmnet(x_all_int.12, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.12 = cv.glmnet(x_all_int.12, y_all_int, alpha=0, family="gaussian")

plot(la.model_all_int.1)
plot(ri.model_all_int.1)

plot(la.model_all_int.2)
plot(ri.model_all_int.2)

plot(la.model_all_int.12)
plot(ri.model_all_int.12)

la.all_int_best_lambda.1 = la.model_all_int.1$lambda.min
la.all_int_best_lambda.1
ri.all_int_best_lambda.1 = ri.model_all_int.1$lambda.min
ri.all_int_best_lambda.1

la.all_int_best_lambda.2 = la.model_all_int.2$lambda.min
la.all_int_best_lambda.2
ri.all_int_best_lambda.2 = ri.model_all_int.2$lambda.min
ri.all_int_best_lambda.2

la.all_int_best_lambda.12 = la.model_all_int.12$lambda.min
la.all_int_best_lambda.12
ri.all_int_best_lambda.12 = ri.model_all_int.12$lambda.min
ri.all_int_best_lambda.12

# Predict values for train dataset:
la.model_all_int.train_pred.1 = predict(la.model_all_int.1, newx = x_all_int.1, s=la.all_int_best_lambda.1)
ri.model_all_int.train_pred.1 = predict(ri.model_all_int.1, newx = x_all_int.1, s=ri.all_int_best_lambda.1)

la.model_all_int.train_pred.2 = predict(la.model_all_int.2, newx = x_all_int.2, s=la.all_int_best_lambda.2)
ri.model_all_int.train_pred.2 = predict(ri.model_all_int.2, newx = x_all_int.2, s=ri.all_int_best_lambda.2)

la.model_all_int.train_pred.12 = predict(la.model_all_int.12, newx = x_all_int.12, s=la.all_int_best_lambda.12)
ri.model_all_int.train_pred.12 = predict(ri.model_all_int.12, newx = x_all_int.12, s=ri.all_int_best_lambda.12)

# Prediction accuracy with train data:
cor(la.model_all_int.train_pred.1, dat_no_int$C0, method="pearson")
# 0.6970752
cor(ri.model_all_int.train_pred.1, dat_no_int$C0, method="pearson")
# 0.7022866
cor(la.model_all_int.train_pred.2, dat_no_int$C0, method="pearson")
# 0.7645926
cor(ri.model_all_int.train_pred.2, dat_no_int$C0, method="pearson")
# 0.9038084
cor(la.model_all_int.train_pred.12, dat_no_int$C0, method="pearson")
# 0.7915266
cor(ri.model_all_int.train_pred.12, dat_no_int$C0, method="pearson")
# 0.9055516


# Prepare test data:
X_all_int_test.1 = dat_test %>%
  select(all_of(model_cols_all_int.1))

X_all_int_test.2 = dat_test %>%
  select(all_of(model_cols_all_int.2))

x_all_int_test.1 = sparse.model.matrix(~.^2, data=X_all_int_test.1)

x_all_int_test.2 = sparse.model.matrix(~.^2, data=X_all_int_test.2)

x_all_int_test.12.dimnames<-list(union(rownames(x_all_int_test.1),rownames(x_all_int_test.2)),
                                 c(colnames(x_all_int_test.1),colnames(x_all_int_test.2)))
x_all_int_test.12<- Matrix(0,nrow=length(x_all_int_test.12.dimnames[[1]]),
                           ncol=length(x_all_int_test.12.dimnames[[2]]),
                           dimnames=x_all_int_test.12.dimnames)
x_all_int_test.12[rownames(x_all_int_test.1),colnames(x_all_int_test.1)]<-x_all_int_test.1
x_all_int_test.12[rownames(x_all_int_test.2),colnames(x_all_int_test.2)]<-x_all_int_test.2

# Predict values for test dataset:
la.model_all_int.pred.1 = predict(la.model_all_int.1, newx = x_all_int_test.1, s=la.all_int_best_lambda.1)
ri.model_all_int.pred.1 = predict(ri.model_all_int.1, newx = x_all_int_test.1, s=ri.all_int_best_lambda.1)

la.model_all_int.pred.2 = predict(la.model_all_int.2, newx = x_all_int_test.2, s=la.all_int_best_lambda.2)
ri.model_all_int.pred.2 = predict(ri.model_all_int.2, newx = x_all_int_test.2, s=ri.all_int_best_lambda.2)

la.model_all_int.pred.12 = predict(la.model_all_int.12, newx = x_all_int_test.12, s=la.all_int_best_lambda.12)
ri.model_all_int.pred.12 = predict(ri.model_all_int.12, newx = x_all_int_test.12, s=ri.all_int_best_lambda.12)

# Prediction accuracy with test data:
cor(la.model_all_int.pred.1, dat_test$C0, method="pearson")
# 0.5913363

cor(ri.model_all_int.pred.1, dat_test$C0, method="pearson")
# 0.5830442

cor(la.model_all_int.pred.2, dat_test$C0, method="pearson")
# 0.4583323

cor(ri.model_all_int.pred.2, dat_test$C0, method="pearson")
# 0.5013422

cor(la.model_all_int.pred.12, dat_test$C0, method="pearson")
# 0.4849626

cor(ri.model_all_int.pred.12, dat_test$C0, method="pearson")
# 0.5166425










# Set up models that include the Markov feature:
model_cols_all_int_inc_ratio.1 <- c(ps1, "ratio_score_2nd")

model_cols_all_int_inc_ratio.2 <- c(ps2, "ratio_score_2nd")

X_all_int_inc_ratio.1 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio.1))

X_all_int_inc_ratio.2 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio.2))

x_all_int_inc_ratio.1 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio.1)
x_all_int_inc_ratio.2 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio.2)

y_all_int_inc_ratio = dat$C0

# Train models:
set.seed(50)
la.model_all_int_inc_ratio.1 = cv.glmnet(x_all_int_inc_ratio.1, y_all_int_inc_ratio, alpha=1, family="gaussian")
ri.model_all_int_inc_ratio.1 = cv.glmnet(x_all_int_inc_ratio.1, y_all_int_inc_ratio, alpha=0, family="gaussian")

set.seed(50)
la.model_all_int_inc_ratio.2 = cv.glmnet(x_all_int_inc_ratio.2, y_all_int_inc_ratio, alpha=1, family="gaussian")
ri.model_all_int_inc_ratio.2 = cv.glmnet(x_all_int_inc_ratio.2, y_all_int_inc_ratio, alpha=0, family="gaussian")

plot(la.model_all_int_inc_ratio.1)
plot(ri.model_all_int_inc_ratio.1)

plot(la.model_all_int_inc_ratio.2)
plot(ri.model_all_int_inc_ratio.2)

la.all_int_inc_ratio_best_lambda.1 = la.model_all_int_inc_ratio.1$lambda.min
la.all_int_inc_ratio_best_lambda.1
ri.all_int_inc_ratio_best_lambda.1 = ri.model_all_int_inc_ratio.1$lambda.min
ri.all_int_inc_ratio_best_lambda.1

la.all_int_inc_ratio_best_lambda.2 = la.model_all_int_inc_ratio.2$lambda.min
la.all_int_inc_ratio_best_lambda.2
ri.all_int_inc_ratio_best_lambda.2 = ri.model_all_int_inc_ratio.2$lambda.min
ri.all_int_inc_ratio_best_lambda.2

# Predict values for train dataset:
la.model_all_int_inc_ratio.train_pred.1 = predict(la.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio.1, s=la.all_int_inc_ratio_best_lambda.1)
ri.model_all_int_inc_ratio.train_pred.1 = predict(ri.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio.1, s=ri.all_int_inc_ratio_best_lambda.1)

la.model_all_int_inc_ratio.train_pred.2 = predict(la.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio.2, s=la.all_int_inc_ratio_best_lambda.2)
ri.model_all_int_inc_ratio.train_pred.2 = predict(ri.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio.2, s=ri.all_int_inc_ratio_best_lambda.2)

# Prediction accuracy with train data:
cor(la.model_all_int_inc_ratio.train_pred.1, dat$C0, method="pearson")
# 0.7020738
cor(ri.model_all_int_inc_ratio.train_pred.1, dat$C0, method="pearson")
# 0.7059261
cor(la.model_all_int_inc_ratio.train_pred.2, dat$C0, method="pearson")
# 0.7215706
cor(ri.model_all_int_inc_ratio.train_pred.2, dat$C0, method="pearson")
# 0.8939136

# Prepare test data:
X_all_int_inc_ratio_test.1 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio.1))

X_all_int_inc_ratio_test.2 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio.2))

x_all_int_inc_ratio_test.1 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio_test.1)
x_all_int_inc_ratio_test.2 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio_test.2)

# Predict values for test dataset:
la.model_all_int_inc_ratio.pred.1 = predict(la.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio_test.1, s=la.all_int_inc_ratio_best_lambda.1)
ri.model_all_int_inc_ratio.pred.1 = predict(ri.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio_test.1, s=ri.all_int_inc_ratio_best_lambda.1)

la.model_all_int_inc_ratio.pred.2 = predict(la.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio_test.2, s=la.all_int_inc_ratio_best_lambda.2)
ri.model_all_int_inc_ratio.pred.2 = predict(ri.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio_test.2, s=ri.all_int_inc_ratio_best_lambda.2)

# Prediction accuracy with test data:
cor(la.model_all_int_inc_ratio.pred.1, dat_test$C0, method="pearson")
# 0.6498516
cor(ri.model_all_int_inc_ratio.pred.1, dat_test$C0, method="pearson")
# 0.6499146
cor(la.model_all_int_inc_ratio.pred.2, dat_test$C0, method="pearson")
# 0.5746391
cor(ri.model_all_int_inc_ratio.pred.2, dat_test$C0, method="pearson")
# 0.5138724










# Set up models with larger feature sets:
model_cols_all_int_inc_ratio_all_ps.1 <- c(ps1, ps3, ps4, ps5, nucleotides, dinucleotides, 
                                           trinucleotides, "ratio_score_2nd")
model_cols_all_int_inc_ratio_all_ps.2 <- c(ps1, ps2, ps3, ps5, nucleotides, dinucleotides, 
                                           trinucleotides, "ratio_score_2nd")

X_all_int_inc_ratio_all_ps.1 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.1))
X_all_int_inc_ratio_all_ps.2 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.2))


f_all_int_inc_ratio_all_ps.1 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                                  ")^2 + (. - ", paste(ps1, collapse = " - "), ")"))
x_all_int_inc_ratio_all_ps.1 = sparse.model.matrix(f_all_int_inc_ratio_all_ps.1, 
                                                   data=X_all_int_inc_ratio_all_ps.1)

y_all_int_inc_ratio_all_ps = dat$C0

# Train models:
set.seed(50)
la.model_all_int_inc_ratio_all_ps.1 = cv.glmnet(x_all_int_inc_ratio_all_ps.1, y_all_int_inc_ratio_all_ps,
                                                alpha=1, family="gaussian")
ri.model_all_int_inc_ratio_all_ps.1 = cv.glmnet(x_all_int_inc_ratio_all_ps.1, y_all_int_inc_ratio_all_ps,
                                                alpha=0, family="gaussian")

plot(la.model_all_int_inc_ratio_all_ps.1)
plot(ri.model_all_int_inc_ratio_all_ps.1)

la.all_int_inc_ratio_all_ps_best_lambda.1 = la.model_all_int_inc_ratio_all_ps.1$lambda.min
la.all_int_inc_ratio_all_ps_best_lambda.1
ri.all_int_inc_ratio_all_ps_best_lambda.1 = ri.model_all_int_inc_ratio_all_ps.1$lambda.min
ri.all_int_inc_ratio_all_ps_best_lambda.1

# Predict values for train dataset:
la.model_all_int_inc_ratio_all_ps.train_pred.1 = predict(la.model_all_int_inc_ratio_all_ps.1, 
                                                         newx = x_all_int_inc_ratio_all_ps.1, 
                                                         s=la.all_int_inc_ratio_all_ps_best_lambda.1)
ri.model_all_int_inc_ratio_all_ps.train_pred.1 = predict(ri.model_all_int_inc_ratio_all_ps.1,
                                                         newx = x_all_int_inc_ratio_all_ps.1, 
                                                         s=ri.all_int_inc_ratio_all_ps_best_lambda.1)

# Prediction accuracy with train data:
cor(la.model_all_int_inc_ratio_all_ps.train_pred.1, dat$C0, method="pearson")
# 0.692659
cor(ri.model_all_int_inc_ratio_all_ps.train_pred.1, dat$C0, method="pearson")
# 0.7009625

# Prepare test data:
X_all_int_inc_ratio_all_ps_test.1 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.1))

x_all_int_inc_ratio_all_ps_test.1 = sparse.model.matrix(f_all_int_inc_ratio_all_ps.1, data=X_all_int_inc_ratio_all_ps_test.1)

# Predict values for test dataset:
la.model_all_int_inc_ratio_all_ps.pred.1 = predict(la.model_all_int_inc_ratio_all_ps.1, 
                                                   newx = x_all_int_inc_ratio_all_ps_test.1,
                                                   s=la.all_int_inc_ratio_all_ps_best_lambda.1)
ri.model_all_int_inc_ratio_all_ps.pred.1 = predict(ri.model_all_int_inc_ratio_all_ps.1,
                                                   newx = x_all_int_inc_ratio_all_ps_test.1,
                                                   s=ri.all_int_inc_ratio_all_ps_best_lambda.1)

# Prediction accuracy with test data:
cor(la.model_all_int_inc_ratio_all_ps.pred.1, dat_test$C0, method="pearson")
# 0.5943428
cor(ri.model_all_int_inc_ratio_all_ps.pred.1, dat_test$C0, method="pearson")
# 0.4854617










# Set up model with binary A or T variables:
model_cols_all_int_inc_ratio_comb.1 <- c(ps1, "ratio_score_2nd", ps1_AT)

X_all_int_inc_ratio_comb.1 = dat_comb %>%
  select(all_of(model_cols_all_int_inc_ratio_comb.1))

f_all_int_inc_ratio_comb.1 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                                ")^2 + (", paste(ps1_AT, collapse = " + "),
                                                ")^2 + ratio_score_2nd"))
x_all_int_inc_ratio_comb.1 = sparse.model.matrix(f_all_int_inc_ratio_comb.1, data=X_all_int_inc_ratio_comb.1)

y_all_int_inc_ratio_comb = dat$C0

# Train models:
set.seed(50)
la.model_all_int_inc_ratio_comb.1 = cv.glmnet(x_all_int_inc_ratio_comb.1, y_all_int_inc_ratio_comb, alpha=1, family="gaussian")
ri.model_all_int_inc_ratio_comb.1 = cv.glmnet(x_all_int_inc_ratio_comb.1, y_all_int_inc_ratio_comb, alpha=0, family="gaussian")

plot(la.model_all_int_inc_ratio_comb.1)
plot(ri.model_all_int_inc_ratio_comb.1)

la.all_int_inc_ratio_comb_best_lambda.1 = la.model_all_int_inc_ratio_comb.1$lambda.min
la.all_int_inc_ratio_comb_best_lambda.1
ri.all_int_inc_ratio_comb_best_lambda.1 = ri.model_all_int_inc_ratio_comb.1$lambda.min
ri.all_int_inc_ratio_comb_best_lambda.1

# Predict values for train dataset:
la.model_all_int_inc_ratio_comb.train_pred.1 = predict(la.model_all_int_inc_ratio_comb.1, 
                                                       newx = x_all_int_inc_ratio_comb.1, 
                                                       s=la.all_int_inc_ratio_comb_best_lambda.1)
ri.model_all_int_inc_ratio_comb.train_pred.1 = predict(ri.model_all_int_inc_ratio_comb.1, 
                                                       newx = x_all_int_inc_ratio_comb.1, 
                                                       s=ri.all_int_inc_ratio_comb_best_lambda.1)

# Prediction accuracy with train data:
cor(la.model_all_int_inc_ratio_comb.train_pred.1, dat$C0, method="pearson")
# 0.6934904
cor(ri.model_all_int_inc_ratio_comb.train_pred.1, dat$C0, method="pearson")
# 0.706607

# Prepare test data:
X_all_int_inc_ratio_comb_test.1 = dat_test_comb %>%
  select(all_of(model_cols_all_int_inc_ratio_comb.1))

x_all_int_inc_ratio_comb_test.1 = sparse.model.matrix(f_all_int_inc_ratio_comb.1, data=X_all_int_inc_ratio_comb_test.1)

# Predict values for test dataset:
la.model_all_int_inc_ratio_comb.pred.1 = predict(la.model_all_int_inc_ratio_comb.1, 
                                                 newx = x_all_int_inc_ratio_comb_test.1, 
                                                 s=la.all_int_inc_ratio_comb_best_lambda.1)
ri.model_all_int_inc_ratio_comb.pred.1 = predict(ri.model_all_int_inc_ratio_comb.1, 
                                                 newx = x_all_int_inc_ratio_comb_test.1, 
                                                 s=ri.all_int_inc_ratio_comb_best_lambda.1)

# Prediction accuracy with test data:
cor(la.model_all_int_inc_ratio_comb.pred.1, dat_test$C0, method="pearson")
# 0.6541115
cor(ri.model_all_int_inc_ratio_comb.pred.1, dat_test$C0, method="pearson")
# 0.655878


