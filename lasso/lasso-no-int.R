library(tidyverse)
library(glmnet)

dat_no_int = readRDS("data/Created/processed_105_int.rds")
dat = readRDS("data/Created/processed_ratio.rds")

dat_AT = readRDS("data/Created/processed_AT_int.rds")

dat_AT = dat_AT %>% setNames(paste0(names(.), "_AT"))

dat_comb = bind_cols(dat, dat_AT)

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

la.model_no_int = cv.glmnet(x_no_int, y_no_int, alpha=1, family="gaussian")
ri.model_no_int = cv.glmnet(x_no_int, y_no_int, alpha=0, family="gaussian")

plot(la.model_no_int)
plot(ri.model_no_int)

la.no_int_best_lambda = la.model_no_int$lambda.min
la.no_int_best_lambda

ri.no_int_best_lambda = ri.model_no_int$lambda.min
ri.no_int_best_lambda

la.model_no_int.train_pred = predict(la.model_no_int, newx = x_no_int, s=la.no_int_best_lambda)
ri.model_no_int.train_pred = predict(ri.model_no_int, newx = x_no_int, s=ri.no_int_best_lambda)

cor(la.model_no_int.train_pred, dat_no_int$C0, method="pearson")
# 0.4984689

cor(ri.model_no_int.train_pred, dat_no_int$C0, method="pearson")
# 0.6206618

dat_test_no_int = readRDS("data/Created/processed_105_int_test.rds")

dat_test_no_int = dat_test_no_int %>%
  rename(gc_count_no_int = gc_count_105)

X_no_int_test = dat_test_no_int %>%
  select(all_of(model_cols_no_int))

x_no_int_test = sparse.model.matrix(~., data=X_no_int_test)

la.model_no_int.pred = predict(la.model_no_int, newx = x_no_int_test, s=la.no_int_best_lambda)
ri.model_no_int.pred = predict(ri.model_no_int, newx = x_no_int_test, s=ri.no_int_best_lambda)

cor(la.model_no_int.pred, dat_test_no_int$C0, method="pearson")
# 0.3968694

cor(ri.model_no_int.pred, dat_test_no_int$C0, method="pearson")
# 0.3833119

saveRDS(la.model_no_int, "model/lasso-no-int")
saveRDS(ri.model_no_int, "model/ridge-no-int")



# model_cols_all_int <- c(ps1, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count_no_int")
# model_cols_all_int <- c(ps1, nucleotides, dinucleotides, trinucleotides, thermo_feat, "gc_count_no_int")
model_cols_all_int.1 <- c(ps1)

X_all_int.1 = dat_no_int %>%
  select(all_of(model_cols_all_int.1))

x_all_int.1 = sparse.model.matrix(~.^2, data=X_all_int.1)
y_all_int = dat_no_int$C0


model_cols_all_int.2 <- c(ps2)

X_all_int.2 = dat_no_int %>%
  select(all_of(model_cols_all_int.2))

x_all_int.2 = sparse.model.matrix(~.^2, data=X_all_int.2)

## Option: construct multiple sparse model matrices and combine them using the following:
# m12.dimnames<-list(union(rownames(sparse_m1),rownames(sparse_m2)),c(colnames(sparse_m1),colnames(sparse_m2)))
# m12<- Matrix(0,nrow=length(m12.dimnames[[1]]),ncol=length(m12.dimnames[[2]]),dimnames=m12.dimnames)
# m12[rownames(sparse_m2),colnames(sparse_m2)]<-sparse_m1
# m12[rownames(sparse_m2),colnames(sparse_m2)]<-sparse_m2

x_all_int.12.dimnames<-list(union(rownames(x_all_int.1),rownames(x_all_int.2)),c(colnames(x_all_int.1),colnames(x_all_int.2)))
x_all_int.12<- Matrix(0,nrow=length(x_all_int.12.dimnames[[1]]),ncol=length(x_all_int.12.dimnames[[2]]),dimnames=x_all_int.12.dimnames)
x_all_int.12[rownames(x_all_int.1),colnames(x_all_int.1)]<-x_all_int.1
x_all_int.12[rownames(x_all_int.2),colnames(x_all_int.2)]<-x_all_int.2



x_all_int.3 = sparse.model.matrix(~.^3, data=X_all_int.1)


set.seed(50)
la.model_all_int.1 = cv.glmnet(x_all_int.1, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.1 = cv.glmnet(x_all_int.1, y_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_all_int.2 = cv.glmnet(x_all_int.2, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.2 = cv.glmnet(x_all_int.2, y_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_all_int.12 = cv.glmnet(x_all_int.12, y_all_int, alpha=1, family="gaussian")
ri.model_all_int.12 = cv.glmnet(x_all_int.12, y_all_int, alpha=0, family="gaussian")

# la.model_all_int.3 = cv.glmnet(x_all_int.3, y_all_int, alpha=1, family="gaussian")
# ri.model_all_int.3 = cv.glmnet(x_all_int.3, y_all_int, alpha=0, family="gaussian")

plot(la.model_all_int.1)
plot(ri.model_all_int.1)

plot(la.model_all_int.2)
plot(ri.model_all_int.2)

plot(la.model_all_int.12)
plot(ri.model_all_int.12)

# plot(la.model_all_int.3)
# plot(ri.model_all_int.3)

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

# la.all_int_best_lambda.3 = la.model_all_int.3$lambda.min
# la.all_int_best_lambda.3
# 
# ri.all_int_best_lambda.3 = ri.model_all_int.3$lambda.min
# ri.all_int_best_lambda.3

la.model_all_int.train_pred.1 = predict(la.model_all_int.1, newx = x_all_int.1, s=la.all_int_best_lambda.1)
ri.model_all_int.train_pred.1 = predict(ri.model_all_int.1, newx = x_all_int.1, s=ri.all_int_best_lambda.1)

la.model_all_int.train_pred.2 = predict(la.model_all_int.2, newx = x_all_int.2, s=la.all_int_best_lambda.2)
ri.model_all_int.train_pred.2 = predict(ri.model_all_int.2, newx = x_all_int.2, s=ri.all_int_best_lambda.2)

la.model_all_int.train_pred.12 = predict(la.model_all_int.12, newx = x_all_int.12, s=la.all_int_best_lambda.12)
ri.model_all_int.train_pred.12 = predict(ri.model_all_int.12, newx = x_all_int.12, s=ri.all_int_best_lambda.12)

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

dat_test_all_int = readRDS("data/Created/processed_105_int_test.rds")

dat_test_all_int = dat_test_all_int %>%
  rename(gc_count_no_int = gc_count_105)

X_all_int_test.1 = dat_test_all_int %>%
  select(all_of(model_cols_all_int.1))

X_all_int_test.2 = dat_test_all_int %>%
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

X_all_int_test.3 = dat_test_all_int %>%
  select(all_of(model_cols_all_int.3))

la.model_all_int.pred.1 = predict(la.model_all_int.1, newx = x_all_int_test.1, s=la.all_int_best_lambda.1)
ri.model_all_int.pred.1 = predict(ri.model_all_int.1, newx = x_all_int_test.1, s=ri.all_int_best_lambda.1)

la.model_all_int.pred.2 = predict(la.model_all_int.2, newx = x_all_int_test.2, s=la.all_int_best_lambda.2)
ri.model_all_int.pred.2 = predict(ri.model_all_int.2, newx = x_all_int_test.2, s=ri.all_int_best_lambda.2)

la.model_all_int.pred.12 = predict(la.model_all_int.12, newx = x_all_int_test.12, s=la.all_int_best_lambda.12)
ri.model_all_int.pred.12 = predict(ri.model_all_int.12, newx = x_all_int_test.12, s=ri.all_int_best_lambda.12)

# la.model_all_int.pred.3 = predict(la.model_all_int.3, newx = X_all_int_test.3, s=la.all_int_best_lambda.3)
# ri.model_all_int.pred.3 = predict(ri.model_all_int.3, newx = x_all_int_test.3, s=ri.all_int_best_lambda.3)

cor(la.model_all_int.pred.1, dat_test_all_int$C0, method="pearson")
# 0.5913363 with all pairwise interactions of ps1 (".1")

cor(ri.model_all_int.pred.1, dat_test_all_int$C0, method="pearson")
# 0.5830442 with all pairwise interactions of ps1 (".1")

cor(la.model_all_int.pred.2, dat_test_all_int$C0, method="pearson")
# 0.4583323 with all pairwise interactions of ps2 (".2")

cor(ri.model_all_int.pred.2, dat_test_all_int$C0, method="pearson")
# 0.5013422

cor(la.model_all_int.pred.12, dat_test_all_int$C0, method="pearson")
# 0.4849626

cor(ri.model_all_int.pred.12, dat_test_all_int$C0, method="pearson")
# 0.5166425

# cor(la.model_all_int.pred.3, dat_test_all_int$C0, method="pearson")
# # 
# 
# cor(ri.model_all_int.pred.3, dat_test_all_int$C0, method="pearson")
# # 

saveRDS(la.model_all_int.1, "model/lasso-all-int-1")
saveRDS(ri.model_all_int.1, "model/ridge-all-int-1")

saveRDS(la.model_all_int.2, "model/lasso-all-int-2")
saveRDS(ri.model_all_int.2, "model/ridge-all-int-2")

saveRDS(la.model_all_int.12, "model/lasso-all-int-12")
saveRDS(ri.model_all_int.12, "model/ridge-all-int-12")










dat_test = readRDS("data/Created/processed_test_ratio.rds")

model_cols_all_int_inc_ratio.1 <- c(ps1, "ratio_score_2nd")

model_cols_all_int_inc_ratio.2 <- c(ps2, "ratio_score_2nd")

X_all_int_inc_ratio.1 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio.1))

X_all_int_inc_ratio.2 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio.2))

x_all_int_inc_ratio.1 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio.1)
x_all_int_inc_ratio.2 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio.2)

y_all_int_inc_ratio = dat$C0

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

la.model_all_int_inc_ratio.train_pred.1 = predict(la.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio.1, s=la.all_int_inc_ratio_best_lambda.1)
ri.model_all_int_inc_ratio.train_pred.1 = predict(ri.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio.1, s=ri.all_int_inc_ratio_best_lambda.1)

la.model_all_int_inc_ratio.train_pred.2 = predict(la.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio.2, s=la.all_int_inc_ratio_best_lambda.2)
ri.model_all_int_inc_ratio.train_pred.2 = predict(ri.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio.2, s=ri.all_int_inc_ratio_best_lambda.2)

cor(la.model_all_int_inc_ratio.train_pred.1, dat$C0, method="pearson")
# 0.7020738

cor(ri.model_all_int_inc_ratio.train_pred.1, dat$C0, method="pearson")
# 0.7059261

cor(la.model_all_int_inc_ratio.train_pred.2, dat$C0, method="pearson")
# 0.7215706

cor(ri.model_all_int_inc_ratio.train_pred.2, dat$C0, method="pearson")
# 0.8939136

X_all_int_inc_ratio_test.1 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio.1))

X_all_int_inc_ratio_test.2 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio.2))

x_all_int_inc_ratio_test.1 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio_test.1)
x_all_int_inc_ratio_test.2 = sparse.model.matrix(~(.-ratio_score_2nd)^2 + ratio_score_2nd, data=X_all_int_inc_ratio_test.2)

la.model_all_int_inc_ratio.pred.1 = predict(la.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio_test.1, s=la.all_int_inc_ratio_best_lambda.1)
ri.model_all_int_inc_ratio.pred.1 = predict(ri.model_all_int_inc_ratio.1, newx = x_all_int_inc_ratio_test.1, s=ri.all_int_inc_ratio_best_lambda.1)

la.model_all_int_inc_ratio.pred.2 = predict(la.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio_test.2, s=la.all_int_inc_ratio_best_lambda.2)
ri.model_all_int_inc_ratio.pred.2 = predict(ri.model_all_int_inc_ratio.2, newx = x_all_int_inc_ratio_test.2, s=ri.all_int_inc_ratio_best_lambda.2)

cor(la.model_all_int_inc_ratio.pred.1, dat_test$C0, method="pearson")
# 0.6498516

cor(ri.model_all_int_inc_ratio.pred.1, dat_test$C0, method="pearson")
# 0.6499146

cor(la.model_all_int_inc_ratio.pred.2, dat_test$C0, method="pearson")
# 0.5746391

cor(ri.model_all_int_inc_ratio.pred.2, dat_test$C0, method="pearson")
# 0.5138724


coef(la.model_all_int_inc_ratio.1, la.all_int_inc_ratio_best_lambda.1)[which(coef(la.model_all_int_inc_ratio.1, la.all_int_inc_ratio_best_lambda.1) != 0),]









model_cols_all_int_inc_ratio_all_ps.1 <- c(ps1, ps3, ps4, ps5, nucleotides, dinucleotides, trinucleotides, 
                                           thermo_feat, "gc_count", "ratio_score_2nd")
model_cols_all_int_inc_ratio_all_ps.2 <- c(ps1, ps2, ps3, ps5, nucleotides, dinucleotides, trinucleotides, 
                                           thermo_feat, "gc_count", "ratio_score_2nd")

X_all_int_inc_ratio_all_ps.1 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.1))
X_all_int_inc_ratio_all_ps.2 = dat %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.2))


f_all_int_inc_ratio_all_ps.1 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                                  ")^2 + (. - ", paste(ps1, collapse = " - "), ")"))
f_all_int_inc_ratio_all_ps.2 <- as.formula(paste0("~(", paste(ps2, collapse = " + "), 
                                                  ")^2 + (. - ", paste(ps2, collapse = " - "), ")"))

x_all_int_inc_ratio_all_ps.1 = sparse.model.matrix(f_all_int_inc_ratio_all_ps.1, 
                                                   data=X_all_int_inc_ratio_all_ps.1)

y_all_int_inc_ratio_all_ps = dat$C0

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

la.model_all_int_inc_ratio_all_ps.train_pred.1 = predict(la.model_all_int_inc_ratio_all_ps.1, 
                                                         newx = x_all_int_inc_ratio_all_ps.1, 
                                                         s=la.all_int_inc_ratio_all_ps_best_lambda.1)
ri.model_all_int_inc_ratio_all_ps.train_pred.1 = predict(ri.model_all_int_inc_ratio_all_ps.1,
                                                         newx = x_all_int_inc_ratio_all_ps.1, 
                                                         s=ri.all_int_inc_ratio_all_ps_best_lambda.1)

la.model_all_int_inc_ratio_all_ps.train_pred.2 = predict(la.model_all_int_inc_ratio_all_ps.2,
                                                         newx = x_all_int_inc_ratio_all_ps.2,
                                                         s=la.all_int_inc_ratio_all_ps_best_lambda.2)
ri.model_all_int_inc_ratio_all_ps.train_pred.2 = predict(ri.model_all_int_inc_ratio_all_ps.2,
                                                         newx = x_all_int_inc_ratio_all_ps.2,
                                                         s=ri.all_int_inc_ratio_all_ps_best_lambda.2)

cor(la.model_all_int_inc_ratio_all_ps.train_pred.1, dat$C0, method="pearson")
# 0.692659

cor(ri.model_all_int_inc_ratio_all_ps.train_pred.1, dat$C0, method="pearson")
# 0.7009758

X_all_int_inc_ratio_all_ps_test.1 = dat_test %>%
  select(all_of(model_cols_all_int_inc_ratio_all_ps.1))

x_all_int_inc_ratio_all_ps_test.1 = sparse.model.matrix(f_all_int_inc_ratio_all_ps.1, data=X_all_int_inc_ratio_all_ps_test.1)

la.model_all_int_inc_ratio_all_ps.pred.1 = predict(la.model_all_int_inc_ratio_all_ps.1, 
                                                   newx = x_all_int_inc_ratio_all_ps_test.1,
                                                   s=la.all_int_inc_ratio_all_ps_best_lambda.1)
ri.model_all_int_inc_ratio_all_ps.pred.1 = predict(ri.model_all_int_inc_ratio_all_ps.1,
                                                   newx = x_all_int_inc_ratio_all_ps_test.1,
                                                   s=ri.all_int_inc_ratio_all_ps_best_lambda.1)

cor(la.model_all_int_inc_ratio_all_ps.pred.1, dat_test$C0, method="pearson")
# 0.5943428

cor(ri.model_all_int_inc_ratio_all_ps.pred.1, dat_test$C0, method="pearson")
# 0.4854651







dat_test_AT = readRDS("data/Created/processed_AT_int_test.rds")
dat_test_AT = dat_test_AT %>% setNames(paste0(names(.), "_AT"))
dat_test_comb = bind_cols(dat_test, dat_test_AT)

ps1_AT <- paste0("X", 1:50, "mono_AT")

model_cols_all_int_inc_ratio_comb.1 <- c(ps1, "ratio_score_2nd", ps1_AT)

X_all_int_inc_ratio_comb.1 = dat_comb %>%
  select(all_of(model_cols_all_int_inc_ratio_comb.1))

f_all_int_inc_ratio_comb.1 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                                ")^2 + (", paste(ps1_AT, collapse = " + "),
                                                ")^2 + ratio_score_2nd"))

x_all_int_inc_ratio_comb.1 = sparse.model.matrix(f_all_int_inc_ratio_comb.1, data=X_all_int_inc_ratio_comb.1)

y_all_int_inc_ratio_comb = dat$C0

set.seed(50)
la.model_all_int_inc_ratio_comb.1 = cv.glmnet(x_all_int_inc_ratio_comb.1, y_all_int_inc_ratio_comb, alpha=1, family="gaussian")
ri.model_all_int_inc_ratio_comb.1 = cv.glmnet(x_all_int_inc_ratio_comb.1, y_all_int_inc_ratio_comb, alpha=0, family="gaussian")

plot(la.model_all_int_inc_ratio_comb.1)
plot(ri.model_all_int_inc_ratio_comb.1)

la.all_int_inc_ratio_comb_best_lambda.1 = la.model_all_int_inc_ratio_comb.1$lambda.min
la.all_int_inc_ratio_comb_best_lambda.1

ri.all_int_inc_ratio_comb_best_lambda.1 = ri.model_all_int_inc_ratio_comb.1$lambda.min
ri.all_int_inc_ratio_comb_best_lambda.1

la.model_all_int_inc_ratio_comb.train_pred.1 = predict(la.model_all_int_inc_ratio_comb.1, 
                                                       newx = x_all_int_inc_ratio_comb.1, 
                                                       s=la.all_int_inc_ratio_comb_best_lambda.1)
ri.model_all_int_inc_ratio_comb.train_pred.1 = predict(ri.model_all_int_inc_ratio_comb.1, 
                                                       newx = x_all_int_inc_ratio_comb.1, 
                                                       s=ri.all_int_inc_ratio_comb_best_lambda.1)

cor(la.model_all_int_inc_ratio_comb.train_pred.1, dat$C0, method="pearson")
# 0.6934904

cor(ri.model_all_int_inc_ratio_comb.train_pred.1, dat$C0, method="pearson")
# 0.706607

X_all_int_inc_ratio_comb_test.1 = dat_test_comb %>%
  select(all_of(model_cols_all_int_inc_ratio_comb.1))

x_all_int_inc_ratio_comb_test.1 = sparse.model.matrix(f_all_int_inc_ratio_comb.1, data=X_all_int_inc_ratio_comb_test.1)

la.model_all_int_inc_ratio_comb.pred.1 = predict(la.model_all_int_inc_ratio_comb.1, 
                                                 newx = x_all_int_inc_ratio_comb_test.1, 
                                                 s=la.all_int_inc_ratio_comb_best_lambda.1)
ri.model_all_int_inc_ratio_comb.pred.1 = predict(ri.model_all_int_inc_ratio_comb.1, 
                                                 newx = x_all_int_inc_ratio_comb_test.1, 
                                                 s=ri.all_int_inc_ratio_comb_best_lambda.1)

cor(la.model_all_int_inc_ratio_comb.pred.1, dat_test$C0, method="pearson")
# 0.6541115

cor(ri.model_all_int_inc_ratio_comb.pred.1, dat_test$C0, method="pearson")
# 0.655878

coef(la.model_all_int_inc_ratio_comb.1, la.all_int_inc_ratio_comb_best_lambda.1)[which(coef(la.model_all_int_inc_ratio_comb.1, la.all_int_inc_ratio_comb_best_lambda.1) != 0),]










# Add 1-D principal components:
dat_pca_ps1_only = data.frame(x_no_int_ps1_only.pca$x)
dat_pca_ps1_only = dat_pca_ps1_only %>% setNames(paste0(names(.), "_ps1"))
dat_pca_ps2_only = data.frame(x_no_int_ps2_only.pca$x)
dat_pca_ps2_only = dat_pca_ps2_only %>% setNames(paste0(names(.), "_ps2"))
dat_pca = dat %>%
  bind_cols(dat_pca_ps1_only,
            dat_pca_ps2_only)

dat_test_ps1_only = dat_test %>%
  select(all_of(ps1))

x_test_no_int_ps1_only <- sparse.model.matrix(~ ., data = dat_test_ps1_only,
                                         contrasts.arg = lapply(dat_test_ps1_only[, sapply(dat_test_ps1_only, 
                                                                                           is.factor), 
                                                                                  drop = FALSE], 
                                                                contrasts, contrasts = FALSE))

dat_test_ps2_only = dat_test %>%
  select(all_of(ps2))

x_test_no_int_ps2_only <- sparse.model.matrix(~ ., data = dat_test_ps2_only,
                                              contrasts.arg = lapply(dat_test_ps2_only[, sapply(dat_test_ps2_only, 
                                                                                                is.factor), 
                                                                                       drop = FALSE], 
                                                                     contrasts, contrasts = FALSE))

dat_test_pca_ps1_only = x_test_no_int_ps1_only %*% x_no_int_ps1_only.pca$rotation
dat_test_pca_ps1_only = data.frame(matrix(dat_test_pca_ps1_only@x, nrow = nrow(x_test_no_int_ps1_only), 
                                          ncol = ncol(x_test_no_int_ps1_only)))
colnames(dat_test_pca_ps1_only) = paste0("PC", 1:ncol(x_test_no_int_ps1_only), "_ps1")
dat_test_pca_ps2_only = x_test_no_int_ps2_only %*% x_no_int_ps2_only.pca$rotation
dat_test_pca_ps2_only = data.frame(matrix(dat_test_pca_ps2_only@x, nrow = nrow(x_test_no_int_ps2_only), 
                                          ncol = ncol(x_test_no_int_ps2_only)))
colnames(dat_test_pca_ps2_only) = paste0("PC", 1:ncol(x_test_no_int_ps2_only), "_ps2")

dat_test_pca = dat_test %>%
  bind_cols(dat_test_pca_ps1_only,
            dat_test_pca_ps2_only)


PC_ps2.1 <- paste0("PC", 4:10, "_ps2")
PC_ps2.2 <- paste0("PC", 1:200, "_ps2")



model_cols_pca_ps2.1 <- c(ps1, "ratio_score_2nd", PC_ps2.1)
model_cols_pca_ps2.2 <- c(ps1, "ratio_score_2nd", PC_ps2.2)


X_pca_ps2.1 = dat_pca %>%
  select(all_of(model_cols_pca_ps2.1))
X_pca_ps2.2 = dat_pca %>%
  select(all_of(model_cols_pca_ps2.2))

f_pca_ps2.1 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                 ")^2 + ", paste(PC_ps2.1, collapse = " + "),
                                 " + ratio_score_2nd"))
f_pca_ps2.2 <- as.formula(paste0("~(", paste(ps1, collapse = " + "), 
                                 ")^2 + ", paste(PC_ps2.2, collapse = " + "),
                                 " + ratio_score_2nd"))

x_pca_ps2.1 = sparse.model.matrix(f_pca_ps2.1, data=X_pca_ps2.1)
x_pca_ps2.2 = sparse.model.matrix(f_pca_ps2.2, data=X_pca_ps2.2)

y_pca_ps2 = dat$C0

set.seed(50)
la.model_pca_ps2.1 = cv.glmnet(x_pca_ps2.1, y_pca_ps2, alpha=1, family="gaussian")
ri.model_pca_ps2.1 = cv.glmnet(x_pca_ps2.1, y_pca_ps2, alpha=0, family="gaussian")

set.seed(50)
la.model_pca_ps2.2 = cv.glmnet(x_pca_ps2.2, y_pca_ps2, alpha=1, family="gaussian")
ri.model_pca_ps2.2 = cv.glmnet(x_pca_ps2.2, y_pca_ps2, alpha=0, family="gaussian")

plot(la.model_pca_ps2.1)
plot(ri.model_pca_ps2.1)

plot(la.model_pca_ps2.2)
plot(ri.model_pca_ps2.2)

la.pca_ps2_best_lambda.1 = la.model_pca_ps2.1$lambda.min
la.pca_ps2_best_lambda.1

ri.pca_ps2_best_lambda.1 = ri.model_pca_ps2.1$lambda.min
ri.pca_ps2_best_lambda.1

la.pca_ps2_best_lambda.2 = la.model_pca_ps2.2$lambda.min
la.pca_ps2_best_lambda.2

ri.pca_ps2_best_lambda.2 = ri.model_pca_ps2.2$lambda.min
ri.pca_ps2_best_lambda.2

la.model_pca_ps2.train_pred.1 = predict(la.model_pca_ps2.1, 
                                        newx = x_pca_ps2.1, 
                                        s=la.pca_ps2_best_lambda.1)
ri.model_pca_ps2.train_pred.1 = predict(ri.model_pca_ps2.1, 
                                        newx = x_pca_ps2.1, 
                                        s=ri.pca_ps2_best_lambda.1)

la.model_pca_ps2.train_pred.2 = predict(la.model_pca_ps2.2, 
                                        newx = x_pca_ps2.2, 
                                        s=la.pca_ps2_best_lambda.2)
ri.model_pca_ps2.train_pred.2 = predict(ri.model_pca_ps2.2, 
                                        newx = x_pca_ps2.2, 
                                        s=ri.pca_ps2_best_lambda.2)

cor(la.model_pca_ps2.train_pred.1, dat$C0, method="pearson")
# 0.7021148

cor(ri.model_pca_ps2.train_pred.1, dat$C0, method="pearson")
# 0.7060015

cor(la.model_pca_ps2.train_pred.2, dat$C0, method="pearson")
# 0.7026752

cor(ri.model_pca_ps2.train_pred.2, dat$C0, method="pearson")
# 0.7066613

X_pca_ps2_test.1 = dat_test_pca %>%
  select(all_of(model_cols_pca_ps2.1))
X_pca_ps2_test.2 = dat_test_pca %>%
  select(all_of(model_cols_pca_ps2.2))

x_pca_ps2_test.1 = sparse.model.matrix(f_pca_ps2.1, data=X_pca_ps2_test.1)
x_pca_ps2_test.2 = sparse.model.matrix(f_pca_ps2.2, data=X_pca_ps2_test.2)


la.model_pca_ps2.pred.1 = predict(la.model_pca_ps2.1, 
                                  newx = x_pca_ps2_test.1, 
                                  s=la.pca_ps2_best_lambda.1)
ri.model_pca_ps2.pred.1 = predict(ri.model_pca_ps2.1, 
                                  newx = x_pca_ps2_test.1, 
                                  s=ri.pca_ps2_best_lambda.1)

la.model_pca_ps2.pred.2 = predict(la.model_pca_ps2.2, 
                                  newx = x_pca_ps2_test.2, 
                                  s=la.pca_ps2_best_lambda.2)
ri.model_pca_ps2.pred.2 = predict(ri.model_pca_ps2.2, 
                                  newx = x_pca_ps2_test.2, 
                                  s=ri.pca_ps2_best_lambda.2)

cor(la.model_pca_ps2.pred.1, dat_test$C0, method="pearson")
# 0.6500292

cor(ri.model_pca_ps2.pred.1, dat_test$C0, method="pearson")
# 0.6500057

cor(la.model_pca_ps2.pred.2, dat_test$C0, method="pearson")
# 0.6506983

cor(ri.model_pca_ps2.pred.2, dat_test$C0, method="pearson")
# 0.6506983

coef(la.model_pca_ps2.1, la.pca_ps2_best_lambda.1)[which(coef(la.model_pca_ps2.1, la.pca_ps2_best_lambda.1) != 0),]
coef(la.model_pca_ps2.2, la.pca_ps2_best_lambda.2)[which(coef(la.model_pca_ps2.2, la.pca_ps2_best_lambda.2) != 0),]

