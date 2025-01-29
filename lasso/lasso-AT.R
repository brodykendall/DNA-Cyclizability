library(tidyverse)
library(glmnet)

dat_AT_no_int = readRDS("data/Created/processed_AT_int.rds")

nucleotides <- c("A", "C", "G", "T")
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")
ps4 <- paste0("X", 1:47, "tetra")
ps5 <- paste0("X", 1:46, "penta")

nc_AT = c("nc1_AorT", "nc1_no_AorT", "nc2_AorT", "nc2_no_AorT", "nc3_AorT", "nc3_no_AorT",
          "nc_1_int9_AT", "nc_1_int9_no_AorT", "nc_1_int10_AT", "nc_1_int10_no_AorT",
          "nc_1_int11_AT", "nc_1_int11_no_AorT")

thermo_feat <-c("tm1", "tm2", "tm3", "tm4", "grna_energy", "grna_scaffold_energy")

model_cols_AT_no_int <- c(ps1, ps2, ps3, ps4, ps5, nc_AT, thermo_feat)

X_AT_no_int = dat_AT_no_int %>%
  select(all_of(model_cols_AT_no_int))

x_AT_no_int = sparse.model.matrix(~., data=X_AT_no_int)
y_AT_no_int = dat_AT_no_int$C0

set.seed(50)

la.model_AT_no_int = cv.glmnet(x_AT_no_int, y_AT_no_int, alpha=1, family="gaussian")
ri.model_AT_no_int = cv.glmnet(x_AT_no_int, y_AT_no_int, alpha=0, family="gaussian")

plot(la.model_AT_no_int)
plot(ri.model_AT_no_int)

la.AT_no_int_best_lambda = la.model_AT_no_int$lambda.min
la.AT_no_int_best_lambda

ri.AT_no_int_best_lambda = ri.model_AT_no_int$lambda.min
ri.AT_no_int_best_lambda

la.model_AT_no_int.train_pred = predict(la.model_AT_no_int, newx = x_AT_no_int, s=la.AT_no_int_best_lambda)
ri.model_AT_no_int.train_pred = predict(ri.model_AT_no_int, newx = x_AT_no_int, s=ri.AT_no_int_best_lambda)

cor(la.model_AT_no_int.train_pred, dat_AT_no_int$C0, method="pearson")
# 0.2358255

cor(ri.model_AT_no_int.train_pred, dat_AT_no_int$C0, method="pearson")
# 

dat_test_AT_no_int = readRDS("data/Created/processed_AT_int_test.rds")

X_AT_no_int_test = dat_test_AT_no_int %>%
  select(all_of(model_cols_AT_no_int))

x_AT_no_int_test = sparse.model.matrix(~., data=X_AT_no_int_test)

la.model_AT_no_int.pred = predict(la.model_AT_no_int, newx = x_AT_no_int_test, s=la.AT_no_int_best_lambda)
ri.model_AT_no_int.pred = predict(ri.model_AT_no_int, newx = x_AT_no_int_test, s=ri.AT_no_int_best_lambda)

cor(la.model_AT_no_int.pred, dat_test_AT_no_int$C0, method="pearson")
# 0.2288743

cor(ri.model_AT_no_int.pred, dat_test_AT_no_int$C0, method="pearson")
# 

saveRDS(la.model_AT_no_int, "model/lasso-AT-no-int")
saveRDS(ri.model_AT_no_int, "model/ridge-AT-no-int")





model_cols_AT_all_int.1 <- c(ps1)

X_AT_all_int.1 = dat_AT_no_int %>%
  select(all_of(model_cols_AT_all_int.1))

x_AT_all_int.1 = sparse.model.matrix(~.^2, data=X_AT_all_int.1)
y_AT_all_int = dat_AT_no_int$C0


model_cols_AT_all_int.2 <- c(ps2)

X_AT_all_int.2 = dat_AT_no_int %>%
  select(all_of(model_cols_AT_all_int.2))

x_AT_all_int.2 = sparse.model.matrix(~.^2, data=X_AT_all_int.2)

x_AT_all_int.12.dimnames<-list(union(rownames(x_AT_all_int.1),rownames(x_AT_all_int.2)),c(colnames(x_AT_all_int.1),colnames(x_AT_all_int.2)))
x_AT_all_int.12<- Matrix(0,nrow=length(x_AT_all_int.12.dimnames[[1]]),ncol=length(x_AT_all_int.12.dimnames[[2]]),dimnames=x_AT_all_int.12.dimnames)
x_AT_all_int.12[rownames(x_AT_all_int.1),colnames(x_AT_all_int.1)]<-x_AT_all_int.1
x_AT_all_int.12[rownames(x_AT_all_int.2),colnames(x_AT_all_int.2)]<-x_AT_all_int.2



x_AT_all_int.3 = sparse.model.matrix(~.^3, data=X_AT_all_int.1)


set.seed(50)
la.model_AT_all_int.1 = cv.glmnet(x_AT_all_int.1, y_AT_all_int, alpha=1, family="gaussian")
ri.model_AT_all_int.1 = cv.glmnet(x_AT_all_int.1, y_AT_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_AT_all_int.2 = cv.glmnet(x_AT_all_int.2, y_AT_all_int, alpha=1, family="gaussian")
ri.model_AT_all_int.2 = cv.glmnet(x_AT_all_int.2, y_AT_all_int, alpha=0, family="gaussian")

set.seed(50)
la.model_AT_all_int.12 = cv.glmnet(x_AT_all_int.12, y_AT_all_int, alpha=1, family="gaussian")
ri.model_AT_all_int.12 = cv.glmnet(x_AT_all_int.12, y_AT_all_int, alpha=0, family="gaussian")

# la.model_AT_all_int.3 = cv.glmnet(x_AT_all_int.3, y_AT_all_int, alpha=1, family="gaussian")
# ri.model_AT_all_int.3 = cv.glmnet(x_AT_all_int.3, y_AT_all_int, alpha=0, family="gaussian")

plot(la.model_AT_all_int.1)
plot(ri.model_AT_all_int.1)

plot(la.model_AT_all_int.2)
plot(ri.model_AT_all_int.2)

plot(la.model_AT_all_int.12)
plot(ri.model_AT_all_int.12)

# plot(la.model_AT_all_int.3)
# plot(ri.model_AT_all_int.3)

la.AT_all_int_best_lambda.1 = la.model_AT_all_int.1$lambda.min
la.AT_all_int_best_lambda.1

ri.AT_all_int_best_lambda.1 = ri.model_AT_all_int.1$lambda.min
ri.AT_all_int_best_lambda.1

la.AT_all_int_best_lambda.2 = la.model_AT_all_int.2$lambda.min
la.AT_all_int_best_lambda.2

ri.AT_all_int_best_lambda.2 = ri.model_AT_all_int.2$lambda.min
ri.AT_all_int_best_lambda.2

la.AT_all_int_best_lambda.12 = la.model_AT_all_int.12$lambda.min
la.AT_all_int_best_lambda.12

ri.AT_all_int_best_lambda.12 = ri.model_AT_all_int.12$lambda.min
ri.AT_all_int_best_lambda.12

# la.AT_all_int_best_lambda.3 = la.model_AT_all_int.3$lambda.min
# la.AT_all_int_best_lambda.3
# 
# ri.AT_all_int_best_lambda.3 = ri.model_AT_all_int.3$lambda.min
# ri.AT_all_int_best_lambda.3

la.model_AT_all_int.train_pred.1 = predict(la.model_AT_all_int.1, newx = x_AT_all_int.1, s=la.AT_all_int_best_lambda.1)
ri.model_AT_all_int.train_pred.1 = predict(ri.model_AT_all_int.1, newx = x_AT_all_int.1, s=ri.AT_all_int_best_lambda.1)

la.model_AT_all_int.train_pred.2 = predict(la.model_AT_all_int.2, newx = x_AT_all_int.2, s=la.AT_all_int_best_lambda.2)
ri.model_AT_all_int.train_pred.2 = predict(ri.model_AT_all_int.2, newx = x_AT_all_int.2, s=ri.AT_all_int_best_lambda.2)

la.model_AT_all_int.train_pred.12 = predict(la.model_AT_all_int.12, newx = x_AT_all_int.12, s=la.AT_all_int_best_lambda.12)
ri.model_AT_all_int.train_pred.12 = predict(ri.model_AT_all_int.12, newx = x_AT_all_int.12, s=ri.AT_all_int_best_lambda.12)

cor(la.model_AT_all_int.train_pred.1, dat_AT_no_int$C0, method="pearson")
# 0.5209508

cor(ri.model_AT_all_int.train_pred.1, dat_AT_no_int$C0, method="pearson")
# 

cor(la.model_AT_all_int.train_pred.2, dat_AT_no_int$C0, method="pearson")
# 0.2831633

cor(ri.model_AT_all_int.train_pred.2, dat_AT_no_int$C0, method="pearson")
# 

cor(la.model_AT_all_int.train_pred.12, dat_AT_no_int$C0, method="pearson")
# 

cor(ri.model_AT_all_int.train_pred.12, dat_AT_no_int$C0, method="pearson")
# 

dat_test_AT_all_int = readRDS("data/Created/processed_AT_int_test.rds")

X_AT_all_int_test.1 = dat_test_AT_all_int %>%
  select(all_of(model_cols_AT_all_int.1))

X_AT_all_int_test.2 = dat_test_AT_all_int %>%
  select(all_of(model_cols_AT_all_int.2))

x_AT_all_int_test.1 = sparse.model.matrix(~.^2, data=X_AT_all_int_test.1)

x_AT_all_int_test.2 = sparse.model.matrix(~.^2, data=X_AT_all_int_test.2)

x_AT_all_int_test.12.dimnames<-list(union(rownames(x_AT_all_int_test.1),rownames(x_AT_all_int_test.2)),
                                 c(colnames(x_AT_all_int_test.1),colnames(x_AT_all_int_test.2)))
x_AT_all_int_test.12<- Matrix(0,nrow=length(x_AT_all_int_test.12.dimnames[[1]]),
                           ncol=length(x_AT_all_int_test.12.dimnames[[2]]),
                           dimnames=x_AT_all_int_test.12.dimnames)
x_AT_all_int_test.12[rownames(x_AT_all_int_test.1),colnames(x_AT_all_int_test.1)]<-x_AT_all_int_test.1
x_AT_all_int_test.12[rownames(x_AT_all_int_test.2),colnames(x_AT_all_int_test.2)]<-x_AT_all_int_test.2

X_AT_all_int_test.3 = dat_test_AT_all_int %>%
  select(all_of(model_cols_AT_all_int.3))

la.model_AT_all_int.pred.1 = predict(la.model_AT_all_int.1, newx = x_AT_all_int_test.1, s=la.AT_all_int_best_lambda.1)
ri.model_AT_all_int.pred.1 = predict(ri.model_AT_all_int.1, newx = x_AT_all_int_test.1, s=ri.AT_all_int_best_lambda.1)

la.model_AT_all_int.pred.2 = predict(la.model_AT_all_int.2, newx = x_AT_all_int_test.2, s=la.AT_all_int_best_lambda.2)
ri.model_AT_all_int.pred.2 = predict(ri.model_AT_all_int.2, newx = x_AT_all_int_test.2, s=ri.AT_all_int_best_lambda.2)

la.model_AT_all_int.pred.12 = predict(la.model_AT_all_int.12, newx = x_AT_all_int_test.12, s=la.AT_all_int_best_lambda.12)
ri.model_AT_all_int.pred.12 = predict(ri.model_AT_all_int.12, newx = x_AT_all_int_test.12, s=ri.AT_all_int_best_lambda.12)

# la.model_AT_all_int.pred.3 = predict(la.model_AT_all_int.3, newx = X_AT_all_int_test.3, s=la.AT_all_int_best_lambda.3)
# ri.model_AT_all_int.pred.3 = predict(ri.model_AT_all_int.3, newx = x_AT_all_int_test.3, s=ri.AT_all_int_best_lambda.3)

cor(la.model_AT_all_int.pred.1, dat_test_AT_all_int$C0, method="pearson")
# 0.4911253

cor(ri.model_AT_all_int.pred.1, dat_test_AT_all_int$C0, method="pearson")
# 

cor(la.model_AT_all_int.pred.2, dat_test_AT_all_int$C0, method="pearson")
# 0.225942

cor(ri.model_AT_all_int.pred.2, dat_test_AT_all_int$C0, method="pearson")
# 

cor(la.model_AT_all_int.pred.12, dat_test_AT_all_int$C0, method="pearson")
# 

cor(ri.model_AT_all_int.pred.12, dat_test_AT_all_int$C0, method="pearson")
# 

# cor(la.model_AT_all_int.pred.3, dat_test_AT_all_int$C0, method="pearson")
# # 
# 
# cor(ri.model_AT_all_int.pred.3, dat_test_AT_all_int$C0, method="pearson")
# # 

saveRDS(la.model_AT_all_int.1, "model/lasso-AT-all-int-1")
saveRDS(ri.model_AT_all_int.1, "model/ridge-AT-all-int-1")

saveRDS(la.model_AT_all_int.2, "model/lasso-AT-all-int-2")
saveRDS(ri.model_AT_all_int.2, "model/ridge-AT-all-int-2")

saveRDS(la.model_AT_all_int.12, "model/lasso-AT-all-int-12")
saveRDS(ri.model_AT_all_int.12, "model/ridge-AT-all-int-12")
