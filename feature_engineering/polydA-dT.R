library(tidyverse)
library(stringi)
library(mgcv)

dat = readRDS("data/Created/processed_tiling_newC0.rds")
y = dat$C0_new
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
y_test = dat_test$C0_new
dat_tiling_all = rbind(dat, dat_test)
y_tiling_all = c(y, y_test)
dat_random = readRDS("data/Created/processed_random_newC0.rds")
y_random = dat_random$C0_new
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")
y_random_test = dat_random_test$C0_new
dat_random_all = rbind(dat_random, dat_random_test)
y_random_all = c(y_random, y_random_test)
dat_chrV = readRDS("data/Created/processed_chrV_newC0.rds")
y_chrV = dat_chrV$C0_new
dat_chrV_test = readRDS("data/Created/processed_chrV_test_newC0.rds")
y_chrV_test = dat_chrV_test$C0_new
dat_chrV_all = rbind(dat_chrV, dat_chrV_test)
y_chrV_all = c(y_chrV, y_chrV_test)
dat_yeast = readRDS("data/Created/processed_yeast_newC0.rds")
y_yeast = dat_yeast$C0_new
dat_yeast_test = readRDS("data/Created/processed_yeast_test_newC0.rds")
y_yeast_test = dat_yeast_test$C0_new
dat_yeast_all = rbind(dat_yeast, dat_yeast_test)
y_yeast_all = c(y_yeast, y_yeast_test)



# Regex for 3 or more in a row without double counting 6, 9, etc.: [AT]{3}(?=[^AT])
# Regex for exactly 3 in a row: (^|[^AT])[AT]{3}(?=[^AT])

# lengths = 4:19
lengths = 4:10

npolyAorT = map(lengths,
                ~stri_count(dat$x50mer,
                            regex=paste0("(^|[^AT])[AT]{", .x, "}(?=[^AT])"))) %>%
  bind_cols(stri_count(dat$x50mer,
                       regex=paste0("(^|[^AT])[AT]{20,}(?=[^AT])")))

colnames(npolyAorT) = paste0("npolyAorT_", c(lengths, paste0(max(lengths)+1, "plus")))

npolyAorT_test = map(lengths,
                     ~stri_count(dat_test$x50mer,
                                 regex=paste0("(^|[^AT])[AT]{", .x, "}(?=[^AT])"))) %>%
  bind_cols(stri_count(dat_test$x50mer,
                       regex=paste0("(^|[^AT])[AT]{20,}(?=[^AT])")))

colnames(npolyAorT_test) = paste0("npolyAorT_", c(lengths, paste0(max(lengths)+1, "plus")))



npolyAorT_random_all = map(lengths,
                           ~stri_count(dat_random_all$x50mer,
                                       regex=paste0("(^|[^AT])[AT]{", .x, "}(?=[^AT])"))) %>%
  bind_cols(stri_count(dat_random_all$x50mer,
                       regex=paste0("(^|[^AT])[AT]{20,}(?=[^AT])")))

colnames(npolyAorT_random_all) = paste0("npolyAorT_", c(lengths, paste0(max(lengths)+1, "plus")))



dat_npolyAorT = cbind(npolyAorT, dat %>% select(C0_new))

npolyAorT_lm = lm(C0_new~., dat_npolyAorT)
cor(npolyAorT_lm$fitted.values, dat$C0_new)
# 0.1772687

dat_npolyAorT_test = npolyAorT_test
npolyAorT_lm_pred = predict(npolyAorT_lm, dat_npolyAorT_test)
cor(npolyAorT_lm_pred, dat_test$C0_new)
# 0.18956

dat_npolyAorT_random_all = npolyAorT_random_all
npolyAorT_lm_pred_random = predict(npolyAorT_lm, dat_npolyAorT_random_all)
cor(npolyAorT_lm_pred_random, dat_random_all$C0_new)
# 0.1268386





npolyA = map(lengths, 
             ~stri_count(dat$x50mer,
                         regex=paste0("(^|[^A])[A]{", .x, "}(?=[^A])"))) %>%
  bind_cols(stri_count(dat$x50mer,
                       regex=paste0("(^|[^A])[A]{20,}(?=[^A])")))

colnames(npolyA) = paste0("npolyA_", c(lengths, paste0(max(lengths)+1, "plus")))

npolyT = map(lengths, 
             ~stri_count(dat$x50mer,
                         regex=paste0("(^|[^T])[T]{", .x, "}(?=[^T])"))) %>%
  bind_cols(stri_count(dat$x50mer,
                       regex=paste0("(^|[^T])[T]{20,}(?=[^T])")))

colnames(npolyT) = paste0("npolyT_", c(lengths, paste0(max(lengths)+1, "plus")))


npolyA_test = map(lengths, 
                  ~stri_count(dat_test$x50mer,
                              regex=paste0("(^|[^A])[A]{", .x, "}(?=[^A])"))) %>%
  bind_cols(stri_count(dat_test$x50mer,
                       regex=paste0("(^|[^A])[A]{20,}(?=[^A])")))

colnames(npolyA_test) = paste0("npolyA_", c(lengths, paste0(max(lengths)+1, "plus")))

npolyT_test = map(lengths, 
                  ~stri_count(dat_test$x50mer,
                              regex=paste0("(^|[^T])[T]{", .x, "}(?=[^T])"))) %>%
  bind_cols(stri_count(dat_test$x50mer,
                       regex=paste0("(^|[^T])[T]{20,}(?=[^T])")))

colnames(npolyT_test) = paste0("npolyT_", c(lengths, paste0(max(lengths)+1, "plus")))


npolyA_random_all = map(lengths, 
                        ~stri_count(dat_random_all$x50mer,
                                    regex=paste0("(^|[^A])[A]{", .x, "}(?=[^A])"))) %>%
  bind_cols(stri_count(dat_random_all$x50mer,
                       regex=paste0("(^|[^A])[A]{20,}(?=[^A])")))

colnames(npolyA_random_all) = paste0("npolyA_", c(lengths, paste0(max(lengths)+1, "plus")))

npolyT_random_all = map(lengths, 
                        ~stri_count(dat_random_all$x50mer,
                                    regex=paste0("(^|[^T])[T]{", .x, "}(?=[^T])"))) %>%
  bind_cols(stri_count(dat_random_all$x50mer,
                       regex=paste0("(^|[^T])[T]{20,}(?=[^T])")))

colnames(npolyT_random_all) = paste0("npolyT_", c(lengths, paste0(max(lengths)+1, "plus")))


dat_npoly = cbind(npolyA, npolyT, npolyAorT, dat%>%select(C0_new))

npoly_lm = lm(C0_new~., dat_npoly)
cor(npoly_lm$fitted.values, dat$C0_new)
# 0.2566713

dat_npoly_test = cbind(npolyA_test, npolyT_test, npolyAorT_test)
npoly_lm_pred = predict(npoly_lm, dat_npoly_test)
cor(npoly_lm_pred, dat_test$C0_new)
# 0.2582381

dat_npoly_random_all = cbind(npolyA_random_all, npolyT_random_all, npolyAorT_random_all)
npoly_lm_pred_random = predict(npoly_lm, dat_npoly_random_all)
cor(npoly_lm_pred_random, dat_random_all$C0_new)
# 0.1251304






# Inputs: polyA/T/AorT + All Distance Frequency Between AA/TT/AT/TA and CC/GG/CG/GC Dinucleotides

dist_freq_AorTdi_AorTdi_bid2 = readRDS("data/Created/tiling_dist_freq_AorTdi_AorTdi_bid2.rds")
dist_freq_AorTdi_CorGdi_bid2 = readRDS("data/Created/tiling_dist_freq_AorTdi_CorGdi_bid2.rds")
dist_freq_CorGdi_CorGdi_bid2 = readRDS("data/Created/tiling_dist_freq_CorGdi_CorGdi_bid2.rds")
dist_freq_AorTdi_AorTdi_bid2_test = readRDS("data/Created/tiling_dist_freq_AorTdi_AorTdi_bid2_test.rds")
dist_freq_AorTdi_CorGdi_bid2_test = readRDS("data/Created/tiling_dist_freq_AorTdi_CorGdi_bid2_test.rds")
dist_freq_CorGdi_CorGdi_bid2_test = readRDS("data/Created/tiling_dist_freq_CorGdi_CorGdi_bid2_test.rds")
dist_freq_AorTdi_AorTdi_bid2_random_all = readRDS("data/Created/random_all_dist_freq_AorTdi_AorTdi_bid2.rds")
dist_freq_AorTdi_CorGdi_bid2_random_all = readRDS("data/Created/random_all_dist_freq_AorTdi_CorGdi_bid2.rds")
dist_freq_CorGdi_CorGdi_bid2_random_all = readRDS("data/Created/random_all_dist_freq_CorGdi_CorGdi_bid2.rds")


dat_npoly2 = cbind(npolyA, npolyT, npolyAorT, 
                   dist_freq_AorTdi_AorTdi_bid2,
                   dist_freq_AorTdi_CorGdi_bid2,
                   dist_freq_CorGdi_CorGdi_bid2,
                   dat %>% select(C0_new))

npoly2_lm = lm(C0_new~., data=dat_npoly2)
cor(npoly2_lm$fitted.values, y)
# 0.5414124

dat_test_npoly2 = cbind(npolyA_test, npolyT_test, npolyAorT_test,
                        dist_freq_AorTdi_AorTdi_bid2_test,
                        dist_freq_AorTdi_CorGdi_bid2_test,
                        dist_freq_CorGdi_CorGdi_bid2_test,
                        dat_test %>% select(C0_new))

npoly2_pred = predict(npoly2_lm, 
                      dat_test_npoly2)

cor(npoly2_pred, y_test)
# 0.5398178

mean((npoly2_pred - y_test)^2)
# 0.1707931


dat_random_all_npoly2 = cbind(npolyA_random_all, npolyT_random_all, npolyAorT_random_all,
                              dist_freq_AorTdi_AorTdi_bid2_random_all,
                              dist_freq_AorTdi_CorGdi_bid2_random_all,
                              dist_freq_CorGdi_CorGdi_bid2_random_all,
                              dat_random_all %>% select(C0_new))

npoly2_pred_random_all = predict(npoly2_lm, 
                                 dat_random_all_npoly2)

cor(npoly2_pred_random_all, y_random_all)
# 0.4172967

mean((npoly2_pred_random_all - y_random_all)^2)
# 0.1238838








# Inputs: polyA/T/AorT + Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides

dist_freq_di_AorTdi_bid2 = readRDS("data/Created/tiling_dist_freq_di_AorTdi_bid2.rds")
dist_freq_di_CorGdi_bid2 = readRDS("data/Created/tiling_dist_freq_di_CorGdi_bid2.rds")
dist_freq_di_AorTdi_bid2_test = readRDS("data/Created/tiling_dist_freq_di_AorTdi_bid2_test.rds")
dist_freq_di_CorGdi_bid2_test = readRDS("data/Created/tiling_dist_freq_di_CorGdi_bid2_test.rds")

dist_freq_di_AorTdi_bid2_random_all = readRDS("data/Created/random_all_dist_freq_di_AorTdi_bid2.rds")
dist_freq_di_CorGdi_bid2_random_all = readRDS("data/Created/random_all_dist_freq_di_CorGdi_bid2.rds")



dat_npoly3 = cbind(npolyA, npolyT, npolyAorT,
                   dist_freq_di_AorTdi_bid2,
                   dist_freq_di_CorGdi_bid2,
                   dat %>% select(C0_new))

npoly3_lm = lm(C0_new~., data=dat_npoly3)
cor(npoly3_lm$fitted.values, y)
# 0.7407449

dat_test_npoly3 = cbind(npolyA_test, npolyT_test, npolyAorT_test,
                        dist_freq_di_AorTdi_bid2_test,
                        dist_freq_di_CorGdi_bid2_test,
                        dat_test %>% select(C0_new))

npoly3_pred = predict(npoly3_lm, dat_test_npoly3)

cor(npoly3_pred, y_test)
# 0.740785

mean((npoly3_pred - y_test)^2)
# 0.1087789

dat_random_npoly3 = cbind(npolyA_random_all, npolyT_random_all, npolyAorT_random_all,
                          dist_freq_di_AorTdi_bid2_random_all,
                          dist_freq_di_CorGdi_bid2_random_all,
                          dat_random_all %>% select(C0_new))

random_npoly3_pred = predict(npoly3_lm, dat_random_npoly3)

cor(random_npoly3_pred, y_random_all)
# 0.6533972





















# Count the number of times there is a contiguous sequence of A's, T's, and A's or T's of at least 3:
# npolyA = stri_count(dat$x50mer, 
#                     regex="[A]{3,}")
# cor(npolyA, dat$C0_new)
# # -0.04198643
# 
# npolyA_test = stri_count(dat_test$x50mer, 
#                          regex="[A]{3,}")
# 
# npolyT = stri_count(dat$x50mer,
#                     regex="[T]{3,}")
# cor(npolyT, dat$C0_new)
# # -0.04092385
# 
# npolyT_test = stri_count(dat_test$x50mer,
#                          regex="[T]{3,}")
# 
# npolyAorT = stri_count(dat$x50mer,
#                        regex="[AT]{3,}")
# cor(npolyAorT, dat$C0_new)
# # 0.05311018
# 
# npolyAorT_test = stri_count(dat_test$x50mer,
#                             regex="[AT]{3,}")
# 
# dat_npoly = cbind(dat %>% select(C0_new),
#                   npolyA, npolyT, npolyAorT)
# npoly_lm = lm(C0_new~., dat_npoly)
# cor(npoly_lm$fitted.values, dat$C0_new)
# # 0.1201335
# 
# dat_npoly_test = data.frame("npolyA"=npolyA_test, "npolyT"=npolyT_test, "npolyAorT"=npolyAorT_test)
# npoly_lm_pred = predict(npoly_lm, dat_npoly_test)
# cor(npoly_lm_pred, dat_test$C0_new)
# # 0.1321701

# Find the length of the longest sequence of A's, T's, and A's or T's:
max_polyA = map(stri_extract(dat$x50mer,
                             regex="[A]{1,}",
                             mode="all"), 
                ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyA[which(is.na(max_polyA))] = 0
cor(max_polyA, dat$C0_new)
# -0.09055186

max_polyA_test = map(stri_extract(dat_test$x50mer,
                                  regex="[A]{1,}",
                                  mode="all"), 
                     ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyA_test[which(is.na(max_polyA_test))] = 0

max_polyT = map(stri_extract(dat$x50mer,
                             regex="[T]{1,}",
                             mode="all"), 
                ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyT[which(is.na(max_polyT))] = 0
cor(max_polyT, dat$C0_new)
# -0.09000319

max_polyT_test = map(stri_extract(dat_test$x50mer,
                                  regex="[T]{1,}",
                                  mode="all"), 
                     ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyT_test[which(is.na(max_polyT_test))] = 0

max_polyAorT = map(stri_extract(dat$x50mer,
                                regex="[AT]{1,}",
                                mode="all"), 
                   ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyAorT[which(is.na(max_polyAorT))] = 0
cor(max_polyAorT, dat$C0_new)
# -0.05613836

max_polyAorT_test = map(stri_extract(dat_test$x50mer,
                                     regex="[AT]{1,}",
                                     mode="all"), 
                        ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyAorT_test[which(is.na(max_polyAorT_test))] = 0

dat_max_poly = cbind(dat %>% select(C0_new),
                     max_polyA, max_polyT, max_polyAorT)

dat_max_poly_test = cbind(dat_test %>% select(C0_new),
                          max_polyA_test, max_polyT_test, max_polyAorT_test)





# RANDOM:

# Count the number of times there is a contiguous sequence of A's, T's, and A's or T's of at least 3:
# npolyA_random = stri_count(dat_random$x50mer, 
#                            regex="[A]{3,}")
# cor(npolyA_random, dat_random$C0_new)
# # 0.01988172
# 
# npolyA_random_test = stri_count(dat_random_test$x50mer, 
#                                 regex="[A]{3,}")
# 
# npolyT_random = stri_count(dat_random$x50mer,
#                            regex="[T]{3,}")
# cor(npolyT_random, dat_random$C0_new)
# # 0.01674767
# 
# npolyT_random_test = stri_count(dat_random_test$x50mer,
#                                 regex="[T]{3,}")
# 
# npolyAorT_random = stri_count(dat_random$x50mer,
#                               regex="[AT]{3,}")
# cor(npolyAorT_random, dat_random$C0_new)
# # 0.06074041
# 
# npolyAorT_random_test = stri_count(dat_random_test$x50mer,
#                                    regex="[AT]{3,}")
# 
dat_npoly_random_all = cbind(dat_random_all %>% select(C0_new),
                             npolyA_random_all, npolyT_random_all, npolyAorT_random_all)

# Find the length of the longest sequence of A's, T's, and A's or T's:
max_polyA_random_all = map(stri_extract(dat_random_all$x50mer,
                                        regex="[A]{1,}",
                                        mode="all"), 
                           ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyA_random_all[which(is.na(max_polyA_random_all))] = 0
cor(max_polyA_random_all, dat_random_all$C0_new)
# 0.02204305

# max_polyA_random_test = map(stri_extract(dat_random_test$x50mer,
#                                          regex="[A]{1,}",
#                                          mode="all"), 
#                             ~max(sapply(.x, nchar))) %>%
#   unlist()
# max_polyA_random_test[which(is.na(max_polyA_random_test))] = 0

max_polyT_random_all = map(stri_extract(dat_random_all$x50mer,
                                        regex="[T]{1,}",
                                        mode="all"), 
                           ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyT_random_all[which(is.na(max_polyT_random_all))] = 0
cor(max_polyT_random_all, dat_random_all$C0_new)
# 0.00534089

# max_polyT_random_test = map(stri_extract(dat_random_test$x50mer,
#                                          regex="[T]{1,}",
#                                          mode="all"), 
#                             ~max(sapply(.x, nchar))) %>%
#   unlist()
# max_polyT_random_test[which(is.na(max_polyT_random_test))] = 0

max_polyAorT_random_all = map(stri_extract(dat_random_all$x50mer,
                                           regex="[AT]{1,}",
                                           mode="all"), 
                              ~max(sapply(.x, nchar))) %>%
  unlist()
max_polyAorT_random_all[which(is.na(max_polyAorT_random_all))] = 0
cor(max_polyAorT_random_all, dat_random_all$C0_new)
# 0.06000565

# max_polyAorT_random_test = map(stri_extract(dat_random_test$x50mer,
#                                             regex="[AT]{1,}",
#                                             mode="all"), 
#                                ~max(sapply(.x, nchar))) %>%
#   unlist()
# max_polyAorT_random_test[which(is.na(max_polyAorT_random_test))] = 0

dat_max_poly_random_all = cbind(dat_random_all %>% select(C0_new),
                                "max_polyA"=max_polyA_random_all, 
                                "max_polyT"=max_polyT_random_all, 
                                "max_polyAorT"=max_polyAorT_random_all)


dat_npoly = cbind(dat %>% select(C0_new),
                  npolyA, npolyT, npolyAorT)
npoly_lm = lm(C0_new~., dat_npoly)
cor(npoly_lm$fitted.values, dat$C0_new)
# 0.2566713

dat_npoly_test = data.frame(npolyA_test, npolyT_test, npolyAorT_test)
npoly_lm_pred = predict(npoly_lm, dat_npoly_test)
cor(npoly_lm_pred, dat_test$C0_new)
# 0.2582381

# dat_npoly_random_all = rbind(dat_npoly_random %>% select(-C0_new), dat_npoly_random_test)

npoly_lm_pred_random_all = predict(npoly_lm, dat_npoly_random_all)
cor(npoly_lm_pred_random_all, dat_random_all$C0_new)
# 0.1251304

max_poly_lm = lm(C0_new~., dat_max_poly)
cor(max_poly_lm$fitted.values, dat$C0_new)
# 0.13421

dat_max_poly_test = data.frame("max_polyA"=max_polyA_test, "max_polyT"=max_polyT_test, "max_polyAorT"=max_polyAorT_test)
max_poly_lm_pred = predict(max_poly_lm, dat_max_poly_test)
cor(max_poly_lm_pred, dat_test$C0_new)
# 0.1182864

# dat_max_poly_random_test = data.frame("max_polyA"=max_polyA_random_test, "max_polyT"=max_polyT_random_test, "max_polyAorT"=max_polyAorT_random_test)

# dat_max_poly_random_all = rbind(dat_max_poly_random %>% select(-C0_new), dat_max_poly_random_test)

max_poly_lm_pred_random_all = predict(max_poly_lm, dat_max_poly_random_all)
cor(max_poly_lm_pred_random_all, dat_random_all$C0_new)
# -0.00829418

dat_poly = cbind(dat_npoly, dat_max_poly) %>%
  select(-C0_new)
dat_poly_test = cbind(dat_npoly_test, dat_max_poly_test)
dat_poly_random_all = cbind(dat_npoly_random_all, dat_max_poly_random_all) %>%
  select(-C0_new)
# saveRDS(dat_poly_random, "data/Created/random_poly_AorT.rds")


poly_lm = lm(C0_new~., dat_poly)
cor(poly_lm$fitted.values, dat$C0_new)
# 0.2592245

dat_poly_test = cbind(dat_npoly_test, dat_max_poly_test)
# saveRDS(dat_poly_test, "data/Created/tiling_poly_AorT_test.rds")
poly_lm_pred = predict(poly_lm, dat_poly_test)
cor(poly_lm_pred, dat_test$C0_new)
# 0.2627554

dat_poly_random_all = cbind(dat_npoly_random_all, dat_max_poly_random_all)
# saveRDS(dat_poly_random_all, "data/Created/tiling_poly_AorT_random_all.rds")
poly_lm_pred_random_all = predict(poly_lm, dat_poly_random_all)
cor(poly_lm_pred_random_all, dat_random_all$C0_new)
# 0.1249751

# set.seed(50)
# poly_gam = gam(C0_new~s(npolyA, k=8)+s(npolyT, k=10)+s(npolyAorT, k=10)+s(max_polyA, k=27)+s(max_polyT,k=30)+s(max_polyAorT,k=43), data=dat_poly)
# 
# cor(poly_gam$fitted.values, dat$C0_new)
# # 0.2485865
# 
# poly_gam_pred = predict(poly_gam, dat_poly_test)
# cor(poly_gam_pred, dat_test$C0_new)
# # 0.2425483


par(mfrow=c(2,2))
plot(max_polyA, dat$C0_new, col=alpha("blue", 0.1))
plot(max_polyT, dat$C0_new, col=alpha("blue", 0.1))

plot(max_polyA_random_all, dat_random_all$C0_new, col=alpha("blue", 0.1))
plot(max_polyT_random_all, dat_random_all$C0_new, col=alpha("blue", 0.1))

# Notice that there is a negative trend from ~2-9 but a positive trend from 10+ 
# in the tiling library, but there are no sequences with 10+ in the random library,
# so we should break this feature down more. Try counting each individual polyA or polyT
# subsequence for lengths 3-30

# npolyA10 = stri_count(dat$x50mer, 
#                       regex="[A]{10,}")
# cor(npolyA10, dat$C0_new)
# # 
# 
# npolyA_test = stri_count(dat_test$x50mer, 
#                          regex="[A]{3,}")
# 
# npolyT = stri_count(dat$x50mer,
#                     regex="[T]{3,}")
# cor(npolyT, dat$C0_new)
# # 
# 
# npolyT_test = stri_count(dat_test$x50mer,
#                          regex="[T]{3,}")
# 
# npolyAorT = stri_count(dat$x50mer,
#                        regex="[AT]{3,}")
# cor(npolyAorT, dat$C0_new)
# # 
# 
# npolyAorT_test = stri_count(dat_test$x50mer,
#                             regex="[AT]{3,}")
# 
# dat_npoly = cbind(dat %>% select(C0_new),
#                   npolyA, npolyT, npolyAorT)
# npoly_lm = lm(C0_new~., dat_npoly)
# cor(npoly_lm$fitted.values, dat$C0_new)
# # 
# 
# dat_npoly_test = data.frame("npolyA"=npolyA_test, "npolyT"=npolyT_test, "npolyAorT"=npolyAorT_test)
# npoly_lm_pred = predict(npoly_lm, dat_npoly_test)
# cor(npoly_lm_pred, dat_test$C0_new)
# # 






