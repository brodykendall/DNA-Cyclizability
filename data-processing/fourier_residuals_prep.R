############ fourier_di3:
######### Fourier features based on all dinucleotides calculated using fft


fourier_di_cos = readRDS("data/Created/tiling_fourier_di_cos.rds")
fourier_di_sin = readRDS("data/Created/tiling_fourier_di_sin.rds")
fourier_di_cos_test = readRDS("data/Created/tiling_fourier_di_cos_test.rds")
fourier_di_sin_test = readRDS("data/Created/tiling_fourier_di_sin_test.rds")

dat_fourier_di3 = cbind((fourier_di_cos)^2, (fourier_di_sin)^2, 
                        dat %>% select(C0_new, x50mer))

dat_test_fourier_di3 = cbind((fourier_di_cos_test)^2, (fourier_di_sin_test)^2, 
                             dat_test %>% select(C0_new, x50mer))

dat_tiling_all_fourier_di3 = rbind(dat_fourier_di3, dat_test_fourier_di3)

all_fourier_di3_lm = lm(C0_new~., data=(dat_tiling_all_fourier_di3 %>% select(-x50mer)))
cor(all_fourier_di3_lm$fitted.values, dat_tiling_all_fourier_di3$C0_new)
# 0.4547124

dat_tiling_all_fourier_di3.seq_and_resids = cbind(dat_tiling_all_fourier_di3 %>% select(x50mer),
                                                  "residuals"=all_fourier_di3_lm$residuals)

write.csv(dat_tiling_all_fourier_di3.seq_and_resids, "data/Created/tiling_all_fourier_di3_resids.csv")





##### Random:

fourier_di_cos_random = readRDS("data/Created/random_fourier_di_cos.rds")
fourier_di_sin_random = readRDS("data/Created/random_fourier_di_sin.rds")
fourier_di_cos_random_test = readRDS("data/Created/random_fourier_di_cos_test.rds")
fourier_di_sin_random_test = readRDS("data/Created/random_fourier_di_sin_test.rds")

dat_random_fourier_di3 = cbind((fourier_di_cos_random)^2, (fourier_di_sin_random)^2, 
                               dat_random %>% select(C0_new, x50mer))

dat_random_test_fourier_di3 = cbind((fourier_di_cos_random_test)^2, (fourier_di_sin_random_test)^2, 
                                    dat_random_test %>% select(C0_new, x50mer))

dat_random_all_fourier_di3 = rbind(dat_random_fourier_di3, dat_random_test_fourier_di3)


all_fourier_di3.pred_random_all = predict(all_fourier_di3_lm, dat_random_all_fourier_di3)
cor(all_fourier_di3.pred_random_all, dat_random_all_fourier_di3$C0_new)
# 0.4274853

dat_random_all_fourier_di3.seq_and_resids = cbind(dat_random_all_fourier_di3 %>% select(x50mer),
                                                  "residuals"=(dat_random_all_fourier_di3$C0_new - all_fourier_di3.pred_random_all))

write.csv(dat_random_all_fourier_di3.seq_and_resids, "data/Created/random_all_fourier_di3_resids.csv")




##### chrV:

fourier_di_cos_chrV = readRDS("data/Created/chrV_fourier_di_cos.rds")
fourier_di_sin_chrV = readRDS("data/Created/chrV_fourier_di_sin.rds")
fourier_di_cos_chrV_test = readRDS("data/Created/chrV_fourier_di_cos_test.rds")
fourier_di_sin_chrV_test = readRDS("data/Created/chrV_fourier_di_sin_test.rds")

dat_chrV_fourier_di3 = cbind((fourier_di_cos_chrV)^2, (fourier_di_sin_chrV)^2, 
                             dat_chrV %>% select(C0_new, x50mer))

dat_chrV_test_fourier_di3 = cbind((fourier_di_cos_chrV_test)^2, (fourier_di_sin_chrV_test)^2, 
                                  dat_chrV_test %>% select(C0_new, x50mer))

dat_chrV_all_fourier_di3 = rbind(dat_chrV_fourier_di3, dat_chrV_test_fourier_di3)


all_fourier_di3.pred_chrV_all = predict(all_fourier_di3_lm, dat_chrV_all_fourier_di3)
cor(all_fourier_di3.pred_chrV_all, dat_chrV_all_fourier_di3$C0_new)
# 0.331892

dat_chrV_all_fourier_di3.seq_and_resids = cbind(dat_chrV_all_fourier_di3 %>% select(x50mer),
                                                "residuals"=(dat_chrV_all_fourier_di3$C0_new - all_fourier_di3.pred_chrV_all))

write.csv(dat_chrV_all_fourier_di3.seq_and_resids, "data/Created/chrV_all_fourier_di3_resids.csv")





##### Yeast:

fourier_di_cos_yeast = readRDS("data/Created/yeast_fourier_di_cos.rds")
fourier_di_sin_yeast = readRDS("data/Created/yeast_fourier_di_sin.rds")
fourier_di_cos_yeast_test = readRDS("data/Created/yeast_fourier_di_cos_test.rds")
fourier_di_sin_yeast_test = readRDS("data/Created/yeast_fourier_di_sin_test.rds")

dat_yeast_fourier_di3 = cbind((fourier_di_cos_yeast)^2, (fourier_di_sin_yeast)^2, 
                              dat_yeast %>% select(C0_new, x50mer))

dat_yeast_test_fourier_di3 = cbind((fourier_di_cos_yeast_test)^2, (fourier_di_sin_yeast_test)^2, 
                                   dat_yeast_test %>% select(C0_new, x50mer))

dat_yeast_all_fourier_di3 = rbind(dat_yeast_fourier_di3, dat_yeast_test_fourier_di3)


all_fourier_di3.pred_yeast_all = predict(all_fourier_di3_lm, dat_yeast_all_fourier_di3)
cor(all_fourier_di3.pred_yeast_all, dat_yeast_all_fourier_di3$C0_new)
# 0.3810148

dat_yeast_all_fourier_di3.seq_and_resids = cbind(dat_yeast_all_fourier_di3 %>% select(x50mer),
                                                 "residuals"=(dat_yeast_all_fourier_di3$C0_new - all_fourier_di3.pred_yeast_all))

write.csv(dat_yeast_all_fourier_di3.seq_and_resids, "data/Created/yeast_all_fourier_di3_resids.csv")







############ fourier_6_mod3:
######### Fourier Features Based on all Dinucleotides + Modulus of Difference of 
######### Fourier Features Based on all Dinucleotides with h=6 (period 9.8) Calculated using fft


fourier_di_cos = readRDS("data/Created/tiling_fourier_di_cos.rds")
fourier_di_sin = readRDS("data/Created/tiling_fourier_di_sin.rds")
fourier_di_cos_test = readRDS("data/Created/tiling_fourier_di_cos_test.rds")
fourier_di_sin_test = readRDS("data/Created/tiling_fourier_di_sin_test.rds")

fourier_6_mod = readRDS("data/Created/tiling_fourier_6_mod.rds")
fourier_6_mod_test = readRDS("data/Created/tiling_fourier_6_mod_test.rds")

dat_fourier_6_mod3 = cbind((fourier_di_cos)^2, (fourier_di_sin)^2, 
                           (fourier_6_mod)^2, dat %>% select(C0_new, x50mer))

dat_test_fourier_6_mod3 = cbind((fourier_di_cos_test)^2, (fourier_di_sin_test)^2, 
                                (fourier_6_mod_test)^2, dat_test %>% select(C0_new, x50mer))

dat_tiling_all_fourier_6_mod3 = rbind(dat_fourier_6_mod3, dat_test_fourier_6_mod3)

all_fourier_6_mod3_lm = lm(C0_new~., data=(dat_tiling_all_fourier_6_mod3 %>% select(-x50mer)))
cor(all_fourier_6_mod3_lm$fitted.values, dat_tiling_all_fourier_6_mod3$C0_new)
# 0.5939086

dat_tiling_all_fourier_6_mod3.seq_and_resids = cbind(dat_tiling_all_fourier_6_mod3 %>% select(x50mer),
                                                     "residuals"=all_fourier_6_mod3_lm$residuals)

write.csv(dat_tiling_all_fourier_6_mod3.seq_and_resids, "data/Created/tiling_all_fourier_6_mod3_resids.csv")





##### Random:

fourier_di_cos_random = readRDS("data/Created/random_fourier_di_cos.rds")
fourier_di_sin_random = readRDS("data/Created/random_fourier_di_sin.rds")
fourier_di_cos_random_test = readRDS("data/Created/random_fourier_di_cos_test.rds")
fourier_di_sin_random_test = readRDS("data/Created/random_fourier_di_sin_test.rds")

fourier_6_mod_random = readRDS("data/Created/random_fourier_6_mod.rds")
fourier_6_mod_random_test = readRDS("data/Created/random_fourier_6_mod_test.rds")


dat_random_fourier_6_mod3 = cbind((fourier_di_cos_random)^2, (fourier_di_sin_random)^2, 
                                  (fourier_6_mod_random)^2, dat_random %>% select(C0_new, x50mer))

dat_random_test_fourier_6_mod3 = cbind((fourier_di_cos_random_test)^2, (fourier_di_sin_random_test)^2, 
                                       (fourier_6_mod_random_test)^2, dat_random_test %>% select(C0_new, x50mer))

dat_random_all_fourier_6_mod3 = rbind(dat_random_fourier_6_mod3, dat_random_test_fourier_6_mod3)


all_fourier_6_mod3.pred_random_all = predict(all_fourier_6_mod3_lm, dat_random_all_fourier_6_mod3)
cor(all_fourier_6_mod3.pred_random_all, dat_random_all_fourier_6_mod3$C0_new)
# 0.5554498

dat_random_all_fourier_6_mod3.seq_and_resids = cbind(dat_random_all_fourier_6_mod3 %>% select(x50mer),
                                                     "residuals"=(dat_random_all_fourier_6_mod3$C0_new - all_fourier_6_mod3.pred_random_all))

write.csv(dat_random_all_fourier_6_mod3.seq_and_resids, "data/Created/random_all_fourier_6_mod3_resids.csv")




##### chrV:

fourier_di_cos_chrV = readRDS("data/Created/chrV_fourier_di_cos.rds")
fourier_di_sin_chrV = readRDS("data/Created/chrV_fourier_di_sin.rds")
fourier_di_cos_chrV_test = readRDS("data/Created/chrV_fourier_di_cos_test.rds")
fourier_di_sin_chrV_test = readRDS("data/Created/chrV_fourier_di_sin_test.rds")

fourier_6_mod_chrV = readRDS("data/Created/chrV_fourier_6_mod.rds")
fourier_6_mod_chrV_test = readRDS("data/Created/chrV_fourier_6_mod_test.rds")


dat_chrV_fourier_6_mod3 = cbind((fourier_di_cos_chrV)^2, (fourier_di_sin_chrV)^2, 
                                  (fourier_6_mod_chrV)^2, dat_chrV %>% select(C0_new, x50mer))

dat_chrV_test_fourier_6_mod3 = cbind((fourier_di_cos_chrV_test)^2, (fourier_di_sin_chrV_test)^2, 
                                       (fourier_6_mod_chrV_test)^2, dat_chrV_test %>% select(C0_new, x50mer))

dat_chrV_all_fourier_6_mod3 = rbind(dat_chrV_fourier_6_mod3, dat_chrV_test_fourier_6_mod3)


all_fourier_6_mod3.pred_chrV_all = predict(all_fourier_6_mod3_lm, dat_chrV_all_fourier_6_mod3)
cor(all_fourier_6_mod3.pred_chrV_all, dat_chrV_all_fourier_6_mod3$C0_new)
# 0.4642375

dat_chrV_all_fourier_6_mod3.seq_and_resids = cbind(dat_chrV_all_fourier_6_mod3 %>% select(x50mer),
                                                     "residuals"=(dat_chrV_all_fourier_6_mod3$C0_new - all_fourier_6_mod3.pred_chrV_all))

write.csv(dat_chrV_all_fourier_6_mod3.seq_and_resids, "data/Created/chrV_all_fourier_6_mod3_resids.csv")




##### Yeast:

fourier_di_cos_yeast = readRDS("data/Created/yeast_fourier_di_cos.rds")
fourier_di_sin_yeast = readRDS("data/Created/yeast_fourier_di_sin.rds")
fourier_di_cos_yeast_test = readRDS("data/Created/yeast_fourier_di_cos_test.rds")
fourier_di_sin_yeast_test = readRDS("data/Created/yeast_fourier_di_sin_test.rds")

fourier_6_mod_yeast = readRDS("data/Created/yeast_fourier_6_mod.rds")
fourier_6_mod_yeast_test = readRDS("data/Created/yeast_fourier_6_mod_test.rds")


dat_yeast_fourier_6_mod3 = cbind((fourier_di_cos_yeast)^2, (fourier_di_sin_yeast)^2, 
                                  (fourier_6_mod_yeast)^2, dat_yeast %>% select(C0_new, x50mer))

dat_yeast_test_fourier_6_mod3 = cbind((fourier_di_cos_yeast_test)^2, (fourier_di_sin_yeast_test)^2, 
                                       (fourier_6_mod_yeast_test)^2, dat_yeast_test %>% select(C0_new, x50mer))

dat_yeast_all_fourier_6_mod3 = rbind(dat_yeast_fourier_6_mod3, dat_yeast_test_fourier_6_mod3)


all_fourier_6_mod3.pred_yeast_all = predict(all_fourier_6_mod3_lm, dat_yeast_all_fourier_6_mod3)
cor(all_fourier_6_mod3.pred_yeast_all, dat_yeast_all_fourier_6_mod3$C0_new)
# 0.5316907

dat_yeast_all_fourier_6_mod3.seq_and_resids = cbind(dat_yeast_all_fourier_6_mod3 %>% select(x50mer),
                                                     "residuals"=(dat_yeast_all_fourier_6_mod3$C0_new - all_fourier_6_mod3.pred_yeast_all))

write.csv(dat_yeast_all_fourier_6_mod3.seq_and_resids, "data/Created/yeast_all_fourier_6_mod3_resids.csv")




