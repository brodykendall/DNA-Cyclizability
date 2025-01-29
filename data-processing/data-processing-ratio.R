library(tidyverse)
library(markovchain)

source("scripts/functions/markov-functions.R")
source("scripts/functions/sequence-functions.R")

dat <- readRDS("data/Created/processed.rds")

#Set seed for reproducibility
set.seed(50)

cutoffs <- quantile(dat$C0, c(0.2, 0.8))


dat_efficient <- dat[dat$C0 > cutoffs[2],]
dat_inefficient <- dat[dat$C0 < cutoffs[1],]

eff_mat_2nd <- position_specific_score_matrices(dat_efficient, order = 2)
ineff_mat_2nd <- position_specific_score_matrices(dat_inefficient, order = 2)
ratio_mat_2nd <- map2(eff_mat_2nd, ineff_mat_2nd, `-`)

dat$ratio_score_2nd <- pssm_sum_score(dat, ratio_mat_2nd)

saveRDS(ratio_mat_2nd, "model/ratio-matrices-order-2.rds")
saveRDS(dat, "data/Created/processed_ratio.rds")


dat_test <- readRDS("data/Created/processed_test.rds")

# dat_efficient_test <- dat_test[dat_test$C0 > cutoffs_test[2],]
# dat_inefficient_test <- dat_test[dat_test$C0 < cutoffs_test[1],]
# 
# eff_mat_2nd_test <- position_specific_score_matrices(dat_efficient_test, order = 2)
# ineff_mat_2nd_test <- position_specific_score_matrices(dat_inefficient_test, order = 2)
# ratio_mat_2nd_test <- map2(eff_mat_2nd_test, ineff_mat_2nd_test, `-`)

dat_test$ratio_score_2nd <- pssm_sum_score(dat_test, ratio_mat_2nd)

# saveRDS(ratio_mat_2nd_test, "model/ratio-matrices-order-2-test.rds")
saveRDS(dat_test, "data/Created/processed_test_ratio.rds")

