library(tidyverse)
dat = readRDS("data/Created/processed_chrV_post_smooth_C0.rds")
y = dat$post_smooth_5bp
dat_test = readRDS("data/Created/processed_chrV_post_smooth_C0_test.rds")
y_test = dat_test$post_smooth_5bp

ps2 <- paste0("X", 1:49, "di")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

max_freq_at_dist = function(row, dist, di1, di2) {
  return(sum(row[ps2][which(row[ps2] == di1) [which(row[ps2] == di1) < (50-dist)] + dist] == di2))
}

temp_freq_at_dist10_AA_AA = apply(dat, 1, max_freq_at_dist, dist=10, di1="AA", di2="AA")
temp_freq_at_dist10_AA_AT2 = apply(dat, 1, max_freq_at_dist, dist=10, di1="AA", di2="AT")
temp_freq_at_dist10_TT_TT2 = apply(dat, 1, max_freq_at_dist, dist=10, di1="TT", di2="TT")

find_HRI = function(row, di1, di2) {
  # First distance at 10, 20, 30:
  max_freq_10 = max(max_freq_at_dist(row, 9, di1, di2),
                    max_freq_at_dist(row, 10, di1, di2),
                    max_freq_at_dist(row, 11, di1, di2))
  max_freq_20 = max(max_freq_at_dist(row, 19, di1, di2),
                    max_freq_at_dist(row, 20, di1, di2),
                    max_freq_at_dist(row, 21, di1, di2))
  max_freq_30 = max(max_freq_at_dist(row, 29, di1, di2),
                    max_freq_at_dist(row, 30, di1, di2),
                    max_freq_at_dist(row, 31, di1, di2))
  
  # Now for distance at 5, 15, 25:
  max_freq_5 = max(max_freq_at_dist(row, 4, di1, di2),
                   max_freq_at_dist(row, 5, di1, di2),
                   max_freq_at_dist(row, 6, di1, di2))
  max_freq_15 = max(max_freq_at_dist(row, 14, di1, di2),
                    max_freq_at_dist(row, 15, di1, di2),
                    max_freq_at_dist(row, 16, di1, di2))
  max_freq_25 = max(max_freq_at_dist(row, 24, di1, di2),
                    max_freq_at_dist(row, 25, di1, di2),
                    max_freq_at_dist(row, 26, di1, di2))
  
  return(max_freq_10 + max_freq_20 + max_freq_30 - max_freq_5 - max_freq_15 - max_freq_25)
}

HRI_AA_AA = apply(dat, 1, find_HRI, di1="AA", di2="AA")
HRI_AA_AC = apply(dat, 1, find_HRI, di1="AA", di2="AC")
HRI_AA_AG = apply(dat, 1, find_HRI, di1="AA", di2="AG")
HRI_AA_AT = apply(dat, 1, find_HRI, di1="AA", di2="AT")
HRI_AA_CA = apply(dat, 1, find_HRI, di1="AA", di2="CA")
HRI_AA_CC = apply(dat, 1, find_HRI, di1="AA", di2="CC")
HRI_AA_CG = apply(dat, 1, find_HRI, di1="AA", di2="CG")
HRI_AA_CT = apply(dat, 1, find_HRI, di1="AA", di2="CT")
HRI_AA_GA = apply(dat, 1, find_HRI, di1="AA", di2="GA")
HRI_AA_GC = apply(dat, 1, find_HRI, di1="AA", di2="GC")
HRI_AA_GG = apply(dat, 1, find_HRI, di1="AA", di2="GG")
HRI_AA_GT = apply(dat, 1, find_HRI, di1="AA", di2="GT")
HRI_AA_TA = apply(dat, 1, find_HRI, di1="AA", di2="TA")
HRI_AA_TC = apply(dat, 1, find_HRI, di1="AA", di2="TC")
HRI_AA_TG = apply(dat, 1, find_HRI, di1="AA", di2="TG")
HRI_AA_TT = apply(dat, 1, find_HRI, di1="AA", di2="TT")

HRI_AA = paste0("HRI_AA_", dinucleotides)
HRI_AA_chrV_post_smooth_C0 = data.frame(map(HRI_AA, get))
colnames(HRI_AA_chrV_post_smooth_C0) = HRI_AA
saveRDS(HRI_AA_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_AA.rds")

HRI_AC_AC = apply(dat, 1, find_HRI, di1="AC", di2="AC")
HRI_AC_AG = apply(dat, 1, find_HRI, di1="AC", di2="AG")
HRI_AC_AT = apply(dat, 1, find_HRI, di1="AC", di2="AT")
HRI_AC_CA = apply(dat, 1, find_HRI, di1="AC", di2="CA")
HRI_AC_CC = apply(dat, 1, find_HRI, di1="AC", di2="CC")
HRI_AC_CG = apply(dat, 1, find_HRI, di1="AC", di2="CG")
HRI_AC_CT = apply(dat, 1, find_HRI, di1="AC", di2="CT")
HRI_AC_GA = apply(dat, 1, find_HRI, di1="AC", di2="GA")
HRI_AC_GC = apply(dat, 1, find_HRI, di1="AC", di2="GC")
HRI_AC_GG = apply(dat, 1, find_HRI, di1="AC", di2="GG")
HRI_AC_GT = apply(dat, 1, find_HRI, di1="AC", di2="GT")
HRI_AC_TA = apply(dat, 1, find_HRI, di1="AC", di2="TA")
HRI_AC_TC = apply(dat, 1, find_HRI, di1="AC", di2="TC")
HRI_AC_TG = apply(dat, 1, find_HRI, di1="AC", di2="TG")
HRI_AC_TT = apply(dat, 1, find_HRI, di1="AC", di2="TT")

HRI_AC = paste0("HRI_AC_", dinucleotides[-1])
HRI_AC_chrV_post_smooth_C0 = data.frame(map(HRI_AC, get))
colnames(HRI_AC_chrV_post_smooth_C0) = HRI_AC
saveRDS(HRI_AC_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_AC.rds")

HRI_AG_AG = apply(dat, 1, find_HRI, di1="AG", di2="AG")
HRI_AG_AT = apply(dat, 1, find_HRI, di1="AG", di2="AT")
HRI_AG_CA = apply(dat, 1, find_HRI, di1="AG", di2="CA")
HRI_AG_CC = apply(dat, 1, find_HRI, di1="AG", di2="CC")
HRI_AG_CG = apply(dat, 1, find_HRI, di1="AG", di2="CG")
HRI_AG_CT = apply(dat, 1, find_HRI, di1="AG", di2="CT")
HRI_AG_GA = apply(dat, 1, find_HRI, di1="AG", di2="GA")
HRI_AG_GC = apply(dat, 1, find_HRI, di1="AG", di2="GC")
HRI_AG_GG = apply(dat, 1, find_HRI, di1="AG", di2="GG")
HRI_AG_GT = apply(dat, 1, find_HRI, di1="AG", di2="GT")
HRI_AG_TA = apply(dat, 1, find_HRI, di1="AG", di2="TA")
HRI_AG_TC = apply(dat, 1, find_HRI, di1="AG", di2="TC")
HRI_AG_TG = apply(dat, 1, find_HRI, di1="AG", di2="TG")
HRI_AG_TT = apply(dat, 1, find_HRI, di1="AG", di2="TT")

HRI_AG = paste0("HRI_AG_", dinucleotides[-(1:2)])
HRI_AG_chrV_post_smooth_C0 = data.frame(map(HRI_AG, get))
colnames(HRI_AG_chrV_post_smooth_C0) = HRI_AG
saveRDS(HRI_AG_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_AG.rds")

HRI_AT_AT = apply(dat, 1, find_HRI, di1="AT", di2="AT")
HRI_AT_CA = apply(dat, 1, find_HRI, di1="AT", di2="CA")
HRI_AT_CC = apply(dat, 1, find_HRI, di1="AT", di2="CC")
HRI_AT_CG = apply(dat, 1, find_HRI, di1="AT", di2="CG")
HRI_AT_CT = apply(dat, 1, find_HRI, di1="AT", di2="CT")
HRI_AT_GA = apply(dat, 1, find_HRI, di1="AT", di2="GA")
HRI_AT_GC = apply(dat, 1, find_HRI, di1="AT", di2="GC")
HRI_AT_GG = apply(dat, 1, find_HRI, di1="AT", di2="GG")
HRI_AT_GT = apply(dat, 1, find_HRI, di1="AT", di2="GT")
HRI_AT_TA = apply(dat, 1, find_HRI, di1="AT", di2="TA")
HRI_AT_TC = apply(dat, 1, find_HRI, di1="AT", di2="TC")
HRI_AT_TG = apply(dat, 1, find_HRI, di1="AT", di2="TG")
HRI_AT_TT = apply(dat, 1, find_HRI, di1="AT", di2="TT")

HRI_AT = paste0("HRI_AT_", dinucleotides[-(1:3)])
HRI_AT_chrV_post_smooth_C0 = data.frame(map(HRI_AT, get))
colnames(HRI_AT_chrV_post_smooth_C0) = HRI_AT
saveRDS(HRI_AT_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_AT.rds")

HRI_CA_CA = apply(dat, 1, find_HRI, di1="CA", di2="CA")
HRI_CA_CC = apply(dat, 1, find_HRI, di1="CA", di2="CC")
HRI_CA_CG = apply(dat, 1, find_HRI, di1="CA", di2="CG")
HRI_CA_CT = apply(dat, 1, find_HRI, di1="CA", di2="CT")
HRI_CA_GA = apply(dat, 1, find_HRI, di1="CA", di2="GA")
HRI_CA_GC = apply(dat, 1, find_HRI, di1="CA", di2="GC")
HRI_CA_GG = apply(dat, 1, find_HRI, di1="CA", di2="GG")
HRI_CA_GT = apply(dat, 1, find_HRI, di1="CA", di2="GT")
HRI_CA_TA = apply(dat, 1, find_HRI, di1="CA", di2="TA")
HRI_CA_TC = apply(dat, 1, find_HRI, di1="CA", di2="TC")
HRI_CA_TG = apply(dat, 1, find_HRI, di1="CA", di2="TG")
HRI_CA_TT = apply(dat, 1, find_HRI, di1="CA", di2="TT")

HRI_CA = paste0("HRI_CA_", dinucleotides[-(1:4)])
HRI_CA_chrV_post_smooth_C0 = data.frame(map(HRI_CA, get))
colnames(HRI_CA_chrV_post_smooth_C0) = HRI_CA
saveRDS(HRI_CA_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_CA.rds")

HRI_CC_CC = apply(dat, 1, find_HRI, di1="CC", di2="CC")
HRI_CC_CG = apply(dat, 1, find_HRI, di1="CC", di2="CG")
HRI_CC_CT = apply(dat, 1, find_HRI, di1="CC", di2="CT")
HRI_CC_GA = apply(dat, 1, find_HRI, di1="CC", di2="GA")
HRI_CC_GC = apply(dat, 1, find_HRI, di1="CC", di2="GC")
HRI_CC_GG = apply(dat, 1, find_HRI, di1="CC", di2="GG")
HRI_CC_GT = apply(dat, 1, find_HRI, di1="CC", di2="GT")
HRI_CC_TA = apply(dat, 1, find_HRI, di1="CC", di2="TA")
HRI_CC_TC = apply(dat, 1, find_HRI, di1="CC", di2="TC")
HRI_CC_TG = apply(dat, 1, find_HRI, di1="CC", di2="TG")
HRI_CC_TT = apply(dat, 1, find_HRI, di1="CC", di2="TT")

HRI_CC = paste0("HRI_CC_", dinucleotides[-(1:5)])
HRI_CC_chrV_post_smooth_C0 = data.frame(map(HRI_CC, get))
colnames(HRI_CC_chrV_post_smooth_C0) = HRI_CC
saveRDS(HRI_CC_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_CC.rds")

HRI_CG_CG = apply(dat, 1, find_HRI, di1="CG", di2="CG")
HRI_CG_CT = apply(dat, 1, find_HRI, di1="CG", di2="CT")
HRI_CG_GA = apply(dat, 1, find_HRI, di1="CG", di2="GA")
HRI_CG_GC = apply(dat, 1, find_HRI, di1="CG", di2="GC")
HRI_CG_GG = apply(dat, 1, find_HRI, di1="CG", di2="GG")
HRI_CG_GT = apply(dat, 1, find_HRI, di1="CG", di2="GT")
HRI_CG_TA = apply(dat, 1, find_HRI, di1="CG", di2="TA")
HRI_CG_TC = apply(dat, 1, find_HRI, di1="CG", di2="TC")
HRI_CG_TG = apply(dat, 1, find_HRI, di1="CG", di2="TG")
HRI_CG_TT = apply(dat, 1, find_HRI, di1="CG", di2="TT")

HRI_CG = paste0("HRI_CG_", dinucleotides[-(1:6)])
HRI_CG_chrV_post_smooth_C0 = data.frame(map(HRI_CG, get))
colnames(HRI_CG_chrV_post_smooth_C0) = HRI_CG
saveRDS(HRI_CG_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_CG.rds")

HRI_CT_CT = apply(dat, 1, find_HRI, di1="CT", di2="CT")
HRI_CT_GA = apply(dat, 1, find_HRI, di1="CT", di2="GA")
HRI_CT_GC = apply(dat, 1, find_HRI, di1="CT", di2="GC")
HRI_CT_GG = apply(dat, 1, find_HRI, di1="CT", di2="GG")
HRI_CT_GT = apply(dat, 1, find_HRI, di1="CT", di2="GT")
HRI_CT_TA = apply(dat, 1, find_HRI, di1="CT", di2="TA")
HRI_CT_TC = apply(dat, 1, find_HRI, di1="CT", di2="TC")
HRI_CT_TG = apply(dat, 1, find_HRI, di1="CT", di2="TG")
HRI_CT_TT = apply(dat, 1, find_HRI, di1="CT", di2="TT")

HRI_CT = paste0("HRI_CT_", dinucleotides[-(1:7)])
HRI_CT_chrV_post_smooth_C0 = data.frame(map(HRI_CT, get))
colnames(HRI_CT_chrV_post_smooth_C0) = HRI_CT
saveRDS(HRI_CT_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_CT.rds")

HRI_GA_GA = apply(dat, 1, find_HRI, di1="GA", di2="GA")
HRI_GA_GC = apply(dat, 1, find_HRI, di1="GA", di2="GC")
HRI_GA_GG = apply(dat, 1, find_HRI, di1="GA", di2="GG")
HRI_GA_GT = apply(dat, 1, find_HRI, di1="GA", di2="GT")
HRI_GA_TA = apply(dat, 1, find_HRI, di1="GA", di2="TA")
HRI_GA_TC = apply(dat, 1, find_HRI, di1="GA", di2="TC")
HRI_GA_TG = apply(dat, 1, find_HRI, di1="GA", di2="TG")
HRI_GA_TT = apply(dat, 1, find_HRI, di1="GA", di2="TT")

HRI_GA = paste0("HRI_GA_", dinucleotides[-(1:8)])
HRI_GA_chrV_post_smooth_C0 = data.frame(map(HRI_GA, get))
colnames(HRI_GA_chrV_post_smooth_C0) = HRI_GA
saveRDS(HRI_GA_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_GA.rds")

HRI_GC_GC = apply(dat, 1, find_HRI, di1="GC", di2="GC")
HRI_GC_GG = apply(dat, 1, find_HRI, di1="GC", di2="GG")
HRI_GC_GT = apply(dat, 1, find_HRI, di1="GC", di2="GT")
HRI_GC_TA = apply(dat, 1, find_HRI, di1="GC", di2="TA")
HRI_GC_TC = apply(dat, 1, find_HRI, di1="GC", di2="TC")
HRI_GC_TG = apply(dat, 1, find_HRI, di1="GC", di2="TG")
HRI_GC_TT = apply(dat, 1, find_HRI, di1="GC", di2="TT")

HRI_GC = paste0("HRI_GC_", dinucleotides[-(1:9)])
HRI_GC_chrV_post_smooth_C0 = data.frame(map(HRI_GC, get))
colnames(HRI_GC_chrV_post_smooth_C0) = HRI_GC
saveRDS(HRI_GC_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_GC.rds")

HRI_GG_GG = apply(dat, 1, find_HRI, di1="GG", di2="GG")
HRI_GG_GT = apply(dat, 1, find_HRI, di1="GG", di2="GT")
HRI_GG_TA = apply(dat, 1, find_HRI, di1="GG", di2="TA")
HRI_GG_TC = apply(dat, 1, find_HRI, di1="GG", di2="TC")
HRI_GG_TG = apply(dat, 1, find_HRI, di1="GG", di2="TG")
HRI_GG_TT = apply(dat, 1, find_HRI, di1="GG", di2="TT")

HRI_GG = paste0("HRI_GG_", dinucleotides[-(1:10)])
HRI_GG_chrV_post_smooth_C0 = data.frame(map(HRI_GG, get))
colnames(HRI_GG_chrV_post_smooth_C0) = HRI_GG
saveRDS(HRI_GG_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_GG.rds")

HRI_GT_GT = apply(dat, 1, find_HRI, di1="GT", di2="GT")
HRI_GT_TA = apply(dat, 1, find_HRI, di1="GT", di2="TA")
HRI_GT_TC = apply(dat, 1, find_HRI, di1="GT", di2="TC")
HRI_GT_TG = apply(dat, 1, find_HRI, di1="GT", di2="TG")
HRI_GT_TT = apply(dat, 1, find_HRI, di1="GT", di2="TT")

HRI_GT = paste0("HRI_GT_", dinucleotides[-(1:11)])
HRI_GT_chrV_post_smooth_C0 = data.frame(map(HRI_GT, get))
colnames(HRI_GT_chrV_post_smooth_C0) = HRI_GT
saveRDS(HRI_GT_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_GT.rds")

HRI_TA_TA = apply(dat, 1, find_HRI, di1="TA", di2="TA")
HRI_TA_TC = apply(dat, 1, find_HRI, di1="TA", di2="TC")
HRI_TA_TG = apply(dat, 1, find_HRI, di1="TA", di2="TG")
HRI_TA_TT = apply(dat, 1, find_HRI, di1="TA", di2="TT")

HRI_TA = paste0("HRI_TA_", dinucleotides[-(1:12)])
HRI_TA_chrV_post_smooth_C0 = data.frame(map(HRI_TA, get))
colnames(HRI_TA_chrV_post_smooth_C0) = HRI_TA
saveRDS(HRI_TA_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_TA.rds")

HRI_TC_TC = apply(dat, 1, find_HRI, di1="TC", di2="TC")
HRI_TC_TG = apply(dat, 1, find_HRI, di1="TC", di2="TG")
HRI_TC_TT = apply(dat, 1, find_HRI, di1="TC", di2="TT")

HRI_TC = paste0("HRI_TC_", dinucleotides[-(1:13)])
HRI_TC_chrV_post_smooth_C0 = data.frame(map(HRI_TC, get))
colnames(HRI_TC_chrV_post_smooth_C0) = HRI_TC
saveRDS(HRI_TC_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_TC.rds")

HRI_TG_TG = apply(dat, 1, find_HRI, di1="TG", di2="TG")
HRI_TG_TT = apply(dat, 1, find_HRI, di1="TG", di2="TT")

HRI_TG = paste0("HRI_TG_", dinucleotides[-(1:14)])
HRI_TG_chrV_post_smooth_C0 = data.frame(map(HRI_TG, get))
colnames(HRI_TG_chrV_post_smooth_C0) = HRI_TG
saveRDS(HRI_TG_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_TG.rds")

HRI_TT_TT = apply(dat, 1, find_HRI, di1="TT", di2="TT")

HRI_TT_chrV_post_smooth_C0 = data.frame(HRI_TT_TT)
colnames(HRI_TT_chrV_post_smooth_C0) = "HRI_TT_TT"
saveRDS(HRI_TT_chrV_post_smooth_C0, "data/Created/chrV_post_smooth_C0_HRI_TT.rds")

HRI = data.frame(map(paste0("HRI_", dinucleotides, "_chrV_post_smooth_C0"), get))
saveRDS(HRI, "data/Created/chrV_post_smooth_C0_HRI.rds")

HRI_AA_AA_test = apply(dat_test, 1, find_HRI, di1="AA", di2="AA")
HRI_AA_AC_test = apply(dat_test, 1, find_HRI, di1="AA", di2="AC")
HRI_AA_AG_test = apply(dat_test, 1, find_HRI, di1="AA", di2="AG")
HRI_AA_AT_test = apply(dat_test, 1, find_HRI, di1="AA", di2="AT")
HRI_AA_CA_test = apply(dat_test, 1, find_HRI, di1="AA", di2="CA")
HRI_AA_CC_test = apply(dat_test, 1, find_HRI, di1="AA", di2="CC")
HRI_AA_CG_test = apply(dat_test, 1, find_HRI, di1="AA", di2="CG")
HRI_AA_CT_test = apply(dat_test, 1, find_HRI, di1="AA", di2="CT")
HRI_AA_GA_test = apply(dat_test, 1, find_HRI, di1="AA", di2="GA")
HRI_AA_GC_test = apply(dat_test, 1, find_HRI, di1="AA", di2="GC")
HRI_AA_GG_test = apply(dat_test, 1, find_HRI, di1="AA", di2="GG")
HRI_AA_GT_test = apply(dat_test, 1, find_HRI, di1="AA", di2="GT")
HRI_AA_TA_test = apply(dat_test, 1, find_HRI, di1="AA", di2="TA")
HRI_AA_TC_test = apply(dat_test, 1, find_HRI, di1="AA", di2="TC")
HRI_AA_TG_test = apply(dat_test, 1, find_HRI, di1="AA", di2="TG")
HRI_AA_TT_test = apply(dat_test, 1, find_HRI, di1="AA", di2="TT")

HRI_AA_test = paste0("HRI_AA_", dinucleotides, "_test")
HRI_AA_chrV_post_smooth_C0_test = data.frame(map(HRI_AA_test, get))
colnames(HRI_AA_chrV_post_smooth_C0_test) = paste0("HRI_AA_", dinucleotides)
saveRDS(HRI_AA_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_AA_test.rds")

HRI_AC_AC_test = apply(dat_test, 1, find_HRI, di1="AC", di2="AC")
HRI_AC_AG_test = apply(dat_test, 1, find_HRI, di1="AC", di2="AG")
HRI_AC_AT_test = apply(dat_test, 1, find_HRI, di1="AC", di2="AT")
HRI_AC_CA_test = apply(dat_test, 1, find_HRI, di1="AC", di2="CA")
HRI_AC_CC_test = apply(dat_test, 1, find_HRI, di1="AC", di2="CC")
HRI_AC_CG_test = apply(dat_test, 1, find_HRI, di1="AC", di2="CG")
HRI_AC_CT_test = apply(dat_test, 1, find_HRI, di1="AC", di2="CT")
HRI_AC_GA_test = apply(dat_test, 1, find_HRI, di1="AC", di2="GA")
HRI_AC_GC_test = apply(dat_test, 1, find_HRI, di1="AC", di2="GC")
HRI_AC_GG_test = apply(dat_test, 1, find_HRI, di1="AC", di2="GG")
HRI_AC_GT_test = apply(dat_test, 1, find_HRI, di1="AC", di2="GT")
HRI_AC_TA_test = apply(dat_test, 1, find_HRI, di1="AC", di2="TA")
HRI_AC_TC_test = apply(dat_test, 1, find_HRI, di1="AC", di2="TC")
HRI_AC_TG_test = apply(dat_test, 1, find_HRI, di1="AC", di2="TG")
HRI_AC_TT_test = apply(dat_test, 1, find_HRI, di1="AC", di2="TT")

HRI_AC_test = paste0("HRI_AC_", dinucleotides[-1], "_test")
HRI_AC_chrV_post_smooth_C0_test = data.frame(map(HRI_AC_test, get))
colnames(HRI_AC_chrV_post_smooth_C0_test) = paste0("HRI_AC_", dinucleotides[-1])
saveRDS(HRI_AC_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_AC_test.rds")

HRI_AG_AG_test = apply(dat_test, 1, find_HRI, di1="AG", di2="AG")
HRI_AG_AT_test = apply(dat_test, 1, find_HRI, di1="AG", di2="AT")
HRI_AG_CA_test = apply(dat_test, 1, find_HRI, di1="AG", di2="CA")
HRI_AG_CC_test = apply(dat_test, 1, find_HRI, di1="AG", di2="CC")
HRI_AG_CG_test = apply(dat_test, 1, find_HRI, di1="AG", di2="CG")
HRI_AG_CT_test = apply(dat_test, 1, find_HRI, di1="AG", di2="CT")
HRI_AG_GA_test = apply(dat_test, 1, find_HRI, di1="AG", di2="GA")
HRI_AG_GC_test = apply(dat_test, 1, find_HRI, di1="AG", di2="GC")
HRI_AG_GG_test = apply(dat_test, 1, find_HRI, di1="AG", di2="GG")
HRI_AG_GT_test = apply(dat_test, 1, find_HRI, di1="AG", di2="GT")
HRI_AG_TA_test = apply(dat_test, 1, find_HRI, di1="AG", di2="TA")
HRI_AG_TC_test = apply(dat_test, 1, find_HRI, di1="AG", di2="TC")
HRI_AG_TG_test = apply(dat_test, 1, find_HRI, di1="AG", di2="TG")
HRI_AG_TT_test = apply(dat_test, 1, find_HRI, di1="AG", di2="TT")

HRI_AG_test = paste0("HRI_AG_", dinucleotides[-(1:2)], "_test")
HRI_AG_chrV_post_smooth_C0_test = data.frame(map(HRI_AG_test, get))
colnames(HRI_AG_chrV_post_smooth_C0_test) = paste0("HRI_AG_", dinucleotides[-(1:2)])
saveRDS(HRI_AG_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_AG_test.rds")

HRI_AT_AT_test = apply(dat_test, 1, find_HRI, di1="AT", di2="AT")
HRI_AT_CA_test = apply(dat_test, 1, find_HRI, di1="AT", di2="CA")
HRI_AT_CC_test = apply(dat_test, 1, find_HRI, di1="AT", di2="CC")
HRI_AT_CG_test = apply(dat_test, 1, find_HRI, di1="AT", di2="CG")
HRI_AT_CT_test = apply(dat_test, 1, find_HRI, di1="AT", di2="CT")
HRI_AT_GA_test = apply(dat_test, 1, find_HRI, di1="AT", di2="GA")
HRI_AT_GC_test = apply(dat_test, 1, find_HRI, di1="AT", di2="GC")
HRI_AT_GG_test = apply(dat_test, 1, find_HRI, di1="AT", di2="GG")
HRI_AT_GT_test = apply(dat_test, 1, find_HRI, di1="AT", di2="GT")
HRI_AT_TA_test = apply(dat_test, 1, find_HRI, di1="AT", di2="TA")
HRI_AT_TC_test = apply(dat_test, 1, find_HRI, di1="AT", di2="TC")
HRI_AT_TG_test = apply(dat_test, 1, find_HRI, di1="AT", di2="TG")
HRI_AT_TT_test = apply(dat_test, 1, find_HRI, di1="AT", di2="TT")

HRI_AT_test = paste0("HRI_AT_", dinucleotides[-(1:3)], "_test")
HRI_AT_chrV_post_smooth_C0_test = data.frame(map(HRI_AT_test, get))
colnames(HRI_AT_chrV_post_smooth_C0_test) = paste0("HRI_AT_", dinucleotides[-(1:3)])
saveRDS(HRI_AT_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_AT_test.rds")

HRI_CA_CA_test = apply(dat_test, 1, find_HRI, di1="CA", di2="CA")
HRI_CA_CC_test = apply(dat_test, 1, find_HRI, di1="CA", di2="CC")
HRI_CA_CG_test = apply(dat_test, 1, find_HRI, di1="CA", di2="CG")
HRI_CA_CT_test = apply(dat_test, 1, find_HRI, di1="CA", di2="CT")
HRI_CA_GA_test = apply(dat_test, 1, find_HRI, di1="CA", di2="GA")
HRI_CA_GC_test = apply(dat_test, 1, find_HRI, di1="CA", di2="GC")
HRI_CA_GG_test = apply(dat_test, 1, find_HRI, di1="CA", di2="GG")
HRI_CA_GT_test = apply(dat_test, 1, find_HRI, di1="CA", di2="GT")
HRI_CA_TA_test = apply(dat_test, 1, find_HRI, di1="CA", di2="TA")
HRI_CA_TC_test = apply(dat_test, 1, find_HRI, di1="CA", di2="TC")
HRI_CA_TG_test = apply(dat_test, 1, find_HRI, di1="CA", di2="TG")
HRI_CA_TT_test = apply(dat_test, 1, find_HRI, di1="CA", di2="TT")

HRI_CA_test = paste0("HRI_CA_", dinucleotides[-(1:4)], "_test")
HRI_CA_chrV_post_smooth_C0_test = data.frame(map(HRI_CA_test, get))
colnames(HRI_CA_chrV_post_smooth_C0_test) = paste0("HRI_CA_", dinucleotides[-(1:4)])
saveRDS(HRI_CA_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_CA_test.rds")

HRI_CC_CC_test = apply(dat_test, 1, find_HRI, di1="CC", di2="CC")
HRI_CC_CG_test = apply(dat_test, 1, find_HRI, di1="CC", di2="CG")
HRI_CC_CT_test = apply(dat_test, 1, find_HRI, di1="CC", di2="CT")
HRI_CC_GA_test = apply(dat_test, 1, find_HRI, di1="CC", di2="GA")
HRI_CC_GC_test = apply(dat_test, 1, find_HRI, di1="CC", di2="GC")
HRI_CC_GG_test = apply(dat_test, 1, find_HRI, di1="CC", di2="GG")
HRI_CC_GT_test = apply(dat_test, 1, find_HRI, di1="CC", di2="GT")
HRI_CC_TA_test = apply(dat_test, 1, find_HRI, di1="CC", di2="TA")
HRI_CC_TC_test = apply(dat_test, 1, find_HRI, di1="CC", di2="TC")
HRI_CC_TG_test = apply(dat_test, 1, find_HRI, di1="CC", di2="TG")
HRI_CC_TT_test = apply(dat_test, 1, find_HRI, di1="CC", di2="TT")

HRI_CC_test = paste0("HRI_CC_", dinucleotides[-(1:5)], "_test")
HRI_CC_chrV_post_smooth_C0_test = data.frame(map(HRI_CC_test, get))
colnames(HRI_CC_chrV_post_smooth_C0_test) = paste0("HRI_CC_", dinucleotides[-(1:5)])
saveRDS(HRI_CC_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_CC_test.rds")

HRI_CG_CG_test = apply(dat_test, 1, find_HRI, di1="CG", di2="CG")
HRI_CG_CT_test = apply(dat_test, 1, find_HRI, di1="CG", di2="CT")
HRI_CG_GA_test = apply(dat_test, 1, find_HRI, di1="CG", di2="GA")
HRI_CG_GC_test = apply(dat_test, 1, find_HRI, di1="CG", di2="GC")
HRI_CG_GG_test = apply(dat_test, 1, find_HRI, di1="CG", di2="GG")
HRI_CG_GT_test = apply(dat_test, 1, find_HRI, di1="CG", di2="GT")
HRI_CG_TA_test = apply(dat_test, 1, find_HRI, di1="CG", di2="TA")
HRI_CG_TC_test = apply(dat_test, 1, find_HRI, di1="CG", di2="TC")
HRI_CG_TG_test = apply(dat_test, 1, find_HRI, di1="CG", di2="TG")
HRI_CG_TT_test = apply(dat_test, 1, find_HRI, di1="CG", di2="TT")

HRI_CG_test = paste0("HRI_CG_", dinucleotides[-(1:6)], "_test")
HRI_CG_chrV_post_smooth_C0_test = data.frame(map(HRI_CG_test, get))
colnames(HRI_CG_chrV_post_smooth_C0_test) = paste0("HRI_CG_", dinucleotides[-(1:6)])
saveRDS(HRI_CG_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_CG_test.rds")

HRI_CT_CT_test = apply(dat_test, 1, find_HRI, di1="CT", di2="CT")
HRI_CT_GA_test = apply(dat_test, 1, find_HRI, di1="CT", di2="GA")
HRI_CT_GC_test = apply(dat_test, 1, find_HRI, di1="CT", di2="GC")
HRI_CT_GG_test = apply(dat_test, 1, find_HRI, di1="CT", di2="GG")
HRI_CT_GT_test = apply(dat_test, 1, find_HRI, di1="CT", di2="GT")
HRI_CT_TA_test = apply(dat_test, 1, find_HRI, di1="CT", di2="TA")
HRI_CT_TC_test = apply(dat_test, 1, find_HRI, di1="CT", di2="TC")
HRI_CT_TG_test = apply(dat_test, 1, find_HRI, di1="CT", di2="TG")
HRI_CT_TT_test = apply(dat_test, 1, find_HRI, di1="CT", di2="TT")

HRI_CT_test = paste0("HRI_CT_", dinucleotides[-(1:7)], "_test")
HRI_CT_chrV_post_smooth_C0_test = data.frame(map(HRI_CT_test, get))
colnames(HRI_CT_chrV_post_smooth_C0_test) = paste0("HRI_CT_", dinucleotides[-(1:7)])
saveRDS(HRI_CT_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_CT_test.rds")

HRI_GA_GA_test = apply(dat_test, 1, find_HRI, di1="GA", di2="GA")
HRI_GA_GC_test = apply(dat_test, 1, find_HRI, di1="GA", di2="GC")
HRI_GA_GG_test = apply(dat_test, 1, find_HRI, di1="GA", di2="GG")
HRI_GA_GT_test = apply(dat_test, 1, find_HRI, di1="GA", di2="GT")
HRI_GA_TA_test = apply(dat_test, 1, find_HRI, di1="GA", di2="TA")
HRI_GA_TC_test = apply(dat_test, 1, find_HRI, di1="GA", di2="TC")
HRI_GA_TG_test = apply(dat_test, 1, find_HRI, di1="GA", di2="TG")
HRI_GA_TT_test = apply(dat_test, 1, find_HRI, di1="GA", di2="TT")

HRI_GA_test = paste0("HRI_GA_", dinucleotides[-(1:8)], "_test")
HRI_GA_chrV_post_smooth_C0_test = data.frame(map(HRI_GA_test, get))
colnames(HRI_GA_chrV_post_smooth_C0_test) = paste0("HRI_GA_", dinucleotides[-(1:8)])
saveRDS(HRI_GA_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_GA_test.rds")

HRI_GC_GC_test = apply(dat_test, 1, find_HRI, di1="GC", di2="GC")
HRI_GC_GG_test = apply(dat_test, 1, find_HRI, di1="GC", di2="GG")
HRI_GC_GT_test = apply(dat_test, 1, find_HRI, di1="GC", di2="GT")
HRI_GC_TA_test = apply(dat_test, 1, find_HRI, di1="GC", di2="TA")
HRI_GC_TC_test = apply(dat_test, 1, find_HRI, di1="GC", di2="TC")
HRI_GC_TG_test = apply(dat_test, 1, find_HRI, di1="GC", di2="TG")
HRI_GC_TT_test = apply(dat_test, 1, find_HRI, di1="GC", di2="TT")

HRI_GC_test = paste0("HRI_GC_", dinucleotides[-(1:9)], "_test")
HRI_GC_chrV_post_smooth_C0_test = data.frame(map(HRI_GC_test, get))
colnames(HRI_GC_chrV_post_smooth_C0_test) = paste0("HRI_GC_", dinucleotides[-(1:9)])
saveRDS(HRI_GC_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_GC_test.rds")

HRI_GG_GG_test = apply(dat_test, 1, find_HRI, di1="GG", di2="GG")
HRI_GG_GT_test = apply(dat_test, 1, find_HRI, di1="GG", di2="GT")
HRI_GG_TA_test = apply(dat_test, 1, find_HRI, di1="GG", di2="TA")
HRI_GG_TC_test = apply(dat_test, 1, find_HRI, di1="GG", di2="TC")
HRI_GG_TG_test = apply(dat_test, 1, find_HRI, di1="GG", di2="TG")
HRI_GG_TT_test = apply(dat_test, 1, find_HRI, di1="GG", di2="TT")

HRI_GG_test = paste0("HRI_GG_", dinucleotides[-(1:10)], "_test")
HRI_GG_chrV_post_smooth_C0_test = data.frame(map(HRI_GG_test, get))
colnames(HRI_GG_chrV_post_smooth_C0_test) = paste0("HRI_GG_", dinucleotides[-(1:10)])
saveRDS(HRI_GG_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_GG_test.rds")

HRI_GT_GT_test = apply(dat_test, 1, find_HRI, di1="GT", di2="GT")
HRI_GT_TA_test = apply(dat_test, 1, find_HRI, di1="GT", di2="TA")
HRI_GT_TC_test = apply(dat_test, 1, find_HRI, di1="GT", di2="TC")
HRI_GT_TG_test = apply(dat_test, 1, find_HRI, di1="GT", di2="TG")
HRI_GT_TT_test = apply(dat_test, 1, find_HRI, di1="GT", di2="TT")

HRI_GT_test = paste0("HRI_GT_", dinucleotides[-(1:11)], "_test")
HRI_GT_chrV_post_smooth_C0_test = data.frame(map(HRI_GT_test, get))
colnames(HRI_GT_chrV_post_smooth_C0_test) = paste0("HRI_GT_", dinucleotides[-(1:11)])
saveRDS(HRI_GT_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_GT_test.rds")

HRI_TA_TA_test = apply(dat_test, 1, find_HRI, di1="TA", di2="TA")
HRI_TA_TC_test = apply(dat_test, 1, find_HRI, di1="TA", di2="TC")
HRI_TA_TG_test = apply(dat_test, 1, find_HRI, di1="TA", di2="TG")
HRI_TA_TT_test = apply(dat_test, 1, find_HRI, di1="TA", di2="TT")

HRI_TA_test = paste0("HRI_TA_", dinucleotides[-(1:12)], "_test")
HRI_TA_chrV_post_smooth_C0_test = data.frame(map(HRI_TA_test, get))
colnames(HRI_TA_chrV_post_smooth_C0_test) = paste0("HRI_TA_", dinucleotides[-(1:12)])
saveRDS(HRI_TA_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_TA_test.rds")

HRI_TC_TC_test = apply(dat_test, 1, find_HRI, di1="TC", di2="TC")
HRI_TC_TG_test = apply(dat_test, 1, find_HRI, di1="TC", di2="TG")
HRI_TC_TT_test = apply(dat_test, 1, find_HRI, di1="TC", di2="TT")

HRI_TC_test = paste0("HRI_TC_", dinucleotides[-(1:13)], "_test")
HRI_TC_chrV_post_smooth_C0_test = data.frame(map(HRI_TC_test, get))
colnames(HRI_TC_chrV_post_smooth_C0_test) = paste0("HRI_TC_", dinucleotides[-(1:13)])
saveRDS(HRI_TC_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_TC_test.rds")

HRI_TG_TG_test = apply(dat_test, 1, find_HRI, di1="TG", di2="TG")
HRI_TG_TT_test = apply(dat_test, 1, find_HRI, di1="TG", di2="TT")

HRI_TG_test = paste0("HRI_TG_", dinucleotides[-(1:14)], "_test")
HRI_TG_chrV_post_smooth_C0_test = data.frame(map(HRI_TG_test, get))
colnames(HRI_TG_chrV_post_smooth_C0_test) = paste0("HRI_TG_", dinucleotides[-(1:14)])
saveRDS(HRI_TG_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_TG_test.rds")

HRI_TT_TT_test = apply(dat_test, 1, find_HRI, di1="TT", di2="TT")

HRI_TT_chrV_post_smooth_C0_test = data.frame(HRI_TT_TT_test)
colnames(HRI_TT_chrV_post_smooth_C0_test) = "HRI_TT_TT"
saveRDS(HRI_TT_chrV_post_smooth_C0_test, "data/Created/chrV_post_smooth_C0_HRI_TT_test.rds")

HRI_test = data.frame(map(paste0("HRI_", dinucleotides, "_chrV_post_smooth_C0_test"), get))
saveRDS(HRI_test, "data/Created/chrV_post_smooth_C0_HRI_test.rds")



