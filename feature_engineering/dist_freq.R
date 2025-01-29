library(tidyverse)
dat = readRDS("data/Created/processed_tiling_newC0.rds")
y = dat$C0_new
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
y_test = dat_test$C0_new

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

ps2 <- paste0("X", 1:49, "di")

contains_di = function(di_sub, di_sup) {di_sub %in% di_sup}
max_freq_at_dist = function(row, dist, di1, di2) {
  to_sum = sapply(row[ps2][which(sapply(row[ps2], contains_di, di_sup = di1))[which(sapply(row[ps2], contains_di, di_sup = di1)) < (50-dist)] + dist], 
                  contains_di, di_sup = di2)
  if (length(to_sum) == 0) {return(0)}
  return(sum(to_sum))
}

max_freq_at_dist_bidirectional = function(row, dist, di1, di2) {
  to_sum1 = sapply(row[ps2][which(sapply(row[ps2], contains_di, di_sup = di1))[which(sapply(row[ps2], contains_di, di_sup = di1)) < (50-dist)] + dist], 
                  contains_di, di_sup = di2)
  to_sum2 = sapply(row[ps2][which(sapply(row[ps2], contains_di, di_sup = di1))[which(sapply(row[ps2], contains_di, di_sup = di1)) > dist] - dist], 
                   contains_di, di_sup = di2)
  if ((length(to_sum1) == 0) & (length(to_sum2) == 0)) {return(0)}
  else if (length(to_sum1) == 0) {return(sum(to_sum2))}
  else if (length(to_sum2) == 0) {return(sum(to_sum1))}
  return(sum(to_sum1) + sum(to_sum2))
}



## AA/AT/TA/TT vs CC/CG/GC/GG

dist5_freq_AorT_AorTdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist5_freq_AorT_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist5_freq_CorG_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist5_freq_AorTdi = data.frame(dist5_freq_AorT_AorTdi, dist5_freq_AorT_CorGdi, 
                               dist5_freq_CorG_CorGdi)

dist10_freq_AorT_AorTdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist10_freq_AorT_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist10_freq_CorG_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist10_freq_AorTdi = data.frame(dist10_freq_AorT_AorTdi, dist10_freq_AorT_CorGdi, 
                                dist10_freq_CorG_CorGdi)

dist15_freq_AorT_AorTdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist15_freq_AorT_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist15_freq_CorG_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist15_freq_AorTdi = data.frame(dist15_freq_AorT_AorTdi, dist15_freq_AorT_CorGdi, 
                                dist15_freq_CorG_CorGdi)

dist20_freq_AorT_AorTdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist20_freq_AorT_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist20_freq_CorG_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist20_freq_AorTdi = data.frame(dist20_freq_AorT_AorTdi, dist20_freq_AorT_CorGdi, 
                                dist20_freq_CorG_CorGdi)

dist25_freq_AorT_AorTdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist25_freq_AorT_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist25_freq_CorG_CorGdi = apply(dat, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist25_freq_AorTdi = data.frame(dist25_freq_AorT_AorTdi, dist25_freq_AorT_CorGdi, 
                                dist25_freq_CorG_CorGdi)

dist_freq_AorTdi = data.frame(map(paste0("dist", c(5, 10, 15, 20, 25), "_freq_AorTdi"), get))
saveRDS(dist_freq_AorTdi, "data/Created/tiling_dist_freq_AorTdi.rds")






dist5_freq_AorT_AorTdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist5_freq_AorT_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist5_freq_CorG_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=5, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist5_freq_AorTdi_test = data.frame(dist5_freq_AorT_AorTdi_test, dist5_freq_AorT_CorGdi_test, 
                                    dist5_freq_CorG_CorGdi_test)

dist10_freq_AorT_AorTdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist10_freq_AorT_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist10_freq_CorG_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=10, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist10_freq_AorTdi_test = data.frame(dist10_freq_AorT_AorTdi_test, dist10_freq_AorT_CorGdi_test, 
                                     dist10_freq_CorG_CorGdi_test)

dist15_freq_AorT_AorTdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist15_freq_AorT_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist15_freq_CorG_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=15, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist15_freq_AorTdi_test = data.frame(dist15_freq_AorT_AorTdi_test, dist15_freq_AorT_CorGdi_test, 
                                     dist15_freq_CorG_CorGdi_test)

dist20_freq_AorT_AorTdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist20_freq_AorT_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist20_freq_CorG_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=20, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist20_freq_AorTdi_test = data.frame(dist20_freq_AorT_AorTdi_test, dist20_freq_AorT_CorGdi_test, 
                                     dist20_freq_CorG_CorGdi_test)

dist25_freq_AorT_AorTdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("AA", "AT", "TA", "TT"), di2=c("AA", "AT", "TA", "TT"))
dist25_freq_AorT_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("AA", "AT", "TA", "TT"), di2=c("CC", "CG", "GC", "GG"))
dist25_freq_CorG_CorGdi_test = apply(dat_test, 1, max_freq_at_dist_bidirectional, dist=25, di1=c("CC", "CG", "GC", "GG"), di2=c("CC", "CG", "GC", "GG"))

dist25_freq_AorTdi_test = data.frame(dist25_freq_AorT_AorTdi_test, dist25_freq_AorT_CorGdi_test, 
                                     dist25_freq_CorG_CorGdi_test)

dist_freq_AorTdi_test = data.frame(map(paste0("dist", c(5, 10, 15, 20, 25), "_freq_AorTdi_test"), get))
colnames(dist_freq_AorTdi_test) = colnames(dist_freq_AorTdi)
saveRDS(dist_freq_AorTdi_test, "data/Created/tiling_dist_freq_AorTdi_test.rds")







## Individual dinucleotides vs AA/AT/TA/TT or CC/CG/GC/GG (bidirectional)

max_freq_bidirectional = function(row, di1, di2) {
  di2_string = di2
  if(identical(di2, c("AA", "AT", "TA", "TT"))) {
    di2_string = "AorTdi"
  }
  else if(identical(di2, c("CC", "CG", "GC", "GG"))) {
    di2_string = "CorGdi"
  }
  ret = c(assign(paste0("dist5_freq_", di1, "_", di2_string), max_freq_at_dist_bidirectional(row, 5, di1, di2)),
          assign(paste0("dist10_freq_", di1, "_", di2_string), max_freq_at_dist_bidirectional(row, 10, di1, di2)),
          assign(paste0("dist15_freq_", di1, "_", di2_string), max_freq_at_dist_bidirectional(row, 15, di1, di2)),
          assign(paste0("dist20_freq_", di1, "_", di2_string), max_freq_at_dist_bidirectional(row, 20, di1, di2)),
          assign(paste0("dist25_freq_", di1, "_", di2_string), max_freq_at_dist_bidirectional(row, 25, di1, di2)))
  return(ret)
}

dist_freq_AA_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AA_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AorTdi_bid")
saveRDS(dist_freq_AA_AorTdi_bid, "data/Created/tiling_dist_freq_AA_AorTdi_bid.rds")
dist_freq_AA_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AA_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CorGdi_bid")
saveRDS(dist_freq_AA_CorGdi_bid, "data/Created/tiling_dist_freq_AA_CorGdi_bid.rds")

dist_freq_AC_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AC_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AorTdi_bid")
saveRDS(dist_freq_AC_AorTdi_bid, "data/Created/tiling_dist_freq_AC_AorTdi_bid.rds")
dist_freq_AC_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AC_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CorGdi_bid")
saveRDS(dist_freq_AC_CorGdi_bid, "data/Created/tiling_dist_freq_AC_CorGdi_bid.rds")

dist_freq_AG_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AG_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AorTdi_bid")
saveRDS(dist_freq_AG_AorTdi_bid, "data/Created/tiling_dist_freq_AG_AorTdi_bid.rds")
dist_freq_AG_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AG_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CorGdi_bid")
saveRDS(dist_freq_AG_CorGdi_bid, "data/Created/tiling_dist_freq_AG_CorGdi_bid.rds")

dist_freq_AT_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AT_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_AorTdi_bid")
saveRDS(dist_freq_AT_AorTdi_bid, "data/Created/tiling_dist_freq_AT_AorTdi_bid.rds")
dist_freq_AT_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AT_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CorGdi_bid")
saveRDS(dist_freq_AT_CorGdi_bid, "data/Created/tiling_dist_freq_AT_CorGdi_bid.rds")

dist_freq_CA_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CA_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_AorTdi_bid")
saveRDS(dist_freq_CA_AorTdi_bid, "data/Created/tiling_dist_freq_CA_AorTdi_bid.rds")
dist_freq_CA_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CA_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CorGdi_bid")
saveRDS(dist_freq_CA_CorGdi_bid, "data/Created/tiling_dist_freq_CA_CorGdi_bid.rds")

dist_freq_CC_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CC_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_AorTdi_bid")
saveRDS(dist_freq_CC_AorTdi_bid, "data/Created/tiling_dist_freq_CC_AorTdi_bid.rds")
dist_freq_CC_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CC_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CorGdi_bid")
saveRDS(dist_freq_CC_CorGdi_bid, "data/Created/tiling_dist_freq_CC_CorGdi_bid.rds")

dist_freq_CG_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CG_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_AorTdi_bid")
saveRDS(dist_freq_CG_AorTdi_bid, "data/Created/tiling_dist_freq_CG_AorTdi_bid.rds")
dist_freq_CG_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CG_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CorGdi_bid")
saveRDS(dist_freq_CG_CorGdi_bid, "data/Created/tiling_dist_freq_CG_CorGdi_bid.rds")

dist_freq_CT_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CT_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_AorTdi_bid")
saveRDS(dist_freq_CT_AorTdi_bid, "data/Created/tiling_dist_freq_CT_AorTdi_bid.rds")
dist_freq_CT_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CT_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_CorGdi_bid")
saveRDS(dist_freq_CT_CorGdi_bid, "data/Created/tiling_dist_freq_CT_CorGdi_bid.rds")

dist_freq_GA_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GA_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_AorTdi_bid")
saveRDS(dist_freq_GA_AorTdi_bid, "data/Created/tiling_dist_freq_GA_AorTdi_bid.rds")
dist_freq_GA_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GA_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_CorGdi_bid")
saveRDS(dist_freq_GA_CorGdi_bid, "data/Created/tiling_dist_freq_GA_CorGdi_bid.rds")

dist_freq_GC_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GC_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_AorTdi_bid")
saveRDS(dist_freq_GC_AorTdi_bid, "data/Created/tiling_dist_freq_GC_AorTdi_bid.rds")
dist_freq_GC_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GC_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_CorGdi_bid")
saveRDS(dist_freq_GC_CorGdi_bid, "data/Created/tiling_dist_freq_GC_CorGdi_bid.rds")

dist_freq_GG_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GG_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_AorTdi_bid")
saveRDS(dist_freq_GG_AorTdi_bid, "data/Created/tiling_dist_freq_GG_AorTdi_bid.rds")
dist_freq_GG_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GG_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_CorGdi_bid")
saveRDS(dist_freq_GG_CorGdi_bid, "data/Created/tiling_dist_freq_GG_CorGdi_bid.rds")

dist_freq_GT_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GT_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_AorTdi_bid")
saveRDS(dist_freq_GT_AorTdi_bid, "data/Created/tiling_dist_freq_GT_AorTdi_bid.rds")
dist_freq_GT_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GT_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_CorGdi_bid")
saveRDS(dist_freq_GT_CorGdi_bid, "data/Created/tiling_dist_freq_GT_CorGdi_bid.rds")

dist_freq_TA_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TA_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_AorTdi_bid")
saveRDS(dist_freq_TA_AorTdi_bid, "data/Created/tiling_dist_freq_TA_AorTdi_bid.rds")
dist_freq_TA_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TA_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_CorGdi_bid")
saveRDS(dist_freq_TA_CorGdi_bid, "data/Created/tiling_dist_freq_TA_CorGdi_bid.rds")

dist_freq_TC_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TC_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_AorTdi_bid")
saveRDS(dist_freq_TC_AorTdi_bid, "data/Created/tiling_dist_freq_TC_AorTdi_bid.rds")
dist_freq_TC_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TC_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_CorGdi_bid")
saveRDS(dist_freq_TC_CorGdi_bid, "data/Created/tiling_dist_freq_TC_CorGdi_bid.rds")

dist_freq_TG_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TG_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_AorTdi_bid")
saveRDS(dist_freq_TG_AorTdi_bid, "data/Created/tiling_dist_freq_TG_AorTdi_bid.rds")
dist_freq_TG_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TG_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_CorGdi_bid")
saveRDS(dist_freq_TG_CorGdi_bid, "data/Created/tiling_dist_freq_TG_CorGdi_bid.rds")

dist_freq_TT_AorTdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TT_AorTdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_AorTdi_bid")
saveRDS(dist_freq_TT_AorTdi_bid, "data/Created/tiling_dist_freq_TT_AorTdi_bid.rds")
dist_freq_TT_CorGdi_bid = t(apply(dat, 1, max_freq_bidirectional, di1="TT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TT_CorGdi_bid) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_CorGdi_bid")
saveRDS(dist_freq_TT_CorGdi_bid, "data/Created/tiling_dist_freq_TT_CorGdi_bid.rds")


dist_freq_di_AorTdi_bid = data.frame(map(paste0("dist_freq_", dinucleotides, "_AorTdi_bid"), get),
                                     map(paste0("dist_freq_", dinucleotides, "_CorGdi_bid"), get))
saveRDS(dist_freq_di_AorTdi_bid, "data/Created/tiling_dist_freq_di_AorTdi_bid.rds")




dist_freq_AA_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AA_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AorTdi_bid")
saveRDS(dist_freq_AA_AorTdi_bid_test, "data/Created/tiling_dist_freq_AA_AorTdi_bid_test.rds")
dist_freq_AA_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AA_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CorGdi_bid")
saveRDS(dist_freq_AA_CorGdi_bid_test, "data/Created/tiling_dist_freq_AA_CorGdi_bid_test.rds")

dist_freq_AC_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AC_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AorTdi_bid")
saveRDS(dist_freq_AC_AorTdi_bid_test, "data/Created/tiling_dist_freq_AC_AorTdi_bid_test.rds")
dist_freq_AC_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AC_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CorGdi_bid")
saveRDS(dist_freq_AC_CorGdi_bid_test, "data/Created/tiling_dist_freq_AC_CorGdi_bid_test.rds")

dist_freq_AG_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AG_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AorTdi_bid")
saveRDS(dist_freq_AG_AorTdi_bid_test, "data/Created/tiling_dist_freq_AG_AorTdi_bid_test.rds")
dist_freq_AG_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AG_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CorGdi_bid")
saveRDS(dist_freq_AG_CorGdi_bid_test, "data/Created/tiling_dist_freq_AG_CorGdi_bid_test.rds")

dist_freq_AT_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_AT_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_AorTdi_bid")
saveRDS(dist_freq_AT_AorTdi_bid_test, "data/Created/tiling_dist_freq_AT_AorTdi_bid_test.rds")
dist_freq_AT_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_AT_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CorGdi_bid")
saveRDS(dist_freq_AT_CorGdi_bid_test, "data/Created/tiling_dist_freq_AT_CorGdi_bid_test.rds")

dist_freq_CA_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CA_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_AorTdi_bid")
saveRDS(dist_freq_CA_AorTdi_bid_test, "data/Created/tiling_dist_freq_CA_AorTdi_bid_test.rds")
dist_freq_CA_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CA_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CorGdi_bid")
saveRDS(dist_freq_CA_CorGdi_bid_test, "data/Created/tiling_dist_freq_CA_CorGdi_bid_test.rds")

dist_freq_CC_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CC_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_AorTdi_bid")
saveRDS(dist_freq_CC_AorTdi_bid_test, "data/Created/tiling_dist_freq_CC_AorTdi_bid_test.rds")
dist_freq_CC_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CC_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CorGdi_bid")
saveRDS(dist_freq_CC_CorGdi_bid_test, "data/Created/tiling_dist_freq_CC_CorGdi_bid_test.rds")

dist_freq_CG_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CG_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_AorTdi_bid")
saveRDS(dist_freq_CG_AorTdi_bid_test, "data/Created/tiling_dist_freq_CG_AorTdi_bid_test.rds")
dist_freq_CG_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CG_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CorGdi_bid")
saveRDS(dist_freq_CG_CorGdi_bid_test, "data/Created/tiling_dist_freq_CG_CorGdi_bid_test.rds")

dist_freq_CT_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_CT_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_AorTdi_bid")
saveRDS(dist_freq_CT_AorTdi_bid_test, "data/Created/tiling_dist_freq_CT_AorTdi_bid_test.rds")
dist_freq_CT_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_CT_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_CorGdi_bid")
saveRDS(dist_freq_CT_CorGdi_bid_test, "data/Created/tiling_dist_freq_CT_CorGdi_bid_test.rds")

dist_freq_GA_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GA_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_AorTdi_bid")
saveRDS(dist_freq_GA_AorTdi_bid_test, "data/Created/tiling_dist_freq_GA_AorTdi_bid_test.rds")
dist_freq_GA_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GA_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_CorGdi_bid")
saveRDS(dist_freq_GA_CorGdi_bid_test, "data/Created/tiling_dist_freq_GA_CorGdi_bid_test.rds")

dist_freq_GC_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GC_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_AorTdi_bid")
saveRDS(dist_freq_GC_AorTdi_bid_test, "data/Created/tiling_dist_freq_GC_AorTdi_bid_test.rds")
dist_freq_GC_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GC_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_CorGdi_bid")
saveRDS(dist_freq_GC_CorGdi_bid_test, "data/Created/tiling_dist_freq_GC_CorGdi_bid_test.rds")

dist_freq_GG_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GG_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_AorTdi_bid")
saveRDS(dist_freq_GG_AorTdi_bid_test, "data/Created/tiling_dist_freq_GG_AorTdi_bid_test.rds")
dist_freq_GG_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GG_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_CorGdi_bid")
saveRDS(dist_freq_GG_CorGdi_bid_test, "data/Created/tiling_dist_freq_GG_CorGdi_bid_test.rds")

dist_freq_GT_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_GT_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_AorTdi_bid")
saveRDS(dist_freq_GT_AorTdi_bid_test, "data/Created/tiling_dist_freq_GT_AorTdi_bid_test.rds")
dist_freq_GT_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_GT_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_CorGdi_bid")
saveRDS(dist_freq_GT_CorGdi_bid_test, "data/Created/tiling_dist_freq_GT_CorGdi_bid_test.rds")

dist_freq_TA_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TA_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_AorTdi_bid")
saveRDS(dist_freq_TA_AorTdi_bid_test, "data/Created/tiling_dist_freq_TA_AorTdi_bid_test.rds")
dist_freq_TA_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TA_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_CorGdi_bid")
saveRDS(dist_freq_TA_CorGdi_bid_test, "data/Created/tiling_dist_freq_TA_CorGdi_bid_test.rds")

dist_freq_TC_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TC", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TC_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_AorTdi_bid")
saveRDS(dist_freq_TC_AorTdi_bid_test, "data/Created/tiling_dist_freq_TC_AorTdi_bid_test.rds")
dist_freq_TC_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TC", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TC_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_CorGdi_bid")
saveRDS(dist_freq_TC_CorGdi_bid_test, "data/Created/tiling_dist_freq_TC_CorGdi_bid_test.rds")

dist_freq_TG_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TG", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TG_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_AorTdi_bid")
saveRDS(dist_freq_TG_AorTdi_bid_test, "data/Created/tiling_dist_freq_TG_AorTdi_bid_test.rds")
dist_freq_TG_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TG", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TG_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_CorGdi_bid")
saveRDS(dist_freq_TG_CorGdi_bid_test, "data/Created/tiling_dist_freq_TG_CorGdi_bid_test.rds")

dist_freq_TT_AorTdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TT", di2=c("AA", "AT", "TA", "TT")))
colnames(dist_freq_TT_AorTdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_AorTdi_bid")
saveRDS(dist_freq_TT_AorTdi_bid_test, "data/Created/tiling_dist_freq_TT_AorTdi_bid_test.rds")
dist_freq_TT_CorGdi_bid_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TT", di2=c("CC", "CG", "GC", "GG")))
colnames(dist_freq_TT_CorGdi_bid_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_CorGdi_bid")
saveRDS(dist_freq_TT_CorGdi_bid_test, "data/Created/tiling_dist_freq_TT_CorGdi_bid_test.rds")


dist_freq_di_AorTdi_bid_test = data.frame(map(paste0("dist_freq_", dinucleotides, "_AorTdi_bid_test"), get),
                                          map(paste0("dist_freq_", dinucleotides, "_CorGdi_bid_test"), get))
saveRDS(dist_freq_di_AorTdi_bid_test, "data/Created/tiling_dist_freq_di_AorTdi_bid_test.rds")





## All individual dinucleotides

max_freq = function(row, di1, di2) {
  ret = c(assign(paste0("dist5_freq_", di1, "_", di2), max_freq_at_dist(row, 5, di1, di2)),
          assign(paste0("dist10_freq_", di1, "_", di2), max_freq_at_dist(row, 10, di1, di2)),
          assign(paste0("dist15_freq_", di1, "_", di2), max_freq_at_dist(row, 15, di1, di2)),
          assign(paste0("dist20_freq_", di1, "_", di2), max_freq_at_dist(row, 20, di1, di2)),
          assign(paste0("dist25_freq_", di1, "_", di2), max_freq_at_dist(row, 25, di1, di2)))
  return(ret)
}



dist_freq_AA_AA = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="AA"))
colnames(dist_freq_AA_AA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AA")
dist_freq_AA_AC = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="AC"))
colnames(dist_freq_AA_AC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AC")
dist_freq_AA_AG = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="AG"))
colnames(dist_freq_AA_AG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AG")
dist_freq_AA_AT = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="AT"))
colnames(dist_freq_AA_AT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AT")
dist_freq_AA_CA = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="CA"))
colnames(dist_freq_AA_CA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CA")
dist_freq_AA_CC = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="CC"))
colnames(dist_freq_AA_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CC")
dist_freq_AA_CG = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="CG"))
colnames(dist_freq_AA_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CG")
dist_freq_AA_CT = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="CT"))
colnames(dist_freq_AA_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CT")
dist_freq_AA_GA = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="GA"))
colnames(dist_freq_AA_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GA")
dist_freq_AA_GC = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="GC"))
colnames(dist_freq_AA_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GC")
dist_freq_AA_GG = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="GG"))
colnames(dist_freq_AA_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GG")
dist_freq_AA_GT = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="GT"))
colnames(dist_freq_AA_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GT")
dist_freq_AA_TA = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="TA"))
colnames(dist_freq_AA_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TA")
dist_freq_AA_TC = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="TC"))
colnames(dist_freq_AA_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TC")
dist_freq_AA_TG = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="TG"))
colnames(dist_freq_AA_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TG")
dist_freq_AA_TT = t(apply(dat, 1, max_freq_bidirectional, di1="AA", di2="TT"))
colnames(dist_freq_AA_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TT")

dist_freq_AA = cbind(dist_freq_AA_AA, dist_freq_AA_AC, dist_freq_AA_AG, dist_freq_AA_AT,
                     dist_freq_AA_CA, dist_freq_AA_CC, dist_freq_AA_CG, dist_freq_AA_CT,
                     dist_freq_AA_GA, dist_freq_AA_GC, dist_freq_AA_GG, dist_freq_AA_GT,
                     dist_freq_AA_TA, dist_freq_AA_TC, dist_freq_AA_TG, dist_freq_AA_TT)
saveRDS(dist_freq_AA, "data/Created/tiling_dist_freq_AA.rds")

dist_freq_AC_AC = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="AC"))
colnames(dist_freq_AC_AC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AC")
dist_freq_AC_AG = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="AG"))
colnames(dist_freq_AC_AG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AG")
dist_freq_AC_AT = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="AT"))
colnames(dist_freq_AC_AT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AT")
dist_freq_AC_CA = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="CA"))
colnames(dist_freq_AC_CA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CA")
dist_freq_AC_CC = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="CC"))
colnames(dist_freq_AC_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CC")
dist_freq_AC_CG = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="CG"))
colnames(dist_freq_AC_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CG")
dist_freq_AC_CT = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="CT"))
colnames(dist_freq_AC_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CT")
dist_freq_AC_GA = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="GA"))
colnames(dist_freq_AC_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GA")
dist_freq_AC_GC = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="GC"))
colnames(dist_freq_AC_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GC")
dist_freq_AC_GG = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="GG"))
colnames(dist_freq_AC_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GG")
dist_freq_AC_GT = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="GT"))
colnames(dist_freq_AC_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GT")
dist_freq_AC_TA = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="TA"))
colnames(dist_freq_AC_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TA")
dist_freq_AC_TC = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="TC"))
colnames(dist_freq_AC_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TC")
dist_freq_AC_TG = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="TG"))
colnames(dist_freq_AC_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TG")
dist_freq_AC_TT = t(apply(dat, 1, max_freq_bidirectional, di1="AC", di2="TT"))
colnames(dist_freq_AC_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TT")

dist_freq_AC = cbind(dist_freq_AC_AC, dist_freq_AC_AG, dist_freq_AC_AT,
                     dist_freq_AC_CA, dist_freq_AC_CC, dist_freq_AC_CG, dist_freq_AC_CT,
                     dist_freq_AC_GA, dist_freq_AC_GC, dist_freq_AC_GG, dist_freq_AC_GT,
                     dist_freq_AC_TA, dist_freq_AC_TC, dist_freq_AC_TG, dist_freq_AC_TT)
saveRDS(dist_freq_AC, "data/Created/tiling_dist_freq_AC.rds")


dist_freq_AG_AG = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="AG"))
colnames(dist_freq_AG_AG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AG")
dist_freq_AG_AT = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="AT"))
colnames(dist_freq_AG_AT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AT")
dist_freq_AG_CA = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="CA"))
colnames(dist_freq_AG_CA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CA")
dist_freq_AG_CC = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="CC"))
colnames(dist_freq_AG_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CC")
dist_freq_AG_CG = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="CG"))
colnames(dist_freq_AG_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CG")
dist_freq_AG_CT = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="CT"))
colnames(dist_freq_AG_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CT")
dist_freq_AG_GA = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="GA"))
colnames(dist_freq_AG_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GA")
dist_freq_AG_GC = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="GC"))
colnames(dist_freq_AG_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GC")
dist_freq_AG_GG = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="GG"))
colnames(dist_freq_AG_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GG")
dist_freq_AG_GT = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="GT"))
colnames(dist_freq_AG_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GT")
dist_freq_AG_TA = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="TA"))
colnames(dist_freq_AG_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TA")
dist_freq_AG_TC = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="TC"))
colnames(dist_freq_AG_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TC")
dist_freq_AG_TG = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="TG"))
colnames(dist_freq_AG_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TG")
dist_freq_AG_TT = t(apply(dat, 1, max_freq_bidirectional, di1="AG", di2="TT"))
colnames(dist_freq_AG_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TT")

dist_freq_AG = cbind(dist_freq_AG_AG, dist_freq_AG_AT,
                     dist_freq_AG_CA, dist_freq_AG_CC, dist_freq_AG_CG, dist_freq_AG_CT,
                     dist_freq_AG_GA, dist_freq_AG_GC, dist_freq_AG_GG, dist_freq_AG_GT,
                     dist_freq_AG_TA, dist_freq_AG_TC, dist_freq_AG_TG, dist_freq_AG_TT)
saveRDS(dist_freq_AG, "data/Created/tiling_dist_freq_AG.rds")



dist_freq_AT_AT = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="AT"))
colnames(dist_freq_AT_AT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_AT")
dist_freq_AT_CA = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="CA"))
colnames(dist_freq_AT_CA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CA")
dist_freq_AT_CC = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="CC"))
colnames(dist_freq_AT_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CC")
dist_freq_AT_CG = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="CG"))
colnames(dist_freq_AT_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CG")
dist_freq_AT_CT = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="CT"))
colnames(dist_freq_AT_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CT")
dist_freq_AT_GA = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="GA"))
colnames(dist_freq_AT_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GA")
dist_freq_AT_GC = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="GC"))
colnames(dist_freq_AT_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GC")
dist_freq_AT_GG = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="GG"))
colnames(dist_freq_AT_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GG")
dist_freq_AT_GT = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="GT"))
colnames(dist_freq_AT_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GT")
dist_freq_AT_TA = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="TA"))
colnames(dist_freq_AT_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TA")
dist_freq_AT_TC = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="TC"))
colnames(dist_freq_AT_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TC")
dist_freq_AT_TG = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="TG"))
colnames(dist_freq_AT_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TG")
dist_freq_AT_TT = t(apply(dat, 1, max_freq_bidirectional, di1="AT", di2="TT"))
colnames(dist_freq_AT_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TT")

dist_freq_AT = cbind(dist_freq_AT_AT,
                     dist_freq_AT_CA, dist_freq_AT_CC, dist_freq_AT_CG, dist_freq_AT_CT,
                     dist_freq_AT_GA, dist_freq_AT_GC, dist_freq_AT_GG, dist_freq_AT_GT,
                     dist_freq_AT_TA, dist_freq_AT_TC, dist_freq_AT_TG, dist_freq_AT_TT)
saveRDS(dist_freq_AT, "data/Created/tiling_dist_freq_AT.rds")


dist_freq_CA_CA = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="CA"))
colnames(dist_freq_CA_CA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CA")
dist_freq_CA_CC = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="CC"))
colnames(dist_freq_CA_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CC")
dist_freq_CA_CG = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="CG"))
colnames(dist_freq_CA_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CG")
dist_freq_CA_CT = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="CT"))
colnames(dist_freq_CA_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CT")
dist_freq_CA_GA = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="GA"))
colnames(dist_freq_CA_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GA")
dist_freq_CA_GC = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="GC"))
colnames(dist_freq_CA_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GC")
dist_freq_CA_GG = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="GG"))
colnames(dist_freq_CA_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GG")
dist_freq_CA_GT = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="GT"))
colnames(dist_freq_CA_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GT")
dist_freq_CA_TA = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="TA"))
colnames(dist_freq_CA_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TA")
dist_freq_CA_TC = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="TC"))
colnames(dist_freq_CA_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TC")
dist_freq_CA_TG = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="TG"))
colnames(dist_freq_CA_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TG")
dist_freq_CA_TT = t(apply(dat, 1, max_freq_bidirectional, di1="CA", di2="TT"))
colnames(dist_freq_CA_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TT")

dist_freq_CA = cbind(dist_freq_CA_CA, dist_freq_CA_CC, dist_freq_CA_CG, dist_freq_CA_CT,
                     dist_freq_CA_GA, dist_freq_CA_GC, dist_freq_CA_GG, dist_freq_CA_GT,
                     dist_freq_CA_TA, dist_freq_CA_TC, dist_freq_CA_TG, dist_freq_CA_TT)
saveRDS(dist_freq_CA, "data/Created/tiling_dist_freq_CA.rds")



dist_freq_CC_CC = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="CC"))
colnames(dist_freq_CC_CC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CC")
dist_freq_CC_CG = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="CG"))
colnames(dist_freq_CC_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CG")
dist_freq_CC_CT = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="CT"))
colnames(dist_freq_CC_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CT")
dist_freq_CC_GA = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="GA"))
colnames(dist_freq_CC_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GA")
dist_freq_CC_GC = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="GC"))
colnames(dist_freq_CC_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GC")
dist_freq_CC_GG = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="GG"))
colnames(dist_freq_CC_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GG")
dist_freq_CC_GT = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="GT"))
colnames(dist_freq_CC_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GT")
dist_freq_CC_TA = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="TA"))
colnames(dist_freq_CC_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TA")
dist_freq_CC_TC = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="TC"))
colnames(dist_freq_CC_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TC")
dist_freq_CC_TG = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="TG"))
colnames(dist_freq_CC_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TG")
dist_freq_CC_TT = t(apply(dat, 1, max_freq_bidirectional, di1="CC", di2="TT"))
colnames(dist_freq_CC_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TT")

dist_freq_CC = cbind(dist_freq_CC_CC, dist_freq_CC_CG, dist_freq_CC_CT,
                     dist_freq_CC_GA, dist_freq_CC_GC, dist_freq_CC_GG, dist_freq_CC_GT,
                     dist_freq_CC_TA, dist_freq_CC_TC, dist_freq_CC_TG, dist_freq_CC_TT)
saveRDS(dist_freq_CC, "data/Created/tiling_dist_freq_CC.rds")



dist_freq_CG_CG = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="CG"))
colnames(dist_freq_CG_CG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CG")
dist_freq_CG_CT = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="CT"))
colnames(dist_freq_CG_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CT")
dist_freq_CG_GA = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="GA"))
colnames(dist_freq_CG_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GA")
dist_freq_CG_GC = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="GC"))
colnames(dist_freq_CG_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GC")
dist_freq_CG_GG = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="GG"))
colnames(dist_freq_CG_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GG")
dist_freq_CG_GT = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="GT"))
colnames(dist_freq_CG_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GT")
dist_freq_CG_TA = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="TA"))
colnames(dist_freq_CG_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TA")
dist_freq_CG_TC = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="TC"))
colnames(dist_freq_CG_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TC")
dist_freq_CG_TG = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="TG"))
colnames(dist_freq_CG_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TG")
dist_freq_CG_TT = t(apply(dat, 1, max_freq_bidirectional, di1="CG", di2="TT"))
colnames(dist_freq_CG_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TT")

dist_freq_CG = cbind(dist_freq_CG_CG, dist_freq_CG_CT,
                     dist_freq_CG_GA, dist_freq_CG_GC, dist_freq_CG_GG, dist_freq_CG_GT,
                     dist_freq_CG_TA, dist_freq_CG_TC, dist_freq_CG_TG, dist_freq_CG_TT)
saveRDS(dist_freq_CG, "data/Created/tiling_dist_freq_CG.rds")



dist_freq_CT_CT = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="CT"))
colnames(dist_freq_CT_CT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_CT")
dist_freq_CT_GA = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="GA"))
colnames(dist_freq_CT_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GA")
dist_freq_CT_GC = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="GC"))
colnames(dist_freq_CT_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GC")
dist_freq_CT_GG = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="GG"))
colnames(dist_freq_CT_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GG")
dist_freq_CT_GT = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="GT"))
colnames(dist_freq_CT_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GT")
dist_freq_CT_TA = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="TA"))
colnames(dist_freq_CT_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TA")
dist_freq_CT_TC = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="TC"))
colnames(dist_freq_CT_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TC")
dist_freq_CT_TG = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="TG"))
colnames(dist_freq_CT_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TG")
dist_freq_CT_TT = t(apply(dat, 1, max_freq_bidirectional, di1="CT", di2="TT"))
colnames(dist_freq_CT_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TT")

dist_freq_CT = cbind(dist_freq_CT_CT,
                     dist_freq_CT_GA, dist_freq_CT_GC, dist_freq_CT_GG, dist_freq_CT_GT,
                     dist_freq_CT_TA, dist_freq_CT_TC, dist_freq_CT_TG, dist_freq_CT_TT)
saveRDS(dist_freq_CT, "data/Created/tiling_dist_freq_CT.rds")



dist_freq_GA_GA = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="GA"))
colnames(dist_freq_GA_GA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GA")
dist_freq_GA_GC = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="GC"))
colnames(dist_freq_GA_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GC")
dist_freq_GA_GG = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="GG"))
colnames(dist_freq_GA_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GG")
dist_freq_GA_GT = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="GT"))
colnames(dist_freq_GA_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GT")
dist_freq_GA_TA = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="TA"))
colnames(dist_freq_GA_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TA")
dist_freq_GA_TC = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="TC"))
colnames(dist_freq_GA_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TC")
dist_freq_GA_TG = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="TG"))
colnames(dist_freq_GA_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TG")
dist_freq_GA_TT = t(apply(dat, 1, max_freq_bidirectional, di1="GA", di2="TT"))
colnames(dist_freq_GA_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TT")

dist_freq_GA = cbind(dist_freq_GA_GA, dist_freq_GA_GC, dist_freq_GA_GG, dist_freq_GA_GT,
                     dist_freq_GA_TA, dist_freq_GA_TC, dist_freq_GA_TG, dist_freq_GA_TT)
saveRDS(dist_freq_GA, "data/Created/tiling_dist_freq_GA.rds")



dist_freq_GC_GC = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="GC"))
colnames(dist_freq_GC_GC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GC")
dist_freq_GC_GG = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="GG"))
colnames(dist_freq_GC_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GG")
dist_freq_GC_GT = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="GT"))
colnames(dist_freq_GC_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GT")
dist_freq_GC_TA = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="TA"))
colnames(dist_freq_GC_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TA")
dist_freq_GC_TC = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="TC"))
colnames(dist_freq_GC_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TC")
dist_freq_GC_TG = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="TG"))
colnames(dist_freq_GC_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TG")
dist_freq_GC_TT = t(apply(dat, 1, max_freq_bidirectional, di1="GC", di2="TT"))
colnames(dist_freq_GC_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TT")

dist_freq_GC = cbind(dist_freq_GC_GC, dist_freq_GC_GG, dist_freq_GC_GT,
                     dist_freq_GC_TA, dist_freq_GC_TC, dist_freq_GC_TG, dist_freq_GC_TT)
saveRDS(dist_freq_GC, "data/Created/tiling_dist_freq_GC.rds")



dist_freq_GG_GG = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="GG"))
colnames(dist_freq_GG_GG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_GG")
dist_freq_GG_GT = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="GT"))
colnames(dist_freq_GG_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_GT")
dist_freq_GG_TA = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="TA"))
colnames(dist_freq_GG_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TA")
dist_freq_GG_TC = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="TC"))
colnames(dist_freq_GG_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TC")
dist_freq_GG_TG = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="TG"))
colnames(dist_freq_GG_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TG")
dist_freq_GG_TT = t(apply(dat, 1, max_freq_bidirectional, di1="GG", di2="TT"))
colnames(dist_freq_GG_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TT")

dist_freq_GG = cbind(dist_freq_GG_GG, dist_freq_GG_GT,
                     dist_freq_GG_TA, dist_freq_GG_TC, dist_freq_GG_TG, dist_freq_GG_TT)
saveRDS(dist_freq_GG, "data/Created/tiling_dist_freq_GG.rds")



dist_freq_GT_GT = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2="GT"))
colnames(dist_freq_GT_GT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_GT")
dist_freq_GT_TA = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2="TA"))
colnames(dist_freq_GT_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TA")
dist_freq_GT_TC = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2="TC"))
colnames(dist_freq_GT_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TC")
dist_freq_GT_TG = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2="TG"))
colnames(dist_freq_GT_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TG")
dist_freq_GT_TT = t(apply(dat, 1, max_freq_bidirectional, di1="GT", di2="TT"))
colnames(dist_freq_GT_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TT")

dist_freq_GT = cbind(dist_freq_GT_GT,
                     dist_freq_GT_TA, dist_freq_GT_TC, dist_freq_GT_TG, dist_freq_GT_TT)
saveRDS(dist_freq_GT, "data/Created/tiling_dist_freq_GT.rds")



dist_freq_TA_TA = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2="TA"))
colnames(dist_freq_TA_TA) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TA")
dist_freq_TA_TC = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2="TC"))
colnames(dist_freq_TA_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TC")
dist_freq_TA_TG = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2="TG"))
colnames(dist_freq_TA_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TG")
dist_freq_TA_TT = t(apply(dat, 1, max_freq_bidirectional, di1="TA", di2="TT"))
colnames(dist_freq_TA_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TT")

dist_freq_TA = cbind(dist_freq_TA_TA, dist_freq_TA_TC, dist_freq_TA_TG, dist_freq_TA_TT)
saveRDS(dist_freq_TA, "data/Created/tiling_dist_freq_TA.rds")



dist_freq_TC_TC = t(apply(dat, 1, max_freq_bidirectional, di1="TC", di2="TC"))
colnames(dist_freq_TC_TC) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TC")
dist_freq_TC_TG = t(apply(dat, 1, max_freq_bidirectional, di1="TC", di2="TG"))
colnames(dist_freq_TC_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TG")
dist_freq_TC_TT = t(apply(dat, 1, max_freq_bidirectional, di1="TC", di2="TT"))
colnames(dist_freq_TC_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TT")

dist_freq_TC = cbind(dist_freq_TC_TC, dist_freq_TC_TG, dist_freq_TC_TT)
saveRDS(dist_freq_TC, "data/Created/tiling_dist_freq_TC.rds")




dist_freq_TG_TG = t(apply(dat, 1, max_freq_bidirectional, di1="TG", di2="TG"))
colnames(dist_freq_TG_TG) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_TG")
dist_freq_TG_TT = t(apply(dat, 1, max_freq_bidirectional, di1="TG", di2="TT"))
colnames(dist_freq_TG_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_TT")

dist_freq_TG = cbind(dist_freq_TG_TG, dist_freq_TG_TT)
saveRDS(dist_freq_TG, "data/Created/tiling_dist_freq_TG.rds")



dist_freq_TT_TT = t(apply(dat, 1, max_freq_bidirectional, di1="TT", di2="TT"))
colnames(dist_freq_TT_TT) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_TT")

dist_freq_TT = cbind(dist_freq_TT_TT)
saveRDS(dist_freq_TT, "data/Created/tiling_dist_freq_TT.rds")


dist_freq_di = data.frame(map(paste0("dist_freq_", dinucleotides), get))
saveRDS(dist_freq_di, "data/Created/tiling_dist_freq_di.rds")




dist_freq_AA_AA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="AA"))
colnames(dist_freq_AA_AA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AA")
dist_freq_AA_AC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="AC"))
colnames(dist_freq_AA_AC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AC")
dist_freq_AA_AG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="AG"))
colnames(dist_freq_AA_AG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AG")
dist_freq_AA_AT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="AT"))
colnames(dist_freq_AA_AT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_AT")
dist_freq_AA_CA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="CA"))
colnames(dist_freq_AA_CA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CA")
dist_freq_AA_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="CC"))
colnames(dist_freq_AA_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CC")
dist_freq_AA_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="CG"))
colnames(dist_freq_AA_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CG")
dist_freq_AA_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="CT"))
colnames(dist_freq_AA_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_CT")
dist_freq_AA_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="GA"))
colnames(dist_freq_AA_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GA")
dist_freq_AA_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="GC"))
colnames(dist_freq_AA_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GC")
dist_freq_AA_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="GG"))
colnames(dist_freq_AA_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GG")
dist_freq_AA_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="GT"))
colnames(dist_freq_AA_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_GT")
dist_freq_AA_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="TA"))
colnames(dist_freq_AA_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TA")
dist_freq_AA_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="TC"))
colnames(dist_freq_AA_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TC")
dist_freq_AA_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="TG"))
colnames(dist_freq_AA_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TG")
dist_freq_AA_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AA", di2="TT"))
colnames(dist_freq_AA_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AA_TT")

dist_freq_AA_test = cbind(dist_freq_AA_AA_test, dist_freq_AA_AC_test, dist_freq_AA_AG_test, dist_freq_AA_AT_test,
                          dist_freq_AA_CA_test, dist_freq_AA_CC_test, dist_freq_AA_CG_test, dist_freq_AA_CT_test,
                          dist_freq_AA_GA_test, dist_freq_AA_GC_test, dist_freq_AA_GG_test, dist_freq_AA_GT_test,
                          dist_freq_AA_TA_test, dist_freq_AA_TC_test, dist_freq_AA_TG_test, dist_freq_AA_TT_test)
saveRDS(dist_freq_AA_test, "data/Created/tiling_dist_freq_AA_test.rds")

dist_freq_AC_AC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="AC"))
colnames(dist_freq_AC_AC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AC")
dist_freq_AC_AG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="AG"))
colnames(dist_freq_AC_AG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AG")
dist_freq_AC_AT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="AT"))
colnames(dist_freq_AC_AT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_AT")
dist_freq_AC_CA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="CA"))
colnames(dist_freq_AC_CA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CA")
dist_freq_AC_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="CC"))
colnames(dist_freq_AC_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CC")
dist_freq_AC_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="CG"))
colnames(dist_freq_AC_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CG")
dist_freq_AC_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="CT"))
colnames(dist_freq_AC_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_CT")
dist_freq_AC_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="GA"))
colnames(dist_freq_AC_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GA")
dist_freq_AC_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="GC"))
colnames(dist_freq_AC_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GC")
dist_freq_AC_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="GG"))
colnames(dist_freq_AC_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GG")
dist_freq_AC_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="GT"))
colnames(dist_freq_AC_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_GT")
dist_freq_AC_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="TA"))
colnames(dist_freq_AC_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TA")
dist_freq_AC_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="TC"))
colnames(dist_freq_AC_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TC")
dist_freq_AC_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="TG"))
colnames(dist_freq_AC_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TG")
dist_freq_AC_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AC", di2="TT"))
colnames(dist_freq_AC_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AC_TT")

dist_freq_AC_test = cbind(dist_freq_AC_AC_test, dist_freq_AC_AG_test, dist_freq_AC_AT_test,
                          dist_freq_AC_CA_test, dist_freq_AC_CC_test, dist_freq_AC_CG_test, dist_freq_AC_CT_test,
                          dist_freq_AC_GA_test, dist_freq_AC_GC_test, dist_freq_AC_GG_test, dist_freq_AC_GT_test,
                          dist_freq_AC_TA_test, dist_freq_AC_TC_test, dist_freq_AC_TG_test, dist_freq_AC_TT_test)
saveRDS(dist_freq_AC_test, "data/Created/tiling_dist_freq_AC_test.rds")


dist_freq_AG_AG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="AG"))
colnames(dist_freq_AG_AG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AG")
dist_freq_AG_AT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="AT"))
colnames(dist_freq_AG_AT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_AT")
dist_freq_AG_CA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="CA"))
colnames(dist_freq_AG_CA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CA")
dist_freq_AG_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="CC"))
colnames(dist_freq_AG_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CC")
dist_freq_AG_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="CG"))
colnames(dist_freq_AG_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CG")
dist_freq_AG_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="CT"))
colnames(dist_freq_AG_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_CT")
dist_freq_AG_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="GA"))
colnames(dist_freq_AG_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GA")
dist_freq_AG_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="GC"))
colnames(dist_freq_AG_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GC")
dist_freq_AG_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="GG"))
colnames(dist_freq_AG_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GG")
dist_freq_AG_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="GT"))
colnames(dist_freq_AG_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_GT")
dist_freq_AG_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="TA"))
colnames(dist_freq_AG_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TA")
dist_freq_AG_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="TC"))
colnames(dist_freq_AG_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TC")
dist_freq_AG_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="TG"))
colnames(dist_freq_AG_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TG")
dist_freq_AG_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AG", di2="TT"))
colnames(dist_freq_AG_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AG_TT")

dist_freq_AG_test = cbind(dist_freq_AG_AG_test, dist_freq_AG_AT_test,
                          dist_freq_AG_CA_test, dist_freq_AG_CC_test, dist_freq_AG_CG_test, dist_freq_AG_CT_test,
                          dist_freq_AG_GA_test, dist_freq_AG_GC_test, dist_freq_AG_GG_test, dist_freq_AG_GT_test,
                          dist_freq_AG_TA_test, dist_freq_AG_TC_test, dist_freq_AG_TG_test, dist_freq_AG_TT_test)
saveRDS(dist_freq_AG_test, "data/Created/tiling_dist_freq_AG_test.rds")



dist_freq_AT_AT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="AT"))
colnames(dist_freq_AT_AT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_AT")
dist_freq_AT_CA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="CA"))
colnames(dist_freq_AT_CA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CA")
dist_freq_AT_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="CC"))
colnames(dist_freq_AT_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CC")
dist_freq_AT_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="CG"))
colnames(dist_freq_AT_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CG")
dist_freq_AT_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="CT"))
colnames(dist_freq_AT_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_CT")
dist_freq_AT_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="GA"))
colnames(dist_freq_AT_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GA")
dist_freq_AT_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="GC"))
colnames(dist_freq_AT_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GC")
dist_freq_AT_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="GG"))
colnames(dist_freq_AT_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GG")
dist_freq_AT_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="GT"))
colnames(dist_freq_AT_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_GT")
dist_freq_AT_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="TA"))
colnames(dist_freq_AT_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TA")
dist_freq_AT_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="TC"))
colnames(dist_freq_AT_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TC")
dist_freq_AT_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="TG"))
colnames(dist_freq_AT_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TG")
dist_freq_AT_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="AT", di2="TT"))
colnames(dist_freq_AT_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_AT_TT")

dist_freq_AT_test = cbind(dist_freq_AT_AT_test,
                          dist_freq_AT_CA_test, dist_freq_AT_CC_test, dist_freq_AT_CG_test, dist_freq_AT_CT_test,
                          dist_freq_AT_GA_test, dist_freq_AT_GC_test, dist_freq_AT_GG_test, dist_freq_AT_GT_test,
                          dist_freq_AT_TA_test, dist_freq_AT_TC_test, dist_freq_AT_TG_test, dist_freq_AT_TT_test)
saveRDS(dist_freq_AT_test, "data/Created/tiling_dist_freq_AT_test.rds")


dist_freq_CA_CA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="CA"))
colnames(dist_freq_CA_CA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CA")
dist_freq_CA_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="CC"))
colnames(dist_freq_CA_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CC")
dist_freq_CA_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="CG"))
colnames(dist_freq_CA_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CG")
dist_freq_CA_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="CT"))
colnames(dist_freq_CA_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_CT")
dist_freq_CA_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="GA"))
colnames(dist_freq_CA_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GA")
dist_freq_CA_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="GC"))
colnames(dist_freq_CA_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GC")
dist_freq_CA_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="GG"))
colnames(dist_freq_CA_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GG")
dist_freq_CA_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="GT"))
colnames(dist_freq_CA_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_GT")
dist_freq_CA_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="TA"))
colnames(dist_freq_CA_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TA")
dist_freq_CA_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="TC"))
colnames(dist_freq_CA_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TC")
dist_freq_CA_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="TG"))
colnames(dist_freq_CA_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TG")
dist_freq_CA_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CA", di2="TT"))
colnames(dist_freq_CA_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CA_TT")

dist_freq_CA_test = cbind(dist_freq_CA_CA_test, dist_freq_CA_CC_test, dist_freq_CA_CG_test, dist_freq_CA_CT_test,
                          dist_freq_CA_GA_test, dist_freq_CA_GC_test, dist_freq_CA_GG_test, dist_freq_CA_GT_test,
                          dist_freq_CA_TA_test, dist_freq_CA_TC_test, dist_freq_CA_TG_test, dist_freq_CA_TT_test)
saveRDS(dist_freq_CA_test, "data/Created/tiling_dist_freq_CA_test.rds")



dist_freq_CC_CC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="CC"))
colnames(dist_freq_CC_CC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CC")
dist_freq_CC_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="CG"))
colnames(dist_freq_CC_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CG")
dist_freq_CC_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="CT"))
colnames(dist_freq_CC_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_CT")
dist_freq_CC_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="GA"))
colnames(dist_freq_CC_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GA")
dist_freq_CC_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="GC"))
colnames(dist_freq_CC_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GC")
dist_freq_CC_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="GG"))
colnames(dist_freq_CC_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GG")
dist_freq_CC_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="GT"))
colnames(dist_freq_CC_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_GT")
dist_freq_CC_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="TA"))
colnames(dist_freq_CC_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TA")
dist_freq_CC_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="TC"))
colnames(dist_freq_CC_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TC")
dist_freq_CC_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="TG"))
colnames(dist_freq_CC_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TG")
dist_freq_CC_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CC", di2="TT"))
colnames(dist_freq_CC_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CC_TT")

dist_freq_CC_test = cbind(dist_freq_CC_CC_test, dist_freq_CC_CG_test, dist_freq_CC_CT_test,
                          dist_freq_CC_GA_test, dist_freq_CC_GC_test, dist_freq_CC_GG_test, dist_freq_CC_GT_test,
                          dist_freq_CC_TA_test, dist_freq_CC_TC_test, dist_freq_CC_TG_test, dist_freq_CC_TT_test)
saveRDS(dist_freq_CC_test, "data/Created/tiling_dist_freq_CC_test.rds")



dist_freq_CG_CG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="CG"))
colnames(dist_freq_CG_CG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CG")
dist_freq_CG_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="CT"))
colnames(dist_freq_CG_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_CT")
dist_freq_CG_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="GA"))
colnames(dist_freq_CG_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GA")
dist_freq_CG_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="GC"))
colnames(dist_freq_CG_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GC")
dist_freq_CG_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="GG"))
colnames(dist_freq_CG_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GG")
dist_freq_CG_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="GT"))
colnames(dist_freq_CG_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_GT")
dist_freq_CG_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="TA"))
colnames(dist_freq_CG_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TA")
dist_freq_CG_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="TC"))
colnames(dist_freq_CG_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TC")
dist_freq_CG_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="TG"))
colnames(dist_freq_CG_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TG")
dist_freq_CG_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CG", di2="TT"))
colnames(dist_freq_CG_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CG_TT")

dist_freq_CG_test = cbind(dist_freq_CG_CG_test, dist_freq_CG_CT_test,
                          dist_freq_CG_GA_test, dist_freq_CG_GC_test, dist_freq_CG_GG_test, dist_freq_CG_GT_test,
                          dist_freq_CG_TA_test, dist_freq_CG_TC_test, dist_freq_CG_TG_test, dist_freq_CG_TT_test)
saveRDS(dist_freq_CG_test, "data/Created/tiling_dist_freq_CG_test.rds")



dist_freq_CT_CT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="CT"))
colnames(dist_freq_CT_CT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_CT")
dist_freq_CT_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="GA"))
colnames(dist_freq_CT_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GA")
dist_freq_CT_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="GC"))
colnames(dist_freq_CT_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GC")
dist_freq_CT_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="GG"))
colnames(dist_freq_CT_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GG")
dist_freq_CT_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="GT"))
colnames(dist_freq_CT_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_GT")
dist_freq_CT_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="TA"))
colnames(dist_freq_CT_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TA")
dist_freq_CT_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="TC"))
colnames(dist_freq_CT_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TC")
dist_freq_CT_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="TG"))
colnames(dist_freq_CT_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TG")
dist_freq_CT_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="CT", di2="TT"))
colnames(dist_freq_CT_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_CT_TT")

dist_freq_CT_test = cbind(dist_freq_CT_CT_test,
                          dist_freq_CT_GA_test, dist_freq_CT_GC_test, dist_freq_CT_GG_test, dist_freq_CT_GT_test,
                          dist_freq_CT_TA_test, dist_freq_CT_TC_test, dist_freq_CT_TG_test, dist_freq_CT_TT_test)
saveRDS(dist_freq_CT_test, "data/Created/tiling_dist_freq_CT_test.rds")



dist_freq_GA_GA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="GA"))
colnames(dist_freq_GA_GA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GA")
dist_freq_GA_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="GC"))
colnames(dist_freq_GA_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GC")
dist_freq_GA_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="GG"))
colnames(dist_freq_GA_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GG")
dist_freq_GA_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="GT"))
colnames(dist_freq_GA_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_GT")
dist_freq_GA_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="TA"))
colnames(dist_freq_GA_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TA")
dist_freq_GA_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="TC"))
colnames(dist_freq_GA_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TC")
dist_freq_GA_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="TG"))
colnames(dist_freq_GA_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TG")
dist_freq_GA_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GA", di2="TT"))
colnames(dist_freq_GA_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GA_TT")

dist_freq_GA_test = cbind(dist_freq_GA_GA_test, dist_freq_GA_GC_test, dist_freq_GA_GG_test, dist_freq_GA_GT_test,
                          dist_freq_GA_TA_test, dist_freq_GA_TC_test, dist_freq_GA_TG_test, dist_freq_GA_TT_test)
saveRDS(dist_freq_GA_test, "data/Created/tiling_dist_freq_GA_test.rds")



dist_freq_GC_GC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="GC"))
colnames(dist_freq_GC_GC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GC")
dist_freq_GC_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="GG"))
colnames(dist_freq_GC_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GG")
dist_freq_GC_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="GT"))
colnames(dist_freq_GC_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_GT")
dist_freq_GC_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="TA"))
colnames(dist_freq_GC_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TA")
dist_freq_GC_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="TC"))
colnames(dist_freq_GC_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TC")
dist_freq_GC_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="TG"))
colnames(dist_freq_GC_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TG")
dist_freq_GC_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GC", di2="TT"))
colnames(dist_freq_GC_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GC_TT")

dist_freq_GC_test = cbind(dist_freq_GC_GC_test, dist_freq_GC_GG_test, dist_freq_GC_GT_test,
                          dist_freq_GC_TA_test, dist_freq_GC_TC_test, dist_freq_GC_TG_test, dist_freq_GC_TT_test)
saveRDS(dist_freq_GC_test, "data/Created/tiling_dist_freq_GC_test.rds")



dist_freq_GG_GG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="GG"))
colnames(dist_freq_GG_GG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_GG")
dist_freq_GG_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="GT"))
colnames(dist_freq_GG_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_GT")
dist_freq_GG_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="TA"))
colnames(dist_freq_GG_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TA")
dist_freq_GG_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="TC"))
colnames(dist_freq_GG_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TC")
dist_freq_GG_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="TG"))
colnames(dist_freq_GG_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TG")
dist_freq_GG_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GG", di2="TT"))
colnames(dist_freq_GG_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GG_TT")

dist_freq_GG_test = cbind(dist_freq_GG_GG_test, dist_freq_GG_GT_test,
                          dist_freq_GG_TA_test, dist_freq_GG_TC_test, dist_freq_GG_TG_test, dist_freq_GG_TT_test)
saveRDS(dist_freq_GG_test, "data/Created/tiling_dist_freq_GG_test.rds")



dist_freq_GT_GT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2="GT"))
colnames(dist_freq_GT_GT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_GT")
dist_freq_GT_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2="TA"))
colnames(dist_freq_GT_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TA")
dist_freq_GT_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2="TC"))
colnames(dist_freq_GT_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TC")
dist_freq_GT_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2="TG"))
colnames(dist_freq_GT_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TG")
dist_freq_GT_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="GT", di2="TT"))
colnames(dist_freq_GT_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_GT_TT")

dist_freq_GT_test = cbind(dist_freq_GT_GT_test,
                          dist_freq_GT_TA_test, dist_freq_GT_TC_test, dist_freq_GT_TG_test, dist_freq_GT_TT_test)
saveRDS(dist_freq_GT_test, "data/Created/tiling_dist_freq_GT_test.rds")



dist_freq_TA_TA_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2="TA"))
colnames(dist_freq_TA_TA_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TA")
dist_freq_TA_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2="TC"))
colnames(dist_freq_TA_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TC")
dist_freq_TA_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2="TG"))
colnames(dist_freq_TA_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TG")
dist_freq_TA_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TA", di2="TT"))
colnames(dist_freq_TA_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TA_TT")

dist_freq_TA_test = cbind(dist_freq_TA_TA_test, dist_freq_TA_TC_test, dist_freq_TA_TG_test, dist_freq_TA_TT_test)
saveRDS(dist_freq_TA_test, "data/Created/tiling_dist_freq_TA_test.rds")



dist_freq_TC_TC_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TC", di2="TC"))
colnames(dist_freq_TC_TC_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TC")
dist_freq_TC_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TC", di2="TG"))
colnames(dist_freq_TC_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TG")
dist_freq_TC_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TC", di2="TT"))
colnames(dist_freq_TC_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TC_TT")

dist_freq_TC_test = cbind(dist_freq_TC_TC_test, dist_freq_TC_TG_test, dist_freq_TC_TT_test)
saveRDS(dist_freq_TC_test, "data/Created/tiling_dist_freq_TC_test.rds")




dist_freq_TG_TG_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TG", di2="TG"))
colnames(dist_freq_TG_TG_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_TG")
dist_freq_TG_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TG", di2="TT"))
colnames(dist_freq_TG_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TG_TT")

dist_freq_TG_test = cbind(dist_freq_TG_TG_test, dist_freq_TG_TT_test)
saveRDS(dist_freq_TG_test, "data/Created/tiling_dist_freq_TG_test.rds")



dist_freq_TT_TT_test = t(apply(dat_test, 1, max_freq_bidirectional, di1="TT", di2="TT"))
colnames(dist_freq_TT_TT_test) = paste0("dist", c(5, 10, 15, 20, 25), "_freq_TT_TT")

dist_freq_TT_test = cbind(dist_freq_TT_TT_test)
saveRDS(dist_freq_TT_test, "data/Created/tiling_dist_freq_TT_test.rds")


dist_freq_di_test = data.frame(map(paste0("dist_freq_", dinucleotides, "_test"), get))
saveRDS(dist_freq_di_test, "data/Created/tiling_dist_freq_di_test.rds")

