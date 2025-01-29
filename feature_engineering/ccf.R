dat = readRDS("data/Created/processed_tiling_newC0.rds")
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")

########################################################################
# Tiling
########################################################################

## Train data:

# Nucleotides to nucleotides:

Xone = dat %>% select(all_of(ps1))
Xone_A = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_C = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_G = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_T = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_A) = colnames(Xone)
colnames(Xone_C) = colnames(Xone)
colnames(Xone_G) = colnames(Xone)
colnames(Xone_T) = colnames(Xone)
Xone_A[] = (Xone == "A") %>% as.matrix() %>% as.numeric()
Xone_C[] = (Xone == "C") %>% as.matrix() %>% as.numeric()
Xone_G[] = (Xone == "G") %>% as.matrix() %>% as.numeric()
Xone_T[] = (Xone == "T") %>% as.matrix() %>% as.numeric()

Xone_A_list = as.list(data.frame(t(Xone_A)))
Xone_C_list = as.list(data.frame(t(Xone_C)))
Xone_G_list = as.list(data.frame(t(Xone_G)))
Xone_T_list = as.list(data.frame(t(Xone_T)))

Xone_A_to_C_ccf = map2(Xone_A_list, Xone_C_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_G_ccf = map2(Xone_A_list, Xone_G_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_T_ccf = map2(Xone_A_list, Xone_T_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_G_ccf = map2(Xone_C_list, Xone_G_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_T_ccf = map2(Xone_C_list, Xone_T_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_G_to_T_ccf = map2(Xone_G_list, Xone_T_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_A_to_C_ccf = do.call(rbind, Xone_A_to_C_ccf)
Xone_A_to_G_ccf = do.call(rbind, Xone_A_to_G_ccf)
Xone_A_to_T_ccf = do.call(rbind, Xone_A_to_T_ccf)
Xone_C_to_G_ccf = do.call(rbind, Xone_C_to_G_ccf)
Xone_C_to_T_ccf = do.call(rbind, Xone_C_to_T_ccf)
Xone_G_to_T_ccf = do.call(rbind, Xone_G_to_T_ccf)

Xone_A_to_C_ccf[is.na(Xone_A_to_C_ccf)] = 0
Xone_A_to_G_ccf[is.na(Xone_A_to_G_ccf)] = 0
Xone_A_to_T_ccf[is.na(Xone_A_to_T_ccf)] = 0
Xone_C_to_G_ccf[is.na(Xone_C_to_G_ccf)] = 0
Xone_C_to_T_ccf[is.na(Xone_C_to_T_ccf)] = 0
Xone_G_to_T_ccf[is.na(Xone_G_to_T_ccf)] = 0

rownames(Xone_A_to_C_ccf) = NULL
rownames(Xone_A_to_G_ccf) = NULL
rownames(Xone_A_to_T_ccf) = NULL
rownames(Xone_C_to_G_ccf) = NULL
rownames(Xone_C_to_T_ccf) = NULL
rownames(Xone_G_to_T_ccf) = NULL

colnames(Xone_A_to_C_ccf) = paste0("ccf_A_to_C_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_G_ccf) = paste0("ccf_A_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_T_ccf) = paste0("ccf_A_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_G_ccf) = paste0("ccf_C_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_T_ccf) = paste0("ccf_C_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_G_to_T_ccf) = paste0("ccf_G_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Pairs of nucleotides to pairs of nucleotides:

Xone_AorC = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_AorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_AorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_GorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AorC) = colnames(Xone)
colnames(Xone_AorG) = colnames(Xone)
colnames(Xone_AorT) = colnames(Xone)
colnames(Xone_CorG) = colnames(Xone)
colnames(Xone_CorT) = colnames(Xone)
colnames(Xone_GorT) = colnames(Xone)
Xone_AorC[] = ((Xone == "A") | (Xone == "C")) %>% as.matrix() %>% as.numeric()
Xone_AorG[] = ((Xone == "A") | (Xone == "G")) %>% as.matrix() %>% as.numeric()
Xone_AorT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()
Xone_CorT[] = ((Xone == "C") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_GorT[] = ((Xone == "G") | (Xone == "T")) %>% as.matrix() %>% as.numeric()

Xone_AorC_list = as.list(data.frame(t(Xone_AorC)))
Xone_AorG_list = as.list(data.frame(t(Xone_AorG)))
Xone_AorT_list = as.list(data.frame(t(Xone_AorT)))
Xone_CorG_list = as.list(data.frame(t(Xone_CorG)))
Xone_CorT_list = as.list(data.frame(t(Xone_CorT)))
Xone_GorT_list = as.list(data.frame(t(Xone_GorT)))

Xone_AorT_to_CorG_ccf = map2(Xone_AorT_list, Xone_CorG_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_AorT_to_CorG_ccf = do.call(rbind, Xone_AorT_to_CorG_ccf)

Xone_AorT_to_CorG_ccf[is.na(Xone_AorT_to_CorG_ccf)] = 0

rownames(Xone_AorT_to_CorG_ccf) = NULL

colnames(Xone_AorT_to_CorG_ccf) = paste0("ccf_AorT_to_CorG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Groups of dinucleotides to groups of dinucleotides:

Xtwo_AAorATorTAorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorCGorGCorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AAorATorTAorTT) = colnames(Xtwo)
colnames(Xtwo_CCorCGorGCorGG) = colnames(Xtwo)
Xtwo_AAorATorTAorTT[] = ((Xtwo == "AA") | (Xtwo == "AT") |
                                (Xtwo == "TA") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG[] = ((Xtwo == "CC") | (Xtwo == "CG") |
                                (Xtwo == "GC") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_list = as.list(data.frame(t(Xtwo_AAorATorTAorTT)))
Xtwo_CCorCGorGCorGG_list = as.list(data.frame(t(Xtwo_CCorCGorGCorGG)))

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf = map2(Xtwo_AAorATorTAorTT_list, Xtwo_CCorCGorGCorGG_list, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf = do.call(rbind, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf)

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf[is.na(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf)] = 0

rownames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf) = NULL

colnames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf) = paste0("ccf_AAorATorTAorTT_to_CCorCGorGCorGG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))



## Test data:

# Nucleotides to nucleotides:

Xone_test = dat_test %>% select(all_of(ps1))
Xone_A_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_C_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_G_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_T_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_A_test) = colnames(Xone_test)
colnames(Xone_C_test) = colnames(Xone_test)
colnames(Xone_G_test) = colnames(Xone_test)
colnames(Xone_T_test) = colnames(Xone_test)
Xone_A_test[] = (Xone_test == "A") %>% as.matrix() %>% as.numeric()
Xone_C_test[] = (Xone_test == "C") %>% as.matrix() %>% as.numeric()
Xone_G_test[] = (Xone_test == "G") %>% as.matrix() %>% as.numeric()
Xone_T_test[] = (Xone_test == "T") %>% as.matrix() %>% as.numeric()

Xone_A_list_test = as.list(data.frame(t(Xone_A_test)))
Xone_C_list_test = as.list(data.frame(t(Xone_C_test)))
Xone_G_list_test = as.list(data.frame(t(Xone_G_test)))
Xone_T_list_test = as.list(data.frame(t(Xone_T_test)))

Xone_A_to_C_ccf_test = map2(Xone_A_list_test, Xone_C_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_G_ccf_test = map2(Xone_A_list_test, Xone_G_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_T_ccf_test = map2(Xone_A_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_G_ccf_test = map2(Xone_C_list_test, Xone_G_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_T_ccf_test = map2(Xone_C_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_G_to_T_ccf_test = map2(Xone_G_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_A_to_C_ccf_test = do.call(rbind, Xone_A_to_C_ccf_test)
Xone_A_to_G_ccf_test = do.call(rbind, Xone_A_to_G_ccf_test)
Xone_A_to_T_ccf_test = do.call(rbind, Xone_A_to_T_ccf_test)
Xone_C_to_G_ccf_test = do.call(rbind, Xone_C_to_G_ccf_test)
Xone_C_to_T_ccf_test = do.call(rbind, Xone_C_to_T_ccf_test)
Xone_G_to_T_ccf_test = do.call(rbind, Xone_G_to_T_ccf_test)

Xone_A_to_C_ccf_test[is.na(Xone_A_to_C_ccf_test)] = 0
Xone_A_to_G_ccf_test[is.na(Xone_A_to_G_ccf_test)] = 0
Xone_A_to_T_ccf_test[is.na(Xone_A_to_T_ccf_test)] = 0
Xone_C_to_G_ccf_test[is.na(Xone_C_to_G_ccf_test)] = 0
Xone_C_to_T_ccf_test[is.na(Xone_C_to_T_ccf_test)] = 0
Xone_G_to_T_ccf_test[is.na(Xone_G_to_T_ccf_test)] = 0

rownames(Xone_A_to_C_ccf_test) = NULL
rownames(Xone_A_to_G_ccf_test) = NULL
rownames(Xone_A_to_T_ccf_test) = NULL
rownames(Xone_C_to_G_ccf_test) = NULL
rownames(Xone_C_to_T_ccf_test) = NULL
rownames(Xone_G_to_T_ccf_test) = NULL

colnames(Xone_A_to_C_ccf_test) = paste0("ccf_A_to_C_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_G_ccf_test) = paste0("ccf_A_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_T_ccf_test) = paste0("ccf_A_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_G_ccf_test) = paste0("ccf_C_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_T_ccf_test) = paste0("ccf_C_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_G_to_T_ccf_test) = paste0("ccf_G_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Pairs of nucleotides to pairs of nucleotides:

Xone_AorC_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_AorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_AorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_GorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AorC_test) = colnames(Xone_test)
colnames(Xone_AorG_test) = colnames(Xone_test)
colnames(Xone_AorT_test) = colnames(Xone_test)
colnames(Xone_CorG_test) = colnames(Xone_test)
colnames(Xone_CorT_test) = colnames(Xone_test)
colnames(Xone_GorT_test) = colnames(Xone_test)
Xone_AorC_test[] = ((Xone_test == "A") | (Xone_test == "C")) %>% as.matrix() %>% as.numeric()
Xone_AorG_test[] = ((Xone_test == "A") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()
Xone_AorT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()
Xone_CorT_test[] = ((Xone_test == "C") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_GorT_test[] = ((Xone_test == "G") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()

Xone_AorC_list_test = as.list(data.frame(t(Xone_AorC_test)))
Xone_AorG_list_test = as.list(data.frame(t(Xone_AorG_test)))
Xone_AorT_list_test = as.list(data.frame(t(Xone_AorT_test)))
Xone_CorG_list_test = as.list(data.frame(t(Xone_CorG_test)))
Xone_CorT_list_test = as.list(data.frame(t(Xone_CorT_test)))
Xone_GorT_list_test = as.list(data.frame(t(Xone_GorT_test)))

Xone_AorT_to_CorG_ccf_test = map2(Xone_AorT_list_test, Xone_CorG_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_AorT_to_CorG_ccf_test = do.call(rbind, Xone_AorT_to_CorG_ccf_test)

Xone_AorT_to_CorG_ccf_test[is.na(Xone_AorT_to_CorG_ccf_test)] = 0

rownames(Xone_AorT_to_CorG_ccf_test) = NULL

colnames(Xone_AorT_to_CorG_ccf_test) = paste0("ccf_AorT_to_CorG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Dinucleotides to dinucleotides:

#FIXME:

Xtwo_test = dat_test %>% select(all_of(ps2))
Xtwo_AA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_AC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_AG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_GA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_GC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_GG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_GT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_TA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_TC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_TG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_TT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AA_test) = colnames(Xtwo_test)
colnames(Xtwo_AC_test) = colnames(Xtwo_test)
colnames(Xtwo_AG_test) = colnames(Xtwo_test)
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
colnames(Xtwo_CA_test) = colnames(Xtwo_test)
colnames(Xtwo_CC_test) = colnames(Xtwo_test)
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
colnames(Xtwo_CT_test) = colnames(Xtwo_test)
colnames(Xtwo_GA_test) = colnames(Xtwo_test)
colnames(Xtwo_GC_test) = colnames(Xtwo_test)
colnames(Xtwo_GG_test) = colnames(Xtwo_test)
colnames(Xtwo_GT_test) = colnames(Xtwo_test)
colnames(Xtwo_TA_test) = colnames(Xtwo_test)
colnames(Xtwo_TC_test) = colnames(Xtwo_test)
colnames(Xtwo_TG_test) = colnames(Xtwo_test)
colnames(Xtwo_TT_test) = colnames(Xtwo_test)
Xtwo_AA_test[] = (Xtwo_test == "AA") %>% as.matrix() %>% as.numeric()
Xtwo_AC_test[] = (Xtwo_test == "AC") %>% as.matrix() %>% as.numeric()
Xtwo_AG_test[] = (Xtwo_test == "AG") %>% as.matrix() %>% as.numeric()
Xtwo_AT_test[] = (Xtwo_test == "AT") %>% as.matrix() %>% as.numeric()
Xtwo_CA_test[] = (Xtwo_test == "CA") %>% as.matrix() %>% as.numeric()
Xtwo_CC_test[] = (Xtwo_test == "CC") %>% as.matrix() %>% as.numeric()
Xtwo_CG_test[] = (Xtwo_test == "CG") %>% as.matrix() %>% as.numeric()
Xtwo_CT_test[] = (Xtwo_test == "CT") %>% as.matrix() %>% as.numeric()
Xtwo_GA_test[] = (Xtwo_test == "GA") %>% as.matrix() %>% as.numeric()
Xtwo_GC_test[] = (Xtwo_test == "GC") %>% as.matrix() %>% as.numeric()
Xtwo_GG_test[] = (Xtwo_test == "GG") %>% as.matrix() %>% as.numeric()
Xtwo_GT_test[] = (Xtwo_test == "GT") %>% as.matrix() %>% as.numeric()
Xtwo_TA_test[] = (Xtwo_test == "TA") %>% as.matrix() %>% as.numeric()
Xtwo_TC_test[] = (Xtwo_test == "TC") %>% as.matrix() %>% as.numeric()
Xtwo_TG_test[] = (Xtwo_test == "TG") %>% as.matrix() %>% as.numeric()
Xtwo_TT_test[] = (Xtwo_test == "TT") %>% as.matrix() %>% as.numeric()

Xtwo_AA_list_test = as.list(data.frame(t(Xtwo_AA_test)))
Xtwo_AC_list_test = as.list(data.frame(t(Xtwo_AC_test)))
Xtwo_AG_list_test = as.list(data.frame(t(Xtwo_AG_test)))
Xtwo_AT_list_test = as.list(data.frame(t(Xtwo_AT_test)))
Xtwo_CA_list_test = as.list(data.frame(t(Xtwo_CA_test)))
Xtwo_CC_list_test = as.list(data.frame(t(Xtwo_CC_test)))
Xtwo_CG_list_test = as.list(data.frame(t(Xtwo_CG_test)))
Xtwo_CT_list_test = as.list(data.frame(t(Xtwo_CT_test)))
Xtwo_GA_list_test = as.list(data.frame(t(Xtwo_GA_test)))
Xtwo_GC_list_test = as.list(data.frame(t(Xtwo_GC_test)))
Xtwo_GG_list_test = as.list(data.frame(t(Xtwo_GG_test)))
Xtwo_GT_list_test = as.list(data.frame(t(Xtwo_GT_test)))
Xtwo_TA_list_test = as.list(data.frame(t(Xtwo_TA_test)))
Xtwo_TC_list_test = as.list(data.frame(t(Xtwo_TC_test)))
Xtwo_TG_list_test = as.list(data.frame(t(Xtwo_TG_test)))
Xtwo_TT_list_test = as.list(data.frame(t(Xtwo_TT_test)))

Xone_AA_to_AT_ccf_test = map2(Xone_A_list_test, Xone_C_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_AA_to_TA_ccf_test = map2(Xone_A_list_test, Xone_G_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_AA_to_TT_ccf_test
Xone_AA_to_CC_ccf_test
Xone_AA_to_CG_ccf_test
Xone_AA_to_GG_ccf_test

Xone_A_to_T_ccf_test = map2(Xone_A_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_G_ccf_test = map2(Xone_C_list_test, Xone_G_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_T_ccf_test = map2(Xone_C_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_G_to_T_ccf_test = map2(Xone_G_list_test, Xone_T_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_A_to_C_ccf_test = do.call(rbind, Xone_A_to_C_ccf_test)
Xone_A_to_G_ccf_test = do.call(rbind, Xone_A_to_G_ccf_test)
Xone_A_to_T_ccf_test = do.call(rbind, Xone_A_to_T_ccf_test)
Xone_C_to_G_ccf_test = do.call(rbind, Xone_C_to_G_ccf_test)
Xone_C_to_T_ccf_test = do.call(rbind, Xone_C_to_T_ccf_test)
Xone_G_to_T_ccf_test = do.call(rbind, Xone_G_to_T_ccf_test)

Xone_A_to_C_ccf_test[is.na(Xone_A_to_C_ccf_test)] = 0
Xone_A_to_G_ccf_test[is.na(Xone_A_to_G_ccf_test)] = 0
Xone_A_to_T_ccf_test[is.na(Xone_A_to_T_ccf_test)] = 0
Xone_C_to_G_ccf_test[is.na(Xone_C_to_G_ccf_test)] = 0
Xone_C_to_T_ccf_test[is.na(Xone_C_to_T_ccf_test)] = 0
Xone_G_to_T_ccf_test[is.na(Xone_G_to_T_ccf_test)] = 0

rownames(Xone_A_to_C_ccf_test) = NULL
rownames(Xone_A_to_G_ccf_test) = NULL
rownames(Xone_A_to_T_ccf_test) = NULL
rownames(Xone_C_to_G_ccf_test) = NULL
rownames(Xone_C_to_T_ccf_test) = NULL
rownames(Xone_G_to_T_ccf_test) = NULL

colnames(Xone_A_to_C_ccf_test) = paste0("ccf_A_to_C_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_G_ccf_test) = paste0("ccf_A_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_T_ccf_test) = paste0("ccf_A_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_G_ccf_test) = paste0("ccf_C_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_T_ccf_test) = paste0("ccf_C_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_G_to_T_ccf_test) = paste0("ccf_G_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Groups of dinucleotides to groups of dinucleotides:

Xtwo_AAorATorTAorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorCGorGCorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AAorATorTAorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorCGorGCorGG_test) = colnames(Xtwo_test)
Xtwo_AAorATorTAorTT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "AT") |
                                (Xtwo_test == "TA") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "CG") |
                                (Xtwo_test == "GC") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_list_test = as.list(data.frame(t(Xtwo_AAorATorTAorTT_test)))
Xtwo_CCorCGorGCorGG_list_test = as.list(data.frame(t(Xtwo_CCorCGorGCorGG_test)))

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test = map2(Xtwo_AAorATorTAorTT_list_test, Xtwo_CCorCGorGCorGG_list_test, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test = do.call(rbind, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test)

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test[is.na(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test)] = 0

rownames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test) = NULL

colnames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test) = paste0("ccf_AAorATorTAorTT_to_CCorCGorGCorGG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))



########################################################################
# Random
########################################################################

# Nucleotides to nucleotides:

Xone_random_all = dat_random_all %>% select(all_of(ps1))
Xone_A_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_C_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_G_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_T_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
colnames(Xone_A_random_all) = colnames(Xone_random_all)
colnames(Xone_C_random_all) = colnames(Xone_random_all)
colnames(Xone_G_random_all) = colnames(Xone_random_all)
colnames(Xone_T_random_all) = colnames(Xone_random_all)
Xone_A_random_all[] = (Xone_random_all == "A") %>% as.matrix() %>% as.numeric()
Xone_C_random_all[] = (Xone_random_all == "C") %>% as.matrix() %>% as.numeric()
Xone_G_random_all[] = (Xone_random_all == "G") %>% as.matrix() %>% as.numeric()
Xone_T_random_all[] = (Xone_random_all == "T") %>% as.matrix() %>% as.numeric()

Xone_A_list_random_all = as.list(data.frame(t(Xone_A_random_all)))
Xone_C_list_random_all = as.list(data.frame(t(Xone_C_random_all)))
Xone_G_list_random_all = as.list(data.frame(t(Xone_G_random_all)))
Xone_T_list_random_all = as.list(data.frame(t(Xone_T_random_all)))

Xone_A_to_C_ccf_random_all = map2(Xone_A_list_random_all, Xone_C_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_G_ccf_random_all = map2(Xone_A_list_random_all, Xone_G_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_A_to_T_ccf_random_all = map2(Xone_A_list_random_all, Xone_T_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_G_ccf_random_all = map2(Xone_C_list_random_all, Xone_G_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_C_to_T_ccf_random_all = map2(Xone_C_list_random_all, Xone_T_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])
Xone_G_to_T_ccf_random_all = map2(Xone_G_list_random_all, Xone_T_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_A_to_C_ccf_random_all = do.call(rbind, Xone_A_to_C_ccf_random_all)
Xone_A_to_G_ccf_random_all = do.call(rbind, Xone_A_to_G_ccf_random_all)
Xone_A_to_T_ccf_random_all = do.call(rbind, Xone_A_to_T_ccf_random_all)
Xone_C_to_G_ccf_random_all = do.call(rbind, Xone_C_to_G_ccf_random_all)
Xone_C_to_T_ccf_random_all = do.call(rbind, Xone_C_to_T_ccf_random_all)
Xone_G_to_T_ccf_random_all = do.call(rbind, Xone_G_to_T_ccf_random_all)

Xone_A_to_C_ccf_random_all[is.na(Xone_A_to_C_ccf_random_all)] = 0
Xone_A_to_G_ccf_random_all[is.na(Xone_A_to_G_ccf_random_all)] = 0
Xone_A_to_T_ccf_random_all[is.na(Xone_A_to_T_ccf_random_all)] = 0
Xone_C_to_G_ccf_random_all[is.na(Xone_C_to_G_ccf_random_all)] = 0
Xone_C_to_T_ccf_random_all[is.na(Xone_C_to_T_ccf_random_all)] = 0
Xone_G_to_T_ccf_random_all[is.na(Xone_G_to_T_ccf_random_all)] = 0

rownames(Xone_A_to_C_ccf_random_all) = NULL
rownames(Xone_A_to_G_ccf_random_all) = NULL
rownames(Xone_A_to_T_ccf_random_all) = NULL
rownames(Xone_C_to_G_ccf_random_all) = NULL
rownames(Xone_C_to_T_ccf_random_all) = NULL
rownames(Xone_G_to_T_ccf_random_all) = NULL

colnames(Xone_A_to_C_ccf_random_all) = paste0("ccf_A_to_C_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_G_ccf_random_all) = paste0("ccf_A_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_A_to_T_ccf_random_all) = paste0("ccf_A_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_G_ccf_random_all) = paste0("ccf_C_to_G_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_C_to_T_ccf_random_all) = paste0("ccf_C_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))
colnames(Xone_G_to_T_ccf_random_all) = paste0("ccf_G_to_T_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Pairs of nucleotides to pairs of nucleotides:

Xone_AorC_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_AorG_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_AorT_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_CorG_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_CorT_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
Xone_GorT_random_all = matrix(nrow=nrow(Xone_random_all), ncol=ncol(Xone_random_all))
colnames(Xone_AorC_random_all) = colnames(Xone_random_all)
colnames(Xone_AorG_random_all) = colnames(Xone_random_all)
colnames(Xone_AorT_random_all) = colnames(Xone_random_all)
colnames(Xone_CorG_random_all) = colnames(Xone_random_all)
colnames(Xone_CorT_random_all) = colnames(Xone_random_all)
colnames(Xone_GorT_random_all) = colnames(Xone_random_all)
Xone_AorC_random_all[] = ((Xone_random_all == "A") | (Xone_random_all == "C")) %>% as.matrix() %>% as.numeric()
Xone_AorG_random_all[] = ((Xone_random_all == "A") | (Xone_random_all == "G")) %>% as.matrix() %>% as.numeric()
Xone_AorT_random_all[] = ((Xone_random_all == "A") | (Xone_random_all == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_random_all[] = ((Xone_random_all == "C") | (Xone_random_all == "G")) %>% as.matrix() %>% as.numeric()
Xone_CorT_random_all[] = ((Xone_random_all == "C") | (Xone_random_all == "T")) %>% as.matrix() %>% as.numeric()
Xone_GorT_random_all[] = ((Xone_random_all == "G") | (Xone_random_all == "T")) %>% as.matrix() %>% as.numeric()

Xone_AorC_list_random_all = as.list(data.frame(t(Xone_AorC_random_all)))
Xone_AorG_list_random_all = as.list(data.frame(t(Xone_AorG_random_all)))
Xone_AorT_list_random_all = as.list(data.frame(t(Xone_AorT_random_all)))
Xone_CorG_list_random_all = as.list(data.frame(t(Xone_CorG_random_all)))
Xone_CorT_list_random_all = as.list(data.frame(t(Xone_CorT_random_all)))
Xone_GorT_list_random_all = as.list(data.frame(t(Xone_GorT_random_all)))

Xone_AorT_to_CorG_ccf_random_all = map2(Xone_AorT_list_random_all, Xone_CorG_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xone_AorT_to_CorG_ccf_random_all = do.call(rbind, Xone_AorT_to_CorG_ccf_random_all)

Xone_AorT_to_CorG_ccf_random_all[is.na(Xone_AorT_to_CorG_ccf_random_all)] = 0

rownames(Xone_AorT_to_CorG_ccf_random_all) = NULL

colnames(Xone_AorT_to_CorG_ccf_random_all) = paste0("ccf_AorT_to_CorG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))


# Groups of dinucleotides to groups of dinucleotides:

Xtwo_AAorATorTAorTT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CCorCGorGCorGG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
colnames(Xtwo_AAorATorTAorTT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CCorCGorGCorGG_random_all) = colnames(Xtwo_random_all)
Xtwo_AAorATorTAorTT_random_all[] = ((Xtwo_random_all == "AA") | (Xtwo_random_all == "AT") |
                                (Xtwo_random_all == "TA") | (Xtwo_random_all == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG_random_all[] = ((Xtwo_random_all == "CC") | (Xtwo_random_all == "CG") |
                                (Xtwo_random_all == "GC") | (Xtwo_random_all == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_list_random_all = as.list(data.frame(t(Xtwo_AAorATorTAorTT_random_all)))
Xtwo_CCorCGorGCorGG_list_random_all = as.list(data.frame(t(Xtwo_CCorCGorGCorGG_random_all)))

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all = map2(Xtwo_AAorATorTAorTT_list_random_all, Xtwo_CCorCGorGCorGG_list_random_all, ~ccf(.x,.y, plot=FALSE, lag.max=10)[["acf"]][,1,1])

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all = do.call(rbind, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all)

Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all[is.na(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all)] = 0

rownames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all) = NULL

colnames(Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all) = paste0("ccf_AAorATorTAorTT_to_CCorCGorGCorGG_lag_", c(paste0("neg", 10:1), 0, paste0("pos", 1:10)))








########################################################################
# Modelling
########################################################################

# Nucleotides to nucleotides:
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0),
                 # select(all_of(trinucleotides), C0),
                 # select(C0),
                 Xone_A_to_C_ccf, Xone_A_to_G_ccf, Xone_A_to_T_ccf, 
                 Xone_C_to_G_ccf, Xone_C_to_T_ccf, Xone_G_to_T_ccf)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0),
                      # select(all_of(trinucleotides), C0),
                      # select(C0),
                      Xone_A_to_C_ccf_test, Xone_A_to_G_ccf_test, Xone_A_to_T_ccf_test, 
                      Xone_C_to_G_ccf_test, Xone_C_to_T_ccf_test, Xone_G_to_T_ccf_test)
dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0),
                            # select(all_of(trinucleotides), C0),
                            # select(C0),
                            Xone_A_to_C_ccf_random_all, Xone_A_to_G_ccf_random_all, Xone_A_to_T_ccf_random_all, 
                            Xone_C_to_G_ccf_random_all, Xone_C_to_T_ccf_random_all, Xone_G_to_T_ccf_random_all)

temp_lm = lm(C0~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.4987632
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.493006
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4908533
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


# Nucleotides to nucleotides + Pairs of nucleotides to pairs of nucleotides:
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0),
                 # select(all_of(trinucleotides), C0),
                 # select(C0),
                 Xone_A_to_C_ccf, Xone_A_to_G_ccf, Xone_A_to_T_ccf, 
                 Xone_C_to_G_ccf, Xone_C_to_T_ccf, Xone_G_to_T_ccf,
                 Xone_AorT_to_CorG_ccf)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0),
                      # select(all_of(trinucleotides), C0),
                      # select(C0),
                      Xone_A_to_C_ccf_test, Xone_A_to_G_ccf_test, Xone_A_to_T_ccf_test, 
                      Xone_C_to_G_ccf_test, Xone_C_to_T_ccf_test, Xone_G_to_T_ccf_test,
                      Xone_AorT_to_CorG_ccf_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0),
                            # select(all_of(trinucleotides), C0),
                            # select(C0),
                            Xone_A_to_C_ccf_random_all, Xone_A_to_G_ccf_random_all, Xone_A_to_T_ccf_random_all, 
                            Xone_C_to_G_ccf_random_all, Xone_C_to_T_ccf_random_all, Xone_G_to_T_ccf_random_all,
                            Xone_AorT_to_CorG_ccf_random_all)

temp_lm = lm(C0~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5044582
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5014832
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4932384
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


temp_gam_formula = as.formula(paste0("C0~", paste0(colnames(dat_temp%>%select(all_of(ps2))), 
                                                       collapse="+"),
                                     "+s(", paste0(colnames(dat_temp%>%select(-C0, 
                                                                              -all_of(ps2))), 
                                                   collapse=",k=6)+s("), ",k=6)"))
temp_gam = gam(temp_gam_formula, data=dat_temp)

cor(temp_gam$fitted.values, y)
# 
plot(temp_gam$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_gam_pred = predict(temp_gam, dat_temp_test)
cor(temp_gam_pred, y_test)
# 
plot(temp_gam_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_gam_random_all_pred = predict(temp_gam, dat_temp_random_all)
cor(temp_gam_random_all_pred, y_random_all)
# 
plot(temp_gam_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


# Fourier{Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides} + all polyA/T + N-Gapped Dinucleotides + 
# ccf{Nucleotides to nucleotides}

dat_temp = cbind(dat %>% 
                   # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                   #        all_of(dinucleotides), all_of(trinucleotides), C0),
                   select(all_of(ps2), all_of(trinucleotides), C0),
                   # select(all_of(trinucleotides), C0),
                   # select(C0),
                 Xone_A_fourier_mag^2, Xone_C_fourier_mag^2, Xone_G_fourier_mag^2, Xone_T_fourier_mag^2,
                 Xone_A_fourier_phase^2, Xone_C_fourier_phase^2, Xone_G_fourier_phase^2, Xone_T_fourier_phase^2,
                 Xone_AorT_fourier_mag^2, 
                 Xone_CorG_fourier_mag^2,
                 Xone_AorT_fourier_phase^2, 
                 Xone_CorG_fourier_phase^2, 
                 Xtwo_AA_fourier_mag^2, Xtwo_AT_fourier_mag^2,
                 Xtwo_CC_fourier_mag^2, Xtwo_CG_fourier_mag^2,
                 Xtwo_GC_fourier_mag^2, Xtwo_GG_fourier_mag^2,
                 Xtwo_TA_fourier_mag^2, Xtwo_TT_fourier_mag^2,
                 Xtwo_AA_fourier_phase^2, Xtwo_AT_fourier_phase^2,
                 Xtwo_CC_fourier_phase^2, Xtwo_CG_fourier_phase^2,
                 Xtwo_GC_fourier_phase^2, Xtwo_GG_fourier_phase^2,
                 Xtwo_TA_fourier_phase^2, Xtwo_TT_fourier_phase^2,
                 Xtwo_AAorAT_fourier_mag^2, Xtwo_AAorTA_fourier_mag^2,
                 Xtwo_AAorTT_fourier_mag^2, Xtwo_ATorTA_fourier_mag^2,
                 Xtwo_ATorTT_fourier_mag^2, Xtwo_TAorTT_fourier_mag^2,
                 Xtwo_CCorCG_fourier_mag^2, Xtwo_CCorGC_fourier_mag^2,
                 Xtwo_CCorGG_fourier_mag^2, Xtwo_CGorGC_fourier_mag^2,
                 Xtwo_CGorGG_fourier_mag^2, Xtwo_GCorGG_fourier_mag^2,
                 Xtwo_AAorAT_fourier_phase^2, Xtwo_AAorTA_fourier_phase^2,
                 Xtwo_AAorTT_fourier_phase^2, Xtwo_ATorTA_fourier_phase^2,
                 Xtwo_ATorTT_fourier_phase^2, Xtwo_TAorTT_fourier_phase^2,
                 Xtwo_CCorCG_fourier_phase^2, Xtwo_CCorGC_fourier_phase^2,
                 Xtwo_CCorGG_fourier_phase^2, Xtwo_CGorGC_fourier_phase^2,
                 Xtwo_CGorGG_fourier_phase^2, Xtwo_GCorGG_fourier_phase^2,
                 dat_poly, n_gap_di,
                 Xone_A_to_C_ccf, Xone_A_to_G_ccf, Xone_A_to_T_ccf, 
                 Xone_C_to_G_ccf, Xone_C_to_T_ccf, Xone_G_to_T_ccf)

dat_temp_test = cbind(dat_test %>% 
                        # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                        #        all_of(dinucleotides), all_of(trinucleotides), C0),
                        select(all_of(ps2), all_of(trinucleotides), C0),
                        # select(all_of(trinucleotides), C0),
                        # select(C0),
                      Xone_A_fourier_mag_test^2, Xone_C_fourier_mag_test^2, Xone_G_fourier_mag_test^2, Xone_T_fourier_mag_test^2,
                      Xone_A_fourier_phase_test^2, Xone_C_fourier_phase_test^2, Xone_G_fourier_phase_test^2, Xone_T_fourier_phase_test^2,
                      Xone_AorT_fourier_mag_test^2, 
                      Xone_CorG_fourier_mag_test^2,
                      Xone_AorT_fourier_phase_test^2, 
                      Xone_CorG_fourier_phase_test^2, 
                      Xtwo_AA_fourier_mag_test^2, Xtwo_AT_fourier_mag_test^2,
                      Xtwo_CC_fourier_mag_test^2, Xtwo_CG_fourier_mag_test^2,
                      Xtwo_GC_fourier_mag_test^2, Xtwo_GG_fourier_mag_test^2,
                      Xtwo_TA_fourier_mag_test^2, Xtwo_TT_fourier_mag_test^2,
                      Xtwo_AA_fourier_phase_test^2, Xtwo_AT_fourier_phase_test^2,
                      Xtwo_CC_fourier_phase_test^2, Xtwo_CG_fourier_phase_test^2,
                      Xtwo_GC_fourier_phase_test^2, Xtwo_GG_fourier_phase_test^2,
                      Xtwo_TA_fourier_phase_test^2, Xtwo_TT_fourier_phase_test^2,
                      Xtwo_AAorAT_fourier_mag_test^2, Xtwo_AAorTA_fourier_mag_test^2,
                      Xtwo_AAorTT_fourier_mag_test^2, Xtwo_ATorTA_fourier_mag_test^2,
                      Xtwo_ATorTT_fourier_mag_test^2, Xtwo_TAorTT_fourier_mag_test^2,
                      Xtwo_CCorCG_fourier_mag_test^2, Xtwo_CCorGC_fourier_mag_test^2,
                      Xtwo_CCorGG_fourier_mag_test^2, Xtwo_CGorGC_fourier_mag_test^2,
                      Xtwo_CGorGG_fourier_mag_test^2, Xtwo_GCorGG_fourier_mag_test^2,
                      Xtwo_AAorAT_fourier_phase_test^2, Xtwo_AAorTA_fourier_phase_test^2,
                      Xtwo_AAorTT_fourier_phase_test^2, Xtwo_ATorTA_fourier_phase_test^2,
                      Xtwo_ATorTT_fourier_phase_test^2, Xtwo_TAorTT_fourier_phase_test^2,
                      Xtwo_CCorCG_fourier_phase_test^2, Xtwo_CCorGC_fourier_phase_test^2,
                      Xtwo_CCorGG_fourier_phase_test^2, Xtwo_CGorGC_fourier_phase_test^2,
                      Xtwo_CGorGG_fourier_phase_test^2, Xtwo_GCorGG_fourier_phase_test^2,
                      dat_poly_test, n_gap_di_test,
                      Xone_A_to_C_ccf_test, Xone_A_to_G_ccf_test, Xone_A_to_T_ccf_test, 
                      Xone_C_to_G_ccf_test, Xone_C_to_T_ccf_test, Xone_G_to_T_ccf_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                              #        all_of(dinucleotides), all_of(trinucleotides), C0),
                              select(all_of(ps2), all_of(trinucleotides), C0),
                              # select(all_of(trinucleotides), C0),
                              # select(C0),
                            Xone_A_fourier_mag_random_all^2, Xone_C_fourier_mag_random_all^2, Xone_G_fourier_mag_random_all^2, Xone_T_fourier_mag_random_all^2,
                            Xone_A_fourier_phase_random_all^2, Xone_C_fourier_phase_random_all^2, Xone_G_fourier_phase_random_all^2, Xone_T_fourier_phase_random_all^2,
                            Xone_AorT_fourier_mag_random_all^2, 
                            Xone_CorG_fourier_mag_random_all^2,
                            Xone_AorT_fourier_phase_random_all^2, 
                            Xone_CorG_fourier_phase_random_all^2, 
                            Xtwo_AA_fourier_mag_random_all^2, Xtwo_AT_fourier_mag_random_all^2,
                            Xtwo_CC_fourier_mag_random_all^2, Xtwo_CG_fourier_mag_random_all^2,
                            Xtwo_GC_fourier_mag_random_all^2, Xtwo_GG_fourier_mag_random_all^2,
                            Xtwo_TA_fourier_mag_random_all^2, Xtwo_TT_fourier_mag_random_all^2,
                            Xtwo_AA_fourier_phase_random_all^2, Xtwo_AT_fourier_phase_random_all^2,
                            Xtwo_CC_fourier_phase_random_all^2, Xtwo_CG_fourier_phase_random_all^2,
                            Xtwo_GC_fourier_phase_random_all^2, Xtwo_GG_fourier_phase_random_all^2,
                            Xtwo_TA_fourier_phase_random_all^2, Xtwo_TT_fourier_phase_random_all^2,
                            Xtwo_AAorAT_fourier_mag_random_all^2, Xtwo_AAorTA_fourier_mag_random_all^2,
                            Xtwo_AAorTT_fourier_mag_random_all^2, Xtwo_ATorTA_fourier_mag_random_all^2,
                            Xtwo_ATorTT_fourier_mag_random_all^2, Xtwo_TAorTT_fourier_mag_random_all^2,
                            Xtwo_CCorCG_fourier_mag_random_all^2, Xtwo_CCorGC_fourier_mag_random_all^2,
                            Xtwo_CCorGG_fourier_mag_random_all^2, Xtwo_CGorGC_fourier_mag_random_all^2,
                            Xtwo_CGorGG_fourier_mag_random_all^2, Xtwo_GCorGG_fourier_mag_random_all^2,
                            Xtwo_AAorAT_fourier_phase_random_all^2, Xtwo_AAorTA_fourier_phase_random_all^2,
                            Xtwo_AAorTT_fourier_phase_random_all^2, Xtwo_ATorTA_fourier_phase_random_all^2,
                            Xtwo_ATorTT_fourier_phase_random_all^2, Xtwo_TAorTT_fourier_phase_random_all^2,
                            Xtwo_CCorCG_fourier_phase_random_all^2, Xtwo_CCorGC_fourier_phase_random_all^2,
                            Xtwo_CCorGG_fourier_phase_random_all^2, Xtwo_CGorGC_fourier_phase_random_all^2,
                            Xtwo_CGorGG_fourier_phase_random_all^2, Xtwo_GCorGG_fourier_phase_random_all^2,
                            dat_poly_random_all, n_gap_di_random_all,
                            Xone_A_to_C_ccf_random_all, Xone_A_to_G_ccf_random_all, Xone_A_to_T_ccf_random_all, 
                            Xone_C_to_G_ccf_random_all, Xone_C_to_T_ccf_random_all, Xone_G_to_T_ccf_random_all)

temp_lm = lm(C0~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.6557339
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.6570853
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.6018134
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LM")
abline(0,1)



# Fourier{Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides} + all polyA/T + N-Gapped Dinucleotides + 
# ccf{Nucleotides to nucleotides + Pairs of nucleotides to pairs of nucleotides}

dat_temp = cbind(dat %>% 
                   # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                   #        all_of(dinucleotides), all_of(trinucleotides), C0),
                   select(all_of(ps2), all_of(trinucleotides), C0),
                 # select(all_of(trinucleotides), C0),
                 # select(C0),
                 Xone_A_fourier_mag^2, Xone_C_fourier_mag^2, Xone_G_fourier_mag^2, Xone_T_fourier_mag^2,
                 Xone_A_fourier_phase^2, Xone_C_fourier_phase^2, Xone_G_fourier_phase^2, Xone_T_fourier_phase^2,
                 Xone_AorT_fourier_mag^2, 
                 Xone_CorG_fourier_mag^2,
                 Xone_AorT_fourier_phase^2, 
                 Xone_CorG_fourier_phase^2, 
                 Xtwo_AA_fourier_mag^2, Xtwo_AT_fourier_mag^2,
                 Xtwo_CC_fourier_mag^2, Xtwo_CG_fourier_mag^2,
                 Xtwo_GC_fourier_mag^2, Xtwo_GG_fourier_mag^2,
                 Xtwo_TA_fourier_mag^2, Xtwo_TT_fourier_mag^2,
                 Xtwo_AA_fourier_phase^2, Xtwo_AT_fourier_phase^2,
                 Xtwo_CC_fourier_phase^2, Xtwo_CG_fourier_phase^2,
                 Xtwo_GC_fourier_phase^2, Xtwo_GG_fourier_phase^2,
                 Xtwo_TA_fourier_phase^2, Xtwo_TT_fourier_phase^2,
                 Xtwo_AAorAT_fourier_mag^2, Xtwo_AAorTA_fourier_mag^2,
                 Xtwo_AAorTT_fourier_mag^2, Xtwo_ATorTA_fourier_mag^2,
                 Xtwo_ATorTT_fourier_mag^2, Xtwo_TAorTT_fourier_mag^2,
                 Xtwo_CCorCG_fourier_mag^2, Xtwo_CCorGC_fourier_mag^2,
                 Xtwo_CCorGG_fourier_mag^2, Xtwo_CGorGC_fourier_mag^2,
                 Xtwo_CGorGG_fourier_mag^2, Xtwo_GCorGG_fourier_mag^2,
                 Xtwo_AAorAT_fourier_phase^2, Xtwo_AAorTA_fourier_phase^2,
                 Xtwo_AAorTT_fourier_phase^2, Xtwo_ATorTA_fourier_phase^2,
                 Xtwo_ATorTT_fourier_phase^2, Xtwo_TAorTT_fourier_phase^2,
                 Xtwo_CCorCG_fourier_phase^2, Xtwo_CCorGC_fourier_phase^2,
                 Xtwo_CCorGG_fourier_phase^2, Xtwo_CGorGC_fourier_phase^2,
                 Xtwo_CGorGG_fourier_phase^2, Xtwo_GCorGG_fourier_phase^2,
                 dat_poly, n_gap_di,
                 Xone_A_to_C_ccf, Xone_A_to_G_ccf, Xone_A_to_T_ccf, 
                 Xone_C_to_G_ccf, Xone_C_to_T_ccf, Xone_G_to_T_ccf,
                 Xone_AorT_to_CorG_ccf)

dat_temp_test = cbind(dat_test %>% 
                        # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                        #        all_of(dinucleotides), all_of(trinucleotides), C0),
                        select(all_of(ps2), all_of(trinucleotides), C0),
                      # select(all_of(trinucleotides), C0),
                      # select(C0),
                      Xone_A_fourier_mag_test^2, Xone_C_fourier_mag_test^2, Xone_G_fourier_mag_test^2, Xone_T_fourier_mag_test^2,
                      Xone_A_fourier_phase_test^2, Xone_C_fourier_phase_test^2, Xone_G_fourier_phase_test^2, Xone_T_fourier_phase_test^2,
                      Xone_AorT_fourier_mag_test^2, 
                      Xone_CorG_fourier_mag_test^2,
                      Xone_AorT_fourier_phase_test^2, 
                      Xone_CorG_fourier_phase_test^2, 
                      Xtwo_AA_fourier_mag_test^2, Xtwo_AT_fourier_mag_test^2,
                      Xtwo_CC_fourier_mag_test^2, Xtwo_CG_fourier_mag_test^2,
                      Xtwo_GC_fourier_mag_test^2, Xtwo_GG_fourier_mag_test^2,
                      Xtwo_TA_fourier_mag_test^2, Xtwo_TT_fourier_mag_test^2,
                      Xtwo_AA_fourier_phase_test^2, Xtwo_AT_fourier_phase_test^2,
                      Xtwo_CC_fourier_phase_test^2, Xtwo_CG_fourier_phase_test^2,
                      Xtwo_GC_fourier_phase_test^2, Xtwo_GG_fourier_phase_test^2,
                      Xtwo_TA_fourier_phase_test^2, Xtwo_TT_fourier_phase_test^2,
                      Xtwo_AAorAT_fourier_mag_test^2, Xtwo_AAorTA_fourier_mag_test^2,
                      Xtwo_AAorTT_fourier_mag_test^2, Xtwo_ATorTA_fourier_mag_test^2,
                      Xtwo_ATorTT_fourier_mag_test^2, Xtwo_TAorTT_fourier_mag_test^2,
                      Xtwo_CCorCG_fourier_mag_test^2, Xtwo_CCorGC_fourier_mag_test^2,
                      Xtwo_CCorGG_fourier_mag_test^2, Xtwo_CGorGC_fourier_mag_test^2,
                      Xtwo_CGorGG_fourier_mag_test^2, Xtwo_GCorGG_fourier_mag_test^2,
                      Xtwo_AAorAT_fourier_phase_test^2, Xtwo_AAorTA_fourier_phase_test^2,
                      Xtwo_AAorTT_fourier_phase_test^2, Xtwo_ATorTA_fourier_phase_test^2,
                      Xtwo_ATorTT_fourier_phase_test^2, Xtwo_TAorTT_fourier_phase_test^2,
                      Xtwo_CCorCG_fourier_phase_test^2, Xtwo_CCorGC_fourier_phase_test^2,
                      Xtwo_CCorGG_fourier_phase_test^2, Xtwo_CGorGC_fourier_phase_test^2,
                      Xtwo_CGorGG_fourier_phase_test^2, Xtwo_GCorGG_fourier_phase_test^2,
                      dat_poly_test, n_gap_di_test,
                      Xone_A_to_C_ccf_test, Xone_A_to_G_ccf_test, Xone_A_to_T_ccf_test, 
                      Xone_C_to_G_ccf_test, Xone_C_to_T_ccf_test, Xone_G_to_T_ccf_test,
                      Xone_AorT_to_CorG_ccf_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                              #        all_of(dinucleotides), all_of(trinucleotides), C0),
                              select(all_of(ps2), all_of(trinucleotides), C0),
                            # select(all_of(trinucleotides), C0),
                            # select(C0),
                            Xone_A_fourier_mag_random_all^2, Xone_C_fourier_mag_random_all^2, Xone_G_fourier_mag_random_all^2, Xone_T_fourier_mag_random_all^2,
                            Xone_A_fourier_phase_random_all^2, Xone_C_fourier_phase_random_all^2, Xone_G_fourier_phase_random_all^2, Xone_T_fourier_phase_random_all^2,
                            Xone_AorT_fourier_mag_random_all^2, 
                            Xone_CorG_fourier_mag_random_all^2,
                            Xone_AorT_fourier_phase_random_all^2, 
                            Xone_CorG_fourier_phase_random_all^2, 
                            Xtwo_AA_fourier_mag_random_all^2, Xtwo_AT_fourier_mag_random_all^2,
                            Xtwo_CC_fourier_mag_random_all^2, Xtwo_CG_fourier_mag_random_all^2,
                            Xtwo_GC_fourier_mag_random_all^2, Xtwo_GG_fourier_mag_random_all^2,
                            Xtwo_TA_fourier_mag_random_all^2, Xtwo_TT_fourier_mag_random_all^2,
                            Xtwo_AA_fourier_phase_random_all^2, Xtwo_AT_fourier_phase_random_all^2,
                            Xtwo_CC_fourier_phase_random_all^2, Xtwo_CG_fourier_phase_random_all^2,
                            Xtwo_GC_fourier_phase_random_all^2, Xtwo_GG_fourier_phase_random_all^2,
                            Xtwo_TA_fourier_phase_random_all^2, Xtwo_TT_fourier_phase_random_all^2,
                            Xtwo_AAorAT_fourier_mag_random_all^2, Xtwo_AAorTA_fourier_mag_random_all^2,
                            Xtwo_AAorTT_fourier_mag_random_all^2, Xtwo_ATorTA_fourier_mag_random_all^2,
                            Xtwo_ATorTT_fourier_mag_random_all^2, Xtwo_TAorTT_fourier_mag_random_all^2,
                            Xtwo_CCorCG_fourier_mag_random_all^2, Xtwo_CCorGC_fourier_mag_random_all^2,
                            Xtwo_CCorGG_fourier_mag_random_all^2, Xtwo_CGorGC_fourier_mag_random_all^2,
                            Xtwo_CGorGG_fourier_mag_random_all^2, Xtwo_GCorGG_fourier_mag_random_all^2,
                            Xtwo_AAorAT_fourier_phase_random_all^2, Xtwo_AAorTA_fourier_phase_random_all^2,
                            Xtwo_AAorTT_fourier_phase_random_all^2, Xtwo_ATorTA_fourier_phase_random_all^2,
                            Xtwo_ATorTT_fourier_phase_random_all^2, Xtwo_TAorTT_fourier_phase_random_all^2,
                            Xtwo_CCorCG_fourier_phase_random_all^2, Xtwo_CCorGC_fourier_phase_random_all^2,
                            Xtwo_CCorGG_fourier_phase_random_all^2, Xtwo_CGorGC_fourier_phase_random_all^2,
                            Xtwo_CGorGG_fourier_phase_random_all^2, Xtwo_GCorGG_fourier_phase_random_all^2,
                            dat_poly_random_all, n_gap_di_random_all,
                            Xone_A_to_C_ccf_random_all, Xone_A_to_G_ccf_random_all, Xone_A_to_T_ccf_random_all, 
                            Xone_C_to_G_ccf_random_all, Xone_C_to_T_ccf_random_all, Xone_G_to_T_ccf_random_all,
                            Xone_AorT_to_CorG_ccf_random_all)

temp_lm = lm(C0~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.6563554
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.6576925
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.6022612
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LM")
abline(0,1)




# Fourier{Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides} + all polyA/T + N-Gapped Dinucleotides + 
# ccf{Nucleotides to nucleotides + Pairs of nucleotides to pairs of nucleotides +
# Groups of dinucleotides to groups of dinucleotides}

dat_temp = cbind(dat %>% 
                   # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                   #        all_of(dinucleotides), all_of(trinucleotides), C0),
                   select(all_of(ps2), all_of(trinucleotides), C0),
                 # select(all_of(trinucleotides), C0),
                 # select(C0),
                 Xone_A_fourier_mag^2, Xone_C_fourier_mag^2, Xone_G_fourier_mag^2, Xone_T_fourier_mag^2,
                 Xone_A_fourier_phase^2, Xone_C_fourier_phase^2, Xone_G_fourier_phase^2, Xone_T_fourier_phase^2,
                 Xone_AorT_fourier_mag^2, 
                 Xone_CorG_fourier_mag^2,
                 Xone_AorT_fourier_phase^2, 
                 Xone_CorG_fourier_phase^2, 
                 Xtwo_AA_fourier_mag^2, Xtwo_AT_fourier_mag^2,
                 Xtwo_CC_fourier_mag^2, Xtwo_CG_fourier_mag^2,
                 Xtwo_GC_fourier_mag^2, Xtwo_GG_fourier_mag^2,
                 Xtwo_TA_fourier_mag^2, Xtwo_TT_fourier_mag^2,
                 Xtwo_AA_fourier_phase^2, Xtwo_AT_fourier_phase^2,
                 Xtwo_CC_fourier_phase^2, Xtwo_CG_fourier_phase^2,
                 Xtwo_GC_fourier_phase^2, Xtwo_GG_fourier_phase^2,
                 Xtwo_TA_fourier_phase^2, Xtwo_TT_fourier_phase^2,
                 Xtwo_AAorAT_fourier_mag^2, Xtwo_AAorTA_fourier_mag^2,
                 Xtwo_AAorTT_fourier_mag^2, Xtwo_ATorTA_fourier_mag^2,
                 Xtwo_ATorTT_fourier_mag^2, Xtwo_TAorTT_fourier_mag^2,
                 Xtwo_CCorCG_fourier_mag^2, Xtwo_CCorGC_fourier_mag^2,
                 Xtwo_CCorGG_fourier_mag^2, Xtwo_CGorGC_fourier_mag^2,
                 Xtwo_CGorGG_fourier_mag^2, Xtwo_GCorGG_fourier_mag^2,
                 Xtwo_AAorAT_fourier_phase^2, Xtwo_AAorTA_fourier_phase^2,
                 Xtwo_AAorTT_fourier_phase^2, Xtwo_ATorTA_fourier_phase^2,
                 Xtwo_ATorTT_fourier_phase^2, Xtwo_TAorTT_fourier_phase^2,
                 Xtwo_CCorCG_fourier_phase^2, Xtwo_CCorGC_fourier_phase^2,
                 Xtwo_CCorGG_fourier_phase^2, Xtwo_CGorGC_fourier_phase^2,
                 Xtwo_CGorGG_fourier_phase^2, Xtwo_GCorGG_fourier_phase^2,
                 dat_poly, n_gap_di,
                 Xone_A_to_C_ccf, Xone_A_to_G_ccf, Xone_A_to_T_ccf, 
                 Xone_C_to_G_ccf, Xone_C_to_T_ccf, Xone_G_to_T_ccf,
                 Xone_AorT_to_CorG_ccf, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf)

dat_temp_test = cbind(dat_test %>% 
                        # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                        #        all_of(dinucleotides), all_of(trinucleotides), C0),
                        select(all_of(ps2), all_of(trinucleotides), C0),
                      # select(all_of(trinucleotides), C0),
                      # select(C0),
                      Xone_A_fourier_mag_test^2, Xone_C_fourier_mag_test^2, Xone_G_fourier_mag_test^2, Xone_T_fourier_mag_test^2,
                      Xone_A_fourier_phase_test^2, Xone_C_fourier_phase_test^2, Xone_G_fourier_phase_test^2, Xone_T_fourier_phase_test^2,
                      Xone_AorT_fourier_mag_test^2, 
                      Xone_CorG_fourier_mag_test^2,
                      Xone_AorT_fourier_phase_test^2, 
                      Xone_CorG_fourier_phase_test^2, 
                      Xtwo_AA_fourier_mag_test^2, Xtwo_AT_fourier_mag_test^2,
                      Xtwo_CC_fourier_mag_test^2, Xtwo_CG_fourier_mag_test^2,
                      Xtwo_GC_fourier_mag_test^2, Xtwo_GG_fourier_mag_test^2,
                      Xtwo_TA_fourier_mag_test^2, Xtwo_TT_fourier_mag_test^2,
                      Xtwo_AA_fourier_phase_test^2, Xtwo_AT_fourier_phase_test^2,
                      Xtwo_CC_fourier_phase_test^2, Xtwo_CG_fourier_phase_test^2,
                      Xtwo_GC_fourier_phase_test^2, Xtwo_GG_fourier_phase_test^2,
                      Xtwo_TA_fourier_phase_test^2, Xtwo_TT_fourier_phase_test^2,
                      Xtwo_AAorAT_fourier_mag_test^2, Xtwo_AAorTA_fourier_mag_test^2,
                      Xtwo_AAorTT_fourier_mag_test^2, Xtwo_ATorTA_fourier_mag_test^2,
                      Xtwo_ATorTT_fourier_mag_test^2, Xtwo_TAorTT_fourier_mag_test^2,
                      Xtwo_CCorCG_fourier_mag_test^2, Xtwo_CCorGC_fourier_mag_test^2,
                      Xtwo_CCorGG_fourier_mag_test^2, Xtwo_CGorGC_fourier_mag_test^2,
                      Xtwo_CGorGG_fourier_mag_test^2, Xtwo_GCorGG_fourier_mag_test^2,
                      Xtwo_AAorAT_fourier_phase_test^2, Xtwo_AAorTA_fourier_phase_test^2,
                      Xtwo_AAorTT_fourier_phase_test^2, Xtwo_ATorTA_fourier_phase_test^2,
                      Xtwo_ATorTT_fourier_phase_test^2, Xtwo_TAorTT_fourier_phase_test^2,
                      Xtwo_CCorCG_fourier_phase_test^2, Xtwo_CCorGC_fourier_phase_test^2,
                      Xtwo_CCorGG_fourier_phase_test^2, Xtwo_CGorGC_fourier_phase_test^2,
                      Xtwo_CGorGG_fourier_phase_test^2, Xtwo_GCorGG_fourier_phase_test^2,
                      dat_poly_test, n_gap_di_test,
                      Xone_A_to_C_ccf_test, Xone_A_to_G_ccf_test, Xone_A_to_T_ccf_test, 
                      Xone_C_to_G_ccf_test, Xone_C_to_T_ccf_test, Xone_G_to_T_ccf_test,
                      Xone_AorT_to_CorG_ccf_test, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                              #        all_of(dinucleotides), all_of(trinucleotides), C0),
                              select(all_of(ps2), all_of(trinucleotides), C0),
                            # select(all_of(trinucleotides), C0),
                            # select(C0),
                            Xone_A_fourier_mag_random_all^2, Xone_C_fourier_mag_random_all^2, Xone_G_fourier_mag_random_all^2, Xone_T_fourier_mag_random_all^2,
                            Xone_A_fourier_phase_random_all^2, Xone_C_fourier_phase_random_all^2, Xone_G_fourier_phase_random_all^2, Xone_T_fourier_phase_random_all^2,
                            Xone_AorT_fourier_mag_random_all^2, 
                            Xone_CorG_fourier_mag_random_all^2,
                            Xone_AorT_fourier_phase_random_all^2, 
                            Xone_CorG_fourier_phase_random_all^2, 
                            Xtwo_AA_fourier_mag_random_all^2, Xtwo_AT_fourier_mag_random_all^2,
                            Xtwo_CC_fourier_mag_random_all^2, Xtwo_CG_fourier_mag_random_all^2,
                            Xtwo_GC_fourier_mag_random_all^2, Xtwo_GG_fourier_mag_random_all^2,
                            Xtwo_TA_fourier_mag_random_all^2, Xtwo_TT_fourier_mag_random_all^2,
                            Xtwo_AA_fourier_phase_random_all^2, Xtwo_AT_fourier_phase_random_all^2,
                            Xtwo_CC_fourier_phase_random_all^2, Xtwo_CG_fourier_phase_random_all^2,
                            Xtwo_GC_fourier_phase_random_all^2, Xtwo_GG_fourier_phase_random_all^2,
                            Xtwo_TA_fourier_phase_random_all^2, Xtwo_TT_fourier_phase_random_all^2,
                            Xtwo_AAorAT_fourier_mag_random_all^2, Xtwo_AAorTA_fourier_mag_random_all^2,
                            Xtwo_AAorTT_fourier_mag_random_all^2, Xtwo_ATorTA_fourier_mag_random_all^2,
                            Xtwo_ATorTT_fourier_mag_random_all^2, Xtwo_TAorTT_fourier_mag_random_all^2,
                            Xtwo_CCorCG_fourier_mag_random_all^2, Xtwo_CCorGC_fourier_mag_random_all^2,
                            Xtwo_CCorGG_fourier_mag_random_all^2, Xtwo_CGorGC_fourier_mag_random_all^2,
                            Xtwo_CGorGG_fourier_mag_random_all^2, Xtwo_GCorGG_fourier_mag_random_all^2,
                            Xtwo_AAorAT_fourier_phase_random_all^2, Xtwo_AAorTA_fourier_phase_random_all^2,
                            Xtwo_AAorTT_fourier_phase_random_all^2, Xtwo_ATorTA_fourier_phase_random_all^2,
                            Xtwo_ATorTT_fourier_phase_random_all^2, Xtwo_TAorTT_fourier_phase_random_all^2,
                            Xtwo_CCorCG_fourier_phase_random_all^2, Xtwo_CCorGC_fourier_phase_random_all^2,
                            Xtwo_CCorGG_fourier_phase_random_all^2, Xtwo_CGorGC_fourier_phase_random_all^2,
                            Xtwo_CGorGG_fourier_phase_random_all^2, Xtwo_GCorGG_fourier_phase_random_all^2,
                            dat_poly_random_all, n_gap_di_random_all,
                            Xone_A_to_C_ccf_random_all, Xone_A_to_G_ccf_random_all, Xone_A_to_T_ccf_random_all, 
                            Xone_C_to_G_ccf_random_all, Xone_C_to_T_ccf_random_all, Xone_G_to_T_ccf_random_all,
                            Xone_AorT_to_CorG_ccf_random_all, Xtwo_AAorATorTAorTT_to_CCorCGorGCorGG_ccf_random_all)

temp_lm = lm(C0~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.6578441
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.6592657
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.6035103
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LM")
abline(0,1)

