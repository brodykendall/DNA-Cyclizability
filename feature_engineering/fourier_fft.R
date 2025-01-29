dat = readRDS("data/Created/processed_tiling_newC0.rds")
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
dat_random = readRDS("data/Created/processed_random_newC0.rds")
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")
dat_chrV = readRDS("data/Created/processed_chrV_newC0.rds")
dat_chrV_test = readRDS("data/Created/processed_chrV_test_newC0.rds")
dat_yeast = readRDS("data/Created/processed_yeast_newC0.rds")
dat_yeast_test = readRDS("data/Created/processed_yeast_test_newC0.rds")




########################################################################
# Tiling
########################################################################

# Train data
ps1 <- paste0("X", 1:50, "mono")

Xone = dat %>% select(all_of(ps1))
Xone_AorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AorT) = colnames(Xone)
colnames(Xone_CorG) = colnames(Xone)
Xone_AorT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier = t(apply(Xone_AorT, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos = Re(Xone_AorT_fourier)
Xone_AorT_fourier_sin = Im(Xone_AorT_fourier)

colnames(Xone_AorT_fourier_cos) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin) = paste0("fouriersinAorT", c(1:25, 50))

colnames(Xone_CorG_fourier_cos) = paste0("fouriercosCorG", c(1:25, 50))
colnames(Xone_CorG_fourier_sin) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos, "data/Created/tiling_Xone_AorT_fourier_cos.rds")
saveRDS(Xone_AorT_fourier_sin, "data/Created/tiling_Xone_AorT_fourier_sin.rds")

saveRDS(Xone_CorG_fourier_cos, "data/Created/tiling_Xone_CorG_fourier_cos.rds")
saveRDS(Xone_CorG_fourier_sin, "data/Created/tiling_Xone_CorG_fourier_sin.rds")

ps2 <- paste0("X", 1:49, "di")

Xtwo = dat %>% select(all_of(ps2))
Xtwo_AorT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CorG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AorT) = colnames(Xtwo)
colnames(Xtwo_CorG) = colnames(Xtwo)
Xtwo_AorT[] = ((Xtwo == "AA") | (Xtwo == "TT") | (Xtwo == "AT") | (Xtwo == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CorG[] = ((Xtwo == "CC") | (Xtwo == "GG") | (Xtwo == "CG") | (Xtwo == "GC")) %>%
  as.matrix() %>% as.numeric()

Xtwo_AorT_fourier = t(apply(Xtwo_AorT, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos = Re(Xtwo_AorT_fourier)
Xtwo_AorT_fourier_sin = Im(Xtwo_AorT_fourier)

colnames(Xtwo_AorT_fourier_cos) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin) = paste0("fouriersinAorT2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos, "data/Created/tiling_Xtwo_AorT_fourier_cos.rds")
saveRDS(Xtwo_AorT_fourier_sin, "data/Created/tiling_Xtwo_AorT_fourier_sin.rds")



Xtwo_AA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AA) = colnames(Xtwo)
Xtwo_AA[] = (Xtwo == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier = t(apply(Xtwo_AA, 1, fft))[,1:25]
colnames(AA_fourier) = paste0("fourierAA_", 1:25)
AA_fourier_cos = Re(AA_fourier)
AA_fourier_sin = Im(AA_fourier)
colnames(AA_fourier_cos) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos, "data/Created/tiling_AA_fourier_cos.rds")
saveRDS(AA_fourier_sin, "data/Created/tiling_AA_fourier_sin.rds")

Xtwo_AC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AC) = colnames(Xtwo)
Xtwo_AC[] = (Xtwo == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier = t(apply(Xtwo_AC, 1, fft))[,1:25]
colnames(AC_fourier) = paste0("fourierAC_", 1:25)
AC_fourier_cos = Re(AC_fourier)
AC_fourier_sin = Im(AC_fourier)
colnames(AC_fourier_cos) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos, "data/Created/tiling_AC_fourier_cos.rds")
saveRDS(AC_fourier_sin, "data/Created/tiling_AC_fourier_sin.rds")

Xtwo_AG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AG) = colnames(Xtwo)
Xtwo_AG[] = (Xtwo == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier = t(apply(Xtwo_AG, 1, fft))[,1:25]
colnames(AG_fourier) = paste0("fourierAG_", 1:25)
AG_fourier_cos = Re(AG_fourier)
AG_fourier_sin = Im(AG_fourier)
colnames(AG_fourier_cos) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos, "data/Created/tiling_AG_fourier_cos.rds")
saveRDS(AG_fourier_sin, "data/Created/tiling_AG_fourier_sin.rds")

Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AT) = colnames(Xtwo)
Xtwo_AT[] = (Xtwo == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier = t(apply(Xtwo_AT, 1, fft))[,1:25]
colnames(AT_fourier) = paste0("fourierAT_", 1:25)
AT_fourier_cos = Re(AT_fourier)
AT_fourier_sin = Im(AT_fourier)
colnames(AT_fourier_cos) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos, "data/Created/tiling_AT_fourier_cos.rds")
saveRDS(AT_fourier_sin, "data/Created/tiling_AT_fourier_sin.rds")

Xtwo_CA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CA) = colnames(Xtwo)
Xtwo_CA[] = (Xtwo == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier = t(apply(Xtwo_CA, 1, fft))[,1:25]
colnames(CA_fourier) = paste0("fourierCA_", 1:25)
CA_fourier_cos = Re(CA_fourier)
CA_fourier_sin = Im(CA_fourier)
colnames(CA_fourier_cos) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos, "data/Created/tiling_CA_fourier_cos.rds")
saveRDS(CA_fourier_sin, "data/Created/tiling_CA_fourier_sin.rds")

Xtwo_CC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CC) = colnames(Xtwo)
Xtwo_CC[] = (Xtwo == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier = t(apply(Xtwo_CC, 1, fft))[,1:25]
colnames(CC_fourier) = paste0("fourierCC_", 1:25)
CC_fourier_cos = Re(CC_fourier)
CC_fourier_sin = Im(CC_fourier)
colnames(CC_fourier_cos) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos, "data/Created/tiling_CC_fourier_cos.rds")
saveRDS(CC_fourier_sin, "data/Created/tiling_CC_fourier_sin.rds")

Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CG) = colnames(Xtwo)
Xtwo_CG[] = (Xtwo == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier = t(apply(Xtwo_CG, 1, fft))[,1:25]
colnames(CG_fourier) = paste0("fourierCG_", 1:25)
CG_fourier_cos = Re(CG_fourier)
CG_fourier_sin = Im(CG_fourier)
colnames(CG_fourier_cos) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos, "data/Created/tiling_CG_fourier_cos.rds")
saveRDS(CG_fourier_sin, "data/Created/tiling_CG_fourier_sin.rds")

Xtwo_CT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CT) = colnames(Xtwo)
Xtwo_CT[] = (Xtwo == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier = t(apply(Xtwo_CT, 1, fft))[,1:25]
colnames(CT_fourier) = paste0("fourierCT_", 1:25)
CT_fourier_cos = Re(CT_fourier)
CT_fourier_sin = Im(CT_fourier)
colnames(CT_fourier_cos) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos, "data/Created/tiling_CT_fourier_cos.rds")
saveRDS(CT_fourier_sin, "data/Created/tiling_CT_fourier_sin.rds")

Xtwo_GA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GA) = colnames(Xtwo)
Xtwo_GA[] = (Xtwo == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier = t(apply(Xtwo_GA, 1, fft))[,1:25]
colnames(GA_fourier) = paste0("fourierGA_", 1:25)
GA_fourier_cos = Re(GA_fourier)
GA_fourier_sin = Im(GA_fourier)
colnames(GA_fourier_cos) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos, "data/Created/tiling_GA_fourier_cos.rds")
saveRDS(GA_fourier_sin, "data/Created/tiling_GA_fourier_sin.rds")

Xtwo_GC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GC) = colnames(Xtwo)
Xtwo_GC[] = (Xtwo == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier = t(apply(Xtwo_GC, 1, fft))[,1:25]
colnames(GC_fourier) = paste0("fourierGC_", 1:25)
GC_fourier_cos = Re(GC_fourier)
GC_fourier_sin = Im(GC_fourier)
colnames(GC_fourier_cos) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos, "data/Created/tiling_GC_fourier_cos.rds")
saveRDS(GC_fourier_sin, "data/Created/tiling_GC_fourier_sin.rds")

Xtwo_GG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GG) = colnames(Xtwo)
Xtwo_GG[] = (Xtwo == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier = t(apply(Xtwo_GG, 1, fft))[,1:25]
colnames(GG_fourier) = paste0("fourierGG_", 1:25)
GG_fourier_cos = Re(GG_fourier)
GG_fourier_sin = Im(GG_fourier)
colnames(GG_fourier_cos) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos, "data/Created/tiling_GG_fourier_cos.rds")
saveRDS(GG_fourier_sin, "data/Created/tiling_GG_fourier_sin.rds")

Xtwo_GT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GT) = colnames(Xtwo)
Xtwo_GT[] = (Xtwo == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier = t(apply(Xtwo_GT, 1, fft))[,1:25]
colnames(GT_fourier) = paste0("fourierGT_", 1:25)
GT_fourier_cos = Re(GT_fourier)
GT_fourier_sin = Im(GT_fourier)
colnames(GT_fourier_cos) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos, "data/Created/tiling_GT_fourier_cos.rds")
saveRDS(GT_fourier_sin, "data/Created/tiling_GT_fourier_sin.rds")

Xtwo_TA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TA) = colnames(Xtwo)
Xtwo_TA[] = (Xtwo == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier = t(apply(Xtwo_TA, 1, fft))[,1:25]
colnames(TA_fourier) = paste0("fourierTA_", 1:25)
TA_fourier_cos = Re(TA_fourier)
TA_fourier_sin = Im(TA_fourier)
colnames(TA_fourier_cos) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos, "data/Created/tiling_TA_fourier_cos.rds")
saveRDS(TA_fourier_sin, "data/Created/tiling_TA_fourier_sin.rds")

Xtwo_TC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TC) = colnames(Xtwo)
Xtwo_TC[] = (Xtwo == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier = t(apply(Xtwo_TC, 1, fft))[,1:25]
colnames(TC_fourier) = paste0("fourierTC_", 1:25)
TC_fourier_cos = Re(TC_fourier)
TC_fourier_sin = Im(TC_fourier)
colnames(TC_fourier_cos) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos, "data/Created/tiling_TC_fourier_cos.rds")
saveRDS(TC_fourier_sin, "data/Created/tiling_TC_fourier_sin.rds")

Xtwo_TG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TG) = colnames(Xtwo)
Xtwo_TG[] = (Xtwo == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier = t(apply(Xtwo_TG, 1, fft))[,1:25]
colnames(TG_fourier) = paste0("fourierTG_", 1:25)
TG_fourier_cos = Re(TG_fourier)
TG_fourier_sin = Im(TG_fourier)
colnames(TG_fourier_cos) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos, "data/Created/tiling_TG_fourier_cos.rds")
saveRDS(TG_fourier_sin, "data/Created/tiling_TG_fourier_sin.rds")

Xtwo_TT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TT) = colnames(Xtwo)
Xtwo_TT[] = (Xtwo == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier = t(apply(Xtwo_TT, 1, fft))[,1:25]
colnames(TT_fourier) = paste0("fourierTT_", 1:25)
TT_fourier_cos = Re(TT_fourier)
TT_fourier_sin = Im(TT_fourier)
colnames(TT_fourier_cos) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos, "data/Created/tiling_TT_fourier_cos.rds")
saveRDS(TT_fourier_sin, "data/Created/tiling_TT_fourier_sin.rds")

fourier_di_cos = cbind(AA_fourier_cos, AC_fourier_cos, AG_fourier_cos, AT_fourier_cos, 
                       CA_fourier_cos, CC_fourier_cos, CG_fourier_cos, CT_fourier_cos,
                       GA_fourier_cos, GC_fourier_cos, GG_fourier_cos, GT_fourier_cos,
                       TA_fourier_cos, TC_fourier_cos, TG_fourier_cos, TT_fourier_cos)
fourier_di_sin = cbind(AA_fourier_sin, AC_fourier_sin, AG_fourier_sin, AT_fourier_sin, 
                       CA_fourier_sin, CC_fourier_sin, CG_fourier_sin, CT_fourier_sin,
                       GA_fourier_sin, GC_fourier_sin, GG_fourier_sin, GT_fourier_sin,
                       TA_fourier_sin, TC_fourier_sin, TG_fourier_sin, TT_fourier_sin)
saveRDS(fourier_di_cos, "data/Created/tiling_fourier_di_cos.rds")
saveRDS(fourier_di_sin, "data/Created/tiling_fourier_di_sin.rds")

fourier_AA_AC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AC_fourier[,"fourierAC_6"])
fourier_AA_AG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AG_fourier[,"fourierAG_6"])
fourier_AA_AT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AT_fourier[,"fourierAT_6"])
fourier_AA_CA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CA_fourier[,"fourierCA_6"])
fourier_AA_CC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CC_fourier[,"fourierCC_6"])
fourier_AA_CG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CG_fourier[,"fourierCG_6"])
fourier_AA_CT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CT_fourier[,"fourierCT_6"])
fourier_AA_GA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GA_fourier[,"fourierGA_6"])
fourier_AA_GC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GC_fourier[,"fourierGC_6"])
fourier_AA_GG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GG_fourier[,"fourierGG_6"])
fourier_AA_GT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GT_fourier[,"fourierGT_6"])
fourier_AA_TA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TA_fourier[,"fourierTA_6"])
fourier_AA_TC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TC_fourier[,"fourierTC_6"])
fourier_AA_TG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TG_fourier[,"fourierTG_6"])
fourier_AA_TT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TT_fourier[,"fourierTT_6"])

fourier_AA_6_mod = paste0("fourier_AA_", dinucleotides[-1], "_6_mod")
fourier_AA_6_mod_tiling = data.frame(map(fourier_AA_6_mod, get))
colnames(fourier_AA_6_mod_tiling) = fourier_AA_6_mod

fourier_AC_AG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AG_fourier[,"fourierAG_6"])
fourier_AC_AT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AT_fourier[,"fourierAT_6"])
fourier_AC_CA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CA_fourier[,"fourierCA_6"])
fourier_AC_CC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CC_fourier[,"fourierCC_6"])
fourier_AC_CG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CG_fourier[,"fourierCG_6"])
fourier_AC_CT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CT_fourier[,"fourierCT_6"])
fourier_AC_GA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GA_fourier[,"fourierGA_6"])
fourier_AC_GC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GC_fourier[,"fourierGC_6"])
fourier_AC_GG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GG_fourier[,"fourierGG_6"])
fourier_AC_GT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GT_fourier[,"fourierGT_6"])
fourier_AC_TA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TA_fourier[,"fourierTA_6"])
fourier_AC_TC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TC_fourier[,"fourierTC_6"])
fourier_AC_TG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TG_fourier[,"fourierTG_6"])
fourier_AC_TT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TT_fourier[,"fourierTT_6"])

fourier_AC_6_mod = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod")
fourier_AC_6_mod_tiling = data.frame(map(fourier_AC_6_mod, get))
colnames(fourier_AC_6_mod_tiling) = fourier_AC_6_mod

fourier_AG_AT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - AT_fourier[,"fourierAT_6"])
fourier_AG_CA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CA_fourier[,"fourierCA_6"])
fourier_AG_CC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CC_fourier[,"fourierCC_6"])
fourier_AG_CG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CG_fourier[,"fourierCG_6"])
fourier_AG_CT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CT_fourier[,"fourierCT_6"])
fourier_AG_GA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GA_fourier[,"fourierGA_6"])
fourier_AG_GC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GC_fourier[,"fourierGC_6"])
fourier_AG_GG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GG_fourier[,"fourierGG_6"])
fourier_AG_GT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GT_fourier[,"fourierGT_6"])
fourier_AG_TA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TA_fourier[,"fourierTA_6"])
fourier_AG_TC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TC_fourier[,"fourierTC_6"])
fourier_AG_TG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TG_fourier[,"fourierTG_6"])
fourier_AG_TT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TT_fourier[,"fourierTT_6"])

fourier_AG_6_mod = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod")
fourier_AG_6_mod_tiling = data.frame(map(fourier_AG_6_mod, get))
colnames(fourier_AG_6_mod_tiling) = fourier_AG_6_mod

fourier_AT_CA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CA_fourier[,"fourierCA_6"])
fourier_AT_CC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CC_fourier[,"fourierCC_6"])
fourier_AT_CG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CG_fourier[,"fourierCG_6"])
fourier_AT_CT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CT_fourier[,"fourierCT_6"])
fourier_AT_GA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GA_fourier[,"fourierGA_6"])
fourier_AT_GC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GC_fourier[,"fourierGC_6"])
fourier_AT_GG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GG_fourier[,"fourierGG_6"])
fourier_AT_GT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GT_fourier[,"fourierGT_6"])
fourier_AT_TA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TA_fourier[,"fourierTA_6"])
fourier_AT_TC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TC_fourier[,"fourierTC_6"])
fourier_AT_TG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TG_fourier[,"fourierTG_6"])
fourier_AT_TT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TT_fourier[,"fourierTT_6"])

fourier_AT_6_mod = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod")
fourier_AT_6_mod_tiling = data.frame(map(fourier_AT_6_mod, get))
colnames(fourier_AT_6_mod_tiling) = fourier_AT_6_mod

fourier_CA_CC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CC_fourier[,"fourierCC_6"])
fourier_CA_CG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CG_fourier[,"fourierCG_6"])
fourier_CA_CT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CT_fourier[,"fourierCT_6"])
fourier_CA_GA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GA_fourier[,"fourierGA_6"])
fourier_CA_GC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GC_fourier[,"fourierGC_6"])
fourier_CA_GG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GG_fourier[,"fourierGG_6"])
fourier_CA_GT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GT_fourier[,"fourierGT_6"])
fourier_CA_TA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TA_fourier[,"fourierTA_6"])
fourier_CA_TC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TC_fourier[,"fourierTC_6"])
fourier_CA_TG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TG_fourier[,"fourierTG_6"])
fourier_CA_TT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TT_fourier[,"fourierTT_6"])

fourier_CA_6_mod = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod")
fourier_CA_6_mod_tiling = data.frame(map(fourier_CA_6_mod, get))
colnames(fourier_CA_6_mod_tiling) = fourier_CA_6_mod

fourier_CC_CG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CG_fourier[,"fourierCG_6"])
fourier_CC_CT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CT_fourier[,"fourierCT_6"])
fourier_CC_GA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GA_fourier[,"fourierGA_6"])
fourier_CC_GC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GC_fourier[,"fourierGC_6"])
fourier_CC_GG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GG_fourier[,"fourierGG_6"])
fourier_CC_GT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GT_fourier[,"fourierGT_6"])
fourier_CC_TA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TA_fourier[,"fourierTA_6"])
fourier_CC_TC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TC_fourier[,"fourierTC_6"])
fourier_CC_TG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TG_fourier[,"fourierTG_6"])
fourier_CC_TT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TT_fourier[,"fourierTT_6"])

fourier_CC_6_mod = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod")
fourier_CC_6_mod_tiling = data.frame(map(fourier_CC_6_mod, get))
colnames(fourier_CC_6_mod_tiling) = fourier_CC_6_mod

fourier_CG_CT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - CT_fourier[,"fourierCT_6"])
fourier_CG_GA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GA_fourier[,"fourierGA_6"])
fourier_CG_GC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GC_fourier[,"fourierGC_6"])
fourier_CG_GG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GG_fourier[,"fourierGG_6"])
fourier_CG_GT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GT_fourier[,"fourierGT_6"])
fourier_CG_TA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TA_fourier[,"fourierTA_6"])
fourier_CG_TC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TC_fourier[,"fourierTC_6"])
fourier_CG_TG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TG_fourier[,"fourierTG_6"])
fourier_CG_TT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TT_fourier[,"fourierTT_6"])

fourier_CG_6_mod = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod")
fourier_CG_6_mod_tiling = data.frame(map(fourier_CG_6_mod, get))
colnames(fourier_CG_6_mod_tiling) = fourier_CG_6_mod

fourier_CT_GA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GA_fourier[,"fourierGA_6"])
fourier_CT_GC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GC_fourier[,"fourierGC_6"])
fourier_CT_GG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GG_fourier[,"fourierGG_6"])
fourier_CT_GT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GT_fourier[,"fourierGT_6"])
fourier_CT_TA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TA_fourier[,"fourierTA_6"])
fourier_CT_TC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TC_fourier[,"fourierTC_6"])
fourier_CT_TG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TG_fourier[,"fourierTG_6"])
fourier_CT_TT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TT_fourier[,"fourierTT_6"])

fourier_CT_6_mod = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod")
fourier_CT_6_mod_tiling = data.frame(map(fourier_CT_6_mod, get))
colnames(fourier_CT_6_mod_tiling) = fourier_CT_6_mod

fourier_GA_GC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GC_fourier[,"fourierGC_6"])
fourier_GA_GG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GG_fourier[,"fourierGG_6"])
fourier_GA_GT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GT_fourier[,"fourierGT_6"])
fourier_GA_TA_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TA_fourier[,"fourierTA_6"])
fourier_GA_TC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TC_fourier[,"fourierTC_6"])
fourier_GA_TG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TG_fourier[,"fourierTG_6"])
fourier_GA_TT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TT_fourier[,"fourierTT_6"])

fourier_GA_6_mod = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod")
fourier_GA_6_mod_tiling = data.frame(map(fourier_GA_6_mod, get))
colnames(fourier_GA_6_mod_tiling) = fourier_GA_6_mod

fourier_GC_GG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GG_fourier[,"fourierGG_6"])
fourier_GC_GT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GT_fourier[,"fourierGT_6"])
fourier_GC_TA_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TA_fourier[,"fourierTA_6"])
fourier_GC_TC_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TC_fourier[,"fourierTC_6"])
fourier_GC_TG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TG_fourier[,"fourierTG_6"])
fourier_GC_TT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TT_fourier[,"fourierTT_6"])

fourier_GC_6_mod = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod")
fourier_GC_6_mod_tiling = data.frame(map(fourier_GC_6_mod, get))
colnames(fourier_GC_6_mod_tiling) = fourier_GC_6_mod

fourier_GG_GT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - GT_fourier[,"fourierGT_6"])
fourier_GG_TA_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TA_fourier[,"fourierTA_6"])
fourier_GG_TC_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TC_fourier[,"fourierTC_6"])
fourier_GG_TG_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TG_fourier[,"fourierTG_6"])
fourier_GG_TT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TT_fourier[,"fourierTT_6"])

fourier_GG_6_mod = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod")
fourier_GG_6_mod_tiling = data.frame(map(fourier_GG_6_mod, get))
colnames(fourier_GG_6_mod_tiling) = fourier_GG_6_mod

fourier_GT_TA_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TA_fourier[,"fourierTA_6"])
fourier_GT_TC_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TC_fourier[,"fourierTC_6"])
fourier_GT_TG_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TG_fourier[,"fourierTG_6"])
fourier_GT_TT_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TT_fourier[,"fourierTT_6"])

fourier_GT_6_mod = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod")
fourier_GT_6_mod_tiling = data.frame(map(fourier_GT_6_mod, get))
colnames(fourier_GT_6_mod_tiling) = fourier_GT_6_mod

fourier_TA_TC_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TC_fourier[,"fourierTC_6"])
fourier_TA_TG_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TG_fourier[,"fourierTG_6"])
fourier_TA_TT_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TT_fourier[,"fourierTT_6"])

fourier_TA_6_mod = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod")
fourier_TA_6_mod_tiling = data.frame(map(fourier_TA_6_mod, get))
colnames(fourier_TA_6_mod_tiling) = fourier_TA_6_mod

fourier_TC_TG_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TG_fourier[,"fourierTG_6"])
fourier_TC_TT_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TT_fourier[,"fourierTT_6"])

fourier_TC_6_mod = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod")
fourier_TC_6_mod_tiling = data.frame(map(fourier_TC_6_mod, get))
colnames(fourier_TC_6_mod_tiling) = fourier_TC_6_mod

fourier_TG_TT_6_mod = Mod(TG_fourier[,"fourierTG_6"] - TT_fourier[,"fourierTT_6"])

fourier_TG_6_mod = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod")
fourier_TG_6_mod_tiling = data.frame(map(fourier_TG_6_mod, get))
colnames(fourier_TG_6_mod_tiling) = fourier_TG_6_mod

fourier_6_mod = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_tiling"), get))
saveRDS(fourier_6_mod, "data/Created/tiling_fourier_6_mod.rds")




# Test data:
Xone_test = dat_test %>% select(all_of(ps1))
Xone_AorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AorT_test) = colnames(Xone_test)
colnames(Xone_CorG_test) = colnames(Xone_test)
Xone_AorT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_test = Re(Xone_AorT_fourier_test)
Xone_AorT_fourier_sin_test = Im(Xone_AorT_fourier_test)

Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,c(1:25, 50)]
Xone_CorG_fourier_cos_test = Re(Xone_CorG_fourier_test)
Xone_CorG_fourier_sin_test = Im(Xone_CorG_fourier_test)

colnames(Xone_AorT_fourier_cos_test) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_test) = paste0("fouriersinAorT", c(1:25, 50))

colnames(Xone_CorG_fourier_cos_test) = paste0("fouriercosCorG", c(1:25, 50))
colnames(Xone_CorG_fourier_sin_test) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_test, "data/Created/tiling_Xone_AorT_fourier_cos_test.rds")
saveRDS(Xone_AorT_fourier_sin_test, "data/Created/tiling_Xone_AorT_fourier_sin_test.rds")

saveRDS(Xone_CorG_fourier_cos_test, "data/Created/tiling_Xone_CorG_fourier_cos_test.rds")
saveRDS(Xone_CorG_fourier_sin_test, "data/Created/tiling_Xone_CorG_fourier_sin_test.rds")


Xtwo_test = dat_test %>% select(all_of(ps2))
Xtwo_AorT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
# Xtwo_CorG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AorT_test) = colnames(Xtwo_test)
# colnames(Xtwo_CorG_test) = colnames(Xtwo_test)
Xtwo_AorT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT") | (Xtwo_test == "AT") | (Xtwo_test == "TA")) %>% 
  as.matrix() %>% as.numeric()
# Xtwo_CorG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG") | (Xtwo_test == "CG") | (Xtwo_test == "GC")) %>%
#   as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_test = t(apply(Xtwo_AorT_test, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_test = Re(Xtwo_AorT_fourier_test)
Xtwo_AorT_fourier_sin_test = Im(Xtwo_AorT_fourier_test)

# Xtwo_CorG_fourier = t(apply(Xtwo_CorG_test, 1, fft))[,c(1:25, 50)]
# Xtwo_CorG_fourier_cos_test = Re(Xtwo_CorG_fourier)
# Xtwo_CorG_fourier_sin_test = Im(Xtwo_CorG_fourier)

colnames(Xtwo_AorT_fourier_cos_test) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_test) = paste0("fouriersinAorT2_", 1:25)

# colnames(Xtwo_CorG_fourier_cos_test) = paste0("fouriercosCorG2_", 1:25)
# colnames(Xtwo_CorG_fourier_sin_test) = paste0("fouriersinCorG2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_test, "data/Created/tiling_Xtwo_AorT_fourier_cos_test.rds")
saveRDS(Xtwo_AorT_fourier_sin_test, "data/Created/tiling_Xtwo_AorT_fourier_sin_test.rds")

# saveRDS(Xtwo_CorG_fourier_cos_test, "data/Created/tiling_Xtwo_CorG_fourier_cos_test.rds")
# saveRDS(Xtwo_CorG_fourier_sin_test, "data/Created/tiling_Xtwo_CorG_fourier_sin_test.rds")


Xtwo_AA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AA_test) = colnames(Xtwo_test)
Xtwo_AA_test[] = (Xtwo_test == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,1:25]
colnames(AA_fourier_test) = paste0("fourierAA_", 1:25)
AA_fourier_cos_test = Re(AA_fourier_test)
AA_fourier_sin_test = Im(AA_fourier_test)
colnames(AA_fourier_cos_test) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin_test) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos_test, "data/Created/tiling_AA_fourier_cos_test.rds")
saveRDS(AA_fourier_sin_test, "data/Created/tiling_AA_fourier_sin_test.rds")

Xtwo_AC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AC_test) = colnames(Xtwo_test)
Xtwo_AC_test[] = (Xtwo_test == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,1:25]
colnames(AC_fourier_test) = paste0("fourierAC_", 1:25)
AC_fourier_cos_test = Re(AC_fourier_test)
AC_fourier_sin_test = Im(AC_fourier_test)
colnames(AC_fourier_cos_test) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin_test) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos_test, "data/Created/tiling_AC_fourier_cos_test.rds")
saveRDS(AC_fourier_sin_test, "data/Created/tiling_AC_fourier_sin_test.rds")

Xtwo_AG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AG_test) = colnames(Xtwo_test)
Xtwo_AG_test[] = (Xtwo_test == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,1:25]
colnames(AG_fourier_test) = paste0("fourierAG_", 1:25)
AG_fourier_cos_test = Re(AG_fourier_test)
AG_fourier_sin_test = Im(AG_fourier_test)
colnames(AG_fourier_cos_test) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin_test) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos_test, "data/Created/tiling_AG_fourier_cos_test.rds")
saveRDS(AG_fourier_sin_test, "data/Created/tiling_AG_fourier_sin_test.rds")

Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
Xtwo_AT_test[] = (Xtwo_test == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,1:25]
colnames(AT_fourier_test) = paste0("fourierAT_", 1:25)
AT_fourier_cos_test = Re(AT_fourier_test)
AT_fourier_sin_test = Im(AT_fourier_test)
colnames(AT_fourier_cos_test) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin_test) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos_test, "data/Created/tiling_AT_fourier_cos_test.rds")
saveRDS(AT_fourier_sin_test, "data/Created/tiling_AT_fourier_sin_test.rds")

Xtwo_CA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CA_test) = colnames(Xtwo_test)
Xtwo_CA_test[] = (Xtwo_test == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,1:25]
colnames(CA_fourier_test) = paste0("fourierCA_", 1:25)
CA_fourier_cos_test = Re(CA_fourier_test)
CA_fourier_sin_test = Im(CA_fourier_test)
colnames(CA_fourier_cos_test) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin_test) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos_test, "data/Created/tiling_CA_fourier_cos_test.rds")
saveRDS(CA_fourier_sin_test, "data/Created/tiling_CA_fourier_sin_test.rds")

Xtwo_CC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CC_test) = colnames(Xtwo_test)
Xtwo_CC_test[] = (Xtwo_test == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,1:25]
colnames(CC_fourier_test) = paste0("fourierCC_", 1:25)
CC_fourier_cos_test = Re(CC_fourier_test)
CC_fourier_sin_test = Im(CC_fourier_test)
colnames(CC_fourier_cos_test) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin_test) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos_test, "data/Created/tiling_CC_fourier_cos_test.rds")
saveRDS(CC_fourier_sin_test, "data/Created/tiling_CC_fourier_sin_test.rds")

Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
Xtwo_CG_test[] = (Xtwo_test == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,1:25]
colnames(CG_fourier_test) = paste0("fourierCG_", 1:25)
CG_fourier_cos_test = Re(CG_fourier_test)
CG_fourier_sin_test = Im(CG_fourier_test)
colnames(CG_fourier_cos_test) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin_test) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos_test, "data/Created/tiling_CG_fourier_cos_test.rds")
saveRDS(CG_fourier_sin_test, "data/Created/tiling_CG_fourier_sin_test.rds")

Xtwo_CT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CT_test) = colnames(Xtwo_test)
Xtwo_CT_test[] = (Xtwo_test == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,1:25]
colnames(CT_fourier_test) = paste0("fourierCT_", 1:25)
CT_fourier_cos_test = Re(CT_fourier_test)
CT_fourier_sin_test = Im(CT_fourier_test)
colnames(CT_fourier_cos_test) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin_test) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos_test, "data/Created/tiling_CT_fourier_cos_test.rds")
saveRDS(CT_fourier_sin_test, "data/Created/tiling_CT_fourier_sin_test.rds")

Xtwo_GA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GA_test) = colnames(Xtwo_test)
Xtwo_GA_test[] = (Xtwo_test == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,1:25]
colnames(GA_fourier_test) = paste0("fourierGA_", 1:25)
GA_fourier_cos_test = Re(GA_fourier_test)
GA_fourier_sin_test = Im(GA_fourier_test)
colnames(GA_fourier_cos_test) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin_test) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos_test, "data/Created/tiling_GA_fourier_cos_test.rds")
saveRDS(GA_fourier_sin_test, "data/Created/tiling_GA_fourier_sin_test.rds")

Xtwo_GC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GC_test) = colnames(Xtwo_test)
Xtwo_GC_test[] = (Xtwo_test == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,1:25]
colnames(GC_fourier_test) = paste0("fourierGC_", 1:25)
GC_fourier_cos_test = Re(GC_fourier_test)
GC_fourier_sin_test = Im(GC_fourier_test)
colnames(GC_fourier_cos_test) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin_test) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos_test, "data/Created/tiling_GC_fourier_cos_test.rds")
saveRDS(GC_fourier_sin_test, "data/Created/tiling_GC_fourier_sin_test.rds")

Xtwo_GG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GG_test) = colnames(Xtwo_test)
Xtwo_GG_test[] = (Xtwo_test == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,1:25]
colnames(GG_fourier_test) = paste0("fourierGG_", 1:25)
GG_fourier_cos_test = Re(GG_fourier_test)
GG_fourier_sin_test = Im(GG_fourier_test)
colnames(GG_fourier_cos_test) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin_test) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos_test, "data/Created/tiling_GG_fourier_cos_test.rds")
saveRDS(GG_fourier_sin_test, "data/Created/tiling_GG_fourier_sin_test.rds")

Xtwo_GT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GT_test) = colnames(Xtwo_test)
Xtwo_GT_test[] = (Xtwo_test == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,1:25]
colnames(GT_fourier_test) = paste0("fourierGT_", 1:25)
GT_fourier_cos_test = Re(GT_fourier_test)
GT_fourier_sin_test = Im(GT_fourier_test)
colnames(GT_fourier_cos_test) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin_test) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos_test, "data/Created/tiling_GT_fourier_cos_test.rds")
saveRDS(GT_fourier_sin_test, "data/Created/tiling_GT_fourier_sin_test.rds")

Xtwo_TA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TA_test) = colnames(Xtwo_test)
Xtwo_TA_test[] = (Xtwo_test == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,1:25]
colnames(TA_fourier_test) = paste0("fourierTA_", 1:25)
TA_fourier_cos_test = Re(TA_fourier_test)
TA_fourier_sin_test = Im(TA_fourier_test)
colnames(TA_fourier_cos_test) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin_test) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos_test, "data/Created/tiling_TA_fourier_cos_test.rds")
saveRDS(TA_fourier_sin_test, "data/Created/tiling_TA_fourier_sin_test.rds")

Xtwo_TC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TC_test) = colnames(Xtwo_test)
Xtwo_TC_test[] = (Xtwo_test == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,1:25]
colnames(TC_fourier_test) = paste0("fourierTC_", 1:25)
TC_fourier_cos_test = Re(TC_fourier_test)
TC_fourier_sin_test = Im(TC_fourier_test)
colnames(TC_fourier_cos_test) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin_test) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos_test, "data/Created/tiling_TC_fourier_cos_test.rds")
saveRDS(TC_fourier_sin_test, "data/Created/tiling_TC_fourier_sin_test.rds")

Xtwo_TG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TG_test) = colnames(Xtwo_test)
Xtwo_TG_test[] = (Xtwo_test == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,1:25]
colnames(TG_fourier_test) = paste0("fourierTG_", 1:25)
TG_fourier_cos_test = Re(TG_fourier_test)
TG_fourier_sin_test = Im(TG_fourier_test)
colnames(TG_fourier_cos_test) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin_test) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos_test, "data/Created/tiling_TG_fourier_cos_test.rds")
saveRDS(TG_fourier_sin_test, "data/Created/tiling_TG_fourier_sin_test.rds")

Xtwo_TT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TT_test) = colnames(Xtwo_test)
Xtwo_TT_test[] = (Xtwo_test == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,1:25]
colnames(TT_fourier_test) = paste0("fourierTT_", 1:25)
TT_fourier_cos_test = Re(TT_fourier_test)
TT_fourier_sin_test = Im(TT_fourier_test)
colnames(TT_fourier_cos_test) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin_test) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos_test, "data/Created/tiling_TT_fourier_cos_test.rds")
saveRDS(TT_fourier_sin_test, "data/Created/tiling_TT_fourier_sin_test.rds")

fourier_di_cos_test = cbind(AA_fourier_cos_test, AC_fourier_cos_test, AG_fourier_cos_test, AT_fourier_cos_test, 
                            CA_fourier_cos_test, CC_fourier_cos_test, CG_fourier_cos_test, CT_fourier_cos_test,
                            GA_fourier_cos_test, GC_fourier_cos_test, GG_fourier_cos_test, GT_fourier_cos_test,
                            TA_fourier_cos_test, TC_fourier_cos_test, TG_fourier_cos_test, TT_fourier_cos_test)
fourier_di_sin_test = cbind(AA_fourier_sin_test, AC_fourier_sin_test, AG_fourier_sin_test, AT_fourier_sin_test, 
                            CA_fourier_sin_test, CC_fourier_sin_test, CG_fourier_sin_test, CT_fourier_sin_test,
                            GA_fourier_sin_test, GC_fourier_sin_test, GG_fourier_sin_test, GT_fourier_sin_test,
                            TA_fourier_sin_test, TC_fourier_sin_test, TG_fourier_sin_test, TT_fourier_sin_test)
saveRDS(fourier_di_cos_test, "data/Created/tiling_fourier_di_cos_test.rds")
saveRDS(fourier_di_sin_test, "data/Created/tiling_fourier_di_sin_test.rds")



fourier_AA_AC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AC_fourier_test[,"fourierAC_6"])
fourier_AA_AG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AA_AT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AA_CA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AA_CC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AA_CG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AA_CT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AA_GA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AA_GC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AA_GG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AA_GT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AA_TA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AA_TC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AA_TG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AA_TT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AA_6_mod_test = paste0("fourier_AA_", dinucleotides[-1], "_6_mod_test")
fourier_AA_6_mod_tiling_test = data.frame(map(fourier_AA_6_mod_test, get))
colnames(fourier_AA_6_mod_tiling_test) = fourier_AA_6_mod

fourier_AC_AG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AC_AT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AC_CA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AC_CC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AC_CG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AC_CT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AC_GA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AC_GC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AC_GG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AC_GT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AC_TA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AC_TC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AC_TG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AC_TT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AC_6_mod_test = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod_test")
fourier_AC_6_mod_tiling_test = data.frame(map(fourier_AC_6_mod_test, get))
colnames(fourier_AC_6_mod_tiling_test) = fourier_AC_6_mod

fourier_AG_AT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AG_CA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AG_CC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AG_CG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AG_CT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AG_GA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AG_GC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AG_GG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AG_GT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AG_TA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AG_TC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AG_TG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AG_TT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AG_6_mod_test = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod_test")
fourier_AG_6_mod_tiling_test = data.frame(map(fourier_AG_6_mod_test, get))
colnames(fourier_AG_6_mod_tiling_test) = fourier_AG_6_mod

fourier_AT_CA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AT_CC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AT_CG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AT_CT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AT_GA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AT_GC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AT_GG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AT_GT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AT_TA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AT_TC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AT_TG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AT_TT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AT_6_mod_test = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod_test")
fourier_AT_6_mod_tiling_test = data.frame(map(fourier_AT_6_mod_test, get))
colnames(fourier_AT_6_mod_tiling_test) = fourier_AT_6_mod

fourier_CA_CC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_CA_CG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CA_CT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CA_GA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CA_GC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CA_GG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CA_GT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CA_TA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CA_TC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CA_TG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CA_TT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CA_6_mod_test = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod_test")
fourier_CA_6_mod_tiling_test = data.frame(map(fourier_CA_6_mod_test, get))
colnames(fourier_CA_6_mod_tiling_test) = fourier_CA_6_mod

fourier_CC_CG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CC_CT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CC_GA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CC_GC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CC_GG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CC_GT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CC_TA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CC_TC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CC_TG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CC_TT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CC_6_mod_test = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod_test")
fourier_CC_6_mod_tiling_test = data.frame(map(fourier_CC_6_mod_test, get))
colnames(fourier_CC_6_mod_tiling_test) = fourier_CC_6_mod

fourier_CG_CT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CG_GA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CG_GC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CG_GG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CG_GT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CG_TA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CG_TC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CG_TG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CG_TT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CG_6_mod_test = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod_test")
fourier_CG_6_mod_tiling_test = data.frame(map(fourier_CG_6_mod_test, get))
colnames(fourier_CG_6_mod_tiling_test) = fourier_CG_6_mod

fourier_CT_GA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CT_GC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CT_GG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CT_GT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CT_TA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CT_TC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CT_TG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CT_TT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CT_6_mod_test = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod_test")
fourier_CT_6_mod_tiling_test = data.frame(map(fourier_CT_6_mod_test, get))
colnames(fourier_CT_6_mod_tiling_test) = fourier_CT_6_mod

fourier_GA_GC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_GA_GG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GA_GT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GA_TA_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GA_TC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GA_TG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GA_TT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GA_6_mod_test = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod_test")
fourier_GA_6_mod_tiling_test = data.frame(map(fourier_GA_6_mod_test, get))
colnames(fourier_GA_6_mod_tiling_test) = fourier_GA_6_mod

fourier_GC_GG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GC_GT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GC_TA_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GC_TC_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GC_TG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GC_TT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GC_6_mod_test = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod_test")
fourier_GC_6_mod_tiling_test = data.frame(map(fourier_GC_6_mod_test, get))
colnames(fourier_GC_6_mod_tiling_test) = fourier_GC_6_mod

fourier_GG_GT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GG_TA_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GG_TC_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GG_TG_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GG_TT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GG_6_mod_test = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod_test")
fourier_GG_6_mod_tiling_test = data.frame(map(fourier_GG_6_mod_test, get))
colnames(fourier_GG_6_mod_tiling_test) = fourier_GG_6_mod

fourier_GT_TA_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GT_TC_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GT_TG_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GT_TT_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GT_6_mod_test = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod_test")
fourier_GT_6_mod_tiling_test = data.frame(map(fourier_GT_6_mod_test, get))
colnames(fourier_GT_6_mod_tiling_test) = fourier_GT_6_mod

fourier_TA_TC_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_TA_TG_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TA_TT_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TA_6_mod_test = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod_test")
fourier_TA_6_mod_tiling_test = data.frame(map(fourier_TA_6_mod_test, get))
colnames(fourier_TA_6_mod_tiling_test) = fourier_TA_6_mod

fourier_TC_TG_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TC_TT_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TC_6_mod_test = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod_test")
fourier_TC_6_mod_tiling_test = data.frame(map(fourier_TC_6_mod_test, get))
colnames(fourier_TC_6_mod_tiling_test) = fourier_TC_6_mod

fourier_TG_TT_6_mod_test = Mod(TG_fourier_test[,"fourierTG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TG_6_mod_test = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod_test")
fourier_TG_6_mod_tiling_test = data.frame(map(fourier_TG_6_mod_test, get))
colnames(fourier_TG_6_mod_tiling_test) = fourier_TG_6_mod

fourier_6_mod_test = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_tiling_test"), get))
saveRDS(fourier_6_mod_test, "data/Created/tiling_fourier_6_mod_test.rds")











########################################################################
# Random
########################################################################

# Train data
ps1 <- paste0("X", 1:50, "mono")

Xone = dat_random %>% select(all_of(ps1))
Xone_AorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AorT) = colnames(Xone)
colnames(Xone_CorG) = colnames(Xone)
Xone_AorT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_random = t(apply(Xone_AorT, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_random = Re(Xone_AorT_fourier_random)
Xone_AorT_fourier_sin_random = Im(Xone_AorT_fourier_random)

colnames(Xone_AorT_fourier_cos_random) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_random) = paste0("fouriersinAorT", c(1:25, 50))

# colnames(Xone_CorG_fourier_cos_random) = paste0("fouriercosCorG", c(1:25, 50))
# colnames(Xone_CorG_fourier_sin_random) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_random, "data/Created/random_Xone_AorT_fourier_cos.rds")
saveRDS(Xone_AorT_fourier_sin_random, "data/Created/random_Xone_AorT_fourier_sin.rds")

# saveRDS(Xone_CorG_fourier_cos, "data/Created/tiling_Xone_CorG_fourier_cos.rds")
# saveRDS(Xone_CorG_fourier_sin, "data/Created/tiling_Xone_CorG_fourier_sin.rds")

ps2 <- paste0("X", 1:49, "di")

Xtwo = dat_random %>% select(all_of(ps2))
Xtwo_AorT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CorG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AorT) = colnames(Xtwo)
colnames(Xtwo_CorG) = colnames(Xtwo)
Xtwo_AorT[] = ((Xtwo == "AA") | (Xtwo == "TT") | (Xtwo == "AT") | (Xtwo == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CorG[] = ((Xtwo == "CC") | (Xtwo == "GG") | (Xtwo == "CG") | (Xtwo == "GC")) %>%
  as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_random = t(apply(Xtwo_AorT, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_random = Re(Xtwo_AorT_fourier_random)
Xtwo_AorT_fourier_sin_random = Im(Xtwo_AorT_fourier_random)

colnames(Xtwo_AorT_fourier_cos_random) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_random) = paste0("fouriersinAorT2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_random, "data/Created/random_Xtwo_AorT_fourier_cos.rds")
saveRDS(Xtwo_AorT_fourier_sin_random, "data/Created/random_Xtwo_AorT_fourier_sin.rds")



Xtwo_AA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AA) = colnames(Xtwo)
Xtwo_AA[] = (Xtwo == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier = t(apply(Xtwo_AA, 1, fft))[,1:25]
colnames(AA_fourier) = paste0("fourierAA_", 1:25)
AA_fourier_cos = Re(AA_fourier)
AA_fourier_sin = Im(AA_fourier)
colnames(AA_fourier_cos) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos, "data/Created/random_AA_fourier_cos.rds")
saveRDS(AA_fourier_sin, "data/Created/random_AA_fourier_sin.rds")

Xtwo_AC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AC) = colnames(Xtwo)
Xtwo_AC[] = (Xtwo == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier = t(apply(Xtwo_AC, 1, fft))[,1:25]
colnames(AC_fourier) = paste0("fourierAC_", 1:25)
AC_fourier_cos = Re(AC_fourier)
AC_fourier_sin = Im(AC_fourier)
colnames(AC_fourier_cos) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos, "data/Created/random_AC_fourier_cos.rds")
saveRDS(AC_fourier_sin, "data/Created/random_AC_fourier_sin.rds")

Xtwo_AG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AG) = colnames(Xtwo)
Xtwo_AG[] = (Xtwo == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier = t(apply(Xtwo_AG, 1, fft))[,1:25]
colnames(AG_fourier) = paste0("fourierAG_", 1:25)
AG_fourier_cos = Re(AG_fourier)
AG_fourier_sin = Im(AG_fourier)
colnames(AG_fourier_cos) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos, "data/Created/random_AG_fourier_cos.rds")
saveRDS(AG_fourier_sin, "data/Created/random_AG_fourier_sin.rds")

Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AT) = colnames(Xtwo)
Xtwo_AT[] = (Xtwo == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier = t(apply(Xtwo_AT, 1, fft))[,1:25]
colnames(AT_fourier) = paste0("fourierAT_", 1:25)
AT_fourier_cos = Re(AT_fourier)
AT_fourier_sin = Im(AT_fourier)
colnames(AT_fourier_cos) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos, "data/Created/random_AT_fourier_cos.rds")
saveRDS(AT_fourier_sin, "data/Created/random_AT_fourier_sin.rds")

Xtwo_CA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CA) = colnames(Xtwo)
Xtwo_CA[] = (Xtwo == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier = t(apply(Xtwo_CA, 1, fft))[,1:25]
colnames(CA_fourier) = paste0("fourierCA_", 1:25)
CA_fourier_cos = Re(CA_fourier)
CA_fourier_sin = Im(CA_fourier)
colnames(CA_fourier_cos) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos, "data/Created/random_CA_fourier_cos.rds")
saveRDS(CA_fourier_sin, "data/Created/random_CA_fourier_sin.rds")

Xtwo_CC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CC) = colnames(Xtwo)
Xtwo_CC[] = (Xtwo == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier = t(apply(Xtwo_CC, 1, fft))[,1:25]
colnames(CC_fourier) = paste0("fourierCC_", 1:25)
CC_fourier_cos = Re(CC_fourier)
CC_fourier_sin = Im(CC_fourier)
colnames(CC_fourier_cos) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos, "data/Created/random_CC_fourier_cos.rds")
saveRDS(CC_fourier_sin, "data/Created/random_CC_fourier_sin.rds")

Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CG) = colnames(Xtwo)
Xtwo_CG[] = (Xtwo == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier = t(apply(Xtwo_CG, 1, fft))[,1:25]
colnames(CG_fourier) = paste0("fourierCG_", 1:25)
CG_fourier_cos = Re(CG_fourier)
CG_fourier_sin = Im(CG_fourier)
colnames(CG_fourier_cos) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos, "data/Created/random_CG_fourier_cos.rds")
saveRDS(CG_fourier_sin, "data/Created/random_CG_fourier_sin.rds")

Xtwo_CT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CT) = colnames(Xtwo)
Xtwo_CT[] = (Xtwo == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier = t(apply(Xtwo_CT, 1, fft))[,1:25]
colnames(CT_fourier) = paste0("fourierCT_", 1:25)
CT_fourier_cos = Re(CT_fourier)
CT_fourier_sin = Im(CT_fourier)
colnames(CT_fourier_cos) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos, "data/Created/random_CT_fourier_cos.rds")
saveRDS(CT_fourier_sin, "data/Created/random_CT_fourier_sin.rds")

Xtwo_GA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GA) = colnames(Xtwo)
Xtwo_GA[] = (Xtwo == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier = t(apply(Xtwo_GA, 1, fft))[,1:25]
colnames(GA_fourier) = paste0("fourierGA_", 1:25)
GA_fourier_cos = Re(GA_fourier)
GA_fourier_sin = Im(GA_fourier)
colnames(GA_fourier_cos) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos, "data/Created/random_GA_fourier_cos.rds")
saveRDS(GA_fourier_sin, "data/Created/random_GA_fourier_sin.rds")

Xtwo_GC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GC) = colnames(Xtwo)
Xtwo_GC[] = (Xtwo == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier = t(apply(Xtwo_GC, 1, fft))[,1:25]
colnames(GC_fourier) = paste0("fourierGC_", 1:25)
GC_fourier_cos = Re(GC_fourier)
GC_fourier_sin = Im(GC_fourier)
colnames(GC_fourier_cos) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos, "data/Created/random_GC_fourier_cos.rds")
saveRDS(GC_fourier_sin, "data/Created/random_GC_fourier_sin.rds")

Xtwo_GG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GG) = colnames(Xtwo)
Xtwo_GG[] = (Xtwo == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier = t(apply(Xtwo_GG, 1, fft))[,1:25]
colnames(GG_fourier) = paste0("fourierGG_", 1:25)
GG_fourier_cos = Re(GG_fourier)
GG_fourier_sin = Im(GG_fourier)
colnames(GG_fourier_cos) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos, "data/Created/random_GG_fourier_cos.rds")
saveRDS(GG_fourier_sin, "data/Created/random_GG_fourier_sin.rds")

Xtwo_GT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GT) = colnames(Xtwo)
Xtwo_GT[] = (Xtwo == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier = t(apply(Xtwo_GT, 1, fft))[,1:25]
colnames(GT_fourier) = paste0("fourierGT_", 1:25)
GT_fourier_cos = Re(GT_fourier)
GT_fourier_sin = Im(GT_fourier)
colnames(GT_fourier_cos) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos, "data/Created/random_GT_fourier_cos.rds")
saveRDS(GT_fourier_sin, "data/Created/random_GT_fourier_sin.rds")

Xtwo_TA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TA) = colnames(Xtwo)
Xtwo_TA[] = (Xtwo == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier = t(apply(Xtwo_TA, 1, fft))[,1:25]
colnames(TA_fourier) = paste0("fourierTA_", 1:25)
TA_fourier_cos = Re(TA_fourier)
TA_fourier_sin = Im(TA_fourier)
colnames(TA_fourier_cos) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos, "data/Created/random_TA_fourier_cos.rds")
saveRDS(TA_fourier_sin, "data/Created/random_TA_fourier_sin.rds")

Xtwo_TC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TC) = colnames(Xtwo)
Xtwo_TC[] = (Xtwo == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier = t(apply(Xtwo_TC, 1, fft))[,1:25]
colnames(TC_fourier) = paste0("fourierTC_", 1:25)
TC_fourier_cos = Re(TC_fourier)
TC_fourier_sin = Im(TC_fourier)
colnames(TC_fourier_cos) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos, "data/Created/random_TC_fourier_cos.rds")
saveRDS(TC_fourier_sin, "data/Created/random_TC_fourier_sin.rds")

Xtwo_TG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TG) = colnames(Xtwo)
Xtwo_TG[] = (Xtwo == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier = t(apply(Xtwo_TG, 1, fft))[,1:25]
colnames(TG_fourier) = paste0("fourierTG_", 1:25)
TG_fourier_cos = Re(TG_fourier)
TG_fourier_sin = Im(TG_fourier)
colnames(TG_fourier_cos) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos, "data/Created/random_TG_fourier_cos.rds")
saveRDS(TG_fourier_sin, "data/Created/random_TG_fourier_sin.rds")

Xtwo_TT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TT) = colnames(Xtwo)
Xtwo_TT[] = (Xtwo == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier = t(apply(Xtwo_TT, 1, fft))[,1:25]
colnames(TT_fourier) = paste0("fourierTT_", 1:25)
TT_fourier_cos = Re(TT_fourier)
TT_fourier_sin = Im(TT_fourier)
colnames(TT_fourier_cos) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos, "data/Created/random_TT_fourier_cos.rds")
saveRDS(TT_fourier_sin, "data/Created/random_TT_fourier_sin.rds")

fourier_di_cos = cbind(AA_fourier_cos, AC_fourier_cos, AG_fourier_cos, AT_fourier_cos, 
                       CA_fourier_cos, CC_fourier_cos, CG_fourier_cos, CT_fourier_cos,
                       GA_fourier_cos, GC_fourier_cos, GG_fourier_cos, GT_fourier_cos,
                       TA_fourier_cos, TC_fourier_cos, TG_fourier_cos, TT_fourier_cos)
fourier_di_sin = cbind(AA_fourier_sin, AC_fourier_sin, AG_fourier_sin, AT_fourier_sin, 
                       CA_fourier_sin, CC_fourier_sin, CG_fourier_sin, CT_fourier_sin,
                       GA_fourier_sin, GC_fourier_sin, GG_fourier_sin, GT_fourier_sin,
                       TA_fourier_sin, TC_fourier_sin, TG_fourier_sin, TT_fourier_sin)
saveRDS(fourier_di_cos, "data/Created/random_fourier_di_cos.rds")
saveRDS(fourier_di_sin, "data/Created/random_fourier_di_sin.rds")

fourier_AA_AC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AC_fourier[,"fourierAC_6"])
fourier_AA_AG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AG_fourier[,"fourierAG_6"])
fourier_AA_AT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AT_fourier[,"fourierAT_6"])
fourier_AA_CA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CA_fourier[,"fourierCA_6"])
fourier_AA_CC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CC_fourier[,"fourierCC_6"])
fourier_AA_CG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CG_fourier[,"fourierCG_6"])
fourier_AA_CT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CT_fourier[,"fourierCT_6"])
fourier_AA_GA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GA_fourier[,"fourierGA_6"])
fourier_AA_GC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GC_fourier[,"fourierGC_6"])
fourier_AA_GG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GG_fourier[,"fourierGG_6"])
fourier_AA_GT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GT_fourier[,"fourierGT_6"])
fourier_AA_TA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TA_fourier[,"fourierTA_6"])
fourier_AA_TC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TC_fourier[,"fourierTC_6"])
fourier_AA_TG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TG_fourier[,"fourierTG_6"])
fourier_AA_TT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TT_fourier[,"fourierTT_6"])

fourier_AA_6_mod = paste0("fourier_AA_", dinucleotides[-1], "_6_mod")
fourier_AA_6_mod_random = data.frame(map(fourier_AA_6_mod, get))
colnames(fourier_AA_6_mod_random) = fourier_AA_6_mod

fourier_AC_AG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AG_fourier[,"fourierAG_6"])
fourier_AC_AT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AT_fourier[,"fourierAT_6"])
fourier_AC_CA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CA_fourier[,"fourierCA_6"])
fourier_AC_CC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CC_fourier[,"fourierCC_6"])
fourier_AC_CG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CG_fourier[,"fourierCG_6"])
fourier_AC_CT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CT_fourier[,"fourierCT_6"])
fourier_AC_GA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GA_fourier[,"fourierGA_6"])
fourier_AC_GC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GC_fourier[,"fourierGC_6"])
fourier_AC_GG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GG_fourier[,"fourierGG_6"])
fourier_AC_GT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GT_fourier[,"fourierGT_6"])
fourier_AC_TA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TA_fourier[,"fourierTA_6"])
fourier_AC_TC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TC_fourier[,"fourierTC_6"])
fourier_AC_TG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TG_fourier[,"fourierTG_6"])
fourier_AC_TT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TT_fourier[,"fourierTT_6"])

fourier_AC_6_mod = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod")
fourier_AC_6_mod_random = data.frame(map(fourier_AC_6_mod, get))
colnames(fourier_AC_6_mod_random) = fourier_AC_6_mod

fourier_AG_AT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - AT_fourier[,"fourierAT_6"])
fourier_AG_CA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CA_fourier[,"fourierCA_6"])
fourier_AG_CC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CC_fourier[,"fourierCC_6"])
fourier_AG_CG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CG_fourier[,"fourierCG_6"])
fourier_AG_CT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CT_fourier[,"fourierCT_6"])
fourier_AG_GA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GA_fourier[,"fourierGA_6"])
fourier_AG_GC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GC_fourier[,"fourierGC_6"])
fourier_AG_GG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GG_fourier[,"fourierGG_6"])
fourier_AG_GT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GT_fourier[,"fourierGT_6"])
fourier_AG_TA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TA_fourier[,"fourierTA_6"])
fourier_AG_TC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TC_fourier[,"fourierTC_6"])
fourier_AG_TG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TG_fourier[,"fourierTG_6"])
fourier_AG_TT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TT_fourier[,"fourierTT_6"])

fourier_AG_6_mod = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod")
fourier_AG_6_mod_random = data.frame(map(fourier_AG_6_mod, get))
colnames(fourier_AG_6_mod_random) = fourier_AG_6_mod

fourier_AT_CA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CA_fourier[,"fourierCA_6"])
fourier_AT_CC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CC_fourier[,"fourierCC_6"])
fourier_AT_CG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CG_fourier[,"fourierCG_6"])
fourier_AT_CT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CT_fourier[,"fourierCT_6"])
fourier_AT_GA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GA_fourier[,"fourierGA_6"])
fourier_AT_GC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GC_fourier[,"fourierGC_6"])
fourier_AT_GG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GG_fourier[,"fourierGG_6"])
fourier_AT_GT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GT_fourier[,"fourierGT_6"])
fourier_AT_TA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TA_fourier[,"fourierTA_6"])
fourier_AT_TC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TC_fourier[,"fourierTC_6"])
fourier_AT_TG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TG_fourier[,"fourierTG_6"])
fourier_AT_TT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TT_fourier[,"fourierTT_6"])

fourier_AT_6_mod = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod")
fourier_AT_6_mod_random = data.frame(map(fourier_AT_6_mod, get))
colnames(fourier_AT_6_mod_random) = fourier_AT_6_mod

fourier_CA_CC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CC_fourier[,"fourierCC_6"])
fourier_CA_CG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CG_fourier[,"fourierCG_6"])
fourier_CA_CT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CT_fourier[,"fourierCT_6"])
fourier_CA_GA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GA_fourier[,"fourierGA_6"])
fourier_CA_GC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GC_fourier[,"fourierGC_6"])
fourier_CA_GG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GG_fourier[,"fourierGG_6"])
fourier_CA_GT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GT_fourier[,"fourierGT_6"])
fourier_CA_TA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TA_fourier[,"fourierTA_6"])
fourier_CA_TC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TC_fourier[,"fourierTC_6"])
fourier_CA_TG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TG_fourier[,"fourierTG_6"])
fourier_CA_TT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TT_fourier[,"fourierTT_6"])

fourier_CA_6_mod = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod")
fourier_CA_6_mod_random = data.frame(map(fourier_CA_6_mod, get))
colnames(fourier_CA_6_mod_random) = fourier_CA_6_mod

fourier_CC_CG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CG_fourier[,"fourierCG_6"])
fourier_CC_CT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CT_fourier[,"fourierCT_6"])
fourier_CC_GA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GA_fourier[,"fourierGA_6"])
fourier_CC_GC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GC_fourier[,"fourierGC_6"])
fourier_CC_GG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GG_fourier[,"fourierGG_6"])
fourier_CC_GT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GT_fourier[,"fourierGT_6"])
fourier_CC_TA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TA_fourier[,"fourierTA_6"])
fourier_CC_TC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TC_fourier[,"fourierTC_6"])
fourier_CC_TG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TG_fourier[,"fourierTG_6"])
fourier_CC_TT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TT_fourier[,"fourierTT_6"])

fourier_CC_6_mod = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod")
fourier_CC_6_mod_random = data.frame(map(fourier_CC_6_mod, get))
colnames(fourier_CC_6_mod_random) = fourier_CC_6_mod

fourier_CG_CT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - CT_fourier[,"fourierCT_6"])
fourier_CG_GA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GA_fourier[,"fourierGA_6"])
fourier_CG_GC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GC_fourier[,"fourierGC_6"])
fourier_CG_GG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GG_fourier[,"fourierGG_6"])
fourier_CG_GT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GT_fourier[,"fourierGT_6"])
fourier_CG_TA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TA_fourier[,"fourierTA_6"])
fourier_CG_TC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TC_fourier[,"fourierTC_6"])
fourier_CG_TG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TG_fourier[,"fourierTG_6"])
fourier_CG_TT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TT_fourier[,"fourierTT_6"])

fourier_CG_6_mod = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod")
fourier_CG_6_mod_random = data.frame(map(fourier_CG_6_mod, get))
colnames(fourier_CG_6_mod_random) = fourier_CG_6_mod

fourier_CT_GA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GA_fourier[,"fourierGA_6"])
fourier_CT_GC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GC_fourier[,"fourierGC_6"])
fourier_CT_GG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GG_fourier[,"fourierGG_6"])
fourier_CT_GT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GT_fourier[,"fourierGT_6"])
fourier_CT_TA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TA_fourier[,"fourierTA_6"])
fourier_CT_TC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TC_fourier[,"fourierTC_6"])
fourier_CT_TG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TG_fourier[,"fourierTG_6"])
fourier_CT_TT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TT_fourier[,"fourierTT_6"])

fourier_CT_6_mod = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod")
fourier_CT_6_mod_random = data.frame(map(fourier_CT_6_mod, get))
colnames(fourier_CT_6_mod_random) = fourier_CT_6_mod

fourier_GA_GC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GC_fourier[,"fourierGC_6"])
fourier_GA_GG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GG_fourier[,"fourierGG_6"])
fourier_GA_GT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GT_fourier[,"fourierGT_6"])
fourier_GA_TA_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TA_fourier[,"fourierTA_6"])
fourier_GA_TC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TC_fourier[,"fourierTC_6"])
fourier_GA_TG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TG_fourier[,"fourierTG_6"])
fourier_GA_TT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TT_fourier[,"fourierTT_6"])

fourier_GA_6_mod = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod")
fourier_GA_6_mod_random = data.frame(map(fourier_GA_6_mod, get))
colnames(fourier_GA_6_mod_random) = fourier_GA_6_mod

fourier_GC_GG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GG_fourier[,"fourierGG_6"])
fourier_GC_GT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GT_fourier[,"fourierGT_6"])
fourier_GC_TA_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TA_fourier[,"fourierTA_6"])
fourier_GC_TC_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TC_fourier[,"fourierTC_6"])
fourier_GC_TG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TG_fourier[,"fourierTG_6"])
fourier_GC_TT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TT_fourier[,"fourierTT_6"])

fourier_GC_6_mod = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod")
fourier_GC_6_mod_random = data.frame(map(fourier_GC_6_mod, get))
colnames(fourier_GC_6_mod_random) = fourier_GC_6_mod

fourier_GG_GT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - GT_fourier[,"fourierGT_6"])
fourier_GG_TA_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TA_fourier[,"fourierTA_6"])
fourier_GG_TC_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TC_fourier[,"fourierTC_6"])
fourier_GG_TG_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TG_fourier[,"fourierTG_6"])
fourier_GG_TT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TT_fourier[,"fourierTT_6"])

fourier_GG_6_mod = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod")
fourier_GG_6_mod_random = data.frame(map(fourier_GG_6_mod, get))
colnames(fourier_GG_6_mod_random) = fourier_GG_6_mod

fourier_GT_TA_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TA_fourier[,"fourierTA_6"])
fourier_GT_TC_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TC_fourier[,"fourierTC_6"])
fourier_GT_TG_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TG_fourier[,"fourierTG_6"])
fourier_GT_TT_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TT_fourier[,"fourierTT_6"])

fourier_GT_6_mod = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod")
fourier_GT_6_mod_random = data.frame(map(fourier_GT_6_mod, get))
colnames(fourier_GT_6_mod_random) = fourier_GT_6_mod

fourier_TA_TC_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TC_fourier[,"fourierTC_6"])
fourier_TA_TG_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TG_fourier[,"fourierTG_6"])
fourier_TA_TT_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TT_fourier[,"fourierTT_6"])

fourier_TA_6_mod = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod")
fourier_TA_6_mod_random = data.frame(map(fourier_TA_6_mod, get))
colnames(fourier_TA_6_mod_random) = fourier_TA_6_mod

fourier_TC_TG_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TG_fourier[,"fourierTG_6"])
fourier_TC_TT_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TT_fourier[,"fourierTT_6"])

fourier_TC_6_mod = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod")
fourier_TC_6_mod_random = data.frame(map(fourier_TC_6_mod, get))
colnames(fourier_TC_6_mod_random) = fourier_TC_6_mod

fourier_TG_TT_6_mod = Mod(TG_fourier[,"fourierTG_6"] - TT_fourier[,"fourierTT_6"])

fourier_TG_6_mod = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod")
fourier_TG_6_mod_random = data.frame(map(fourier_TG_6_mod, get))
colnames(fourier_TG_6_mod_random) = fourier_TG_6_mod

fourier_6_mod = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_random"), get))
saveRDS(fourier_6_mod, "data/Created/random_fourier_6_mod.rds")




# Test data:
Xone_test = dat_random_test %>% select(all_of(ps1))
Xone_AorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AorT_test) = colnames(Xone_test)
colnames(Xone_CorG_test) = colnames(Xone_test)
Xone_AorT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_test = Re(Xone_AorT_fourier_test)
Xone_AorT_fourier_sin_test = Im(Xone_AorT_fourier_test)

Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,c(1:25, 50)]
Xone_CorG_fourier_cos_test = Re(Xone_CorG_fourier_test)
Xone_CorG_fourier_sin_test = Im(Xone_CorG_fourier_test)

colnames(Xone_AorT_fourier_cos_test) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_test) = paste0("fouriersinAorT", c(1:25, 50))

colnames(Xone_CorG_fourier_cos_test) = paste0("fouriercosCorG", c(1:25, 50))
colnames(Xone_CorG_fourier_sin_test) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_test, "data/Created/random_Xone_AorT_fourier_cos_test.rds")
saveRDS(Xone_AorT_fourier_sin_test, "data/Created/random_Xone_AorT_fourier_sin_test.rds")

saveRDS(Xone_CorG_fourier_cos_test, "data/Created/random_Xone_CorG_fourier_cos_test.rds")
saveRDS(Xone_CorG_fourier_sin_test, "data/Created/random_Xone_CorG_fourier_sin_test.rds")


Xtwo_test = dat_random_test %>% select(all_of(ps2))
Xtwo_AorT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
# Xtwo_CorG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AorT_test) = colnames(Xtwo_test)
# colnames(Xtwo_CorG_test) = colnames(Xtwo_test)
Xtwo_AorT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT") | (Xtwo_test == "AT") | (Xtwo_test == "TA")) %>% 
  as.matrix() %>% as.numeric()
# Xtwo_CorG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG") | (Xtwo_test == "CG") | (Xtwo_test == "GC")) %>%
#   as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_test = t(apply(Xtwo_AorT_test, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_test = Re(Xtwo_AorT_fourier_test)
Xtwo_AorT_fourier_sin_test = Im(Xtwo_AorT_fourier_test)

# Xtwo_CorG_fourier = t(apply(Xtwo_CorG_test, 1, fft))[,c(1:25, 50)]
# Xtwo_CorG_fourier_cos_test = Re(Xtwo_CorG_fourier)
# Xtwo_CorG_fourier_sin_test = Im(Xtwo_CorG_fourier)

colnames(Xtwo_AorT_fourier_cos_test) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_test) = paste0("fouriersinAorT2_", 1:25)

# colnames(Xtwo_CorG_fourier_cos_test) = paste0("fouriercosCorG2_", 1:25)
# colnames(Xtwo_CorG_fourier_sin_test) = paste0("fouriersinCorG2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_test, "data/Created/random_Xtwo_AorT_fourier_cos_test.rds")
saveRDS(Xtwo_AorT_fourier_sin_test, "data/Created/random_Xtwo_AorT_fourier_sin_test.rds")

# saveRDS(Xtwo_CorG_fourier_cos_test, "data/Created/random_Xtwo_CorG_fourier_cos_test.rds")
# saveRDS(Xtwo_CorG_fourier_sin_test, "data/Created/random_Xtwo_CorG_fourier_sin_test.rds")


Xtwo_AA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AA_test) = colnames(Xtwo_test)
Xtwo_AA_test[] = (Xtwo_test == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,1:25]
colnames(AA_fourier_test) = paste0("fourierAA_", 1:25)
AA_fourier_cos_test = Re(AA_fourier_test)
AA_fourier_sin_test = Im(AA_fourier_test)
colnames(AA_fourier_cos_test) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin_test) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos_test, "data/Created/random_AA_fourier_cos_test.rds")
saveRDS(AA_fourier_sin_test, "data/Created/random_AA_fourier_sin_test.rds")

Xtwo_AC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AC_test) = colnames(Xtwo_test)
Xtwo_AC_test[] = (Xtwo_test == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,1:25]
colnames(AC_fourier_test) = paste0("fourierAC_", 1:25)
AC_fourier_cos_test = Re(AC_fourier_test)
AC_fourier_sin_test = Im(AC_fourier_test)
colnames(AC_fourier_cos_test) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin_test) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos_test, "data/Created/random_AC_fourier_cos_test.rds")
saveRDS(AC_fourier_sin_test, "data/Created/random_AC_fourier_sin_test.rds")

Xtwo_AG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AG_test) = colnames(Xtwo_test)
Xtwo_AG_test[] = (Xtwo_test == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,1:25]
colnames(AG_fourier_test) = paste0("fourierAG_", 1:25)
AG_fourier_cos_test = Re(AG_fourier_test)
AG_fourier_sin_test = Im(AG_fourier_test)
colnames(AG_fourier_cos_test) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin_test) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos_test, "data/Created/random_AG_fourier_cos_test.rds")
saveRDS(AG_fourier_sin_test, "data/Created/random_AG_fourier_sin_test.rds")

Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
Xtwo_AT_test[] = (Xtwo_test == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,1:25]
colnames(AT_fourier_test) = paste0("fourierAT_", 1:25)
AT_fourier_cos_test = Re(AT_fourier_test)
AT_fourier_sin_test = Im(AT_fourier_test)
colnames(AT_fourier_cos_test) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin_test) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos_test, "data/Created/random_AT_fourier_cos_test.rds")
saveRDS(AT_fourier_sin_test, "data/Created/random_AT_fourier_sin_test.rds")

Xtwo_CA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CA_test) = colnames(Xtwo_test)
Xtwo_CA_test[] = (Xtwo_test == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,1:25]
colnames(CA_fourier_test) = paste0("fourierCA_", 1:25)
CA_fourier_cos_test = Re(CA_fourier_test)
CA_fourier_sin_test = Im(CA_fourier_test)
colnames(CA_fourier_cos_test) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin_test) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos_test, "data/Created/random_CA_fourier_cos_test.rds")
saveRDS(CA_fourier_sin_test, "data/Created/random_CA_fourier_sin_test.rds")

Xtwo_CC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CC_test) = colnames(Xtwo_test)
Xtwo_CC_test[] = (Xtwo_test == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,1:25]
colnames(CC_fourier_test) = paste0("fourierCC_", 1:25)
CC_fourier_cos_test = Re(CC_fourier_test)
CC_fourier_sin_test = Im(CC_fourier_test)
colnames(CC_fourier_cos_test) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin_test) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos_test, "data/Created/random_CC_fourier_cos_test.rds")
saveRDS(CC_fourier_sin_test, "data/Created/random_CC_fourier_sin_test.rds")

Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
Xtwo_CG_test[] = (Xtwo_test == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,1:25]
colnames(CG_fourier_test) = paste0("fourierCG_", 1:25)
CG_fourier_cos_test = Re(CG_fourier_test)
CG_fourier_sin_test = Im(CG_fourier_test)
colnames(CG_fourier_cos_test) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin_test) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos_test, "data/Created/random_CG_fourier_cos_test.rds")
saveRDS(CG_fourier_sin_test, "data/Created/random_CG_fourier_sin_test.rds")

Xtwo_CT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CT_test) = colnames(Xtwo_test)
Xtwo_CT_test[] = (Xtwo_test == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,1:25]
colnames(CT_fourier_test) = paste0("fourierCT_", 1:25)
CT_fourier_cos_test = Re(CT_fourier_test)
CT_fourier_sin_test = Im(CT_fourier_test)
colnames(CT_fourier_cos_test) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin_test) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos_test, "data/Created/random_CT_fourier_cos_test.rds")
saveRDS(CT_fourier_sin_test, "data/Created/random_CT_fourier_sin_test.rds")

Xtwo_GA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GA_test) = colnames(Xtwo_test)
Xtwo_GA_test[] = (Xtwo_test == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,1:25]
colnames(GA_fourier_test) = paste0("fourierGA_", 1:25)
GA_fourier_cos_test = Re(GA_fourier_test)
GA_fourier_sin_test = Im(GA_fourier_test)
colnames(GA_fourier_cos_test) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin_test) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos_test, "data/Created/random_GA_fourier_cos_test.rds")
saveRDS(GA_fourier_sin_test, "data/Created/random_GA_fourier_sin_test.rds")

Xtwo_GC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GC_test) = colnames(Xtwo_test)
Xtwo_GC_test[] = (Xtwo_test == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,1:25]
colnames(GC_fourier_test) = paste0("fourierGC_", 1:25)
GC_fourier_cos_test = Re(GC_fourier_test)
GC_fourier_sin_test = Im(GC_fourier_test)
colnames(GC_fourier_cos_test) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin_test) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos_test, "data/Created/random_GC_fourier_cos_test.rds")
saveRDS(GC_fourier_sin_test, "data/Created/random_GC_fourier_sin_test.rds")

Xtwo_GG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GG_test) = colnames(Xtwo_test)
Xtwo_GG_test[] = (Xtwo_test == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,1:25]
colnames(GG_fourier_test) = paste0("fourierGG_", 1:25)
GG_fourier_cos_test = Re(GG_fourier_test)
GG_fourier_sin_test = Im(GG_fourier_test)
colnames(GG_fourier_cos_test) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin_test) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos_test, "data/Created/random_GG_fourier_cos_test.rds")
saveRDS(GG_fourier_sin_test, "data/Created/random_GG_fourier_sin_test.rds")

Xtwo_GT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GT_test) = colnames(Xtwo_test)
Xtwo_GT_test[] = (Xtwo_test == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,1:25]
colnames(GT_fourier_test) = paste0("fourierGT_", 1:25)
GT_fourier_cos_test = Re(GT_fourier_test)
GT_fourier_sin_test = Im(GT_fourier_test)
colnames(GT_fourier_cos_test) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin_test) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos_test, "data/Created/random_GT_fourier_cos_test.rds")
saveRDS(GT_fourier_sin_test, "data/Created/random_GT_fourier_sin_test.rds")

Xtwo_TA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TA_test) = colnames(Xtwo_test)
Xtwo_TA_test[] = (Xtwo_test == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,1:25]
colnames(TA_fourier_test) = paste0("fourierTA_", 1:25)
TA_fourier_cos_test = Re(TA_fourier_test)
TA_fourier_sin_test = Im(TA_fourier_test)
colnames(TA_fourier_cos_test) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin_test) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos_test, "data/Created/random_TA_fourier_cos_test.rds")
saveRDS(TA_fourier_sin_test, "data/Created/random_TA_fourier_sin_test.rds")

Xtwo_TC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TC_test) = colnames(Xtwo_test)
Xtwo_TC_test[] = (Xtwo_test == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,1:25]
colnames(TC_fourier_test) = paste0("fourierTC_", 1:25)
TC_fourier_cos_test = Re(TC_fourier_test)
TC_fourier_sin_test = Im(TC_fourier_test)
colnames(TC_fourier_cos_test) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin_test) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos_test, "data/Created/random_TC_fourier_cos_test.rds")
saveRDS(TC_fourier_sin_test, "data/Created/random_TC_fourier_sin_test.rds")

Xtwo_TG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TG_test) = colnames(Xtwo_test)
Xtwo_TG_test[] = (Xtwo_test == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,1:25]
colnames(TG_fourier_test) = paste0("fourierTG_", 1:25)
TG_fourier_cos_test = Re(TG_fourier_test)
TG_fourier_sin_test = Im(TG_fourier_test)
colnames(TG_fourier_cos_test) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin_test) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos_test, "data/Created/random_TG_fourier_cos_test.rds")
saveRDS(TG_fourier_sin_test, "data/Created/random_TG_fourier_sin_test.rds")

Xtwo_TT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TT_test) = colnames(Xtwo_test)
Xtwo_TT_test[] = (Xtwo_test == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,1:25]
colnames(TT_fourier_test) = paste0("fourierTT_", 1:25)
TT_fourier_cos_test = Re(TT_fourier_test)
TT_fourier_sin_test = Im(TT_fourier_test)
colnames(TT_fourier_cos_test) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin_test) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos_test, "data/Created/random_TT_fourier_cos_test.rds")
saveRDS(TT_fourier_sin_test, "data/Created/random_TT_fourier_sin_test.rds")

fourier_di_cos_test = cbind(AA_fourier_cos_test, AC_fourier_cos_test, AG_fourier_cos_test, AT_fourier_cos_test, 
                            CA_fourier_cos_test, CC_fourier_cos_test, CG_fourier_cos_test, CT_fourier_cos_test,
                            GA_fourier_cos_test, GC_fourier_cos_test, GG_fourier_cos_test, GT_fourier_cos_test,
                            TA_fourier_cos_test, TC_fourier_cos_test, TG_fourier_cos_test, TT_fourier_cos_test)
fourier_di_sin_test = cbind(AA_fourier_sin_test, AC_fourier_sin_test, AG_fourier_sin_test, AT_fourier_sin_test, 
                            CA_fourier_sin_test, CC_fourier_sin_test, CG_fourier_sin_test, CT_fourier_sin_test,
                            GA_fourier_sin_test, GC_fourier_sin_test, GG_fourier_sin_test, GT_fourier_sin_test,
                            TA_fourier_sin_test, TC_fourier_sin_test, TG_fourier_sin_test, TT_fourier_sin_test)
saveRDS(fourier_di_cos_test, "data/Created/random_fourier_di_cos_test.rds")
saveRDS(fourier_di_sin_test, "data/Created/random_fourier_di_sin_test.rds")



fourier_AA_AC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AC_fourier_test[,"fourierAC_6"])
fourier_AA_AG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AA_AT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AA_CA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AA_CC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AA_CG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AA_CT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AA_GA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AA_GC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AA_GG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AA_GT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AA_TA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AA_TC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AA_TG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AA_TT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AA_6_mod_test = paste0("fourier_AA_", dinucleotides[-1], "_6_mod_test")
fourier_AA_6_mod_random_test = data.frame(map(fourier_AA_6_mod_test, get))
colnames(fourier_AA_6_mod_random_test) = fourier_AA_6_mod

fourier_AC_AG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AC_AT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AC_CA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AC_CC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AC_CG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AC_CT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AC_GA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AC_GC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AC_GG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AC_GT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AC_TA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AC_TC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AC_TG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AC_TT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AC_6_mod_test = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod_test")
fourier_AC_6_mod_random_test = data.frame(map(fourier_AC_6_mod_test, get))
colnames(fourier_AC_6_mod_random_test) = fourier_AC_6_mod

fourier_AG_AT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AG_CA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AG_CC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AG_CG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AG_CT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AG_GA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AG_GC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AG_GG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AG_GT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AG_TA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AG_TC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AG_TG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AG_TT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AG_6_mod_test = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod_test")
fourier_AG_6_mod_random_test = data.frame(map(fourier_AG_6_mod_test, get))
colnames(fourier_AG_6_mod_random_test) = fourier_AG_6_mod

fourier_AT_CA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AT_CC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AT_CG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AT_CT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AT_GA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AT_GC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AT_GG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AT_GT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AT_TA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AT_TC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AT_TG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AT_TT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AT_6_mod_test = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod_test")
fourier_AT_6_mod_random_test = data.frame(map(fourier_AT_6_mod_test, get))
colnames(fourier_AT_6_mod_random_test) = fourier_AT_6_mod

fourier_CA_CC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_CA_CG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CA_CT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CA_GA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CA_GC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CA_GG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CA_GT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CA_TA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CA_TC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CA_TG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CA_TT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CA_6_mod_test = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod_test")
fourier_CA_6_mod_random_test = data.frame(map(fourier_CA_6_mod_test, get))
colnames(fourier_CA_6_mod_random_test) = fourier_CA_6_mod

fourier_CC_CG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CC_CT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CC_GA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CC_GC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CC_GG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CC_GT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CC_TA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CC_TC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CC_TG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CC_TT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CC_6_mod_test = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod_test")
fourier_CC_6_mod_random_test = data.frame(map(fourier_CC_6_mod_test, get))
colnames(fourier_CC_6_mod_random_test) = fourier_CC_6_mod

fourier_CG_CT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CG_GA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CG_GC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CG_GG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CG_GT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CG_TA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CG_TC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CG_TG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CG_TT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CG_6_mod_test = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod_test")
fourier_CG_6_mod_random_test = data.frame(map(fourier_CG_6_mod_test, get))
colnames(fourier_CG_6_mod_random_test) = fourier_CG_6_mod

fourier_CT_GA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CT_GC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CT_GG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CT_GT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CT_TA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CT_TC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CT_TG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CT_TT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CT_6_mod_test = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod_test")
fourier_CT_6_mod_random_test = data.frame(map(fourier_CT_6_mod_test, get))
colnames(fourier_CT_6_mod_random_test) = fourier_CT_6_mod

fourier_GA_GC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_GA_GG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GA_GT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GA_TA_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GA_TC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GA_TG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GA_TT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GA_6_mod_test = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod_test")
fourier_GA_6_mod_random_test = data.frame(map(fourier_GA_6_mod_test, get))
colnames(fourier_GA_6_mod_random_test) = fourier_GA_6_mod

fourier_GC_GG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GC_GT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GC_TA_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GC_TC_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GC_TG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GC_TT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GC_6_mod_test = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod_test")
fourier_GC_6_mod_random_test = data.frame(map(fourier_GC_6_mod_test, get))
colnames(fourier_GC_6_mod_random_test) = fourier_GC_6_mod

fourier_GG_GT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GG_TA_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GG_TC_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GG_TG_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GG_TT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GG_6_mod_test = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod_test")
fourier_GG_6_mod_random_test = data.frame(map(fourier_GG_6_mod_test, get))
colnames(fourier_GG_6_mod_random_test) = fourier_GG_6_mod

fourier_GT_TA_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GT_TC_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GT_TG_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GT_TT_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GT_6_mod_test = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod_test")
fourier_GT_6_mod_random_test = data.frame(map(fourier_GT_6_mod_test, get))
colnames(fourier_GT_6_mod_random_test) = fourier_GT_6_mod

fourier_TA_TC_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_TA_TG_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TA_TT_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TA_6_mod_test = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod_test")
fourier_TA_6_mod_random_test = data.frame(map(fourier_TA_6_mod_test, get))
colnames(fourier_TA_6_mod_random_test) = fourier_TA_6_mod

fourier_TC_TG_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TC_TT_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TC_6_mod_test = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod_test")
fourier_TC_6_mod_random_test = data.frame(map(fourier_TC_6_mod_test, get))
colnames(fourier_TC_6_mod_random_test) = fourier_TC_6_mod

fourier_TG_TT_6_mod_test = Mod(TG_fourier_test[,"fourierTG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TG_6_mod_test = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod_test")
fourier_TG_6_mod_random_test = data.frame(map(fourier_TG_6_mod_test, get))
colnames(fourier_TG_6_mod_random_test) = fourier_TG_6_mod

fourier_6_mod_test = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_random_test"), get))
saveRDS(fourier_6_mod_test, "data/Created/random_fourier_6_mod_test.rds")






########################################################################
# chrV
########################################################################

# Train data
ps1 <- paste0("X", 1:50, "mono")

Xone = dat_chrV %>% select(all_of(ps1))
Xone_AorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AorT) = colnames(Xone)
colnames(Xone_CorG) = colnames(Xone)
Xone_AorT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier = t(apply(Xone_AorT, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos = Re(Xone_AorT_fourier)
Xone_AorT_fourier_sin = Im(Xone_AorT_fourier)

colnames(Xone_AorT_fourier_cos) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin) = paste0("fouriersinAorT", c(1:25, 50))

# colnames(Xone_CorG_fourier_cos) = paste0("fouriercosCorG", c(1:25, 50))
# colnames(Xone_CorG_fourier_sin) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos, "data/Created/chrV_Xone_AorT_fourier_cos.rds")
saveRDS(Xone_AorT_fourier_sin, "data/Created/chrV_Xone_AorT_fourier_sin.rds")

# saveRDS(Xone_CorG_fourier_cos, "data/Created/chrV_Xone_CorG_fourier_cos.rds")
# saveRDS(Xone_CorG_fourier_sin, "data/Created/chrV_Xone_CorG_fourier_sin.rds")

ps2 <- paste0("X", 1:49, "di")

Xtwo = dat_chrV %>% select(all_of(ps2))
Xtwo_AorT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CorG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AorT) = colnames(Xtwo)
colnames(Xtwo_CorG) = colnames(Xtwo)
Xtwo_AorT[] = ((Xtwo == "AA") | (Xtwo == "TT") | (Xtwo == "AT") | (Xtwo == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CorG[] = ((Xtwo == "CC") | (Xtwo == "GG") | (Xtwo == "CG") | (Xtwo == "GC")) %>%
  as.matrix() %>% as.numeric()

Xtwo_AorT_fourier = t(apply(Xtwo_AorT, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos = Re(Xtwo_AorT_fourier)
Xtwo_AorT_fourier_sin = Im(Xtwo_AorT_fourier)

colnames(Xtwo_AorT_fourier_cos) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin) = paste0("fouriersinAorT2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos, "data/Created/chrV_Xtwo_AorT_fourier_cos.rds")
saveRDS(Xtwo_AorT_fourier_sin, "data/Created/chrV_Xtwo_AorT_fourier_sin.rds")



Xtwo_AA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AA) = colnames(Xtwo)
Xtwo_AA[] = (Xtwo == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier = t(apply(Xtwo_AA, 1, fft))[,1:25]
colnames(AA_fourier) = paste0("fourierAA_", 1:25)
AA_fourier_cos = Re(AA_fourier)
AA_fourier_sin = Im(AA_fourier)
colnames(AA_fourier_cos) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos, "data/Created/chrV_AA_fourier_cos.rds")
saveRDS(AA_fourier_sin, "data/Created/chrV_AA_fourier_sin.rds")

Xtwo_AC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AC) = colnames(Xtwo)
Xtwo_AC[] = (Xtwo == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier = t(apply(Xtwo_AC, 1, fft))[,1:25]
colnames(AC_fourier) = paste0("fourierAC_", 1:25)
AC_fourier_cos = Re(AC_fourier)
AC_fourier_sin = Im(AC_fourier)
colnames(AC_fourier_cos) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos, "data/Created/chrV_AC_fourier_cos.rds")
saveRDS(AC_fourier_sin, "data/Created/chrV_AC_fourier_sin.rds")

Xtwo_AG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AG) = colnames(Xtwo)
Xtwo_AG[] = (Xtwo == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier = t(apply(Xtwo_AG, 1, fft))[,1:25]
colnames(AG_fourier) = paste0("fourierAG_", 1:25)
AG_fourier_cos = Re(AG_fourier)
AG_fourier_sin = Im(AG_fourier)
colnames(AG_fourier_cos) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos, "data/Created/chrV_AG_fourier_cos.rds")
saveRDS(AG_fourier_sin, "data/Created/chrV_AG_fourier_sin.rds")

Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AT) = colnames(Xtwo)
Xtwo_AT[] = (Xtwo == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier = t(apply(Xtwo_AT, 1, fft))[,1:25]
colnames(AT_fourier) = paste0("fourierAT_", 1:25)
AT_fourier_cos = Re(AT_fourier)
AT_fourier_sin = Im(AT_fourier)
colnames(AT_fourier_cos) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos, "data/Created/chrV_AT_fourier_cos.rds")
saveRDS(AT_fourier_sin, "data/Created/chrV_AT_fourier_sin.rds")

Xtwo_CA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CA) = colnames(Xtwo)
Xtwo_CA[] = (Xtwo == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier = t(apply(Xtwo_CA, 1, fft))[,1:25]
colnames(CA_fourier) = paste0("fourierCA_", 1:25)
CA_fourier_cos = Re(CA_fourier)
CA_fourier_sin = Im(CA_fourier)
colnames(CA_fourier_cos) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos, "data/Created/chrV_CA_fourier_cos.rds")
saveRDS(CA_fourier_sin, "data/Created/chrV_CA_fourier_sin.rds")

Xtwo_CC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CC) = colnames(Xtwo)
Xtwo_CC[] = (Xtwo == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier = t(apply(Xtwo_CC, 1, fft))[,1:25]
colnames(CC_fourier) = paste0("fourierCC_", 1:25)
CC_fourier_cos = Re(CC_fourier)
CC_fourier_sin = Im(CC_fourier)
colnames(CC_fourier_cos) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos, "data/Created/chrV_CC_fourier_cos.rds")
saveRDS(CC_fourier_sin, "data/Created/chrV_CC_fourier_sin.rds")

Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CG) = colnames(Xtwo)
Xtwo_CG[] = (Xtwo == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier = t(apply(Xtwo_CG, 1, fft))[,1:25]
colnames(CG_fourier) = paste0("fourierCG_", 1:25)
CG_fourier_cos = Re(CG_fourier)
CG_fourier_sin = Im(CG_fourier)
colnames(CG_fourier_cos) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos, "data/Created/chrV_CG_fourier_cos.rds")
saveRDS(CG_fourier_sin, "data/Created/chrV_CG_fourier_sin.rds")

Xtwo_CT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CT) = colnames(Xtwo)
Xtwo_CT[] = (Xtwo == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier = t(apply(Xtwo_CT, 1, fft))[,1:25]
colnames(CT_fourier) = paste0("fourierCT_", 1:25)
CT_fourier_cos = Re(CT_fourier)
CT_fourier_sin = Im(CT_fourier)
colnames(CT_fourier_cos) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos, "data/Created/chrV_CT_fourier_cos.rds")
saveRDS(CT_fourier_sin, "data/Created/chrV_CT_fourier_sin.rds")

Xtwo_GA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GA) = colnames(Xtwo)
Xtwo_GA[] = (Xtwo == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier = t(apply(Xtwo_GA, 1, fft))[,1:25]
colnames(GA_fourier) = paste0("fourierGA_", 1:25)
GA_fourier_cos = Re(GA_fourier)
GA_fourier_sin = Im(GA_fourier)
colnames(GA_fourier_cos) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos, "data/Created/chrV_GA_fourier_cos.rds")
saveRDS(GA_fourier_sin, "data/Created/chrV_GA_fourier_sin.rds")

Xtwo_GC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GC) = colnames(Xtwo)
Xtwo_GC[] = (Xtwo == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier = t(apply(Xtwo_GC, 1, fft))[,1:25]
colnames(GC_fourier) = paste0("fourierGC_", 1:25)
GC_fourier_cos = Re(GC_fourier)
GC_fourier_sin = Im(GC_fourier)
colnames(GC_fourier_cos) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos, "data/Created/chrV_GC_fourier_cos.rds")
saveRDS(GC_fourier_sin, "data/Created/chrV_GC_fourier_sin.rds")

Xtwo_GG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GG) = colnames(Xtwo)
Xtwo_GG[] = (Xtwo == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier = t(apply(Xtwo_GG, 1, fft))[,1:25]
colnames(GG_fourier) = paste0("fourierGG_", 1:25)
GG_fourier_cos = Re(GG_fourier)
GG_fourier_sin = Im(GG_fourier)
colnames(GG_fourier_cos) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos, "data/Created/chrV_GG_fourier_cos.rds")
saveRDS(GG_fourier_sin, "data/Created/chrV_GG_fourier_sin.rds")

Xtwo_GT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GT) = colnames(Xtwo)
Xtwo_GT[] = (Xtwo == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier = t(apply(Xtwo_GT, 1, fft))[,1:25]
colnames(GT_fourier) = paste0("fourierGT_", 1:25)
GT_fourier_cos = Re(GT_fourier)
GT_fourier_sin = Im(GT_fourier)
colnames(GT_fourier_cos) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos, "data/Created/chrV_GT_fourier_cos.rds")
saveRDS(GT_fourier_sin, "data/Created/chrV_GT_fourier_sin.rds")

Xtwo_TA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TA) = colnames(Xtwo)
Xtwo_TA[] = (Xtwo == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier = t(apply(Xtwo_TA, 1, fft))[,1:25]
colnames(TA_fourier) = paste0("fourierTA_", 1:25)
TA_fourier_cos = Re(TA_fourier)
TA_fourier_sin = Im(TA_fourier)
colnames(TA_fourier_cos) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos, "data/Created/chrV_TA_fourier_cos.rds")
saveRDS(TA_fourier_sin, "data/Created/chrV_TA_fourier_sin.rds")

Xtwo_TC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TC) = colnames(Xtwo)
Xtwo_TC[] = (Xtwo == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier = t(apply(Xtwo_TC, 1, fft))[,1:25]
colnames(TC_fourier) = paste0("fourierTC_", 1:25)
TC_fourier_cos = Re(TC_fourier)
TC_fourier_sin = Im(TC_fourier)
colnames(TC_fourier_cos) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos, "data/Created/chrV_TC_fourier_cos.rds")
saveRDS(TC_fourier_sin, "data/Created/chrV_TC_fourier_sin.rds")

Xtwo_TG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TG) = colnames(Xtwo)
Xtwo_TG[] = (Xtwo == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier = t(apply(Xtwo_TG, 1, fft))[,1:25]
colnames(TG_fourier) = paste0("fourierTG_", 1:25)
TG_fourier_cos = Re(TG_fourier)
TG_fourier_sin = Im(TG_fourier)
colnames(TG_fourier_cos) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos, "data/Created/chrV_TG_fourier_cos.rds")
saveRDS(TG_fourier_sin, "data/Created/chrV_TG_fourier_sin.rds")

Xtwo_TT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TT) = colnames(Xtwo)
Xtwo_TT[] = (Xtwo == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier = t(apply(Xtwo_TT, 1, fft))[,1:25]
colnames(TT_fourier) = paste0("fourierTT_", 1:25)
TT_fourier_cos = Re(TT_fourier)
TT_fourier_sin = Im(TT_fourier)
colnames(TT_fourier_cos) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos, "data/Created/chrV_TT_fourier_cos.rds")
saveRDS(TT_fourier_sin, "data/Created/chrV_TT_fourier_sin.rds")

fourier_di_cos = cbind(AA_fourier_cos, AC_fourier_cos, AG_fourier_cos, AT_fourier_cos, 
                       CA_fourier_cos, CC_fourier_cos, CG_fourier_cos, CT_fourier_cos,
                       GA_fourier_cos, GC_fourier_cos, GG_fourier_cos, GT_fourier_cos,
                       TA_fourier_cos, TC_fourier_cos, TG_fourier_cos, TT_fourier_cos)
fourier_di_sin = cbind(AA_fourier_sin, AC_fourier_sin, AG_fourier_sin, AT_fourier_sin, 
                       CA_fourier_sin, CC_fourier_sin, CG_fourier_sin, CT_fourier_sin,
                       GA_fourier_sin, GC_fourier_sin, GG_fourier_sin, GT_fourier_sin,
                       TA_fourier_sin, TC_fourier_sin, TG_fourier_sin, TT_fourier_sin)
saveRDS(fourier_di_cos, "data/Created/chrV_fourier_di_cos.rds")
saveRDS(fourier_di_sin, "data/Created/chrV_fourier_di_sin.rds")

fourier_AA_AC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AC_fourier[,"fourierAC_6"])
fourier_AA_AG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AG_fourier[,"fourierAG_6"])
fourier_AA_AT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AT_fourier[,"fourierAT_6"])
fourier_AA_CA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CA_fourier[,"fourierCA_6"])
fourier_AA_CC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CC_fourier[,"fourierCC_6"])
fourier_AA_CG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CG_fourier[,"fourierCG_6"])
fourier_AA_CT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CT_fourier[,"fourierCT_6"])
fourier_AA_GA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GA_fourier[,"fourierGA_6"])
fourier_AA_GC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GC_fourier[,"fourierGC_6"])
fourier_AA_GG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GG_fourier[,"fourierGG_6"])
fourier_AA_GT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GT_fourier[,"fourierGT_6"])
fourier_AA_TA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TA_fourier[,"fourierTA_6"])
fourier_AA_TC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TC_fourier[,"fourierTC_6"])
fourier_AA_TG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TG_fourier[,"fourierTG_6"])
fourier_AA_TT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TT_fourier[,"fourierTT_6"])

fourier_AA_6_mod = paste0("fourier_AA_", dinucleotides[-1], "_6_mod")
fourier_AA_6_mod_chrV = data.frame(map(fourier_AA_6_mod, get))
colnames(fourier_AA_6_mod_chrV) = fourier_AA_6_mod

fourier_AC_AG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AG_fourier[,"fourierAG_6"])
fourier_AC_AT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AT_fourier[,"fourierAT_6"])
fourier_AC_CA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CA_fourier[,"fourierCA_6"])
fourier_AC_CC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CC_fourier[,"fourierCC_6"])
fourier_AC_CG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CG_fourier[,"fourierCG_6"])
fourier_AC_CT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CT_fourier[,"fourierCT_6"])
fourier_AC_GA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GA_fourier[,"fourierGA_6"])
fourier_AC_GC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GC_fourier[,"fourierGC_6"])
fourier_AC_GG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GG_fourier[,"fourierGG_6"])
fourier_AC_GT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GT_fourier[,"fourierGT_6"])
fourier_AC_TA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TA_fourier[,"fourierTA_6"])
fourier_AC_TC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TC_fourier[,"fourierTC_6"])
fourier_AC_TG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TG_fourier[,"fourierTG_6"])
fourier_AC_TT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TT_fourier[,"fourierTT_6"])

fourier_AC_6_mod = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod")
fourier_AC_6_mod_chrV = data.frame(map(fourier_AC_6_mod, get))
colnames(fourier_AC_6_mod_chrV) = fourier_AC_6_mod

fourier_AG_AT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - AT_fourier[,"fourierAT_6"])
fourier_AG_CA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CA_fourier[,"fourierCA_6"])
fourier_AG_CC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CC_fourier[,"fourierCC_6"])
fourier_AG_CG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CG_fourier[,"fourierCG_6"])
fourier_AG_CT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CT_fourier[,"fourierCT_6"])
fourier_AG_GA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GA_fourier[,"fourierGA_6"])
fourier_AG_GC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GC_fourier[,"fourierGC_6"])
fourier_AG_GG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GG_fourier[,"fourierGG_6"])
fourier_AG_GT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GT_fourier[,"fourierGT_6"])
fourier_AG_TA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TA_fourier[,"fourierTA_6"])
fourier_AG_TC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TC_fourier[,"fourierTC_6"])
fourier_AG_TG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TG_fourier[,"fourierTG_6"])
fourier_AG_TT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TT_fourier[,"fourierTT_6"])

fourier_AG_6_mod = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod")
fourier_AG_6_mod_chrV = data.frame(map(fourier_AG_6_mod, get))
colnames(fourier_AG_6_mod_chrV) = fourier_AG_6_mod

fourier_AT_CA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CA_fourier[,"fourierCA_6"])
fourier_AT_CC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CC_fourier[,"fourierCC_6"])
fourier_AT_CG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CG_fourier[,"fourierCG_6"])
fourier_AT_CT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CT_fourier[,"fourierCT_6"])
fourier_AT_GA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GA_fourier[,"fourierGA_6"])
fourier_AT_GC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GC_fourier[,"fourierGC_6"])
fourier_AT_GG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GG_fourier[,"fourierGG_6"])
fourier_AT_GT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GT_fourier[,"fourierGT_6"])
fourier_AT_TA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TA_fourier[,"fourierTA_6"])
fourier_AT_TC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TC_fourier[,"fourierTC_6"])
fourier_AT_TG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TG_fourier[,"fourierTG_6"])
fourier_AT_TT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TT_fourier[,"fourierTT_6"])

fourier_AT_6_mod = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod")
fourier_AT_6_mod_chrV = data.frame(map(fourier_AT_6_mod, get))
colnames(fourier_AT_6_mod_chrV) = fourier_AT_6_mod

fourier_CA_CC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CC_fourier[,"fourierCC_6"])
fourier_CA_CG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CG_fourier[,"fourierCG_6"])
fourier_CA_CT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CT_fourier[,"fourierCT_6"])
fourier_CA_GA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GA_fourier[,"fourierGA_6"])
fourier_CA_GC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GC_fourier[,"fourierGC_6"])
fourier_CA_GG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GG_fourier[,"fourierGG_6"])
fourier_CA_GT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GT_fourier[,"fourierGT_6"])
fourier_CA_TA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TA_fourier[,"fourierTA_6"])
fourier_CA_TC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TC_fourier[,"fourierTC_6"])
fourier_CA_TG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TG_fourier[,"fourierTG_6"])
fourier_CA_TT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TT_fourier[,"fourierTT_6"])

fourier_CA_6_mod = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod")
fourier_CA_6_mod_chrV = data.frame(map(fourier_CA_6_mod, get))
colnames(fourier_CA_6_mod_chrV) = fourier_CA_6_mod

fourier_CC_CG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CG_fourier[,"fourierCG_6"])
fourier_CC_CT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CT_fourier[,"fourierCT_6"])
fourier_CC_GA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GA_fourier[,"fourierGA_6"])
fourier_CC_GC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GC_fourier[,"fourierGC_6"])
fourier_CC_GG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GG_fourier[,"fourierGG_6"])
fourier_CC_GT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GT_fourier[,"fourierGT_6"])
fourier_CC_TA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TA_fourier[,"fourierTA_6"])
fourier_CC_TC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TC_fourier[,"fourierTC_6"])
fourier_CC_TG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TG_fourier[,"fourierTG_6"])
fourier_CC_TT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TT_fourier[,"fourierTT_6"])

fourier_CC_6_mod = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod")
fourier_CC_6_mod_chrV = data.frame(map(fourier_CC_6_mod, get))
colnames(fourier_CC_6_mod_chrV) = fourier_CC_6_mod

fourier_CG_CT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - CT_fourier[,"fourierCT_6"])
fourier_CG_GA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GA_fourier[,"fourierGA_6"])
fourier_CG_GC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GC_fourier[,"fourierGC_6"])
fourier_CG_GG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GG_fourier[,"fourierGG_6"])
fourier_CG_GT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GT_fourier[,"fourierGT_6"])
fourier_CG_TA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TA_fourier[,"fourierTA_6"])
fourier_CG_TC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TC_fourier[,"fourierTC_6"])
fourier_CG_TG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TG_fourier[,"fourierTG_6"])
fourier_CG_TT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TT_fourier[,"fourierTT_6"])

fourier_CG_6_mod = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod")
fourier_CG_6_mod_chrV = data.frame(map(fourier_CG_6_mod, get))
colnames(fourier_CG_6_mod_chrV) = fourier_CG_6_mod

fourier_CT_GA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GA_fourier[,"fourierGA_6"])
fourier_CT_GC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GC_fourier[,"fourierGC_6"])
fourier_CT_GG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GG_fourier[,"fourierGG_6"])
fourier_CT_GT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GT_fourier[,"fourierGT_6"])
fourier_CT_TA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TA_fourier[,"fourierTA_6"])
fourier_CT_TC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TC_fourier[,"fourierTC_6"])
fourier_CT_TG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TG_fourier[,"fourierTG_6"])
fourier_CT_TT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TT_fourier[,"fourierTT_6"])

fourier_CT_6_mod = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod")
fourier_CT_6_mod_chrV = data.frame(map(fourier_CT_6_mod, get))
colnames(fourier_CT_6_mod_chrV) = fourier_CT_6_mod

fourier_GA_GC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GC_fourier[,"fourierGC_6"])
fourier_GA_GG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GG_fourier[,"fourierGG_6"])
fourier_GA_GT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GT_fourier[,"fourierGT_6"])
fourier_GA_TA_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TA_fourier[,"fourierTA_6"])
fourier_GA_TC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TC_fourier[,"fourierTC_6"])
fourier_GA_TG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TG_fourier[,"fourierTG_6"])
fourier_GA_TT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TT_fourier[,"fourierTT_6"])

fourier_GA_6_mod = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod")
fourier_GA_6_mod_chrV = data.frame(map(fourier_GA_6_mod, get))
colnames(fourier_GA_6_mod_chrV) = fourier_GA_6_mod

fourier_GC_GG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GG_fourier[,"fourierGG_6"])
fourier_GC_GT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GT_fourier[,"fourierGT_6"])
fourier_GC_TA_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TA_fourier[,"fourierTA_6"])
fourier_GC_TC_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TC_fourier[,"fourierTC_6"])
fourier_GC_TG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TG_fourier[,"fourierTG_6"])
fourier_GC_TT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TT_fourier[,"fourierTT_6"])

fourier_GC_6_mod = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod")
fourier_GC_6_mod_chrV = data.frame(map(fourier_GC_6_mod, get))
colnames(fourier_GC_6_mod_chrV) = fourier_GC_6_mod

fourier_GG_GT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - GT_fourier[,"fourierGT_6"])
fourier_GG_TA_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TA_fourier[,"fourierTA_6"])
fourier_GG_TC_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TC_fourier[,"fourierTC_6"])
fourier_GG_TG_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TG_fourier[,"fourierTG_6"])
fourier_GG_TT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TT_fourier[,"fourierTT_6"])

fourier_GG_6_mod = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod")
fourier_GG_6_mod_chrV = data.frame(map(fourier_GG_6_mod, get))
colnames(fourier_GG_6_mod_chrV) = fourier_GG_6_mod

fourier_GT_TA_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TA_fourier[,"fourierTA_6"])
fourier_GT_TC_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TC_fourier[,"fourierTC_6"])
fourier_GT_TG_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TG_fourier[,"fourierTG_6"])
fourier_GT_TT_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TT_fourier[,"fourierTT_6"])

fourier_GT_6_mod = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod")
fourier_GT_6_mod_chrV = data.frame(map(fourier_GT_6_mod, get))
colnames(fourier_GT_6_mod_chrV) = fourier_GT_6_mod

fourier_TA_TC_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TC_fourier[,"fourierTC_6"])
fourier_TA_TG_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TG_fourier[,"fourierTG_6"])
fourier_TA_TT_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TT_fourier[,"fourierTT_6"])

fourier_TA_6_mod = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod")
fourier_TA_6_mod_chrV = data.frame(map(fourier_TA_6_mod, get))
colnames(fourier_TA_6_mod_chrV) = fourier_TA_6_mod

fourier_TC_TG_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TG_fourier[,"fourierTG_6"])
fourier_TC_TT_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TT_fourier[,"fourierTT_6"])

fourier_TC_6_mod = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod")
fourier_TC_6_mod_chrV = data.frame(map(fourier_TC_6_mod, get))
colnames(fourier_TC_6_mod_chrV) = fourier_TC_6_mod

fourier_TG_TT_6_mod = Mod(TG_fourier[,"fourierTG_6"] - TT_fourier[,"fourierTT_6"])

fourier_TG_6_mod = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod")
fourier_TG_6_mod_chrV = data.frame(map(fourier_TG_6_mod, get))
colnames(fourier_TG_6_mod_chrV) = fourier_TG_6_mod

fourier_6_mod = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_chrV"), get))
saveRDS(fourier_6_mod, "data/Created/chrV_fourier_6_mod.rds")




# Test data:
Xone_test = dat_chrV_test %>% select(all_of(ps1))
Xone_AorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AorT_test) = colnames(Xone_test)
colnames(Xone_CorG_test) = colnames(Xone_test)
Xone_AorT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_test = Re(Xone_AorT_fourier_test)
Xone_AorT_fourier_sin_test = Im(Xone_AorT_fourier_test)

Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,c(1:25, 50)]
Xone_CorG_fourier_cos_test = Re(Xone_CorG_fourier_test)
Xone_CorG_fourier_sin_test = Im(Xone_CorG_fourier_test)

colnames(Xone_AorT_fourier_cos_test) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_test) = paste0("fouriersinAorT", c(1:25, 50))

colnames(Xone_CorG_fourier_cos_test) = paste0("fouriercosCorG", c(1:25, 50))
colnames(Xone_CorG_fourier_sin_test) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_test, "data/Created/chrV_Xone_AorT_fourier_cos_test.rds")
saveRDS(Xone_AorT_fourier_sin_test, "data/Created/chrV_Xone_AorT_fourier_sin_test.rds")

saveRDS(Xone_CorG_fourier_cos_test, "data/Created/chrV_Xone_CorG_fourier_cos_test.rds")
saveRDS(Xone_CorG_fourier_sin_test, "data/Created/chrV_Xone_CorG_fourier_sin_test.rds")


Xtwo_test = dat_chrV_test %>% select(all_of(ps2))
Xtwo_AorT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
# Xtwo_CorG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AorT_test) = colnames(Xtwo_test)
# colnames(Xtwo_CorG_test) = colnames(Xtwo_test)
Xtwo_AorT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT") | (Xtwo_test == "AT") | (Xtwo_test == "TA")) %>% 
  as.matrix() %>% as.numeric()
# Xtwo_CorG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG") | (Xtwo_test == "CG") | (Xtwo_test == "GC")) %>%
#   as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_test = t(apply(Xtwo_AorT_test, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_test = Re(Xtwo_AorT_fourier_test)
Xtwo_AorT_fourier_sin_test = Im(Xtwo_AorT_fourier_test)

# Xtwo_CorG_fourier = t(apply(Xtwo_CorG_test, 1, fft))[,c(1:25, 50)]
# Xtwo_CorG_fourier_cos_test = Re(Xtwo_CorG_fourier)
# Xtwo_CorG_fourier_sin_test = Im(Xtwo_CorG_fourier)

colnames(Xtwo_AorT_fourier_cos_test) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_test) = paste0("fouriersinAorT2_", 1:25)

# colnames(Xtwo_CorG_fourier_cos_test) = paste0("fouriercosCorG2_", 1:25)
# colnames(Xtwo_CorG_fourier_sin_test) = paste0("fouriersinCorG2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_test, "data/Created/chrV_Xtwo_AorT_fourier_cos_test.rds")
saveRDS(Xtwo_AorT_fourier_sin_test, "data/Created/chrV_Xtwo_AorT_fourier_sin_test.rds")

# saveRDS(Xtwo_CorG_fourier_cos_test, "data/Created/chrV_Xtwo_CorG_fourier_cos_test.rds")
# saveRDS(Xtwo_CorG_fourier_sin_test, "data/Created/chrV_Xtwo_CorG_fourier_sin_test.rds")


Xtwo_AA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AA_test) = colnames(Xtwo_test)
Xtwo_AA_test[] = (Xtwo_test == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,1:25]
colnames(AA_fourier_test) = paste0("fourierAA_", 1:25)
AA_fourier_cos_test = Re(AA_fourier_test)
AA_fourier_sin_test = Im(AA_fourier_test)
colnames(AA_fourier_cos_test) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin_test) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos_test, "data/Created/chrV_AA_fourier_cos_test.rds")
saveRDS(AA_fourier_sin_test, "data/Created/chrV_AA_fourier_sin_test.rds")

Xtwo_AC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AC_test) = colnames(Xtwo_test)
Xtwo_AC_test[] = (Xtwo_test == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,1:25]
colnames(AC_fourier_test) = paste0("fourierAC_", 1:25)
AC_fourier_cos_test = Re(AC_fourier_test)
AC_fourier_sin_test = Im(AC_fourier_test)
colnames(AC_fourier_cos_test) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin_test) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos_test, "data/Created/chrV_AC_fourier_cos_test.rds")
saveRDS(AC_fourier_sin_test, "data/Created/chrV_AC_fourier_sin_test.rds")

Xtwo_AG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AG_test) = colnames(Xtwo_test)
Xtwo_AG_test[] = (Xtwo_test == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,1:25]
colnames(AG_fourier_test) = paste0("fourierAG_", 1:25)
AG_fourier_cos_test = Re(AG_fourier_test)
AG_fourier_sin_test = Im(AG_fourier_test)
colnames(AG_fourier_cos_test) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin_test) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos_test, "data/Created/chrV_AG_fourier_cos_test.rds")
saveRDS(AG_fourier_sin_test, "data/Created/chrV_AG_fourier_sin_test.rds")

Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
Xtwo_AT_test[] = (Xtwo_test == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,1:25]
colnames(AT_fourier_test) = paste0("fourierAT_", 1:25)
AT_fourier_cos_test = Re(AT_fourier_test)
AT_fourier_sin_test = Im(AT_fourier_test)
colnames(AT_fourier_cos_test) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin_test) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos_test, "data/Created/chrV_AT_fourier_cos_test.rds")
saveRDS(AT_fourier_sin_test, "data/Created/chrV_AT_fourier_sin_test.rds")

Xtwo_CA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CA_test) = colnames(Xtwo_test)
Xtwo_CA_test[] = (Xtwo_test == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,1:25]
colnames(CA_fourier_test) = paste0("fourierCA_", 1:25)
CA_fourier_cos_test = Re(CA_fourier_test)
CA_fourier_sin_test = Im(CA_fourier_test)
colnames(CA_fourier_cos_test) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin_test) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos_test, "data/Created/chrV_CA_fourier_cos_test.rds")
saveRDS(CA_fourier_sin_test, "data/Created/chrV_CA_fourier_sin_test.rds")

Xtwo_CC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CC_test) = colnames(Xtwo_test)
Xtwo_CC_test[] = (Xtwo_test == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,1:25]
colnames(CC_fourier_test) = paste0("fourierCC_", 1:25)
CC_fourier_cos_test = Re(CC_fourier_test)
CC_fourier_sin_test = Im(CC_fourier_test)
colnames(CC_fourier_cos_test) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin_test) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos_test, "data/Created/chrV_CC_fourier_cos_test.rds")
saveRDS(CC_fourier_sin_test, "data/Created/chrV_CC_fourier_sin_test.rds")

Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
Xtwo_CG_test[] = (Xtwo_test == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,1:25]
colnames(CG_fourier_test) = paste0("fourierCG_", 1:25)
CG_fourier_cos_test = Re(CG_fourier_test)
CG_fourier_sin_test = Im(CG_fourier_test)
colnames(CG_fourier_cos_test) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin_test) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos_test, "data/Created/chrV_CG_fourier_cos_test.rds")
saveRDS(CG_fourier_sin_test, "data/Created/chrV_CG_fourier_sin_test.rds")

Xtwo_CT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CT_test) = colnames(Xtwo_test)
Xtwo_CT_test[] = (Xtwo_test == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,1:25]
colnames(CT_fourier_test) = paste0("fourierCT_", 1:25)
CT_fourier_cos_test = Re(CT_fourier_test)
CT_fourier_sin_test = Im(CT_fourier_test)
colnames(CT_fourier_cos_test) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin_test) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos_test, "data/Created/chrV_CT_fourier_cos_test.rds")
saveRDS(CT_fourier_sin_test, "data/Created/chrV_CT_fourier_sin_test.rds")

Xtwo_GA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GA_test) = colnames(Xtwo_test)
Xtwo_GA_test[] = (Xtwo_test == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,1:25]
colnames(GA_fourier_test) = paste0("fourierGA_", 1:25)
GA_fourier_cos_test = Re(GA_fourier_test)
GA_fourier_sin_test = Im(GA_fourier_test)
colnames(GA_fourier_cos_test) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin_test) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos_test, "data/Created/chrV_GA_fourier_cos_test.rds")
saveRDS(GA_fourier_sin_test, "data/Created/chrV_GA_fourier_sin_test.rds")

Xtwo_GC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GC_test) = colnames(Xtwo_test)
Xtwo_GC_test[] = (Xtwo_test == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,1:25]
colnames(GC_fourier_test) = paste0("fourierGC_", 1:25)
GC_fourier_cos_test = Re(GC_fourier_test)
GC_fourier_sin_test = Im(GC_fourier_test)
colnames(GC_fourier_cos_test) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin_test) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos_test, "data/Created/chrV_GC_fourier_cos_test.rds")
saveRDS(GC_fourier_sin_test, "data/Created/chrV_GC_fourier_sin_test.rds")

Xtwo_GG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GG_test) = colnames(Xtwo_test)
Xtwo_GG_test[] = (Xtwo_test == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,1:25]
colnames(GG_fourier_test) = paste0("fourierGG_", 1:25)
GG_fourier_cos_test = Re(GG_fourier_test)
GG_fourier_sin_test = Im(GG_fourier_test)
colnames(GG_fourier_cos_test) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin_test) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos_test, "data/Created/chrV_GG_fourier_cos_test.rds")
saveRDS(GG_fourier_sin_test, "data/Created/chrV_GG_fourier_sin_test.rds")

Xtwo_GT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GT_test) = colnames(Xtwo_test)
Xtwo_GT_test[] = (Xtwo_test == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,1:25]
colnames(GT_fourier_test) = paste0("fourierGT_", 1:25)
GT_fourier_cos_test = Re(GT_fourier_test)
GT_fourier_sin_test = Im(GT_fourier_test)
colnames(GT_fourier_cos_test) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin_test) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos_test, "data/Created/chrV_GT_fourier_cos_test.rds")
saveRDS(GT_fourier_sin_test, "data/Created/chrV_GT_fourier_sin_test.rds")

Xtwo_TA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TA_test) = colnames(Xtwo_test)
Xtwo_TA_test[] = (Xtwo_test == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,1:25]
colnames(TA_fourier_test) = paste0("fourierTA_", 1:25)
TA_fourier_cos_test = Re(TA_fourier_test)
TA_fourier_sin_test = Im(TA_fourier_test)
colnames(TA_fourier_cos_test) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin_test) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos_test, "data/Created/chrV_TA_fourier_cos_test.rds")
saveRDS(TA_fourier_sin_test, "data/Created/chrV_TA_fourier_sin_test.rds")

Xtwo_TC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TC_test) = colnames(Xtwo_test)
Xtwo_TC_test[] = (Xtwo_test == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,1:25]
colnames(TC_fourier_test) = paste0("fourierTC_", 1:25)
TC_fourier_cos_test = Re(TC_fourier_test)
TC_fourier_sin_test = Im(TC_fourier_test)
colnames(TC_fourier_cos_test) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin_test) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos_test, "data/Created/chrV_TC_fourier_cos_test.rds")
saveRDS(TC_fourier_sin_test, "data/Created/chrV_TC_fourier_sin_test.rds")

Xtwo_TG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TG_test) = colnames(Xtwo_test)
Xtwo_TG_test[] = (Xtwo_test == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,1:25]
colnames(TG_fourier_test) = paste0("fourierTG_", 1:25)
TG_fourier_cos_test = Re(TG_fourier_test)
TG_fourier_sin_test = Im(TG_fourier_test)
colnames(TG_fourier_cos_test) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin_test) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos_test, "data/Created/chrV_TG_fourier_cos_test.rds")
saveRDS(TG_fourier_sin_test, "data/Created/chrV_TG_fourier_sin_test.rds")

Xtwo_TT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TT_test) = colnames(Xtwo_test)
Xtwo_TT_test[] = (Xtwo_test == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,1:25]
colnames(TT_fourier_test) = paste0("fourierTT_", 1:25)
TT_fourier_cos_test = Re(TT_fourier_test)
TT_fourier_sin_test = Im(TT_fourier_test)
colnames(TT_fourier_cos_test) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin_test) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos_test, "data/Created/chrV_TT_fourier_cos_test.rds")
saveRDS(TT_fourier_sin_test, "data/Created/chrV_TT_fourier_sin_test.rds")

fourier_di_cos_test = cbind(AA_fourier_cos_test, AC_fourier_cos_test, AG_fourier_cos_test, AT_fourier_cos_test, 
                            CA_fourier_cos_test, CC_fourier_cos_test, CG_fourier_cos_test, CT_fourier_cos_test,
                            GA_fourier_cos_test, GC_fourier_cos_test, GG_fourier_cos_test, GT_fourier_cos_test,
                            TA_fourier_cos_test, TC_fourier_cos_test, TG_fourier_cos_test, TT_fourier_cos_test)
fourier_di_sin_test = cbind(AA_fourier_sin_test, AC_fourier_sin_test, AG_fourier_sin_test, AT_fourier_sin_test, 
                            CA_fourier_sin_test, CC_fourier_sin_test, CG_fourier_sin_test, CT_fourier_sin_test,
                            GA_fourier_sin_test, GC_fourier_sin_test, GG_fourier_sin_test, GT_fourier_sin_test,
                            TA_fourier_sin_test, TC_fourier_sin_test, TG_fourier_sin_test, TT_fourier_sin_test)
saveRDS(fourier_di_cos_test, "data/Created/chrV_fourier_di_cos_test.rds")
saveRDS(fourier_di_sin_test, "data/Created/chrV_fourier_di_sin_test.rds")



fourier_AA_AC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AC_fourier_test[,"fourierAC_6"])
fourier_AA_AG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AA_AT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AA_CA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AA_CC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AA_CG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AA_CT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AA_GA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AA_GC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AA_GG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AA_GT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AA_TA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AA_TC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AA_TG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AA_TT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AA_6_mod_test = paste0("fourier_AA_", dinucleotides[-1], "_6_mod_test")
fourier_AA_6_mod_chrV_test = data.frame(map(fourier_AA_6_mod_test, get))
colnames(fourier_AA_6_mod_chrV_test) = fourier_AA_6_mod

fourier_AC_AG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AC_AT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AC_CA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AC_CC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AC_CG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AC_CT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AC_GA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AC_GC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AC_GG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AC_GT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AC_TA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AC_TC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AC_TG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AC_TT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AC_6_mod_test = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod_test")
fourier_AC_6_mod_chrV_test = data.frame(map(fourier_AC_6_mod_test, get))
colnames(fourier_AC_6_mod_chrV_test) = fourier_AC_6_mod

fourier_AG_AT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AG_CA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AG_CC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AG_CG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AG_CT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AG_GA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AG_GC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AG_GG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AG_GT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AG_TA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AG_TC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AG_TG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AG_TT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AG_6_mod_test = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod_test")
fourier_AG_6_mod_chrV_test = data.frame(map(fourier_AG_6_mod_test, get))
colnames(fourier_AG_6_mod_chrV_test) = fourier_AG_6_mod

fourier_AT_CA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AT_CC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AT_CG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AT_CT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AT_GA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AT_GC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AT_GG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AT_GT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AT_TA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AT_TC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AT_TG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AT_TT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AT_6_mod_test = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod_test")
fourier_AT_6_mod_chrV_test = data.frame(map(fourier_AT_6_mod_test, get))
colnames(fourier_AT_6_mod_chrV_test) = fourier_AT_6_mod

fourier_CA_CC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_CA_CG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CA_CT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CA_GA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CA_GC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CA_GG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CA_GT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CA_TA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CA_TC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CA_TG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CA_TT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CA_6_mod_test = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod_test")
fourier_CA_6_mod_chrV_test = data.frame(map(fourier_CA_6_mod_test, get))
colnames(fourier_CA_6_mod_chrV_test) = fourier_CA_6_mod

fourier_CC_CG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CC_CT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CC_GA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CC_GC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CC_GG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CC_GT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CC_TA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CC_TC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CC_TG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CC_TT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CC_6_mod_test = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod_test")
fourier_CC_6_mod_chrV_test = data.frame(map(fourier_CC_6_mod_test, get))
colnames(fourier_CC_6_mod_chrV_test) = fourier_CC_6_mod

fourier_CG_CT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CG_GA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CG_GC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CG_GG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CG_GT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CG_TA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CG_TC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CG_TG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CG_TT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CG_6_mod_test = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod_test")
fourier_CG_6_mod_chrV_test = data.frame(map(fourier_CG_6_mod_test, get))
colnames(fourier_CG_6_mod_chrV_test) = fourier_CG_6_mod

fourier_CT_GA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CT_GC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CT_GG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CT_GT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CT_TA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CT_TC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CT_TG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CT_TT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CT_6_mod_test = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod_test")
fourier_CT_6_mod_chrV_test = data.frame(map(fourier_CT_6_mod_test, get))
colnames(fourier_CT_6_mod_chrV_test) = fourier_CT_6_mod

fourier_GA_GC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_GA_GG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GA_GT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GA_TA_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GA_TC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GA_TG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GA_TT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GA_6_mod_test = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod_test")
fourier_GA_6_mod_chrV_test = data.frame(map(fourier_GA_6_mod_test, get))
colnames(fourier_GA_6_mod_chrV_test) = fourier_GA_6_mod

fourier_GC_GG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GC_GT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GC_TA_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GC_TC_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GC_TG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GC_TT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GC_6_mod_test = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod_test")
fourier_GC_6_mod_chrV_test = data.frame(map(fourier_GC_6_mod_test, get))
colnames(fourier_GC_6_mod_chrV_test) = fourier_GC_6_mod

fourier_GG_GT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GG_TA_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GG_TC_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GG_TG_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GG_TT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GG_6_mod_test = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod_test")
fourier_GG_6_mod_chrV_test = data.frame(map(fourier_GG_6_mod_test, get))
colnames(fourier_GG_6_mod_chrV_test) = fourier_GG_6_mod

fourier_GT_TA_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GT_TC_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GT_TG_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GT_TT_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GT_6_mod_test = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod_test")
fourier_GT_6_mod_chrV_test = data.frame(map(fourier_GT_6_mod_test, get))
colnames(fourier_GT_6_mod_chrV_test) = fourier_GT_6_mod

fourier_TA_TC_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_TA_TG_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TA_TT_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TA_6_mod_test = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod_test")
fourier_TA_6_mod_chrV_test = data.frame(map(fourier_TA_6_mod_test, get))
colnames(fourier_TA_6_mod_chrV_test) = fourier_TA_6_mod

fourier_TC_TG_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TC_TT_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TC_6_mod_test = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod_test")
fourier_TC_6_mod_chrV_test = data.frame(map(fourier_TC_6_mod_test, get))
colnames(fourier_TC_6_mod_chrV_test) = fourier_TC_6_mod

fourier_TG_TT_6_mod_test = Mod(TG_fourier_test[,"fourierTG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TG_6_mod_test = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod_test")
fourier_TG_6_mod_chrV_test = data.frame(map(fourier_TG_6_mod_test, get))
colnames(fourier_TG_6_mod_chrV_test) = fourier_TG_6_mod

fourier_6_mod_test = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_chrV_test"), get))
saveRDS(fourier_6_mod_test, "data/Created/chrV_fourier_6_mod_test.rds")











########################################################################
# Yeast
########################################################################

# Train data
ps1 <- paste0("X", 1:50, "mono")

Xone = dat_yeast %>% select(all_of(ps1))
Xone_AorT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CorG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AorT) = colnames(Xone)
colnames(Xone_CorG) = colnames(Xone)
Xone_AorT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_yeast = t(apply(Xone_AorT, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_yeast = Re(Xone_AorT_fourier_yeast)
Xone_AorT_fourier_sin_yeast = Im(Xone_AorT_fourier_yeast)

colnames(Xone_AorT_fourier_cos_yeast) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_yeast) = paste0("fouriersinAorT", c(1:25, 50))

# colnames(Xone_CorG_fourier_cos_yeast) = paste0("fouriercosCorG", c(1:25, 50))
# colnames(Xone_CorG_fourier_sin_yeast) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_yeast, "data/Created/yeast_Xone_AorT_fourier_cos.rds")
saveRDS(Xone_AorT_fourier_sin_yeast, "data/Created/yeast_Xone_AorT_fourier_sin.rds")

# saveRDS(Xone_CorG_fourier_cos, "data/Created/chrV_Xone_CorG_fourier_cos.rds")
# saveRDS(Xone_CorG_fourier_sin, "data/Created/chrV_Xone_CorG_fourier_sin.rds")

ps2 <- paste0("X", 1:49, "di")

Xtwo = dat_yeast %>% select(all_of(ps2))
Xtwo_AorT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CorG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AorT) = colnames(Xtwo)
colnames(Xtwo_CorG) = colnames(Xtwo)
Xtwo_AorT[] = ((Xtwo == "AA") | (Xtwo == "TT") | (Xtwo == "AT") | (Xtwo == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CorG[] = ((Xtwo == "CC") | (Xtwo == "GG") | (Xtwo == "CG") | (Xtwo == "GC")) %>%
  as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_yeast = t(apply(Xtwo_AorT, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_yeast = Re(Xtwo_AorT_fourier_yeast)
Xtwo_AorT_fourier_sin_yeast = Im(Xtwo_AorT_fourier_yeast)

colnames(Xtwo_AorT_fourier_cos_yeast) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_yeast) = paste0("fouriersinAorT2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_yeast, "data/Created/yeast_Xtwo_AorT_fourier_cos.rds")
saveRDS(Xtwo_AorT_fourier_sin_yeast, "data/Created/yeast_Xtwo_AorT_fourier_sin.rds")



Xtwo_AA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AA) = colnames(Xtwo)
Xtwo_AA[] = (Xtwo == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier = t(apply(Xtwo_AA, 1, fft))[,1:25]
colnames(AA_fourier) = paste0("fourierAA_", 1:25)
AA_fourier_cos = Re(AA_fourier)
AA_fourier_sin = Im(AA_fourier)
colnames(AA_fourier_cos) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos, "data/Created/yeast_AA_fourier_cos.rds")
saveRDS(AA_fourier_sin, "data/Created/yeast_AA_fourier_sin.rds")

Xtwo_AC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AC) = colnames(Xtwo)
Xtwo_AC[] = (Xtwo == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier = t(apply(Xtwo_AC, 1, fft))[,1:25]
colnames(AC_fourier) = paste0("fourierAC_", 1:25)
AC_fourier_cos = Re(AC_fourier)
AC_fourier_sin = Im(AC_fourier)
colnames(AC_fourier_cos) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos, "data/Created/yeast_AC_fourier_cos.rds")
saveRDS(AC_fourier_sin, "data/Created/yeast_AC_fourier_sin.rds")

Xtwo_AG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AG) = colnames(Xtwo)
Xtwo_AG[] = (Xtwo == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier = t(apply(Xtwo_AG, 1, fft))[,1:25]
colnames(AG_fourier) = paste0("fourierAG_", 1:25)
AG_fourier_cos = Re(AG_fourier)
AG_fourier_sin = Im(AG_fourier)
colnames(AG_fourier_cos) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos, "data/Created/yeast_AG_fourier_cos.rds")
saveRDS(AG_fourier_sin, "data/Created/yeast_AG_fourier_sin.rds")

Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AT) = colnames(Xtwo)
Xtwo_AT[] = (Xtwo == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier = t(apply(Xtwo_AT, 1, fft))[,1:25]
colnames(AT_fourier) = paste0("fourierAT_", 1:25)
AT_fourier_cos = Re(AT_fourier)
AT_fourier_sin = Im(AT_fourier)
colnames(AT_fourier_cos) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos, "data/Created/yeast_AT_fourier_cos.rds")
saveRDS(AT_fourier_sin, "data/Created/yeast_AT_fourier_sin.rds")

Xtwo_CA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CA) = colnames(Xtwo)
Xtwo_CA[] = (Xtwo == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier = t(apply(Xtwo_CA, 1, fft))[,1:25]
colnames(CA_fourier) = paste0("fourierCA_", 1:25)
CA_fourier_cos = Re(CA_fourier)
CA_fourier_sin = Im(CA_fourier)
colnames(CA_fourier_cos) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos, "data/Created/yeast_CA_fourier_cos.rds")
saveRDS(CA_fourier_sin, "data/Created/yeast_CA_fourier_sin.rds")

Xtwo_CC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CC) = colnames(Xtwo)
Xtwo_CC[] = (Xtwo == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier = t(apply(Xtwo_CC, 1, fft))[,1:25]
colnames(CC_fourier) = paste0("fourierCC_", 1:25)
CC_fourier_cos = Re(CC_fourier)
CC_fourier_sin = Im(CC_fourier)
colnames(CC_fourier_cos) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos, "data/Created/yeast_CC_fourier_cos.rds")
saveRDS(CC_fourier_sin, "data/Created/yeast_CC_fourier_sin.rds")

Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CG) = colnames(Xtwo)
Xtwo_CG[] = (Xtwo == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier = t(apply(Xtwo_CG, 1, fft))[,1:25]
colnames(CG_fourier) = paste0("fourierCG_", 1:25)
CG_fourier_cos = Re(CG_fourier)
CG_fourier_sin = Im(CG_fourier)
colnames(CG_fourier_cos) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos, "data/Created/yeast_CG_fourier_cos.rds")
saveRDS(CG_fourier_sin, "data/Created/yeast_CG_fourier_sin.rds")

Xtwo_CT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_CT) = colnames(Xtwo)
Xtwo_CT[] = (Xtwo == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier = t(apply(Xtwo_CT, 1, fft))[,1:25]
colnames(CT_fourier) = paste0("fourierCT_", 1:25)
CT_fourier_cos = Re(CT_fourier)
CT_fourier_sin = Im(CT_fourier)
colnames(CT_fourier_cos) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos, "data/Created/yeast_CT_fourier_cos.rds")
saveRDS(CT_fourier_sin, "data/Created/yeast_CT_fourier_sin.rds")

Xtwo_GA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GA) = colnames(Xtwo)
Xtwo_GA[] = (Xtwo == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier = t(apply(Xtwo_GA, 1, fft))[,1:25]
colnames(GA_fourier) = paste0("fourierGA_", 1:25)
GA_fourier_cos = Re(GA_fourier)
GA_fourier_sin = Im(GA_fourier)
colnames(GA_fourier_cos) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos, "data/Created/yeast_GA_fourier_cos.rds")
saveRDS(GA_fourier_sin, "data/Created/yeast_GA_fourier_sin.rds")

Xtwo_GC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GC) = colnames(Xtwo)
Xtwo_GC[] = (Xtwo == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier = t(apply(Xtwo_GC, 1, fft))[,1:25]
colnames(GC_fourier) = paste0("fourierGC_", 1:25)
GC_fourier_cos = Re(GC_fourier)
GC_fourier_sin = Im(GC_fourier)
colnames(GC_fourier_cos) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos, "data/Created/yeast_GC_fourier_cos.rds")
saveRDS(GC_fourier_sin, "data/Created/yeast_GC_fourier_sin.rds")

Xtwo_GG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GG) = colnames(Xtwo)
Xtwo_GG[] = (Xtwo == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier = t(apply(Xtwo_GG, 1, fft))[,1:25]
colnames(GG_fourier) = paste0("fourierGG_", 1:25)
GG_fourier_cos = Re(GG_fourier)
GG_fourier_sin = Im(GG_fourier)
colnames(GG_fourier_cos) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos, "data/Created/yeast_GG_fourier_cos.rds")
saveRDS(GG_fourier_sin, "data/Created/yeast_GG_fourier_sin.rds")

Xtwo_GT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_GT) = colnames(Xtwo)
Xtwo_GT[] = (Xtwo == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier = t(apply(Xtwo_GT, 1, fft))[,1:25]
colnames(GT_fourier) = paste0("fourierGT_", 1:25)
GT_fourier_cos = Re(GT_fourier)
GT_fourier_sin = Im(GT_fourier)
colnames(GT_fourier_cos) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos, "data/Created/yeast_GT_fourier_cos.rds")
saveRDS(GT_fourier_sin, "data/Created/yeast_GT_fourier_sin.rds")

Xtwo_TA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TA) = colnames(Xtwo)
Xtwo_TA[] = (Xtwo == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier = t(apply(Xtwo_TA, 1, fft))[,1:25]
colnames(TA_fourier) = paste0("fourierTA_", 1:25)
TA_fourier_cos = Re(TA_fourier)
TA_fourier_sin = Im(TA_fourier)
colnames(TA_fourier_cos) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos, "data/Created/yeast_TA_fourier_cos.rds")
saveRDS(TA_fourier_sin, "data/Created/yeast_TA_fourier_sin.rds")

Xtwo_TC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TC) = colnames(Xtwo)
Xtwo_TC[] = (Xtwo == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier = t(apply(Xtwo_TC, 1, fft))[,1:25]
colnames(TC_fourier) = paste0("fourierTC_", 1:25)
TC_fourier_cos = Re(TC_fourier)
TC_fourier_sin = Im(TC_fourier)
colnames(TC_fourier_cos) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos, "data/Created/yeast_TC_fourier_cos.rds")
saveRDS(TC_fourier_sin, "data/Created/yeast_TC_fourier_sin.rds")

Xtwo_TG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TG) = colnames(Xtwo)
Xtwo_TG[] = (Xtwo == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier = t(apply(Xtwo_TG, 1, fft))[,1:25]
colnames(TG_fourier) = paste0("fourierTG_", 1:25)
TG_fourier_cos = Re(TG_fourier)
TG_fourier_sin = Im(TG_fourier)
colnames(TG_fourier_cos) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos, "data/Created/yeast_TG_fourier_cos.rds")
saveRDS(TG_fourier_sin, "data/Created/yeast_TG_fourier_sin.rds")

Xtwo_TT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_TT) = colnames(Xtwo)
Xtwo_TT[] = (Xtwo == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier = t(apply(Xtwo_TT, 1, fft))[,1:25]
colnames(TT_fourier) = paste0("fourierTT_", 1:25)
TT_fourier_cos = Re(TT_fourier)
TT_fourier_sin = Im(TT_fourier)
colnames(TT_fourier_cos) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos, "data/Created/yeast_TT_fourier_cos.rds")
saveRDS(TT_fourier_sin, "data/Created/yeast_TT_fourier_sin.rds")

fourier_di_cos = cbind(AA_fourier_cos, AC_fourier_cos, AG_fourier_cos, AT_fourier_cos, 
                       CA_fourier_cos, CC_fourier_cos, CG_fourier_cos, CT_fourier_cos,
                       GA_fourier_cos, GC_fourier_cos, GG_fourier_cos, GT_fourier_cos,
                       TA_fourier_cos, TC_fourier_cos, TG_fourier_cos, TT_fourier_cos)
fourier_di_sin = cbind(AA_fourier_sin, AC_fourier_sin, AG_fourier_sin, AT_fourier_sin, 
                       CA_fourier_sin, CC_fourier_sin, CG_fourier_sin, CT_fourier_sin,
                       GA_fourier_sin, GC_fourier_sin, GG_fourier_sin, GT_fourier_sin,
                       TA_fourier_sin, TC_fourier_sin, TG_fourier_sin, TT_fourier_sin)
saveRDS(fourier_di_cos, "data/Created/yeast_fourier_di_cos.rds")
saveRDS(fourier_di_sin, "data/Created/yeast_fourier_di_sin.rds")

fourier_AA_AC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AC_fourier[,"fourierAC_6"])
fourier_AA_AG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AG_fourier[,"fourierAG_6"])
fourier_AA_AT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - AT_fourier[,"fourierAT_6"])
fourier_AA_CA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CA_fourier[,"fourierCA_6"])
fourier_AA_CC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CC_fourier[,"fourierCC_6"])
fourier_AA_CG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CG_fourier[,"fourierCG_6"])
fourier_AA_CT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - CT_fourier[,"fourierCT_6"])
fourier_AA_GA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GA_fourier[,"fourierGA_6"])
fourier_AA_GC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GC_fourier[,"fourierGC_6"])
fourier_AA_GG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GG_fourier[,"fourierGG_6"])
fourier_AA_GT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - GT_fourier[,"fourierGT_6"])
fourier_AA_TA_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TA_fourier[,"fourierTA_6"])
fourier_AA_TC_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TC_fourier[,"fourierTC_6"])
fourier_AA_TG_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TG_fourier[,"fourierTG_6"])
fourier_AA_TT_6_mod = Mod(AA_fourier[,"fourierAA_6"] - TT_fourier[,"fourierTT_6"])

fourier_AA_6_mod = paste0("fourier_AA_", dinucleotides[-1], "_6_mod")
fourier_AA_6_mod_yeast = data.frame(map(fourier_AA_6_mod, get))
colnames(fourier_AA_6_mod_yeast) = fourier_AA_6_mod

fourier_AC_AG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AG_fourier[,"fourierAG_6"])
fourier_AC_AT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - AT_fourier[,"fourierAT_6"])
fourier_AC_CA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CA_fourier[,"fourierCA_6"])
fourier_AC_CC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CC_fourier[,"fourierCC_6"])
fourier_AC_CG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CG_fourier[,"fourierCG_6"])
fourier_AC_CT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - CT_fourier[,"fourierCT_6"])
fourier_AC_GA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GA_fourier[,"fourierGA_6"])
fourier_AC_GC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GC_fourier[,"fourierGC_6"])
fourier_AC_GG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GG_fourier[,"fourierGG_6"])
fourier_AC_GT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - GT_fourier[,"fourierGT_6"])
fourier_AC_TA_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TA_fourier[,"fourierTA_6"])
fourier_AC_TC_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TC_fourier[,"fourierTC_6"])
fourier_AC_TG_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TG_fourier[,"fourierTG_6"])
fourier_AC_TT_6_mod = Mod(AC_fourier[,"fourierAC_6"] - TT_fourier[,"fourierTT_6"])

fourier_AC_6_mod = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod")
fourier_AC_6_mod_yeast = data.frame(map(fourier_AC_6_mod, get))
colnames(fourier_AC_6_mod_yeast) = fourier_AC_6_mod

fourier_AG_AT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - AT_fourier[,"fourierAT_6"])
fourier_AG_CA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CA_fourier[,"fourierCA_6"])
fourier_AG_CC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CC_fourier[,"fourierCC_6"])
fourier_AG_CG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CG_fourier[,"fourierCG_6"])
fourier_AG_CT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - CT_fourier[,"fourierCT_6"])
fourier_AG_GA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GA_fourier[,"fourierGA_6"])
fourier_AG_GC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GC_fourier[,"fourierGC_6"])
fourier_AG_GG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GG_fourier[,"fourierGG_6"])
fourier_AG_GT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - GT_fourier[,"fourierGT_6"])
fourier_AG_TA_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TA_fourier[,"fourierTA_6"])
fourier_AG_TC_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TC_fourier[,"fourierTC_6"])
fourier_AG_TG_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TG_fourier[,"fourierTG_6"])
fourier_AG_TT_6_mod = Mod(AG_fourier[,"fourierAG_6"] - TT_fourier[,"fourierTT_6"])

fourier_AG_6_mod = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod")
fourier_AG_6_mod_yeast = data.frame(map(fourier_AG_6_mod, get))
colnames(fourier_AG_6_mod_yeast) = fourier_AG_6_mod

fourier_AT_CA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CA_fourier[,"fourierCA_6"])
fourier_AT_CC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CC_fourier[,"fourierCC_6"])
fourier_AT_CG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CG_fourier[,"fourierCG_6"])
fourier_AT_CT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - CT_fourier[,"fourierCT_6"])
fourier_AT_GA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GA_fourier[,"fourierGA_6"])
fourier_AT_GC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GC_fourier[,"fourierGC_6"])
fourier_AT_GG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GG_fourier[,"fourierGG_6"])
fourier_AT_GT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - GT_fourier[,"fourierGT_6"])
fourier_AT_TA_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TA_fourier[,"fourierTA_6"])
fourier_AT_TC_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TC_fourier[,"fourierTC_6"])
fourier_AT_TG_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TG_fourier[,"fourierTG_6"])
fourier_AT_TT_6_mod = Mod(AT_fourier[,"fourierAT_6"] - TT_fourier[,"fourierTT_6"])

fourier_AT_6_mod = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod")
fourier_AT_6_mod_yeast = data.frame(map(fourier_AT_6_mod, get))
colnames(fourier_AT_6_mod_yeast) = fourier_AT_6_mod

fourier_CA_CC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CC_fourier[,"fourierCC_6"])
fourier_CA_CG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CG_fourier[,"fourierCG_6"])
fourier_CA_CT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - CT_fourier[,"fourierCT_6"])
fourier_CA_GA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GA_fourier[,"fourierGA_6"])
fourier_CA_GC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GC_fourier[,"fourierGC_6"])
fourier_CA_GG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GG_fourier[,"fourierGG_6"])
fourier_CA_GT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - GT_fourier[,"fourierGT_6"])
fourier_CA_TA_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TA_fourier[,"fourierTA_6"])
fourier_CA_TC_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TC_fourier[,"fourierTC_6"])
fourier_CA_TG_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TG_fourier[,"fourierTG_6"])
fourier_CA_TT_6_mod = Mod(CA_fourier[,"fourierCA_6"] - TT_fourier[,"fourierTT_6"])

fourier_CA_6_mod = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod")
fourier_CA_6_mod_yeast = data.frame(map(fourier_CA_6_mod, get))
colnames(fourier_CA_6_mod_yeast) = fourier_CA_6_mod

fourier_CC_CG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CG_fourier[,"fourierCG_6"])
fourier_CC_CT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - CT_fourier[,"fourierCT_6"])
fourier_CC_GA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GA_fourier[,"fourierGA_6"])
fourier_CC_GC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GC_fourier[,"fourierGC_6"])
fourier_CC_GG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GG_fourier[,"fourierGG_6"])
fourier_CC_GT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - GT_fourier[,"fourierGT_6"])
fourier_CC_TA_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TA_fourier[,"fourierTA_6"])
fourier_CC_TC_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TC_fourier[,"fourierTC_6"])
fourier_CC_TG_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TG_fourier[,"fourierTG_6"])
fourier_CC_TT_6_mod = Mod(CC_fourier[,"fourierCC_6"] - TT_fourier[,"fourierTT_6"])

fourier_CC_6_mod = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod")
fourier_CC_6_mod_yeast = data.frame(map(fourier_CC_6_mod, get))
colnames(fourier_CC_6_mod_yeast) = fourier_CC_6_mod

fourier_CG_CT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - CT_fourier[,"fourierCT_6"])
fourier_CG_GA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GA_fourier[,"fourierGA_6"])
fourier_CG_GC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GC_fourier[,"fourierGC_6"])
fourier_CG_GG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GG_fourier[,"fourierGG_6"])
fourier_CG_GT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - GT_fourier[,"fourierGT_6"])
fourier_CG_TA_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TA_fourier[,"fourierTA_6"])
fourier_CG_TC_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TC_fourier[,"fourierTC_6"])
fourier_CG_TG_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TG_fourier[,"fourierTG_6"])
fourier_CG_TT_6_mod = Mod(CG_fourier[,"fourierCG_6"] - TT_fourier[,"fourierTT_6"])

fourier_CG_6_mod = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod")
fourier_CG_6_mod_yeast = data.frame(map(fourier_CG_6_mod, get))
colnames(fourier_CG_6_mod_yeast) = fourier_CG_6_mod

fourier_CT_GA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GA_fourier[,"fourierGA_6"])
fourier_CT_GC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GC_fourier[,"fourierGC_6"])
fourier_CT_GG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GG_fourier[,"fourierGG_6"])
fourier_CT_GT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - GT_fourier[,"fourierGT_6"])
fourier_CT_TA_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TA_fourier[,"fourierTA_6"])
fourier_CT_TC_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TC_fourier[,"fourierTC_6"])
fourier_CT_TG_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TG_fourier[,"fourierTG_6"])
fourier_CT_TT_6_mod = Mod(CT_fourier[,"fourierCT_6"] - TT_fourier[,"fourierTT_6"])

fourier_CT_6_mod = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod")
fourier_CT_6_mod_yeast = data.frame(map(fourier_CT_6_mod, get))
colnames(fourier_CT_6_mod_yeast) = fourier_CT_6_mod

fourier_GA_GC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GC_fourier[,"fourierGC_6"])
fourier_GA_GG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GG_fourier[,"fourierGG_6"])
fourier_GA_GT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - GT_fourier[,"fourierGT_6"])
fourier_GA_TA_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TA_fourier[,"fourierTA_6"])
fourier_GA_TC_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TC_fourier[,"fourierTC_6"])
fourier_GA_TG_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TG_fourier[,"fourierTG_6"])
fourier_GA_TT_6_mod = Mod(GA_fourier[,"fourierGA_6"] - TT_fourier[,"fourierTT_6"])

fourier_GA_6_mod = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod")
fourier_GA_6_mod_yeast = data.frame(map(fourier_GA_6_mod, get))
colnames(fourier_GA_6_mod_yeast) = fourier_GA_6_mod

fourier_GC_GG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GG_fourier[,"fourierGG_6"])
fourier_GC_GT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - GT_fourier[,"fourierGT_6"])
fourier_GC_TA_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TA_fourier[,"fourierTA_6"])
fourier_GC_TC_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TC_fourier[,"fourierTC_6"])
fourier_GC_TG_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TG_fourier[,"fourierTG_6"])
fourier_GC_TT_6_mod = Mod(GC_fourier[,"fourierGC_6"] - TT_fourier[,"fourierTT_6"])

fourier_GC_6_mod = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod")
fourier_GC_6_mod_yeast = data.frame(map(fourier_GC_6_mod, get))
colnames(fourier_GC_6_mod_yeast) = fourier_GC_6_mod

fourier_GG_GT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - GT_fourier[,"fourierGT_6"])
fourier_GG_TA_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TA_fourier[,"fourierTA_6"])
fourier_GG_TC_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TC_fourier[,"fourierTC_6"])
fourier_GG_TG_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TG_fourier[,"fourierTG_6"])
fourier_GG_TT_6_mod = Mod(GG_fourier[,"fourierGG_6"] - TT_fourier[,"fourierTT_6"])

fourier_GG_6_mod = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod")
fourier_GG_6_mod_yeast = data.frame(map(fourier_GG_6_mod, get))
colnames(fourier_GG_6_mod_yeast) = fourier_GG_6_mod

fourier_GT_TA_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TA_fourier[,"fourierTA_6"])
fourier_GT_TC_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TC_fourier[,"fourierTC_6"])
fourier_GT_TG_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TG_fourier[,"fourierTG_6"])
fourier_GT_TT_6_mod = Mod(GT_fourier[,"fourierGT_6"] - TT_fourier[,"fourierTT_6"])

fourier_GT_6_mod = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod")
fourier_GT_6_mod_yeast = data.frame(map(fourier_GT_6_mod, get))
colnames(fourier_GT_6_mod_yeast) = fourier_GT_6_mod

fourier_TA_TC_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TC_fourier[,"fourierTC_6"])
fourier_TA_TG_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TG_fourier[,"fourierTG_6"])
fourier_TA_TT_6_mod = Mod(TA_fourier[,"fourierTA_6"] - TT_fourier[,"fourierTT_6"])

fourier_TA_6_mod = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod")
fourier_TA_6_mod_yeast = data.frame(map(fourier_TA_6_mod, get))
colnames(fourier_TA_6_mod_yeast) = fourier_TA_6_mod

fourier_TC_TG_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TG_fourier[,"fourierTG_6"])
fourier_TC_TT_6_mod = Mod(TC_fourier[,"fourierTC_6"] - TT_fourier[,"fourierTT_6"])

fourier_TC_6_mod = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod")
fourier_TC_6_mod_yeast = data.frame(map(fourier_TC_6_mod, get))
colnames(fourier_TC_6_mod_yeast) = fourier_TC_6_mod

fourier_TG_TT_6_mod = Mod(TG_fourier[,"fourierTG_6"] - TT_fourier[,"fourierTT_6"])

fourier_TG_6_mod = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod")
fourier_TG_6_mod_yeast = data.frame(map(fourier_TG_6_mod, get))
colnames(fourier_TG_6_mod_yeast) = fourier_TG_6_mod

fourier_6_mod = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_yeast"), get))
saveRDS(fourier_6_mod, "data/Created/yeast_fourier_6_mod.rds")




# Test data:
Xone_test = dat_yeast_test %>% select(all_of(ps1))
Xone_AorT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CorG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AorT_test) = colnames(Xone_test)
colnames(Xone_CorG_test) = colnames(Xone_test)
Xone_AorT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CorG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,c(1:25, 50)]
Xone_AorT_fourier_cos_test = Re(Xone_AorT_fourier_test)
Xone_AorT_fourier_sin_test = Im(Xone_AorT_fourier_test)

Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,c(1:25, 50)]
Xone_CorG_fourier_cos_test = Re(Xone_CorG_fourier_test)
Xone_CorG_fourier_sin_test = Im(Xone_CorG_fourier_test)

colnames(Xone_AorT_fourier_cos_test) = paste0("fouriercosAorT", c(1:25, 50))
colnames(Xone_AorT_fourier_sin_test) = paste0("fouriersinAorT", c(1:25, 50))

colnames(Xone_CorG_fourier_cos_test) = paste0("fouriercosCorG", c(1:25, 50))
colnames(Xone_CorG_fourier_sin_test) = paste0("fouriersinCorG", c(1:25, 50))

saveRDS(Xone_AorT_fourier_cos_test, "data/Created/yeast_Xone_AorT_fourier_cos_test.rds")
saveRDS(Xone_AorT_fourier_sin_test, "data/Created/yeast_Xone_AorT_fourier_sin_test.rds")

saveRDS(Xone_CorG_fourier_cos_test, "data/Created/yeast_Xone_CorG_fourier_cos_test.rds")
saveRDS(Xone_CorG_fourier_sin_test, "data/Created/yeast_Xone_CorG_fourier_sin_test.rds")


Xtwo_test = dat_yeast_test %>% select(all_of(ps2))
Xtwo_AorT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
# Xtwo_CorG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AorT_test) = colnames(Xtwo_test)
# colnames(Xtwo_CorG_test) = colnames(Xtwo_test)
Xtwo_AorT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT") | (Xtwo_test == "AT") | (Xtwo_test == "TA")) %>% 
  as.matrix() %>% as.numeric()
# Xtwo_CorG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG") | (Xtwo_test == "CG") | (Xtwo_test == "GC")) %>%
#   as.matrix() %>% as.numeric()

Xtwo_AorT_fourier_test = t(apply(Xtwo_AorT_test, 1, fft))[,1:25]
Xtwo_AorT_fourier_cos_test = Re(Xtwo_AorT_fourier_test)
Xtwo_AorT_fourier_sin_test = Im(Xtwo_AorT_fourier_test)

# Xtwo_CorG_fourier = t(apply(Xtwo_CorG_test, 1, fft))[,c(1:25, 50)]
# Xtwo_CorG_fourier_cos_test = Re(Xtwo_CorG_fourier)
# Xtwo_CorG_fourier_sin_test = Im(Xtwo_CorG_fourier)

colnames(Xtwo_AorT_fourier_cos_test) = paste0("fouriercosAorT2_", 1:25)
colnames(Xtwo_AorT_fourier_sin_test) = paste0("fouriersinAorT2_", 1:25)

# colnames(Xtwo_CorG_fourier_cos_test) = paste0("fouriercosCorG2_", 1:25)
# colnames(Xtwo_CorG_fourier_sin_test) = paste0("fouriersinCorG2_", 1:25)

saveRDS(Xtwo_AorT_fourier_cos_test, "data/Created/yeast_Xtwo_AorT_fourier_cos_test.rds")
saveRDS(Xtwo_AorT_fourier_sin_test, "data/Created/yeast_Xtwo_AorT_fourier_sin_test.rds")

# saveRDS(Xtwo_CorG_fourier_cos_test, "data/Created/yeast_Xtwo_CorG_fourier_cos_test.rds")
# saveRDS(Xtwo_CorG_fourier_sin_test, "data/Created/yeast_Xtwo_CorG_fourier_sin_test.rds")


Xtwo_AA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AA_test) = colnames(Xtwo_test)
Xtwo_AA_test[] = (Xtwo_test == "AA") %>% as.matrix() %>% as.numeric()
AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,1:25]
colnames(AA_fourier_test) = paste0("fourierAA_", 1:25)
AA_fourier_cos_test = Re(AA_fourier_test)
AA_fourier_sin_test = Im(AA_fourier_test)
colnames(AA_fourier_cos_test) = paste0("fouriercosAA_", 1:25)
colnames(AA_fourier_sin_test) = paste0("fouriersinAA_", 1:25)
saveRDS(AA_fourier_cos_test, "data/Created/yeast_AA_fourier_cos_test.rds")
saveRDS(AA_fourier_sin_test, "data/Created/yeast_AA_fourier_sin_test.rds")

Xtwo_AC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AC_test) = colnames(Xtwo_test)
Xtwo_AC_test[] = (Xtwo_test == "AC") %>% as.matrix() %>% as.numeric()
AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,1:25]
colnames(AC_fourier_test) = paste0("fourierAC_", 1:25)
AC_fourier_cos_test = Re(AC_fourier_test)
AC_fourier_sin_test = Im(AC_fourier_test)
colnames(AC_fourier_cos_test) = paste0("fouriercosAC_", 1:25)
colnames(AC_fourier_sin_test) = paste0("fouriersinAC_", 1:25)
saveRDS(AC_fourier_cos_test, "data/Created/yeast_AC_fourier_cos_test.rds")
saveRDS(AC_fourier_sin_test, "data/Created/yeast_AC_fourier_sin_test.rds")

Xtwo_AG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AG_test) = colnames(Xtwo_test)
Xtwo_AG_test[] = (Xtwo_test == "AG") %>% as.matrix() %>% as.numeric()
AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,1:25]
colnames(AG_fourier_test) = paste0("fourierAG_", 1:25)
AG_fourier_cos_test = Re(AG_fourier_test)
AG_fourier_sin_test = Im(AG_fourier_test)
colnames(AG_fourier_cos_test) = paste0("fouriercosAG_", 1:25)
colnames(AG_fourier_sin_test) = paste0("fouriersinAG_", 1:25)
saveRDS(AG_fourier_cos_test, "data/Created/yeast_AG_fourier_cos_test.rds")
saveRDS(AG_fourier_sin_test, "data/Created/yeast_AG_fourier_sin_test.rds")

Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
Xtwo_AT_test[] = (Xtwo_test == "AT") %>% as.matrix() %>% as.numeric()
AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,1:25]
colnames(AT_fourier_test) = paste0("fourierAT_", 1:25)
AT_fourier_cos_test = Re(AT_fourier_test)
AT_fourier_sin_test = Im(AT_fourier_test)
colnames(AT_fourier_cos_test) = paste0("fouriercosAT_", 1:25)
colnames(AT_fourier_sin_test) = paste0("fouriersinAT_", 1:25)
saveRDS(AT_fourier_cos_test, "data/Created/yeast_AT_fourier_cos_test.rds")
saveRDS(AT_fourier_sin_test, "data/Created/yeast_AT_fourier_sin_test.rds")

Xtwo_CA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CA_test) = colnames(Xtwo_test)
Xtwo_CA_test[] = (Xtwo_test == "CA") %>% as.matrix() %>% as.numeric()
CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,1:25]
colnames(CA_fourier_test) = paste0("fourierCA_", 1:25)
CA_fourier_cos_test = Re(CA_fourier_test)
CA_fourier_sin_test = Im(CA_fourier_test)
colnames(CA_fourier_cos_test) = paste0("fouriercosCA_", 1:25)
colnames(CA_fourier_sin_test) = paste0("fouriersinCA_", 1:25)
saveRDS(CA_fourier_cos_test, "data/Created/yeast_CA_fourier_cos_test.rds")
saveRDS(CA_fourier_sin_test, "data/Created/yeast_CA_fourier_sin_test.rds")

Xtwo_CC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CC_test) = colnames(Xtwo_test)
Xtwo_CC_test[] = (Xtwo_test == "CC") %>% as.matrix() %>% as.numeric()
CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,1:25]
colnames(CC_fourier_test) = paste0("fourierCC_", 1:25)
CC_fourier_cos_test = Re(CC_fourier_test)
CC_fourier_sin_test = Im(CC_fourier_test)
colnames(CC_fourier_cos_test) = paste0("fouriercosCC_", 1:25)
colnames(CC_fourier_sin_test) = paste0("fouriersinCC_", 1:25)
saveRDS(CC_fourier_cos_test, "data/Created/yeast_CC_fourier_cos_test.rds")
saveRDS(CC_fourier_sin_test, "data/Created/yeast_CC_fourier_sin_test.rds")

Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
Xtwo_CG_test[] = (Xtwo_test == "CG") %>% as.matrix() %>% as.numeric()
CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,1:25]
colnames(CG_fourier_test) = paste0("fourierCG_", 1:25)
CG_fourier_cos_test = Re(CG_fourier_test)
CG_fourier_sin_test = Im(CG_fourier_test)
colnames(CG_fourier_cos_test) = paste0("fouriercosCG_", 1:25)
colnames(CG_fourier_sin_test) = paste0("fouriersinCG_", 1:25)
saveRDS(CG_fourier_cos_test, "data/Created/yeast_CG_fourier_cos_test.rds")
saveRDS(CG_fourier_sin_test, "data/Created/yeast_CG_fourier_sin_test.rds")

Xtwo_CT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_CT_test) = colnames(Xtwo_test)
Xtwo_CT_test[] = (Xtwo_test == "CT") %>% as.matrix() %>% as.numeric()
CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,1:25]
colnames(CT_fourier_test) = paste0("fourierCT_", 1:25)
CT_fourier_cos_test = Re(CT_fourier_test)
CT_fourier_sin_test = Im(CT_fourier_test)
colnames(CT_fourier_cos_test) = paste0("fouriercosCT_", 1:25)
colnames(CT_fourier_sin_test) = paste0("fouriersinCT_", 1:25)
saveRDS(CT_fourier_cos_test, "data/Created/yeast_CT_fourier_cos_test.rds")
saveRDS(CT_fourier_sin_test, "data/Created/yeast_CT_fourier_sin_test.rds")

Xtwo_GA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GA_test) = colnames(Xtwo_test)
Xtwo_GA_test[] = (Xtwo_test == "GA") %>% as.matrix() %>% as.numeric()
GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,1:25]
colnames(GA_fourier_test) = paste0("fourierGA_", 1:25)
GA_fourier_cos_test = Re(GA_fourier_test)
GA_fourier_sin_test = Im(GA_fourier_test)
colnames(GA_fourier_cos_test) = paste0("fouriercosGA_", 1:25)
colnames(GA_fourier_sin_test) = paste0("fouriersinGA_", 1:25)
saveRDS(GA_fourier_cos_test, "data/Created/yeast_GA_fourier_cos_test.rds")
saveRDS(GA_fourier_sin_test, "data/Created/yeast_GA_fourier_sin_test.rds")

Xtwo_GC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GC_test) = colnames(Xtwo_test)
Xtwo_GC_test[] = (Xtwo_test == "GC") %>% as.matrix() %>% as.numeric()
GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,1:25]
colnames(GC_fourier_test) = paste0("fourierGC_", 1:25)
GC_fourier_cos_test = Re(GC_fourier_test)
GC_fourier_sin_test = Im(GC_fourier_test)
colnames(GC_fourier_cos_test) = paste0("fouriercosGC_", 1:25)
colnames(GC_fourier_sin_test) = paste0("fouriersinGC_", 1:25)
saveRDS(GC_fourier_cos_test, "data/Created/yeast_GC_fourier_cos_test.rds")
saveRDS(GC_fourier_sin_test, "data/Created/yeast_GC_fourier_sin_test.rds")

Xtwo_GG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GG_test) = colnames(Xtwo_test)
Xtwo_GG_test[] = (Xtwo_test == "GG") %>% as.matrix() %>% as.numeric()
GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,1:25]
colnames(GG_fourier_test) = paste0("fourierGG_", 1:25)
GG_fourier_cos_test = Re(GG_fourier_test)
GG_fourier_sin_test = Im(GG_fourier_test)
colnames(GG_fourier_cos_test) = paste0("fouriercosGG_", 1:25)
colnames(GG_fourier_sin_test) = paste0("fouriersinGG_", 1:25)
saveRDS(GG_fourier_cos_test, "data/Created/yeast_GG_fourier_cos_test.rds")
saveRDS(GG_fourier_sin_test, "data/Created/yeast_GG_fourier_sin_test.rds")

Xtwo_GT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_GT_test) = colnames(Xtwo_test)
Xtwo_GT_test[] = (Xtwo_test == "GT") %>% as.matrix() %>% as.numeric()
GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,1:25]
colnames(GT_fourier_test) = paste0("fourierGT_", 1:25)
GT_fourier_cos_test = Re(GT_fourier_test)
GT_fourier_sin_test = Im(GT_fourier_test)
colnames(GT_fourier_cos_test) = paste0("fouriercosGT_", 1:25)
colnames(GT_fourier_sin_test) = paste0("fouriersinGT_", 1:25)
saveRDS(GT_fourier_cos_test, "data/Created/yeast_GT_fourier_cos_test.rds")
saveRDS(GT_fourier_sin_test, "data/Created/yeast_GT_fourier_sin_test.rds")

Xtwo_TA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TA_test) = colnames(Xtwo_test)
Xtwo_TA_test[] = (Xtwo_test == "TA") %>% as.matrix() %>% as.numeric()
TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,1:25]
colnames(TA_fourier_test) = paste0("fourierTA_", 1:25)
TA_fourier_cos_test = Re(TA_fourier_test)
TA_fourier_sin_test = Im(TA_fourier_test)
colnames(TA_fourier_cos_test) = paste0("fouriercosTA_", 1:25)
colnames(TA_fourier_sin_test) = paste0("fouriersinTA_", 1:25)
saveRDS(TA_fourier_cos_test, "data/Created/yeast_TA_fourier_cos_test.rds")
saveRDS(TA_fourier_sin_test, "data/Created/yeast_TA_fourier_sin_test.rds")

Xtwo_TC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TC_test) = colnames(Xtwo_test)
Xtwo_TC_test[] = (Xtwo_test == "TC") %>% as.matrix() %>% as.numeric()
TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,1:25]
colnames(TC_fourier_test) = paste0("fourierTC_", 1:25)
TC_fourier_cos_test = Re(TC_fourier_test)
TC_fourier_sin_test = Im(TC_fourier_test)
colnames(TC_fourier_cos_test) = paste0("fouriercosTC_", 1:25)
colnames(TC_fourier_sin_test) = paste0("fouriersinTC_", 1:25)
saveRDS(TC_fourier_cos_test, "data/Created/yeast_TC_fourier_cos_test.rds")
saveRDS(TC_fourier_sin_test, "data/Created/yeast_TC_fourier_sin_test.rds")

Xtwo_TG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TG_test) = colnames(Xtwo_test)
Xtwo_TG_test[] = (Xtwo_test == "TG") %>% as.matrix() %>% as.numeric()
TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,1:25]
colnames(TG_fourier_test) = paste0("fourierTG_", 1:25)
TG_fourier_cos_test = Re(TG_fourier_test)
TG_fourier_sin_test = Im(TG_fourier_test)
colnames(TG_fourier_cos_test) = paste0("fouriercosTG_", 1:25)
colnames(TG_fourier_sin_test) = paste0("fouriersinTG_", 1:25)
saveRDS(TG_fourier_cos_test, "data/Created/yeast_TG_fourier_cos_test.rds")
saveRDS(TG_fourier_sin_test, "data/Created/yeast_TG_fourier_sin_test.rds")

Xtwo_TT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_TT_test) = colnames(Xtwo_test)
Xtwo_TT_test[] = (Xtwo_test == "TT") %>% as.matrix() %>% as.numeric()
TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,1:25]
colnames(TT_fourier_test) = paste0("fourierTT_", 1:25)
TT_fourier_cos_test = Re(TT_fourier_test)
TT_fourier_sin_test = Im(TT_fourier_test)
colnames(TT_fourier_cos_test) = paste0("fouriercosTT_", 1:25)
colnames(TT_fourier_sin_test) = paste0("fouriersinTT_", 1:25)
saveRDS(TT_fourier_cos_test, "data/Created/yeast_TT_fourier_cos_test.rds")
saveRDS(TT_fourier_sin_test, "data/Created/yeast_TT_fourier_sin_test.rds")

fourier_di_cos_test = cbind(AA_fourier_cos_test, AC_fourier_cos_test, AG_fourier_cos_test, AT_fourier_cos_test, 
                            CA_fourier_cos_test, CC_fourier_cos_test, CG_fourier_cos_test, CT_fourier_cos_test,
                            GA_fourier_cos_test, GC_fourier_cos_test, GG_fourier_cos_test, GT_fourier_cos_test,
                            TA_fourier_cos_test, TC_fourier_cos_test, TG_fourier_cos_test, TT_fourier_cos_test)
fourier_di_sin_test = cbind(AA_fourier_sin_test, AC_fourier_sin_test, AG_fourier_sin_test, AT_fourier_sin_test, 
                            CA_fourier_sin_test, CC_fourier_sin_test, CG_fourier_sin_test, CT_fourier_sin_test,
                            GA_fourier_sin_test, GC_fourier_sin_test, GG_fourier_sin_test, GT_fourier_sin_test,
                            TA_fourier_sin_test, TC_fourier_sin_test, TG_fourier_sin_test, TT_fourier_sin_test)
saveRDS(fourier_di_cos_test, "data/Created/yeast_fourier_di_cos_test.rds")
saveRDS(fourier_di_sin_test, "data/Created/yeast_fourier_di_sin_test.rds")



fourier_AA_AC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AC_fourier_test[,"fourierAC_6"])
fourier_AA_AG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AA_AT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AA_CA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AA_CC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AA_CG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AA_CT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AA_GA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AA_GC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AA_GG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AA_GT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AA_TA_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AA_TC_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AA_TG_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AA_TT_6_mod_test = Mod(AA_fourier_test[,"fourierAA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AA_6_mod_test = paste0("fourier_AA_", dinucleotides[-1], "_6_mod_test")
fourier_AA_6_mod_yeast_test = data.frame(map(fourier_AA_6_mod_test, get))
colnames(fourier_AA_6_mod_yeast_test) = fourier_AA_6_mod

fourier_AC_AG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AG_fourier_test[,"fourierAG_6"])
fourier_AC_AT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AC_CA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AC_CC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AC_CG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AC_CT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AC_GA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AC_GC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AC_GG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AC_GT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AC_TA_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AC_TC_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AC_TG_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AC_TT_6_mod_test = Mod(AC_fourier_test[,"fourierAC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AC_6_mod_test = paste0("fourier_AC_", dinucleotides[-(1:2)], "_6_mod_test")
fourier_AC_6_mod_yeast_test = data.frame(map(fourier_AC_6_mod_test, get))
colnames(fourier_AC_6_mod_yeast_test) = fourier_AC_6_mod

fourier_AG_AT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - AT_fourier_test[,"fourierAT_6"])
fourier_AG_CA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AG_CC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AG_CG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AG_CT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AG_GA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AG_GC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AG_GG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AG_GT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AG_TA_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AG_TC_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AG_TG_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AG_TT_6_mod_test = Mod(AG_fourier_test[,"fourierAG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AG_6_mod_test = paste0("fourier_AG_", dinucleotides[-(1:3)], "_6_mod_test")
fourier_AG_6_mod_yeast_test = data.frame(map(fourier_AG_6_mod_test, get))
colnames(fourier_AG_6_mod_yeast_test) = fourier_AG_6_mod

fourier_AT_CA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CA_fourier_test[,"fourierCA_6"])
fourier_AT_CC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_AT_CG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_AT_CT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_AT_GA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_AT_GC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_AT_GG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_AT_GT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_AT_TA_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_AT_TC_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_AT_TG_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_AT_TT_6_mod_test = Mod(AT_fourier_test[,"fourierAT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_AT_6_mod_test = paste0("fourier_AT_", dinucleotides[-(1:4)], "_6_mod_test")
fourier_AT_6_mod_yeast_test = data.frame(map(fourier_AT_6_mod_test, get))
colnames(fourier_AT_6_mod_yeast_test) = fourier_AT_6_mod

fourier_CA_CC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CC_fourier_test[,"fourierCC_6"])
fourier_CA_CG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CA_CT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CA_GA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CA_GC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CA_GG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CA_GT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CA_TA_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CA_TC_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CA_TG_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CA_TT_6_mod_test = Mod(CA_fourier_test[,"fourierCA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CA_6_mod_test = paste0("fourier_CA_", dinucleotides[-(1:5)], "_6_mod_test")
fourier_CA_6_mod_yeast_test = data.frame(map(fourier_CA_6_mod_test, get))
colnames(fourier_CA_6_mod_yeast_test) = fourier_CA_6_mod

fourier_CC_CG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CG_fourier_test[,"fourierCG_6"])
fourier_CC_CT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CC_GA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CC_GC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CC_GG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CC_GT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CC_TA_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CC_TC_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CC_TG_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CC_TT_6_mod_test = Mod(CC_fourier_test[,"fourierCC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CC_6_mod_test = paste0("fourier_CC_", dinucleotides[-(1:6)], "_6_mod_test")
fourier_CC_6_mod_yeast_test = data.frame(map(fourier_CC_6_mod_test, get))
colnames(fourier_CC_6_mod_yeast_test) = fourier_CC_6_mod

fourier_CG_CT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - CT_fourier_test[,"fourierCT_6"])
fourier_CG_GA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CG_GC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CG_GG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CG_GT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CG_TA_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CG_TC_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CG_TG_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CG_TT_6_mod_test = Mod(CG_fourier_test[,"fourierCG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CG_6_mod_test = paste0("fourier_CG_", dinucleotides[-(1:7)], "_6_mod_test")
fourier_CG_6_mod_yeast_test = data.frame(map(fourier_CG_6_mod_test, get))
colnames(fourier_CG_6_mod_yeast_test) = fourier_CG_6_mod

fourier_CT_GA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GA_fourier_test[,"fourierGA_6"])
fourier_CT_GC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_CT_GG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_CT_GT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_CT_TA_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_CT_TC_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_CT_TG_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_CT_TT_6_mod_test = Mod(CT_fourier_test[,"fourierCT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_CT_6_mod_test = paste0("fourier_CT_", dinucleotides[-(1:8)], "_6_mod_test")
fourier_CT_6_mod_yeast_test = data.frame(map(fourier_CT_6_mod_test, get))
colnames(fourier_CT_6_mod_yeast_test) = fourier_CT_6_mod

fourier_GA_GC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GC_fourier_test[,"fourierGC_6"])
fourier_GA_GG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GA_GT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GA_TA_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GA_TC_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GA_TG_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GA_TT_6_mod_test = Mod(GA_fourier_test[,"fourierGA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GA_6_mod_test = paste0("fourier_GA_", dinucleotides[-(1:9)], "_6_mod_test")
fourier_GA_6_mod_yeast_test = data.frame(map(fourier_GA_6_mod_test, get))
colnames(fourier_GA_6_mod_yeast_test) = fourier_GA_6_mod

fourier_GC_GG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GG_fourier_test[,"fourierGG_6"])
fourier_GC_GT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GC_TA_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GC_TC_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GC_TG_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GC_TT_6_mod_test = Mod(GC_fourier_test[,"fourierGC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GC_6_mod_test = paste0("fourier_GC_", dinucleotides[-(1:10)], "_6_mod_test")
fourier_GC_6_mod_yeast_test = data.frame(map(fourier_GC_6_mod_test, get))
colnames(fourier_GC_6_mod_yeast_test) = fourier_GC_6_mod

fourier_GG_GT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - GT_fourier_test[,"fourierGT_6"])
fourier_GG_TA_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GG_TC_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GG_TG_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GG_TT_6_mod_test = Mod(GG_fourier_test[,"fourierGG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GG_6_mod_test = paste0("fourier_GG_", dinucleotides[-(1:11)], "_6_mod_test")
fourier_GG_6_mod_yeast_test = data.frame(map(fourier_GG_6_mod_test, get))
colnames(fourier_GG_6_mod_yeast_test) = fourier_GG_6_mod

fourier_GT_TA_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TA_fourier_test[,"fourierTA_6"])
fourier_GT_TC_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_GT_TG_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_GT_TT_6_mod_test = Mod(GT_fourier_test[,"fourierGT_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_GT_6_mod_test = paste0("fourier_GT_", dinucleotides[-(1:12)], "_6_mod_test")
fourier_GT_6_mod_yeast_test = data.frame(map(fourier_GT_6_mod_test, get))
colnames(fourier_GT_6_mod_yeast_test) = fourier_GT_6_mod

fourier_TA_TC_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TC_fourier_test[,"fourierTC_6"])
fourier_TA_TG_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TA_TT_6_mod_test = Mod(TA_fourier_test[,"fourierTA_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TA_6_mod_test = paste0("fourier_TA_", dinucleotides[-(1:13)], "_6_mod_test")
fourier_TA_6_mod_yeast_test = data.frame(map(fourier_TA_6_mod_test, get))
colnames(fourier_TA_6_mod_yeast_test) = fourier_TA_6_mod

fourier_TC_TG_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TG_fourier_test[,"fourierTG_6"])
fourier_TC_TT_6_mod_test = Mod(TC_fourier_test[,"fourierTC_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TC_6_mod_test = paste0("fourier_TC_", dinucleotides[-(1:14)], "_6_mod_test")
fourier_TC_6_mod_yeast_test = data.frame(map(fourier_TC_6_mod_test, get))
colnames(fourier_TC_6_mod_yeast_test) = fourier_TC_6_mod

fourier_TG_TT_6_mod_test = Mod(TG_fourier_test[,"fourierTG_6"] - TT_fourier_test[,"fourierTT_6"])

fourier_TG_6_mod_test = paste0("fourier_TG_", dinucleotides[-(1:15)], "_6_mod_test")
fourier_TG_6_mod_yeast_test = data.frame(map(fourier_TG_6_mod_test, get))
colnames(fourier_TG_6_mod_yeast_test) = fourier_TG_6_mod

fourier_6_mod_test = data.frame(map(paste0("fourier_", dinucleotides[-16], "_6_mod_yeast_test"), get))
saveRDS(fourier_6_mod_test, "data/Created/yeast_fourier_6_mod_test.rds")





