dat = readRDS("data/Created/processed_tiling_newC0.rds")
dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")

########################################################################
# Tiling
########################################################################

## Train data:

# Nucleotides:

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

# Xone_A_fourier = t(apply(Xone_A, 1, fft))[,1:26]
# Xone_C_fourier = t(apply(Xone_C, 1, fft))[,1:26]
# Xone_G_fourier = t(apply(Xone_G, 1, fft))[,1:26]
# Xone_T_fourier = t(apply(Xone_T, 1, fft))[,1:26]

Xone_A_fourier = t(apply(Xone_A, 1, fft))[,6]
Xone_C_fourier = t(apply(Xone_C, 1, fft))[,6]
Xone_G_fourier = t(apply(Xone_G, 1, fft))[,6]
Xone_T_fourier = t(apply(Xone_T, 1, fft))[,6]

Xone_A_fourier_mag = as.matrix(Mod(Xone_A_fourier))
Xone_C_fourier_mag = as.matrix(Mod(Xone_C_fourier))
Xone_G_fourier_mag = as.matrix(Mod(Xone_G_fourier))
Xone_T_fourier_mag = as.matrix(Mod(Xone_T_fourier))

Xone_A_fourier_phase = as.matrix(Arg(Xone_A_fourier))
Xone_C_fourier_phase = as.matrix(Arg(Xone_C_fourier))
Xone_G_fourier_phase = as.matrix(Arg(Xone_G_fourier))
Xone_T_fourier_phase = as.matrix(Arg(Xone_T_fourier))

# colnames(Xone_A_fourier_mag) = paste0("fourier_mag_A_", 1:26)
# colnames(Xone_C_fourier_mag) = paste0("fourier_mag_C_", 1:26)
# colnames(Xone_G_fourier_mag) = paste0("fourier_mag_G_", 1:26)
# colnames(Xone_T_fourier_mag) = paste0("fourier_mag_T_", 1:26)

# colnames(Xone_A_fourier_phase) = paste0("fourier_phase_A_", 1:26)
# colnames(Xone_C_fourier_phase) = paste0("fourier_phase_C_", 1:26)
# colnames(Xone_G_fourier_phase) = paste0("fourier_phase_G_", 1:26)
# colnames(Xone_T_fourier_phase) = paste0("fourier_phase_T_", 1:26)

colnames(Xone_A_fourier_mag) = "fourier_mag_A_6"
colnames(Xone_C_fourier_mag) = "fourier_mag_C_6"
colnames(Xone_G_fourier_mag) = "fourier_mag_G_6"
colnames(Xone_T_fourier_mag) = "fourier_mag_T_6"

colnames(Xone_A_fourier_phase) = "fourier_phase_A_6"
colnames(Xone_C_fourier_phase) = "fourier_phase_C_6"
colnames(Xone_G_fourier_phase) = "fourier_phase_G_6"
colnames(Xone_T_fourier_phase) = "fourier_phase_T_6"



# Pairs of nucleotides:

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

# Xone_AorC_fourier = t(apply(Xone_AorC, 1, fft))[,1:26]
# Xone_AorG_fourier = t(apply(Xone_AorG, 1, fft))[,1:26]
# Xone_AorT_fourier = t(apply(Xone_AorT, 1, fft))[,1:26]
# Xone_CorG_fourier = t(apply(Xone_CorG, 1, fft))[,1:26]
# Xone_CorT_fourier = t(apply(Xone_CorT, 1, fft))[,1:26]
# Xone_GorT_fourier = t(apply(Xone_GorT, 1, fft))[,1:26]

Xone_AorC_fourier = t(apply(Xone_AorC, 1, fft))[,6]
Xone_AorG_fourier = t(apply(Xone_AorG, 1, fft))[,6]
Xone_AorT_fourier = t(apply(Xone_AorT, 1, fft))[,6]
Xone_CorG_fourier = t(apply(Xone_CorG, 1, fft))[,6]
Xone_CorT_fourier = t(apply(Xone_CorT, 1, fft))[,6]
Xone_GorT_fourier = t(apply(Xone_GorT, 1, fft))[,6]

Xone_AorC_fourier_mag = as.matrix(Mod(Xone_AorC_fourier))
Xone_AorG_fourier_mag = as.matrix(Mod(Xone_AorG_fourier))
Xone_AorT_fourier_mag = as.matrix(Mod(Xone_AorT_fourier))
Xone_CorG_fourier_mag = as.matrix(Mod(Xone_CorG_fourier))
Xone_CorT_fourier_mag = as.matrix(Mod(Xone_CorT_fourier))
Xone_GorT_fourier_mag = as.matrix(Mod(Xone_GorT_fourier))

Xone_AorC_fourier_phase = as.matrix(Arg(Xone_AorC_fourier))
Xone_AorG_fourier_phase = as.matrix(Arg(Xone_AorG_fourier))
Xone_AorT_fourier_phase = as.matrix(Arg(Xone_AorT_fourier))
Xone_CorG_fourier_phase = as.matrix(Arg(Xone_CorG_fourier))
Xone_CorT_fourier_phase = as.matrix(Arg(Xone_CorT_fourier))
Xone_GorT_fourier_phase = as.matrix(Arg(Xone_GorT_fourier))

colnames(Xone_AorC_fourier_mag) = "fourier_mag_AorC_6"
colnames(Xone_AorG_fourier_mag) = "fourier_mag_AorG_6"
colnames(Xone_AorT_fourier_mag) = "fourier_mag_AorT_6"
colnames(Xone_CorG_fourier_mag) = "fourier_mag_CorG_6"
colnames(Xone_CorT_fourier_mag) = "fourier_mag_CorT_6"
colnames(Xone_GorT_fourier_mag) = "fourier_mag_GorT_6"

colnames(Xone_AorC_fourier_phase) = "fourier_phase_AorC_6"
colnames(Xone_AorG_fourier_phase) = "fourier_phase_AorG_6"
colnames(Xone_AorT_fourier_phase) = "fourier_phase_AorT_6"
colnames(Xone_CorG_fourier_phase) = "fourier_phase_CorG_6"
colnames(Xone_CorT_fourier_phase) = "fourier_phase_CorT_6"
colnames(Xone_GorT_fourier_phase) = "fourier_phase_GorT_6"




# Dinucleotides:

Xtwo = dat %>% select(all_of(ps2))
Xtwo_AA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_AC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_AG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_GA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_GC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_GG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_GT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_TA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_TC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_TG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_TT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AA) = colnames(Xtwo)
colnames(Xtwo_AC) = colnames(Xtwo)
colnames(Xtwo_AG) = colnames(Xtwo)
colnames(Xtwo_AT) = colnames(Xtwo)
colnames(Xtwo_CA) = colnames(Xtwo)
colnames(Xtwo_CC) = colnames(Xtwo)
colnames(Xtwo_CG) = colnames(Xtwo)
colnames(Xtwo_CT) = colnames(Xtwo)
colnames(Xtwo_GA) = colnames(Xtwo)
colnames(Xtwo_GC) = colnames(Xtwo)
colnames(Xtwo_GG) = colnames(Xtwo)
colnames(Xtwo_GT) = colnames(Xtwo)
colnames(Xtwo_TA) = colnames(Xtwo)
colnames(Xtwo_TC) = colnames(Xtwo)
colnames(Xtwo_TG) = colnames(Xtwo)
colnames(Xtwo_TT) = colnames(Xtwo)
Xtwo_AA[] = (Xtwo == "AA") %>% as.matrix() %>% as.numeric()
Xtwo_AC[] = (Xtwo == "AC") %>% as.matrix() %>% as.numeric()
Xtwo_AG[] = (Xtwo == "AG") %>% as.matrix() %>% as.numeric()
Xtwo_AT[] = (Xtwo == "AT") %>% as.matrix() %>% as.numeric()
Xtwo_CA[] = (Xtwo == "CA") %>% as.matrix() %>% as.numeric()
Xtwo_CC[] = (Xtwo == "CC") %>% as.matrix() %>% as.numeric()
Xtwo_CG[] = (Xtwo == "CG") %>% as.matrix() %>% as.numeric()
Xtwo_CT[] = (Xtwo == "CT") %>% as.matrix() %>% as.numeric()
Xtwo_GA[] = (Xtwo == "GA") %>% as.matrix() %>% as.numeric()
Xtwo_GC[] = (Xtwo == "GC") %>% as.matrix() %>% as.numeric()
Xtwo_GG[] = (Xtwo == "GG") %>% as.matrix() %>% as.numeric()
Xtwo_GT[] = (Xtwo == "GT") %>% as.matrix() %>% as.numeric()
Xtwo_TA[] = (Xtwo == "TA") %>% as.matrix() %>% as.numeric()
Xtwo_TC[] = (Xtwo == "TC") %>% as.matrix() %>% as.numeric()
Xtwo_TG[] = (Xtwo == "TG") %>% as.matrix() %>% as.numeric()
Xtwo_TT[] = (Xtwo == "TT") %>% as.matrix() %>% as.numeric()

# Xtwo_AA_fourier = t(apply(Xtwo_AA, 1, fft))[,1:25]
# Xtwo_AC_fourier = t(apply(Xtwo_AC, 1, fft))[,1:25]
# Xtwo_AG_fourier = t(apply(Xtwo_AG, 1, fft))[,1:25]
# Xtwo_AT_fourier = t(apply(Xtwo_AT, 1, fft))[,1:25]
# Xtwo_CA_fourier = t(apply(Xtwo_CA, 1, fft))[,1:25]
# Xtwo_CC_fourier = t(apply(Xtwo_CC, 1, fft))[,1:25]
# Xtwo_CG_fourier = t(apply(Xtwo_CG, 1, fft))[,1:25]
# Xtwo_CT_fourier = t(apply(Xtwo_CT, 1, fft))[,1:25]
# Xtwo_GA_fourier = t(apply(Xtwo_GA, 1, fft))[,1:25]
# Xtwo_GC_fourier = t(apply(Xtwo_GC, 1, fft))[,1:25]
# Xtwo_GG_fourier = t(apply(Xtwo_GG, 1, fft))[,1:25]
# Xtwo_GT_fourier = t(apply(Xtwo_GT, 1, fft))[,1:25]
# Xtwo_TA_fourier = t(apply(Xtwo_TA, 1, fft))[,1:25]
# Xtwo_TC_fourier = t(apply(Xtwo_TC, 1, fft))[,1:25]
# Xtwo_TG_fourier = t(apply(Xtwo_TG, 1, fft))[,1:25]
# Xtwo_TT_fourier = t(apply(Xtwo_TT, 1, fft))[,1:25]

Xtwo_AA_fourier = t(apply(Xtwo_AA, 1, fft))[,6]
Xtwo_AC_fourier = t(apply(Xtwo_AC, 1, fft))[,6]
Xtwo_AG_fourier = t(apply(Xtwo_AG, 1, fft))[,6]
Xtwo_AT_fourier = t(apply(Xtwo_AT, 1, fft))[,6]
Xtwo_CA_fourier = t(apply(Xtwo_CA, 1, fft))[,6]
Xtwo_CC_fourier = t(apply(Xtwo_CC, 1, fft))[,6]
Xtwo_CG_fourier = t(apply(Xtwo_CG, 1, fft))[,6]
Xtwo_CT_fourier = t(apply(Xtwo_CT, 1, fft))[,6]
Xtwo_GA_fourier = t(apply(Xtwo_GA, 1, fft))[,6]
Xtwo_GC_fourier = t(apply(Xtwo_GC, 1, fft))[,6]
Xtwo_GG_fourier = t(apply(Xtwo_GG, 1, fft))[,6]
Xtwo_GT_fourier = t(apply(Xtwo_GT, 1, fft))[,6]
Xtwo_TA_fourier = t(apply(Xtwo_TA, 1, fft))[,6]
Xtwo_TC_fourier = t(apply(Xtwo_TC, 1, fft))[,6]
Xtwo_TG_fourier = t(apply(Xtwo_TG, 1, fft))[,6]
Xtwo_TT_fourier = t(apply(Xtwo_TT, 1, fft))[,6]

Xtwo_AA_fourier_mag = as.matrix(Mod(Xtwo_AA_fourier))
Xtwo_AC_fourier_mag = as.matrix(Mod(Xtwo_AC_fourier))
Xtwo_AG_fourier_mag = as.matrix(Mod(Xtwo_AG_fourier))
Xtwo_AT_fourier_mag = as.matrix(Mod(Xtwo_AT_fourier))
Xtwo_CA_fourier_mag = as.matrix(Mod(Xtwo_CA_fourier))
Xtwo_CC_fourier_mag = as.matrix(Mod(Xtwo_CC_fourier))
Xtwo_CG_fourier_mag = as.matrix(Mod(Xtwo_CG_fourier))
Xtwo_CT_fourier_mag = as.matrix(Mod(Xtwo_CT_fourier))
Xtwo_GA_fourier_mag = as.matrix(Mod(Xtwo_GA_fourier))
Xtwo_GC_fourier_mag = as.matrix(Mod(Xtwo_GC_fourier))
Xtwo_GG_fourier_mag = as.matrix(Mod(Xtwo_GG_fourier))
Xtwo_GT_fourier_mag = as.matrix(Mod(Xtwo_GT_fourier))
Xtwo_TA_fourier_mag = as.matrix(Mod(Xtwo_TA_fourier))
Xtwo_TC_fourier_mag = as.matrix(Mod(Xtwo_TC_fourier))
Xtwo_TG_fourier_mag = as.matrix(Mod(Xtwo_TG_fourier))
Xtwo_TT_fourier_mag = as.matrix(Mod(Xtwo_TT_fourier))

Xtwo_AA_fourier_phase = as.matrix(Arg(Xtwo_AA_fourier))
Xtwo_AC_fourier_phase = as.matrix(Arg(Xtwo_AC_fourier))
Xtwo_AG_fourier_phase = as.matrix(Arg(Xtwo_AG_fourier))
Xtwo_AT_fourier_phase = as.matrix(Arg(Xtwo_AT_fourier))
Xtwo_CA_fourier_phase = as.matrix(Arg(Xtwo_CA_fourier))
Xtwo_CC_fourier_phase = as.matrix(Arg(Xtwo_CC_fourier))
Xtwo_CG_fourier_phase = as.matrix(Arg(Xtwo_CG_fourier))
Xtwo_CT_fourier_phase = as.matrix(Arg(Xtwo_CT_fourier))
Xtwo_GA_fourier_phase = as.matrix(Arg(Xtwo_GA_fourier))
Xtwo_GC_fourier_phase = as.matrix(Arg(Xtwo_GC_fourier))
Xtwo_GG_fourier_phase = as.matrix(Arg(Xtwo_GG_fourier))
Xtwo_GT_fourier_phase = as.matrix(Arg(Xtwo_GT_fourier))
Xtwo_TA_fourier_phase = as.matrix(Arg(Xtwo_TA_fourier))
Xtwo_TC_fourier_phase = as.matrix(Arg(Xtwo_TC_fourier))
Xtwo_TG_fourier_phase = as.matrix(Arg(Xtwo_TG_fourier))
Xtwo_TT_fourier_phase = as.matrix(Arg(Xtwo_TT_fourier))

# colnames(Xtwo_AA_fourier_mag) = paste0("fourier_mag_AA_", 1:25)
# colnames(Xtwo_AC_fourier_mag) = paste0("fourier_mag_AC_", 1:25)
# colnames(Xtwo_AG_fourier_mag) = paste0("fourier_mag_AG_", 1:25)
# colnames(Xtwo_AT_fourier_mag) = paste0("fourier_mag_AT_", 1:25)
# colnames(Xtwo_CA_fourier_mag) = paste0("fourier_mag_CA_", 1:25)
# colnames(Xtwo_CC_fourier_mag) = paste0("fourier_mag_CC_", 1:25)
# colnames(Xtwo_CG_fourier_mag) = paste0("fourier_mag_CG_", 1:25)
# colnames(Xtwo_CT_fourier_mag) = paste0("fourier_mag_CT_", 1:25)
# colnames(Xtwo_GA_fourier_mag) = paste0("fourier_mag_GA_", 1:25)
# colnames(Xtwo_GC_fourier_mag) = paste0("fourier_mag_GC_", 1:25)
# colnames(Xtwo_GG_fourier_mag) = paste0("fourier_mag_GG_", 1:25)
# colnames(Xtwo_GT_fourier_mag) = paste0("fourier_mag_GT_", 1:25)
# colnames(Xtwo_TA_fourier_mag) = paste0("fourier_mag_TA_", 1:25)
# colnames(Xtwo_TC_fourier_mag) = paste0("fourier_mag_TC_", 1:25)
# colnames(Xtwo_TG_fourier_mag) = paste0("fourier_mag_TG_", 1:25)
# colnames(Xtwo_TT_fourier_mag) = paste0("fourier_mag_TT_", 1:25)
# 
# colnames(Xtwo_AA_fourier_phase) = paste0("fourier_phase_AA_", 1:25)
# colnames(Xtwo_AC_fourier_phase) = paste0("fourier_phase_AC_", 1:25)
# colnames(Xtwo_AG_fourier_phase) = paste0("fourier_phase_AG_", 1:25)
# colnames(Xtwo_AT_fourier_phase) = paste0("fourier_phase_AT_", 1:25)
# colnames(Xtwo_CA_fourier_phase) = paste0("fourier_phase_CA_", 1:25)
# colnames(Xtwo_CC_fourier_phase) = paste0("fourier_phase_CC_", 1:25)
# colnames(Xtwo_CG_fourier_phase) = paste0("fourier_phase_CG_", 1:25)
# colnames(Xtwo_CT_fourier_phase) = paste0("fourier_phase_CT_", 1:25)
# colnames(Xtwo_GA_fourier_phase) = paste0("fourier_phase_GA_", 1:25)
# colnames(Xtwo_GC_fourier_phase) = paste0("fourier_phase_GC_", 1:25)
# colnames(Xtwo_GG_fourier_phase) = paste0("fourier_phase_GG_", 1:25)
# colnames(Xtwo_GT_fourier_phase) = paste0("fourier_phase_GT_", 1:25)
# colnames(Xtwo_TA_fourier_phase) = paste0("fourier_phase_TA_", 1:25)
# colnames(Xtwo_TC_fourier_phase) = paste0("fourier_phase_TC_", 1:25)
# colnames(Xtwo_TG_fourier_phase) = paste0("fourier_phase_TG_", 1:25)
# colnames(Xtwo_TT_fourier_phase) = paste0("fourier_phase_TT_", 1:25)

colnames(Xtwo_AA_fourier_mag) = "fourier_mag_AA_6"
colnames(Xtwo_AC_fourier_mag) = "fourier_mag_AC_6"
colnames(Xtwo_AG_fourier_mag) = "fourier_mag_AG_6"
colnames(Xtwo_AT_fourier_mag) = "fourier_mag_AT_6"
colnames(Xtwo_CA_fourier_mag) = "fourier_mag_CA_6"
colnames(Xtwo_CC_fourier_mag) = "fourier_mag_CC_6"
colnames(Xtwo_CG_fourier_mag) = "fourier_mag_CG_6"
colnames(Xtwo_CT_fourier_mag) = "fourier_mag_CT_6"
colnames(Xtwo_GA_fourier_mag) = "fourier_mag_GA_6"
colnames(Xtwo_GC_fourier_mag) = "fourier_mag_GC_6"
colnames(Xtwo_GG_fourier_mag) = "fourier_mag_GG_6"
colnames(Xtwo_GT_fourier_mag) = "fourier_mag_GT_6"
colnames(Xtwo_TA_fourier_mag) = "fourier_mag_TA_6"
colnames(Xtwo_TC_fourier_mag) = "fourier_mag_TC_6"
colnames(Xtwo_TG_fourier_mag) = "fourier_mag_TG_6"
colnames(Xtwo_TT_fourier_mag) = "fourier_mag_TT_6"

colnames(Xtwo_AA_fourier_phase) = "fourier_phase_AA_6"
colnames(Xtwo_AC_fourier_phase) = "fourier_phase_AC_6"
colnames(Xtwo_AG_fourier_phase) = "fourier_phase_AG_6"
colnames(Xtwo_AT_fourier_phase) = "fourier_phase_AT_6"
colnames(Xtwo_CA_fourier_phase) = "fourier_phase_CA_6"
colnames(Xtwo_CC_fourier_phase) = "fourier_phase_CC_6"
colnames(Xtwo_CG_fourier_phase) = "fourier_phase_CG_6"
colnames(Xtwo_CT_fourier_phase) = "fourier_phase_CT_6"
colnames(Xtwo_GA_fourier_phase) = "fourier_phase_GA_6"
colnames(Xtwo_GC_fourier_phase) = "fourier_phase_GC_6"
colnames(Xtwo_GG_fourier_phase) = "fourier_phase_GG_6"
colnames(Xtwo_GT_fourier_phase) = "fourier_phase_GT_6"
colnames(Xtwo_TA_fourier_phase) = "fourier_phase_TA_6"
colnames(Xtwo_TC_fourier_phase) = "fourier_phase_TC_6"
colnames(Xtwo_TG_fourier_phase) = "fourier_phase_TG_6"
colnames(Xtwo_TT_fourier_phase) = "fourier_phase_TT_6"

# Pairs of Dinucleotides:

Xtwo_AAorAT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_AAorTA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_AAorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_ATorTA = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_ATorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_TAorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorCG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorGC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CGorGC = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CGorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_GCorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AAorAT) = colnames(Xtwo)
colnames(Xtwo_AAorTA) = colnames(Xtwo)
colnames(Xtwo_AAorTT) = colnames(Xtwo)
colnames(Xtwo_ATorTA) = colnames(Xtwo)
colnames(Xtwo_ATorTT) = colnames(Xtwo)
colnames(Xtwo_TAorTT) = colnames(Xtwo)
colnames(Xtwo_CCorCG) = colnames(Xtwo)
colnames(Xtwo_CCorGC) = colnames(Xtwo)
colnames(Xtwo_CCorGG) = colnames(Xtwo)
colnames(Xtwo_CGorGC) = colnames(Xtwo)
colnames(Xtwo_CGorGG) = colnames(Xtwo)
colnames(Xtwo_GCorGG) = colnames(Xtwo)
Xtwo_AAorAT[] = ((Xtwo == "AA") | (Xtwo == "AT")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTA[] = ((Xtwo == "AA") | (Xtwo == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTT[] = ((Xtwo == "AA") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTA[] = ((Xtwo == "AT") | (Xtwo == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTT[] = ((Xtwo == "AT") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_TAorTT[] = ((Xtwo == "TA") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCG[] = ((Xtwo == "CC") | (Xtwo == "CG")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGC[] = ((Xtwo == "CC") | (Xtwo == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGG[] = ((Xtwo == "CC") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGC[] = ((Xtwo == "CG") | (Xtwo == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGG[] = ((Xtwo == "CG") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_GCorGG[] = ((Xtwo == "GC") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()

# Xtwo_AAorAT_fourier = t(apply(Xtwo_AAorAT, 1, fft))[,1:25]
# Xtwo_AAorTA_fourier = t(apply(Xtwo_AAorTA, 1, fft))[,1:25]
# Xtwo_AAorTT_fourier = t(apply(Xtwo_AAorTT, 1, fft))[,1:25]
# Xtwo_ATorTA_fourier = t(apply(Xtwo_ATorTA, 1, fft))[,1:25]
# Xtwo_ATorTT_fourier = t(apply(Xtwo_ATorTT, 1, fft))[,1:25]
# Xtwo_TAorTT_fourier = t(apply(Xtwo_TAorTT, 1, fft))[,1:25]
# Xtwo_CCorCG_fourier = t(apply(Xtwo_CCorCG, 1, fft))[,1:25]
# Xtwo_CCorGC_fourier = t(apply(Xtwo_CCorGC, 1, fft))[,1:25]
# Xtwo_CCorGG_fourier = t(apply(Xtwo_CCorGG, 1, fft))[,1:25]
# Xtwo_CGorGC_fourier = t(apply(Xtwo_CGorGC, 1, fft))[,1:25]
# Xtwo_CGorGG_fourier = t(apply(Xtwo_CGorGG, 1, fft))[,1:25]
# Xtwo_GCorGG_fourier = t(apply(Xtwo_GCorGG, 1, fft))[,1:25]

Xtwo_AAorAT_fourier = t(apply(Xtwo_AAorAT, 1, fft))[,6]
Xtwo_AAorTA_fourier = t(apply(Xtwo_AAorTA, 1, fft))[,6]
Xtwo_AAorTT_fourier = t(apply(Xtwo_AAorTT, 1, fft))[,6]
Xtwo_ATorTA_fourier = t(apply(Xtwo_ATorTA, 1, fft))[,6]
Xtwo_ATorTT_fourier = t(apply(Xtwo_ATorTT, 1, fft))[,6]
Xtwo_TAorTT_fourier = t(apply(Xtwo_TAorTT, 1, fft))[,6]
Xtwo_CCorCG_fourier = t(apply(Xtwo_CCorCG, 1, fft))[,6]
Xtwo_CCorGC_fourier = t(apply(Xtwo_CCorGC, 1, fft))[,6]
Xtwo_CCorGG_fourier = t(apply(Xtwo_CCorGG, 1, fft))[,6]
Xtwo_CGorGC_fourier = t(apply(Xtwo_CGorGC, 1, fft))[,6]
Xtwo_CGorGG_fourier = t(apply(Xtwo_CGorGG, 1, fft))[,6]
Xtwo_GCorGG_fourier = t(apply(Xtwo_GCorGG, 1, fft))[,6]

Xtwo_AAorAT_fourier_mag = as.matrix(Mod(Xtwo_AAorAT_fourier))
Xtwo_AAorTA_fourier_mag = as.matrix(Mod(Xtwo_AAorTA_fourier))
Xtwo_AAorTT_fourier_mag = as.matrix(Mod(Xtwo_AAorTT_fourier))
Xtwo_ATorTA_fourier_mag = as.matrix(Mod(Xtwo_ATorTA_fourier))
Xtwo_ATorTT_fourier_mag = as.matrix(Mod(Xtwo_ATorTT_fourier))
Xtwo_TAorTT_fourier_mag = as.matrix(Mod(Xtwo_TAorTT_fourier))
Xtwo_CCorCG_fourier_mag = as.matrix(Mod(Xtwo_CCorCG_fourier))
Xtwo_CCorGC_fourier_mag = as.matrix(Mod(Xtwo_CCorGC_fourier))
Xtwo_CCorGG_fourier_mag = as.matrix(Mod(Xtwo_CCorGG_fourier))
Xtwo_CGorGC_fourier_mag = as.matrix(Mod(Xtwo_CGorGC_fourier))
Xtwo_CGorGG_fourier_mag = as.matrix(Mod(Xtwo_CGorGG_fourier))
Xtwo_GCorGG_fourier_mag = as.matrix(Mod(Xtwo_GCorGG_fourier))

Xtwo_AAorAT_fourier_phase = as.matrix(Arg(Xtwo_AAorAT_fourier))
Xtwo_AAorTA_fourier_phase = as.matrix(Arg(Xtwo_AAorTA_fourier))
Xtwo_AAorTT_fourier_phase = as.matrix(Arg(Xtwo_AAorTT_fourier))
Xtwo_ATorTA_fourier_phase = as.matrix(Arg(Xtwo_ATorTA_fourier))
Xtwo_ATorTT_fourier_phase = as.matrix(Arg(Xtwo_ATorTT_fourier))
Xtwo_TAorTT_fourier_phase = as.matrix(Arg(Xtwo_TAorTT_fourier))
Xtwo_CCorCG_fourier_phase = as.matrix(Arg(Xtwo_CCorCG_fourier))
Xtwo_CCorGC_fourier_phase = as.matrix(Arg(Xtwo_CCorGC_fourier))
Xtwo_CCorGG_fourier_phase = as.matrix(Arg(Xtwo_CCorGG_fourier))
Xtwo_CGorGC_fourier_phase = as.matrix(Arg(Xtwo_CGorGC_fourier))
Xtwo_CGorGG_fourier_phase = as.matrix(Arg(Xtwo_CGorGG_fourier))
Xtwo_GCorGG_fourier_phase = as.matrix(Arg(Xtwo_GCorGG_fourier))


colnames(Xtwo_AAorAT_fourier_mag) = "fourier_mag_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_mag) = "fourier_mag_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_mag) = "fourier_mag_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_mag) = "fourier_mag_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_mag) = "fourier_mag_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_mag) = "fourier_mag_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_mag) = "fourier_mag_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_mag) = "fourier_mag_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_mag) = "fourier_mag_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_mag) = "fourier_mag_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_mag) = "fourier_mag_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_mag) = "fourier_mag_GCorGG_6"

colnames(Xtwo_AAorAT_fourier_phase) = "fourier_phase_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_phase) = "fourier_phase_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_phase) = "fourier_phase_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_phase) = "fourier_phase_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_phase) = "fourier_phase_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_phase) = "fourier_phase_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_phase) = "fourier_phase_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_phase) = "fourier_phase_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_phase) = "fourier_phase_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_phase) = "fourier_phase_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_phase) = "fourier_phase_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_phase) = "fourier_phase_GCorGG_6"


# Groups of Dinucleotides:

Xtwo_AAorATorTAorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorCGorGCorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AAorATorTAorTT) = colnames(Xtwo)
colnames(Xtwo_CCorCGorGCorGG) = colnames(Xtwo)
Xtwo_AAorATorTAorTT[] = ((Xtwo == "AA") | (Xtwo == "AT") | 
                           (Xtwo == "TA") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG[] = ((Xtwo == "CC") | (Xtwo == "CG") |
                           (Xtwo == "GC") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_fourier = t(apply(Xtwo_AAorATorTAorTT, 1, fft))[,6]
Xtwo_CCorCGorGCorGG_fourier = t(apply(Xtwo_CCorCGorGCorGG, 1, fft))[,6]

Xtwo_AAorATorTAorTT_fourier_mag = as.matrix(Mod(Xtwo_AAorATorTAorTT_fourier))
Xtwo_CCorCGorGCorGG_fourier_mag = as.matrix(Mod(Xtwo_CCorCGorGCorGG_fourier))

Xtwo_AAorATorTAorTT_fourier_phase = as.matrix(Arg(Xtwo_AAorATorTAorTT_fourier))
Xtwo_CCorCGorGCorGG_fourier_phase = as.matrix(Arg(Xtwo_CCorCGorGCorGG_fourier))


colnames(Xtwo_AAorATorTAorTT_fourier_mag) = "fourier_mag_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_mag) = "fourier_mag_CCorCGorGCorGG_6"

colnames(Xtwo_AAorATorTAorTT_fourier_phase) = "fourier_phase_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_phase) = "fourier_phase_CCorCGorGCorGG_6"




## Test data:

# Nucleotides:

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

# Xone_A_fourier_test = t(apply(Xone_A_test, 1, fft))[,1:26]
# Xone_C_fourier_test = t(apply(Xone_C_test, 1, fft))[,1:26]
# Xone_G_fourier_test = t(apply(Xone_G_test, 1, fft))[,1:26]
# Xone_T_fourier_test = t(apply(Xone_T_test, 1, fft))[,1:26]

Xone_A_fourier_test = t(apply(Xone_A_test, 1, fft))[,6]
Xone_C_fourier_test = t(apply(Xone_C_test, 1, fft))[,6]
Xone_G_fourier_test = t(apply(Xone_G_test, 1, fft))[,6]
Xone_T_fourier_test = t(apply(Xone_T_test, 1, fft))[,6]

Xone_A_fourier_mag_test = as.matrix(Mod(Xone_A_fourier_test))
Xone_C_fourier_mag_test = as.matrix(Mod(Xone_C_fourier_test))
Xone_G_fourier_mag_test = as.matrix(Mod(Xone_G_fourier_test))
Xone_T_fourier_mag_test = as.matrix(Mod(Xone_T_fourier_test))

Xone_A_fourier_phase_test = as.matrix(Arg(Xone_A_fourier_test))
Xone_C_fourier_phase_test = as.matrix(Arg(Xone_C_fourier_test))
Xone_G_fourier_phase_test = as.matrix(Arg(Xone_G_fourier_test))
Xone_T_fourier_phase_test = as.matrix(Arg(Xone_T_fourier_test))

# colnames(Xone_A_fourier_mag_test) = paste0("fourier_mag_A_", 1:26)
# colnames(Xone_C_fourier_mag_test) = paste0("fourier_mag_C_", 1:26)
# colnames(Xone_G_fourier_mag_test) = paste0("fourier_mag_G_", 1:26)
# colnames(Xone_T_fourier_mag_test) = paste0("fourier_mag_T_", 1:26)

# colnames(Xone_A_fourier_phase_test) = paste0("fourier_phase_A_", 1:26)
# colnames(Xone_C_fourier_phase_test) = paste0("fourier_phase_C_", 1:26)
# colnames(Xone_G_fourier_phase_test) = paste0("fourier_phase_G_", 1:26)
# colnames(Xone_T_fourier_phase_test) = paste0("fourier_phase_T_", 1:26)

colnames(Xone_A_fourier_mag_test) = "fourier_mag_A_6"
colnames(Xone_C_fourier_mag_test) = "fourier_mag_C_6"
colnames(Xone_G_fourier_mag_test) = "fourier_mag_G_6"
colnames(Xone_T_fourier_mag_test) = "fourier_mag_T_6"

colnames(Xone_A_fourier_phase_test) = "fourier_phase_A_6"
colnames(Xone_C_fourier_phase_test) = "fourier_phase_C_6"
colnames(Xone_G_fourier_phase_test) = "fourier_phase_G_6"
colnames(Xone_T_fourier_phase_test) = "fourier_phase_T_6"


# Pairs of nucleotides:

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

# Xone_AorC_fourier_test = t(apply(Xone_AorC_test, 1, fft))[,1:26]
# Xone_AorG_fourier_test = t(apply(Xone_AorG_test, 1, fft))[,1:26]
# Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,1:26]
# Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,1:26]
# Xone_CorT_fourier_test = t(apply(Xone_CorT_test, 1, fft))[,1:26]
# Xone_GorT_fourier_test = t(apply(Xone_GorT_test, 1, fft))[,1:26]

Xone_AorC_fourier_test = t(apply(Xone_AorC_test, 1, fft))[,6]
Xone_AorG_fourier_test = t(apply(Xone_AorG_test, 1, fft))[,6]
Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,6]
Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,6]
Xone_CorT_fourier_test = t(apply(Xone_CorT_test, 1, fft))[,6]
Xone_GorT_fourier_test = t(apply(Xone_GorT_test, 1, fft))[,6]

Xone_AorC_fourier_mag_test = as.matrix(Mod(Xone_AorC_fourier_test))
Xone_AorG_fourier_mag_test = as.matrix(Mod(Xone_AorG_fourier_test))
Xone_AorT_fourier_mag_test = as.matrix(Mod(Xone_AorT_fourier_test))
Xone_CorG_fourier_mag_test = as.matrix(Mod(Xone_CorG_fourier_test))
Xone_CorT_fourier_mag_test = as.matrix(Mod(Xone_CorT_fourier_test))
Xone_GorT_fourier_mag_test = as.matrix(Mod(Xone_GorT_fourier_test))

Xone_AorC_fourier_phase_test = as.matrix(Arg(Xone_AorC_fourier_test))
Xone_AorG_fourier_phase_test = as.matrix(Arg(Xone_AorG_fourier_test))
Xone_AorT_fourier_phase_test = as.matrix(Arg(Xone_AorT_fourier_test))
Xone_CorG_fourier_phase_test = as.matrix(Arg(Xone_CorG_fourier_test))
Xone_CorT_fourier_phase_test = as.matrix(Arg(Xone_CorT_fourier_test))
Xone_GorT_fourier_phase_test = as.matrix(Arg(Xone_GorT_fourier_test))

colnames(Xone_AorC_fourier_mag_test) = "fourier_mag_AorC_6"
colnames(Xone_AorG_fourier_mag_test) = "fourier_mag_AorG_6"
colnames(Xone_AorT_fourier_mag_test) = "fourier_mag_AorT_6"
colnames(Xone_CorG_fourier_mag_test) = "fourier_mag_CorG_6"
colnames(Xone_CorT_fourier_mag_test) = "fourier_mag_CorT_6"
colnames(Xone_GorT_fourier_mag_test) = "fourier_mag_GorT_6"

colnames(Xone_AorC_fourier_phase_test) = "fourier_phase_AorC_6"
colnames(Xone_AorG_fourier_phase_test) = "fourier_phase_AorG_6"
colnames(Xone_AorT_fourier_phase_test) = "fourier_phase_AorT_6"
colnames(Xone_CorG_fourier_phase_test) = "fourier_phase_CorG_6"
colnames(Xone_CorT_fourier_phase_test) = "fourier_phase_CorT_6"
colnames(Xone_GorT_fourier_phase_test) = "fourier_phase_GorT_6"


# Dinucleotides:

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

# Xtwo_AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,1:25]
# Xtwo_AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,1:25]
# Xtwo_AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,1:25]
# Xtwo_AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,1:25]
# Xtwo_CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,1:25]
# Xtwo_CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,1:25]
# Xtwo_CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,1:25]
# Xtwo_CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,1:25]
# Xtwo_GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,1:25]
# Xtwo_GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,1:25]
# Xtwo_GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,1:25]
# Xtwo_GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,1:25]
# Xtwo_TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,1:25]
# Xtwo_TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,1:25]
# Xtwo_TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,1:25]
# Xtwo_TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,1:25]

Xtwo_AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,6]
Xtwo_AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,6]
Xtwo_AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,6]
Xtwo_AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,6]
Xtwo_CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,6]
Xtwo_CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,6]
Xtwo_CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,6]
Xtwo_CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,6]
Xtwo_GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,6]
Xtwo_GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,6]
Xtwo_GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,6]
Xtwo_GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,6]
Xtwo_TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,6]
Xtwo_TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,6]
Xtwo_TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,6]
Xtwo_TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,6]

Xtwo_AA_fourier_mag_test = as.matrix(Mod(Xtwo_AA_fourier_test))
Xtwo_AC_fourier_mag_test = as.matrix(Mod(Xtwo_AC_fourier_test))
Xtwo_AG_fourier_mag_test = as.matrix(Mod(Xtwo_AG_fourier_test))
Xtwo_AT_fourier_mag_test = as.matrix(Mod(Xtwo_AT_fourier_test))
Xtwo_CA_fourier_mag_test = as.matrix(Mod(Xtwo_CA_fourier_test))
Xtwo_CC_fourier_mag_test = as.matrix(Mod(Xtwo_CC_fourier_test))
Xtwo_CG_fourier_mag_test = as.matrix(Mod(Xtwo_CG_fourier_test))
Xtwo_CT_fourier_mag_test = as.matrix(Mod(Xtwo_CT_fourier_test))
Xtwo_GA_fourier_mag_test = as.matrix(Mod(Xtwo_GA_fourier_test))
Xtwo_GC_fourier_mag_test = as.matrix(Mod(Xtwo_GC_fourier_test))
Xtwo_GG_fourier_mag_test = as.matrix(Mod(Xtwo_GG_fourier_test))
Xtwo_GT_fourier_mag_test = as.matrix(Mod(Xtwo_GT_fourier_test))
Xtwo_TA_fourier_mag_test = as.matrix(Mod(Xtwo_TA_fourier_test))
Xtwo_TC_fourier_mag_test = as.matrix(Mod(Xtwo_TC_fourier_test))
Xtwo_TG_fourier_mag_test = as.matrix(Mod(Xtwo_TG_fourier_test))
Xtwo_TT_fourier_mag_test = as.matrix(Mod(Xtwo_TT_fourier_test))

Xtwo_AA_fourier_phase_test = as.matrix(Arg(Xtwo_AA_fourier_test))
Xtwo_AC_fourier_phase_test = as.matrix(Arg(Xtwo_AC_fourier_test))
Xtwo_AG_fourier_phase_test = as.matrix(Arg(Xtwo_AG_fourier_test))
Xtwo_AT_fourier_phase_test = as.matrix(Arg(Xtwo_AT_fourier_test))
Xtwo_CA_fourier_phase_test = as.matrix(Arg(Xtwo_CA_fourier_test))
Xtwo_CC_fourier_phase_test = as.matrix(Arg(Xtwo_CC_fourier_test))
Xtwo_CG_fourier_phase_test = as.matrix(Arg(Xtwo_CG_fourier_test))
Xtwo_CT_fourier_phase_test = as.matrix(Arg(Xtwo_CT_fourier_test))
Xtwo_GA_fourier_phase_test = as.matrix(Arg(Xtwo_GA_fourier_test))
Xtwo_GC_fourier_phase_test = as.matrix(Arg(Xtwo_GC_fourier_test))
Xtwo_GG_fourier_phase_test = as.matrix(Arg(Xtwo_GG_fourier_test))
Xtwo_GT_fourier_phase_test = as.matrix(Arg(Xtwo_GT_fourier_test))
Xtwo_TA_fourier_phase_test = as.matrix(Arg(Xtwo_TA_fourier_test))
Xtwo_TC_fourier_phase_test = as.matrix(Arg(Xtwo_TC_fourier_test))
Xtwo_TG_fourier_phase_test = as.matrix(Arg(Xtwo_TG_fourier_test))
Xtwo_TT_fourier_phase_test = as.matrix(Arg(Xtwo_TT_fourier_test))

# colnames(Xtwo_AA_fourier_mag_test) = paste0("fourier_mag_AA_", 1:25)
# colnames(Xtwo_AC_fourier_mag_test) = paste0("fourier_mag_AC_", 1:25)
# colnames(Xtwo_AG_fourier_mag_test) = paste0("fourier_mag_AG_", 1:25)
# colnames(Xtwo_AT_fourier_mag_test) = paste0("fourier_mag_AT_", 1:25)
# colnames(Xtwo_CA_fourier_mag_test) = paste0("fourier_mag_CA_", 1:25)
# colnames(Xtwo_CC_fourier_mag_test) = paste0("fourier_mag_CC_", 1:25)
# colnames(Xtwo_CG_fourier_mag_test) = paste0("fourier_mag_CG_", 1:25)
# colnames(Xtwo_CT_fourier_mag_test) = paste0("fourier_mag_CT_", 1:25)
# colnames(Xtwo_GA_fourier_mag_test) = paste0("fourier_mag_GA_", 1:25)
# colnames(Xtwo_GC_fourier_mag_test) = paste0("fourier_mag_GC_", 1:25)
# colnames(Xtwo_GG_fourier_mag_test) = paste0("fourier_mag_GG_", 1:25)
# colnames(Xtwo_GT_fourier_mag_test) = paste0("fourier_mag_GT_", 1:25)
# colnames(Xtwo_TA_fourier_mag_test) = paste0("fourier_mag_TA_", 1:25)
# colnames(Xtwo_TC_fourier_mag_test) = paste0("fourier_mag_TC_", 1:25)
# colnames(Xtwo_TG_fourier_mag_test) = paste0("fourier_mag_TG_", 1:25)
# colnames(Xtwo_TT_fourier_mag_test) = paste0("fourier_mag_TT_", 1:25)
# 
# colnames(Xtwo_AA_fourier_phase_test) = paste0("fourier_phase_AA_", 1:25)
# colnames(Xtwo_AC_fourier_phase_test) = paste0("fourier_phase_AC_", 1:25)
# colnames(Xtwo_AG_fourier_phase_test) = paste0("fourier_phase_AG_", 1:25)
# colnames(Xtwo_AT_fourier_phase_test) = paste0("fourier_phase_AT_", 1:25)
# colnames(Xtwo_CA_fourier_phase_test) = paste0("fourier_phase_CA_", 1:25)
# colnames(Xtwo_CC_fourier_phase_test) = paste0("fourier_phase_CC_", 1:25)
# colnames(Xtwo_CG_fourier_phase_test) = paste0("fourier_phase_CG_", 1:25)
# colnames(Xtwo_CT_fourier_phase_test) = paste0("fourier_phase_CT_", 1:25)
# colnames(Xtwo_GA_fourier_phase_test) = paste0("fourier_phase_GA_", 1:25)
# colnames(Xtwo_GC_fourier_phase_test) = paste0("fourier_phase_GC_", 1:25)
# colnames(Xtwo_GG_fourier_phase_test) = paste0("fourier_phase_GG_", 1:25)
# colnames(Xtwo_GT_fourier_phase_test) = paste0("fourier_phase_GT_", 1:25)
# colnames(Xtwo_TA_fourier_phase_test) = paste0("fourier_phase_TA_", 1:25)
# colnames(Xtwo_TC_fourier_phase_test) = paste0("fourier_phase_TC_", 1:25)
# colnames(Xtwo_TG_fourier_phase_test) = paste0("fourier_phase_TG_", 1:25)
# colnames(Xtwo_TT_fourier_phase_test) = paste0("fourier_phase_TT_", 1:25)

colnames(Xtwo_AA_fourier_mag_test) = "fourier_mag_AA_6"
colnames(Xtwo_AC_fourier_mag_test) = "fourier_mag_AC_6"
colnames(Xtwo_AG_fourier_mag_test) = "fourier_mag_AG_6"
colnames(Xtwo_AT_fourier_mag_test) = "fourier_mag_AT_6"
colnames(Xtwo_CA_fourier_mag_test) = "fourier_mag_CA_6"
colnames(Xtwo_CC_fourier_mag_test) = "fourier_mag_CC_6"
colnames(Xtwo_CG_fourier_mag_test) = "fourier_mag_CG_6"
colnames(Xtwo_CT_fourier_mag_test) = "fourier_mag_CT_6"
colnames(Xtwo_GA_fourier_mag_test) = "fourier_mag_GA_6"
colnames(Xtwo_GC_fourier_mag_test) = "fourier_mag_GC_6"
colnames(Xtwo_GG_fourier_mag_test) = "fourier_mag_GG_6"
colnames(Xtwo_GT_fourier_mag_test) = "fourier_mag_GT_6"
colnames(Xtwo_TA_fourier_mag_test) = "fourier_mag_TA_6"
colnames(Xtwo_TC_fourier_mag_test) = "fourier_mag_TC_6"
colnames(Xtwo_TG_fourier_mag_test) = "fourier_mag_TG_6"
colnames(Xtwo_TT_fourier_mag_test) = "fourier_mag_TT_6"

colnames(Xtwo_AA_fourier_phase_test) = "fourier_phase_AA_6"
colnames(Xtwo_AC_fourier_phase_test) = "fourier_phase_AC_6"
colnames(Xtwo_AG_fourier_phase_test) = "fourier_phase_AG_6"
colnames(Xtwo_AT_fourier_phase_test) = "fourier_phase_AT_6"
colnames(Xtwo_CA_fourier_phase_test) = "fourier_phase_CA_6"
colnames(Xtwo_CC_fourier_phase_test) = "fourier_phase_CC_6"
colnames(Xtwo_CG_fourier_phase_test) = "fourier_phase_CG_6"
colnames(Xtwo_CT_fourier_phase_test) = "fourier_phase_CT_6"
colnames(Xtwo_GA_fourier_phase_test) = "fourier_phase_GA_6"
colnames(Xtwo_GC_fourier_phase_test) = "fourier_phase_GC_6"
colnames(Xtwo_GG_fourier_phase_test) = "fourier_phase_GG_6"
colnames(Xtwo_GT_fourier_phase_test) = "fourier_phase_GT_6"
colnames(Xtwo_TA_fourier_phase_test) = "fourier_phase_TA_6"
colnames(Xtwo_TC_fourier_phase_test) = "fourier_phase_TC_6"
colnames(Xtwo_TG_fourier_phase_test) = "fourier_phase_TG_6"
colnames(Xtwo_TT_fourier_phase_test) = "fourier_phase_TT_6"

# Pairs of Dinucleotides:

Xtwo_AAorAT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_AAorTA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_AAorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_ATorTA_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_ATorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_TAorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorCG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorGC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CGorGC_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CGorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_GCorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AAorAT_test) = colnames(Xtwo_test)
colnames(Xtwo_AAorTA_test) = colnames(Xtwo_test)
colnames(Xtwo_AAorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_ATorTA_test) = colnames(Xtwo_test)
colnames(Xtwo_ATorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_TAorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorCG_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorGC_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorGG_test) = colnames(Xtwo_test)
colnames(Xtwo_CGorGC_test) = colnames(Xtwo_test)
colnames(Xtwo_CGorGG_test) = colnames(Xtwo_test)
colnames(Xtwo_GCorGG_test) = colnames(Xtwo_test)
Xtwo_AAorAT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "AT")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTA_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTA_test[] = ((Xtwo_test == "AT") | (Xtwo_test == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTT_test[] = ((Xtwo_test == "AT") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_TAorTT_test[] = ((Xtwo_test == "TA") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "CG")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGC_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGC_test[] = ((Xtwo_test == "CG") | (Xtwo_test == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGG_test[] = ((Xtwo_test == "CG") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_GCorGG_test[] = ((Xtwo_test == "GC") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()

# Xtwo_AAorAT_fourier_test = t(apply(Xtwo_AAorAT_test, 1, fft))[,1:25]
# Xtwo_AAorTA_fourier_test = t(apply(Xtwo_AAorTA_test, 1, fft))[,1:25]
# Xtwo_AAorTT_fourier_test = t(apply(Xtwo_AAorTT_test, 1, fft))[,1:25]
# Xtwo_ATorTA_fourier_test = t(apply(Xtwo_ATorTA_test, 1, fft))[,1:25]
# Xtwo_ATorTT_fourier_test = t(apply(Xtwo_ATorTT_test, 1, fft))[,1:25]
# Xtwo_TAorTT_fourier_test = t(apply(Xtwo_TAorTT_test, 1, fft))[,1:25]
# Xtwo_CCorCG_fourier_test = t(apply(Xtwo_CCorCG_test, 1, fft))[,1:25]
# Xtwo_CCorGC_fourier_test = t(apply(Xtwo_CCorGC_test, 1, fft))[,1:25]
# Xtwo_CCorGG_fourier_test = t(apply(Xtwo_CCorGG_test, 1, fft))[,1:25]
# Xtwo_CGorGC_fourier_test = t(apply(Xtwo_CGorGC_test, 1, fft))[,1:25]
# Xtwo_CGorGG_fourier_test = t(apply(Xtwo_CGorGG_test, 1, fft))[,1:25]
# Xtwo_GCorGG_fourier_test = t(apply(Xtwo_GCorGG_test, 1, fft))[,1:25]

Xtwo_AAorAT_fourier_test = t(apply(Xtwo_AAorAT_test, 1, fft))[,6]
Xtwo_AAorTA_fourier_test = t(apply(Xtwo_AAorTA_test, 1, fft))[,6]
Xtwo_AAorTT_fourier_test = t(apply(Xtwo_AAorTT_test, 1, fft))[,6]
Xtwo_ATorTA_fourier_test = t(apply(Xtwo_ATorTA_test, 1, fft))[,6]
Xtwo_ATorTT_fourier_test = t(apply(Xtwo_ATorTT_test, 1, fft))[,6]
Xtwo_TAorTT_fourier_test = t(apply(Xtwo_TAorTT_test, 1, fft))[,6]
Xtwo_CCorCG_fourier_test = t(apply(Xtwo_CCorCG_test, 1, fft))[,6]
Xtwo_CCorGC_fourier_test = t(apply(Xtwo_CCorGC_test, 1, fft))[,6]
Xtwo_CCorGG_fourier_test = t(apply(Xtwo_CCorGG_test, 1, fft))[,6]
Xtwo_CGorGC_fourier_test = t(apply(Xtwo_CGorGC_test, 1, fft))[,6]
Xtwo_CGorGG_fourier_test = t(apply(Xtwo_CGorGG_test, 1, fft))[,6]
Xtwo_GCorGG_fourier_test = t(apply(Xtwo_GCorGG_test, 1, fft))[,6]

Xtwo_AAorAT_fourier_mag_test = as.matrix(Mod(Xtwo_AAorAT_fourier_test))
Xtwo_AAorTA_fourier_mag_test = as.matrix(Mod(Xtwo_AAorTA_fourier_test))
Xtwo_AAorTT_fourier_mag_test = as.matrix(Mod(Xtwo_AAorTT_fourier_test))
Xtwo_ATorTA_fourier_mag_test = as.matrix(Mod(Xtwo_ATorTA_fourier_test))
Xtwo_ATorTT_fourier_mag_test = as.matrix(Mod(Xtwo_ATorTT_fourier_test))
Xtwo_TAorTT_fourier_mag_test = as.matrix(Mod(Xtwo_TAorTT_fourier_test))
Xtwo_CCorCG_fourier_mag_test = as.matrix(Mod(Xtwo_CCorCG_fourier_test))
Xtwo_CCorGC_fourier_mag_test = as.matrix(Mod(Xtwo_CCorGC_fourier_test))
Xtwo_CCorGG_fourier_mag_test = as.matrix(Mod(Xtwo_CCorGG_fourier_test))
Xtwo_CGorGC_fourier_mag_test = as.matrix(Mod(Xtwo_CGorGC_fourier_test))
Xtwo_CGorGG_fourier_mag_test = as.matrix(Mod(Xtwo_CGorGG_fourier_test))
Xtwo_GCorGG_fourier_mag_test = as.matrix(Mod(Xtwo_GCorGG_fourier_test))

Xtwo_AAorAT_fourier_phase_test = as.matrix(Arg(Xtwo_AAorAT_fourier_test))
Xtwo_AAorTA_fourier_phase_test = as.matrix(Arg(Xtwo_AAorTA_fourier_test))
Xtwo_AAorTT_fourier_phase_test = as.matrix(Arg(Xtwo_AAorTT_fourier_test))
Xtwo_ATorTA_fourier_phase_test = as.matrix(Arg(Xtwo_ATorTA_fourier_test))
Xtwo_ATorTT_fourier_phase_test = as.matrix(Arg(Xtwo_ATorTT_fourier_test))
Xtwo_TAorTT_fourier_phase_test = as.matrix(Arg(Xtwo_TAorTT_fourier_test))
Xtwo_CCorCG_fourier_phase_test = as.matrix(Arg(Xtwo_CCorCG_fourier_test))
Xtwo_CCorGC_fourier_phase_test = as.matrix(Arg(Xtwo_CCorGC_fourier_test))
Xtwo_CCorGG_fourier_phase_test = as.matrix(Arg(Xtwo_CCorGG_fourier_test))
Xtwo_CGorGC_fourier_phase_test = as.matrix(Arg(Xtwo_CGorGC_fourier_test))
Xtwo_CGorGG_fourier_phase_test = as.matrix(Arg(Xtwo_CGorGG_fourier_test))
Xtwo_GCorGG_fourier_phase_test = as.matrix(Arg(Xtwo_GCorGG_fourier_test))


colnames(Xtwo_AAorAT_fourier_mag_test) = "fourier_mag_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_mag_test) = "fourier_mag_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_mag_test) = "fourier_mag_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_mag_test) = "fourier_mag_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_mag_test) = "fourier_mag_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_mag_test) = "fourier_mag_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_mag_test) = "fourier_mag_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_mag_test) = "fourier_mag_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_mag_test) = "fourier_mag_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_mag_test) = "fourier_mag_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_mag_test) = "fourier_mag_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_mag_test) = "fourier_mag_GCorGG_6"

colnames(Xtwo_AAorAT_fourier_phase_test) = "fourier_phase_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_phase_test) = "fourier_phase_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_phase_test) = "fourier_phase_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_phase_test) = "fourier_phase_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_phase_test) = "fourier_phase_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_phase_test) = "fourier_phase_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_phase_test) = "fourier_phase_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_phase_test) = "fourier_phase_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_phase_test) = "fourier_phase_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_phase_test) = "fourier_phase_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_phase_test) = "fourier_phase_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_phase_test) = "fourier_phase_GCorGG_6"


# Groups of Dinucleotides:

Xtwo_AAorATorTAorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorCGorGCorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AAorATorTAorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorCGorGCorGG_test) = colnames(Xtwo_test)
Xtwo_AAorATorTAorTT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "AT") | 
                                (Xtwo_test == "TA") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "CG") |
                                (Xtwo_test == "GC") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_fourier_test = t(apply(Xtwo_AAorATorTAorTT_test, 1, fft))[,6]
Xtwo_CCorCGorGCorGG_fourier_test = t(apply(Xtwo_CCorCGorGCorGG_test, 1, fft))[,6]

Xtwo_AAorATorTAorTT_fourier_mag_test = as.matrix(Mod(Xtwo_AAorATorTAorTT_fourier_test))
Xtwo_CCorCGorGCorGG_fourier_mag_test = as.matrix(Mod(Xtwo_CCorCGorGCorGG_fourier_test))

Xtwo_AAorATorTAorTT_fourier_phase_test = as.matrix(Arg(Xtwo_AAorATorTAorTT_fourier_test))
Xtwo_CCorCGorGCorGG_fourier_phase_test = as.matrix(Arg(Xtwo_CCorCGorGCorGG_fourier_test))


colnames(Xtwo_AAorATorTAorTT_fourier_mag_test) = "fourier_mag_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_mag_test) = "fourier_mag_CCorCGorGCorGG_6"

colnames(Xtwo_AAorATorTAorTT_fourier_phase_test) = "fourier_phase_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_phase_test) = "fourier_phase_CCorCGorGCorGG_6"


########################################################################
# Random
########################################################################

# Nucleotides:

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

# Xone_A_fourier_random_all = t(apply(Xone_A_random_all, 1, fft))[,1:26]
# Xone_C_fourier_random_all = t(apply(Xone_C_random_all, 1, fft))[,1:26]
# Xone_G_fourier_random_all = t(apply(Xone_G_random_all, 1, fft))[,1:26]
# Xone_T_fourier_random_all = t(apply(Xone_T_random_all, 1, fft))[,1:26]

Xone_A_fourier_random_all = t(apply(Xone_A_random_all, 1, fft))[,6]
Xone_C_fourier_random_all = t(apply(Xone_C_random_all, 1, fft))[,6]
Xone_G_fourier_random_all = t(apply(Xone_G_random_all, 1, fft))[,6]
Xone_T_fourier_random_all = t(apply(Xone_T_random_all, 1, fft))[,6]

Xone_A_fourier_mag_random_all = as.matrix(Mod(Xone_A_fourier_random_all))
Xone_C_fourier_mag_random_all = as.matrix(Mod(Xone_C_fourier_random_all))
Xone_G_fourier_mag_random_all = as.matrix(Mod(Xone_G_fourier_random_all))
Xone_T_fourier_mag_random_all = as.matrix(Mod(Xone_T_fourier_random_all))

Xone_A_fourier_phase_random_all = as.matrix(Arg(Xone_A_fourier_random_all))
Xone_C_fourier_phase_random_all = as.matrix(Arg(Xone_C_fourier_random_all))
Xone_G_fourier_phase_random_all = as.matrix(Arg(Xone_G_fourier_random_all))
Xone_T_fourier_phase_random_all = as.matrix(Arg(Xone_T_fourier_random_all))

# colnames(Xone_A_fourier_mag_random_all) = paste0("fourier_mag_A_", 1:26)
# colnames(Xone_C_fourier_mag_random_all) = paste0("fourier_mag_C_", 1:26)
# colnames(Xone_G_fourier_mag_random_all) = paste0("fourier_mag_G_", 1:26)
# colnames(Xone_T_fourier_mag_random_all) = paste0("fourier_mag_T_", 1:26)

# colnames(Xone_A_fourier_phase_random_all) = paste0("fourier_phase_A_", 1:26)
# colnames(Xone_C_fourier_phase_random_all) = paste0("fourier_phase_C_", 1:26)
# colnames(Xone_G_fourier_phase_random_all) = paste0("fourier_phase_G_", 1:26)
# colnames(Xone_T_fourier_phase_random_all) = paste0("fourier_phase_T_", 1:26)

colnames(Xone_A_fourier_mag_random_all) = "fourier_mag_A_6"
colnames(Xone_C_fourier_mag_random_all) = "fourier_mag_C_6"
colnames(Xone_G_fourier_mag_random_all) = "fourier_mag_G_6"
colnames(Xone_T_fourier_mag_random_all) = "fourier_mag_T_6"

colnames(Xone_A_fourier_phase_random_all) = "fourier_phase_A_6"
colnames(Xone_C_fourier_phase_random_all) = "fourier_phase_C_6"
colnames(Xone_G_fourier_phase_random_all) = "fourier_phase_G_6"
colnames(Xone_T_fourier_phase_random_all) = "fourier_phase_T_6"


# Pairs of nucleotides:

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

# Xone_AorC_fourier_random_all = t(apply(Xone_AorC_random_all, 1, fft))[,1:26]
# Xone_AorG_fourier_random_all = t(apply(Xone_AorG_random_all, 1, fft))[,1:26]
# Xone_AorT_fourier_random_all = t(apply(Xone_AorT_random_all, 1, fft))[,1:26]
# Xone_CorG_fourier_random_all = t(apply(Xone_CorG_random_all, 1, fft))[,1:26]
# Xone_CorT_fourier_random_all = t(apply(Xone_CorT_random_all, 1, fft))[,1:26]
# Xone_GorT_fourier_random_all = t(apply(Xone_GorT_random_all, 1, fft))[,1:26]

Xone_AorC_fourier_random_all = t(apply(Xone_AorC_random_all, 1, fft))[,6]
Xone_AorG_fourier_random_all = t(apply(Xone_AorG_random_all, 1, fft))[,6]
Xone_AorT_fourier_random_all = t(apply(Xone_AorT_random_all, 1, fft))[,6]
Xone_CorG_fourier_random_all = t(apply(Xone_CorG_random_all, 1, fft))[,6]
Xone_CorT_fourier_random_all = t(apply(Xone_CorT_random_all, 1, fft))[,6]
Xone_GorT_fourier_random_all = t(apply(Xone_GorT_random_all, 1, fft))[,6]

Xone_AorC_fourier_mag_random_all = as.matrix(Mod(Xone_AorC_fourier_random_all))
Xone_AorG_fourier_mag_random_all = as.matrix(Mod(Xone_AorG_fourier_random_all))
Xone_AorT_fourier_mag_random_all = as.matrix(Mod(Xone_AorT_fourier_random_all))
Xone_CorG_fourier_mag_random_all = as.matrix(Mod(Xone_CorG_fourier_random_all))
Xone_CorT_fourier_mag_random_all = as.matrix(Mod(Xone_CorT_fourier_random_all))
Xone_GorT_fourier_mag_random_all = as.matrix(Mod(Xone_GorT_fourier_random_all))

Xone_AorC_fourier_phase_random_all = as.matrix(Arg(Xone_AorC_fourier_random_all))
Xone_AorG_fourier_phase_random_all = as.matrix(Arg(Xone_AorG_fourier_random_all))
Xone_AorT_fourier_phase_random_all = as.matrix(Arg(Xone_AorT_fourier_random_all))
Xone_CorG_fourier_phase_random_all = as.matrix(Arg(Xone_CorG_fourier_random_all))
Xone_CorT_fourier_phase_random_all = as.matrix(Arg(Xone_CorT_fourier_random_all))
Xone_GorT_fourier_phase_random_all = as.matrix(Arg(Xone_GorT_fourier_random_all))

colnames(Xone_AorC_fourier_mag_random_all) = "fourier_mag_AorC_6"
colnames(Xone_AorG_fourier_mag_random_all) = "fourier_mag_AorG_6"
colnames(Xone_AorT_fourier_mag_random_all) = "fourier_mag_AorT_6"
colnames(Xone_CorG_fourier_mag_random_all) = "fourier_mag_CorG_6"
colnames(Xone_CorT_fourier_mag_random_all) = "fourier_mag_CorT_6"
colnames(Xone_GorT_fourier_mag_random_all) = "fourier_mag_GorT_6"

colnames(Xone_AorC_fourier_phase_random_all) = "fourier_phase_AorC_6"
colnames(Xone_AorG_fourier_phase_random_all) = "fourier_phase_AorG_6"
colnames(Xone_AorT_fourier_phase_random_all) = "fourier_phase_AorT_6"
colnames(Xone_CorG_fourier_phase_random_all) = "fourier_phase_CorG_6"
colnames(Xone_CorT_fourier_phase_random_all) = "fourier_phase_CorT_6"
colnames(Xone_GorT_fourier_phase_random_all) = "fourier_phase_GorT_6"


# Dinucleotides:

Xtwo_random_all = dat_random_all %>% select(all_of(ps2))
Xtwo_AA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_AC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_AG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_AT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_GA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_GC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_GG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_GT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_TA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_TC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_TG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_TT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
colnames(Xtwo_AA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_AC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_AG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_AT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_GA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_GC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_GG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_GT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_TA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_TC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_TG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_TT_random_all) = colnames(Xtwo_random_all)
Xtwo_AA_random_all[] = (Xtwo_random_all == "AA") %>% as.matrix() %>% as.numeric()
Xtwo_AC_random_all[] = (Xtwo_random_all == "AC") %>% as.matrix() %>% as.numeric()
Xtwo_AG_random_all[] = (Xtwo_random_all == "AG") %>% as.matrix() %>% as.numeric()
Xtwo_AT_random_all[] = (Xtwo_random_all == "AT") %>% as.matrix() %>% as.numeric()
Xtwo_CA_random_all[] = (Xtwo_random_all == "CA") %>% as.matrix() %>% as.numeric()
Xtwo_CC_random_all[] = (Xtwo_random_all == "CC") %>% as.matrix() %>% as.numeric()
Xtwo_CG_random_all[] = (Xtwo_random_all == "CG") %>% as.matrix() %>% as.numeric()
Xtwo_CT_random_all[] = (Xtwo_random_all == "CT") %>% as.matrix() %>% as.numeric()
Xtwo_GA_random_all[] = (Xtwo_random_all == "GA") %>% as.matrix() %>% as.numeric()
Xtwo_GC_random_all[] = (Xtwo_random_all == "GC") %>% as.matrix() %>% as.numeric()
Xtwo_GG_random_all[] = (Xtwo_random_all == "GG") %>% as.matrix() %>% as.numeric()
Xtwo_GT_random_all[] = (Xtwo_random_all == "GT") %>% as.matrix() %>% as.numeric()
Xtwo_TA_random_all[] = (Xtwo_random_all == "TA") %>% as.matrix() %>% as.numeric()
Xtwo_TC_random_all[] = (Xtwo_random_all == "TC") %>% as.matrix() %>% as.numeric()
Xtwo_TG_random_all[] = (Xtwo_random_all == "TG") %>% as.matrix() %>% as.numeric()
Xtwo_TT_random_all[] = (Xtwo_random_all == "TT") %>% as.matrix() %>% as.numeric()

# Xtwo_AA_fourier_random_all = t(apply(Xtwo_AA_random_all, 1, fft))[,1:25]
# Xtwo_AC_fourier_random_all = t(apply(Xtwo_AC_random_all, 1, fft))[,1:25]
# Xtwo_AG_fourier_random_all = t(apply(Xtwo_AG_random_all, 1, fft))[,1:25]
# Xtwo_AT_fourier_random_all = t(apply(Xtwo_AT_random_all, 1, fft))[,1:25]
# Xtwo_CA_fourier_random_all = t(apply(Xtwo_CA_random_all, 1, fft))[,1:25]
# Xtwo_CC_fourier_random_all = t(apply(Xtwo_CC_random_all, 1, fft))[,1:25]
# Xtwo_CG_fourier_random_all = t(apply(Xtwo_CG_random_all, 1, fft))[,1:25]
# Xtwo_CT_fourier_random_all = t(apply(Xtwo_CT_random_all, 1, fft))[,1:25]
# Xtwo_GA_fourier_random_all = t(apply(Xtwo_GA_random_all, 1, fft))[,1:25]
# Xtwo_GC_fourier_random_all = t(apply(Xtwo_GC_random_all, 1, fft))[,1:25]
# Xtwo_GG_fourier_random_all = t(apply(Xtwo_GG_random_all, 1, fft))[,1:25]
# Xtwo_GT_fourier_random_all = t(apply(Xtwo_GT_random_all, 1, fft))[,1:25]
# Xtwo_TA_fourier_random_all = t(apply(Xtwo_TA_random_all, 1, fft))[,1:25]
# Xtwo_TC_fourier_random_all = t(apply(Xtwo_TC_random_all, 1, fft))[,1:25]
# Xtwo_TG_fourier_random_all = t(apply(Xtwo_TG_random_all, 1, fft))[,1:25]
# Xtwo_TT_fourier_random_all = t(apply(Xtwo_TT_random_all, 1, fft))[,1:25]

Xtwo_AA_fourier_random_all = t(apply(Xtwo_AA_random_all, 1, fft))[,6]
Xtwo_AC_fourier_random_all = t(apply(Xtwo_AC_random_all, 1, fft))[,6]
Xtwo_AG_fourier_random_all = t(apply(Xtwo_AG_random_all, 1, fft))[,6]
Xtwo_AT_fourier_random_all = t(apply(Xtwo_AT_random_all, 1, fft))[,6]
Xtwo_CA_fourier_random_all = t(apply(Xtwo_CA_random_all, 1, fft))[,6]
Xtwo_CC_fourier_random_all = t(apply(Xtwo_CC_random_all, 1, fft))[,6]
Xtwo_CG_fourier_random_all = t(apply(Xtwo_CG_random_all, 1, fft))[,6]
Xtwo_CT_fourier_random_all = t(apply(Xtwo_CT_random_all, 1, fft))[,6]
Xtwo_GA_fourier_random_all = t(apply(Xtwo_GA_random_all, 1, fft))[,6]
Xtwo_GC_fourier_random_all = t(apply(Xtwo_GC_random_all, 1, fft))[,6]
Xtwo_GG_fourier_random_all = t(apply(Xtwo_GG_random_all, 1, fft))[,6]
Xtwo_GT_fourier_random_all = t(apply(Xtwo_GT_random_all, 1, fft))[,6]
Xtwo_TA_fourier_random_all = t(apply(Xtwo_TA_random_all, 1, fft))[,6]
Xtwo_TC_fourier_random_all = t(apply(Xtwo_TC_random_all, 1, fft))[,6]
Xtwo_TG_fourier_random_all = t(apply(Xtwo_TG_random_all, 1, fft))[,6]
Xtwo_TT_fourier_random_all = t(apply(Xtwo_TT_random_all, 1, fft))[,6]

Xtwo_AA_fourier_mag_random_all = as.matrix(Mod(Xtwo_AA_fourier_random_all))
Xtwo_AC_fourier_mag_random_all = as.matrix(Mod(Xtwo_AC_fourier_random_all))
Xtwo_AG_fourier_mag_random_all = as.matrix(Mod(Xtwo_AG_fourier_random_all))
Xtwo_AT_fourier_mag_random_all = as.matrix(Mod(Xtwo_AT_fourier_random_all))
Xtwo_CA_fourier_mag_random_all = as.matrix(Mod(Xtwo_CA_fourier_random_all))
Xtwo_CC_fourier_mag_random_all = as.matrix(Mod(Xtwo_CC_fourier_random_all))
Xtwo_CG_fourier_mag_random_all = as.matrix(Mod(Xtwo_CG_fourier_random_all))
Xtwo_CT_fourier_mag_random_all = as.matrix(Mod(Xtwo_CT_fourier_random_all))
Xtwo_GA_fourier_mag_random_all = as.matrix(Mod(Xtwo_GA_fourier_random_all))
Xtwo_GC_fourier_mag_random_all = as.matrix(Mod(Xtwo_GC_fourier_random_all))
Xtwo_GG_fourier_mag_random_all = as.matrix(Mod(Xtwo_GG_fourier_random_all))
Xtwo_GT_fourier_mag_random_all = as.matrix(Mod(Xtwo_GT_fourier_random_all))
Xtwo_TA_fourier_mag_random_all = as.matrix(Mod(Xtwo_TA_fourier_random_all))
Xtwo_TC_fourier_mag_random_all = as.matrix(Mod(Xtwo_TC_fourier_random_all))
Xtwo_TG_fourier_mag_random_all = as.matrix(Mod(Xtwo_TG_fourier_random_all))
Xtwo_TT_fourier_mag_random_all = as.matrix(Mod(Xtwo_TT_fourier_random_all))

Xtwo_AA_fourier_phase_random_all = as.matrix(Arg(Xtwo_AA_fourier_random_all))
Xtwo_AC_fourier_phase_random_all = as.matrix(Arg(Xtwo_AC_fourier_random_all))
Xtwo_AG_fourier_phase_random_all = as.matrix(Arg(Xtwo_AG_fourier_random_all))
Xtwo_AT_fourier_phase_random_all = as.matrix(Arg(Xtwo_AT_fourier_random_all))
Xtwo_CA_fourier_phase_random_all = as.matrix(Arg(Xtwo_CA_fourier_random_all))
Xtwo_CC_fourier_phase_random_all = as.matrix(Arg(Xtwo_CC_fourier_random_all))
Xtwo_CG_fourier_phase_random_all = as.matrix(Arg(Xtwo_CG_fourier_random_all))
Xtwo_CT_fourier_phase_random_all = as.matrix(Arg(Xtwo_CT_fourier_random_all))
Xtwo_GA_fourier_phase_random_all = as.matrix(Arg(Xtwo_GA_fourier_random_all))
Xtwo_GC_fourier_phase_random_all = as.matrix(Arg(Xtwo_GC_fourier_random_all))
Xtwo_GG_fourier_phase_random_all = as.matrix(Arg(Xtwo_GG_fourier_random_all))
Xtwo_GT_fourier_phase_random_all = as.matrix(Arg(Xtwo_GT_fourier_random_all))
Xtwo_TA_fourier_phase_random_all = as.matrix(Arg(Xtwo_TA_fourier_random_all))
Xtwo_TC_fourier_phase_random_all = as.matrix(Arg(Xtwo_TC_fourier_random_all))
Xtwo_TG_fourier_phase_random_all = as.matrix(Arg(Xtwo_TG_fourier_random_all))
Xtwo_TT_fourier_phase_random_all = as.matrix(Arg(Xtwo_TT_fourier_random_all))

# colnames(Xtwo_AA_fourier_mag_random_all) = paste0("fourier_mag_AA_", 1:25)
# colnames(Xtwo_AC_fourier_mag_random_all) = paste0("fourier_mag_AC_", 1:25)
# colnames(Xtwo_AG_fourier_mag_random_all) = paste0("fourier_mag_AG_", 1:25)
# colnames(Xtwo_AT_fourier_mag_random_all) = paste0("fourier_mag_AT_", 1:25)
# colnames(Xtwo_CA_fourier_mag_random_all) = paste0("fourier_mag_CA_", 1:25)
# colnames(Xtwo_CC_fourier_mag_random_all) = paste0("fourier_mag_CC_", 1:25)
# colnames(Xtwo_CG_fourier_mag_random_all) = paste0("fourier_mag_CG_", 1:25)
# colnames(Xtwo_CT_fourier_mag_random_all) = paste0("fourier_mag_CT_", 1:25)
# colnames(Xtwo_GA_fourier_mag_random_all) = paste0("fourier_mag_GA_", 1:25)
# colnames(Xtwo_GC_fourier_mag_random_all) = paste0("fourier_mag_GC_", 1:25)
# colnames(Xtwo_GG_fourier_mag_random_all) = paste0("fourier_mag_GG_", 1:25)
# colnames(Xtwo_GT_fourier_mag_random_all) = paste0("fourier_mag_GT_", 1:25)
# colnames(Xtwo_TA_fourier_mag_random_all) = paste0("fourier_mag_TA_", 1:25)
# colnames(Xtwo_TC_fourier_mag_random_all) = paste0("fourier_mag_TC_", 1:25)
# colnames(Xtwo_TG_fourier_mag_random_all) = paste0("fourier_mag_TG_", 1:25)
# colnames(Xtwo_TT_fourier_mag_random_all) = paste0("fourier_mag_TT_", 1:25)
# 
# colnames(Xtwo_AA_fourier_phase_random_all) = paste0("fourier_phase_AA_", 1:25)
# colnames(Xtwo_AC_fourier_phase_random_all) = paste0("fourier_phase_AC_", 1:25)
# colnames(Xtwo_AG_fourier_phase_random_all) = paste0("fourier_phase_AG_", 1:25)
# colnames(Xtwo_AT_fourier_phase_random_all) = paste0("fourier_phase_AT_", 1:25)
# colnames(Xtwo_CA_fourier_phase_random_all) = paste0("fourier_phase_CA_", 1:25)
# colnames(Xtwo_CC_fourier_phase_random_all) = paste0("fourier_phase_CC_", 1:25)
# colnames(Xtwo_CG_fourier_phase_random_all) = paste0("fourier_phase_CG_", 1:25)
# colnames(Xtwo_CT_fourier_phase_random_all) = paste0("fourier_phase_CT_", 1:25)
# colnames(Xtwo_GA_fourier_phase_random_all) = paste0("fourier_phase_GA_", 1:25)
# colnames(Xtwo_GC_fourier_phase_random_all) = paste0("fourier_phase_GC_", 1:25)
# colnames(Xtwo_GG_fourier_phase_random_all) = paste0("fourier_phase_GG_", 1:25)
# colnames(Xtwo_GT_fourier_phase_random_all) = paste0("fourier_phase_GT_", 1:25)
# colnames(Xtwo_TA_fourier_phase_random_all) = paste0("fourier_phase_TA_", 1:25)
# colnames(Xtwo_TC_fourier_phase_random_all) = paste0("fourier_phase_TC_", 1:25)
# colnames(Xtwo_TG_fourier_phase_random_all) = paste0("fourier_phase_TG_", 1:25)
# colnames(Xtwo_TT_fourier_phase_random_all) = paste0("fourier_phase_TT_", 1:25)

colnames(Xtwo_AA_fourier_mag_random_all) = "fourier_mag_AA_6"
colnames(Xtwo_AC_fourier_mag_random_all) = "fourier_mag_AC_6"
colnames(Xtwo_AG_fourier_mag_random_all) = "fourier_mag_AG_6"
colnames(Xtwo_AT_fourier_mag_random_all) = "fourier_mag_AT_6"
colnames(Xtwo_CA_fourier_mag_random_all) = "fourier_mag_CA_6"
colnames(Xtwo_CC_fourier_mag_random_all) = "fourier_mag_CC_6"
colnames(Xtwo_CG_fourier_mag_random_all) = "fourier_mag_CG_6"
colnames(Xtwo_CT_fourier_mag_random_all) = "fourier_mag_CT_6"
colnames(Xtwo_GA_fourier_mag_random_all) = "fourier_mag_GA_6"
colnames(Xtwo_GC_fourier_mag_random_all) = "fourier_mag_GC_6"
colnames(Xtwo_GG_fourier_mag_random_all) = "fourier_mag_GG_6"
colnames(Xtwo_GT_fourier_mag_random_all) = "fourier_mag_GT_6"
colnames(Xtwo_TA_fourier_mag_random_all) = "fourier_mag_TA_6"
colnames(Xtwo_TC_fourier_mag_random_all) = "fourier_mag_TC_6"
colnames(Xtwo_TG_fourier_mag_random_all) = "fourier_mag_TG_6"
colnames(Xtwo_TT_fourier_mag_random_all) = "fourier_mag_TT_6"

colnames(Xtwo_AA_fourier_phase_random_all) = "fourier_phase_AA_6"
colnames(Xtwo_AC_fourier_phase_random_all) = "fourier_phase_AC_6"
colnames(Xtwo_AG_fourier_phase_random_all) = "fourier_phase_AG_6"
colnames(Xtwo_AT_fourier_phase_random_all) = "fourier_phase_AT_6"
colnames(Xtwo_CA_fourier_phase_random_all) = "fourier_phase_CA_6"
colnames(Xtwo_CC_fourier_phase_random_all) = "fourier_phase_CC_6"
colnames(Xtwo_CG_fourier_phase_random_all) = "fourier_phase_CG_6"
colnames(Xtwo_CT_fourier_phase_random_all) = "fourier_phase_CT_6"
colnames(Xtwo_GA_fourier_phase_random_all) = "fourier_phase_GA_6"
colnames(Xtwo_GC_fourier_phase_random_all) = "fourier_phase_GC_6"
colnames(Xtwo_GG_fourier_phase_random_all) = "fourier_phase_GG_6"
colnames(Xtwo_GT_fourier_phase_random_all) = "fourier_phase_GT_6"
colnames(Xtwo_TA_fourier_phase_random_all) = "fourier_phase_TA_6"
colnames(Xtwo_TC_fourier_phase_random_all) = "fourier_phase_TC_6"
colnames(Xtwo_TG_fourier_phase_random_all) = "fourier_phase_TG_6"
colnames(Xtwo_TT_fourier_phase_random_all) = "fourier_phase_TT_6"


# Pairs of Dinucleotides:

Xtwo_AAorAT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_AAorTA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_AAorTT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_ATorTA_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_ATorTT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_TAorTT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CCorCG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CCorGC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CCorGG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CGorGC_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CGorGG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_GCorGG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
colnames(Xtwo_AAorAT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_AAorTA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_AAorTT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_ATorTA_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_ATorTT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_TAorTT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CCorCG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CCorGC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CCorGG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CGorGC_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CGorGG_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_GCorGG_random_all) = colnames(Xtwo_random_all)
Xtwo_AAorAT_random_all[] = ((Xtwo_random_all == "AA") | (Xtwo_random_all == "AT")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTA_random_all[] = ((Xtwo_random_all == "AA") | (Xtwo_random_all == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_AAorTT_random_all[] = ((Xtwo_random_all == "AA") | (Xtwo_random_all == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTA_random_all[] = ((Xtwo_random_all == "AT") | (Xtwo_random_all == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_ATorTT_random_all[] = ((Xtwo_random_all == "AT") | (Xtwo_random_all == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_TAorTT_random_all[] = ((Xtwo_random_all == "TA") | (Xtwo_random_all == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCG_random_all[] = ((Xtwo_random_all == "CC") | (Xtwo_random_all == "CG")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGC_random_all[] = ((Xtwo_random_all == "CC") | (Xtwo_random_all == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorGG_random_all[] = ((Xtwo_random_all == "CC") | (Xtwo_random_all == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGC_random_all[] = ((Xtwo_random_all == "CG") | (Xtwo_random_all == "GC")) %>% as.matrix() %>% as.numeric()
Xtwo_CGorGG_random_all[] = ((Xtwo_random_all == "CG") | (Xtwo_random_all == "GG")) %>% as.matrix() %>% as.numeric()
Xtwo_GCorGG_random_all[] = ((Xtwo_random_all == "GC") | (Xtwo_random_all == "GG")) %>% as.matrix() %>% as.numeric()

# Xtwo_AAorAT_fourier_random_all = t(apply(Xtwo_AAorAT_random_all, 1, fft))[,1:25]
# Xtwo_AAorTA_fourier_random_all = t(apply(Xtwo_AAorTA_random_all, 1, fft))[,1:25]
# Xtwo_AAorTT_fourier_random_all = t(apply(Xtwo_AAorTT_random_all, 1, fft))[,1:25]
# Xtwo_ATorTA_fourier_random_all = t(apply(Xtwo_ATorTA_random_all, 1, fft))[,1:25]
# Xtwo_ATorTT_fourier_random_all = t(apply(Xtwo_ATorTT_random_all, 1, fft))[,1:25]
# Xtwo_TAorTT_fourier_random_all = t(apply(Xtwo_TAorTT_random_all, 1, fft))[,1:25]
# Xtwo_CCorCG_fourier_random_all = t(apply(Xtwo_CCorCG_random_all, 1, fft))[,1:25]
# Xtwo_CCorGC_fourier_random_all = t(apply(Xtwo_CCorGC_random_all, 1, fft))[,1:25]
# Xtwo_CCorGG_fourier_random_all = t(apply(Xtwo_CCorGG_random_all, 1, fft))[,1:25]
# Xtwo_CGorGC_fourier_random_all = t(apply(Xtwo_CGorGC_random_all, 1, fft))[,1:25]
# Xtwo_CGorGG_fourier_random_all = t(apply(Xtwo_CGorGG_random_all, 1, fft))[,1:25]
# Xtwo_GCorGG_fourier_random_all = t(apply(Xtwo_GCorGG_random_all, 1, fft))[,1:25]

Xtwo_AAorAT_fourier_random_all = t(apply(Xtwo_AAorAT_random_all, 1, fft))[,6]
Xtwo_AAorTA_fourier_random_all = t(apply(Xtwo_AAorTA_random_all, 1, fft))[,6]
Xtwo_AAorTT_fourier_random_all = t(apply(Xtwo_AAorTT_random_all, 1, fft))[,6]
Xtwo_ATorTA_fourier_random_all = t(apply(Xtwo_ATorTA_random_all, 1, fft))[,6]
Xtwo_ATorTT_fourier_random_all = t(apply(Xtwo_ATorTT_random_all, 1, fft))[,6]
Xtwo_TAorTT_fourier_random_all = t(apply(Xtwo_TAorTT_random_all, 1, fft))[,6]
Xtwo_CCorCG_fourier_random_all = t(apply(Xtwo_CCorCG_random_all, 1, fft))[,6]
Xtwo_CCorGC_fourier_random_all = t(apply(Xtwo_CCorGC_random_all, 1, fft))[,6]
Xtwo_CCorGG_fourier_random_all = t(apply(Xtwo_CCorGG_random_all, 1, fft))[,6]
Xtwo_CGorGC_fourier_random_all = t(apply(Xtwo_CGorGC_random_all, 1, fft))[,6]
Xtwo_CGorGG_fourier_random_all = t(apply(Xtwo_CGorGG_random_all, 1, fft))[,6]
Xtwo_GCorGG_fourier_random_all = t(apply(Xtwo_GCorGG_random_all, 1, fft))[,6]

Xtwo_AAorAT_fourier_mag_random_all = as.matrix(Mod(Xtwo_AAorAT_fourier_random_all))
Xtwo_AAorTA_fourier_mag_random_all = as.matrix(Mod(Xtwo_AAorTA_fourier_random_all))
Xtwo_AAorTT_fourier_mag_random_all = as.matrix(Mod(Xtwo_AAorTT_fourier_random_all))
Xtwo_ATorTA_fourier_mag_random_all = as.matrix(Mod(Xtwo_ATorTA_fourier_random_all))
Xtwo_ATorTT_fourier_mag_random_all = as.matrix(Mod(Xtwo_ATorTT_fourier_random_all))
Xtwo_TAorTT_fourier_mag_random_all = as.matrix(Mod(Xtwo_TAorTT_fourier_random_all))
Xtwo_CCorCG_fourier_mag_random_all = as.matrix(Mod(Xtwo_CCorCG_fourier_random_all))
Xtwo_CCorGC_fourier_mag_random_all = as.matrix(Mod(Xtwo_CCorGC_fourier_random_all))
Xtwo_CCorGG_fourier_mag_random_all = as.matrix(Mod(Xtwo_CCorGG_fourier_random_all))
Xtwo_CGorGC_fourier_mag_random_all = as.matrix(Mod(Xtwo_CGorGC_fourier_random_all))
Xtwo_CGorGG_fourier_mag_random_all = as.matrix(Mod(Xtwo_CGorGG_fourier_random_all))
Xtwo_GCorGG_fourier_mag_random_all = as.matrix(Mod(Xtwo_GCorGG_fourier_random_all))

Xtwo_AAorAT_fourier_phase_random_all = as.matrix(Arg(Xtwo_AAorAT_fourier_random_all))
Xtwo_AAorTA_fourier_phase_random_all = as.matrix(Arg(Xtwo_AAorTA_fourier_random_all))
Xtwo_AAorTT_fourier_phase_random_all = as.matrix(Arg(Xtwo_AAorTT_fourier_random_all))
Xtwo_ATorTA_fourier_phase_random_all = as.matrix(Arg(Xtwo_ATorTA_fourier_random_all))
Xtwo_ATorTT_fourier_phase_random_all = as.matrix(Arg(Xtwo_ATorTT_fourier_random_all))
Xtwo_TAorTT_fourier_phase_random_all = as.matrix(Arg(Xtwo_TAorTT_fourier_random_all))
Xtwo_CCorCG_fourier_phase_random_all = as.matrix(Arg(Xtwo_CCorCG_fourier_random_all))
Xtwo_CCorGC_fourier_phase_random_all = as.matrix(Arg(Xtwo_CCorGC_fourier_random_all))
Xtwo_CCorGG_fourier_phase_random_all = as.matrix(Arg(Xtwo_CCorGG_fourier_random_all))
Xtwo_CGorGC_fourier_phase_random_all = as.matrix(Arg(Xtwo_CGorGC_fourier_random_all))
Xtwo_CGorGG_fourier_phase_random_all = as.matrix(Arg(Xtwo_CGorGG_fourier_random_all))
Xtwo_GCorGG_fourier_phase_random_all = as.matrix(Arg(Xtwo_GCorGG_fourier_random_all))


colnames(Xtwo_AAorAT_fourier_mag_random_all) = "fourier_mag_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_mag_random_all) = "fourier_mag_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_mag_random_all) = "fourier_mag_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_mag_random_all) = "fourier_mag_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_mag_random_all) = "fourier_mag_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_mag_random_all) = "fourier_mag_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_mag_random_all) = "fourier_mag_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_mag_random_all) = "fourier_mag_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_mag_random_all) = "fourier_mag_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_mag_random_all) = "fourier_mag_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_mag_random_all) = "fourier_mag_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_mag_random_all) = "fourier_mag_GCorGG_6"

colnames(Xtwo_AAorAT_fourier_phase_random_all) = "fourier_phase_AAorAT_6"
colnames(Xtwo_AAorTA_fourier_phase_random_all) = "fourier_phase_AAorTA_6"
colnames(Xtwo_AAorTT_fourier_phase_random_all) = "fourier_phase_AAorTT_6"
colnames(Xtwo_ATorTA_fourier_phase_random_all) = "fourier_phase_ATorTA_6"
colnames(Xtwo_ATorTT_fourier_phase_random_all) = "fourier_phase_ATorTT_6"
colnames(Xtwo_TAorTT_fourier_phase_random_all) = "fourier_phase_TAorTT_6"
colnames(Xtwo_CCorCG_fourier_phase_random_all) = "fourier_phase_CCorCG_6"
colnames(Xtwo_CCorGC_fourier_phase_random_all) = "fourier_phase_CCorGC_6"
colnames(Xtwo_CCorGG_fourier_phase_random_all) = "fourier_phase_CCorGG_6"
colnames(Xtwo_CGorGC_fourier_phase_random_all) = "fourier_phase_CGorGC_6"
colnames(Xtwo_CGorGG_fourier_phase_random_all) = "fourier_phase_CGorGG_6"
colnames(Xtwo_GCorGG_fourier_phase_random_all) = "fourier_phase_GCorGG_6"


# Groups of Dinucleotides:

Xtwo_AAorATorTAorTT_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
Xtwo_CCorCGorGCorGG_random_all = matrix(nrow=nrow(Xtwo_random_all), ncol=ncol(Xtwo_random_all))
colnames(Xtwo_AAorATorTAorTT_random_all) = colnames(Xtwo_random_all)
colnames(Xtwo_CCorCGorGCorGG_random_all) = colnames(Xtwo_random_all)
Xtwo_AAorATorTAorTT_random_all[] = ((Xtwo_random_all == "AA") | (Xtwo_random_all == "AT") | 
                                      (Xtwo_random_all == "TA") | (Xtwo_random_all == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG_random_all[] = ((Xtwo_random_all == "CC") | (Xtwo_random_all == "CG") |
                                      (Xtwo_random_all == "GC") | (Xtwo_random_all == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_fourier_random_all = t(apply(Xtwo_AAorATorTAorTT_random_all, 1, fft))[,6]
Xtwo_CCorCGorGCorGG_fourier_random_all = t(apply(Xtwo_CCorCGorGCorGG_random_all, 1, fft))[,6]

Xtwo_AAorATorTAorTT_fourier_mag_random_all = as.matrix(Mod(Xtwo_AAorATorTAorTT_fourier_random_all))
Xtwo_CCorCGorGCorGG_fourier_mag_random_all = as.matrix(Mod(Xtwo_CCorCGorGCorGG_fourier_random_all))

Xtwo_AAorATorTAorTT_fourier_phase_random_all = as.matrix(Arg(Xtwo_AAorATorTAorTT_fourier_random_all))
Xtwo_CCorCGorGCorGG_fourier_phase_random_all = as.matrix(Arg(Xtwo_CCorCGorGCorGG_fourier_random_all))


colnames(Xtwo_AAorATorTAorTT_fourier_mag_random_all) = "fourier_mag_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_mag_random_all) = "fourier_mag_CCorCGorGCorGG_6"

colnames(Xtwo_AAorATorTAorTT_fourier_phase_random_all) = "fourier_phase_AAorATorTAorTT_6"
colnames(Xtwo_CCorCGorGCorGG_fourier_phase_random_all) = "fourier_phase_CCorCGorGCorGG_6"



########################################################################
# Modelling
########################################################################


# Nucleotides Only
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.471685
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.4555443
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4656087
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9801271
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7113047
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6664056
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()
dat_temp.feature_importance



# Nucleotides + Pairs of Nucleotides
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorC_fourier_mag, Xone_AorG_fourier_mag, Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag, Xone_CorT_fourier_mag, Xone_GorT_fourier_mag,
                 Xone_AorC_fourier_phase, Xone_AorG_fourier_phase, Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, Xone_CorT_fourier_phase, Xone_GorT_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorC_fourier_mag_test, Xone_AorG_fourier_mag_test, Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test, Xone_CorT_fourier_mag_test, Xone_GorT_fourier_mag_test,
                      Xone_AorC_fourier_phase_test, Xone_AorG_fourier_phase_test, Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, Xone_CorT_fourier_phase_test, Xone_GorT_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorC_fourier_mag_random_all, Xone_AorG_fourier_mag_random_all, Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all, Xone_CorT_fourier_mag_random_all, Xone_GorT_fourier_mag_random_all,
                            Xone_AorC_fourier_phase_random_all, Xone_AorG_fourier_phase_random_all, Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, Xone_CorT_fourier_phase_random_all, Xone_GorT_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5594003
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5496307
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4656087
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9297395
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7229777
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6848468
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()
dat_temp.feature_importance


# Dinucleotides only
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xtwo_AA_fourier_mag, Xtwo_AC_fourier_mag, Xtwo_AG_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CA_fourier_mag, Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag, Xtwo_CT_fourier_mag,
                 Xtwo_GA_fourier_mag, Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag, Xtwo_GT_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TC_fourier_mag, Xtwo_TG_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AC_fourier_phase, Xtwo_AG_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CA_fourier_phase, Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase, Xtwo_CT_fourier_phase,
                 Xtwo_GA_fourier_phase, Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase, Xtwo_GT_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TC_fourier_phase, Xtwo_TG_fourier_phase, Xtwo_TT_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xtwo_AA_fourier_mag_test, Xtwo_AC_fourier_mag_test, Xtwo_AG_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CA_fourier_mag_test, Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test, Xtwo_CT_fourier_mag_test,
                      Xtwo_GA_fourier_mag_test, Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test, Xtwo_GT_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TC_fourier_mag_test, Xtwo_TG_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AC_fourier_phase_test, Xtwo_AG_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CA_fourier_phase_test, Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test, Xtwo_CT_fourier_phase_test,
                      Xtwo_GA_fourier_phase_test, Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test, Xtwo_GT_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TC_fourier_phase_test, Xtwo_TG_fourier_phase_test, Xtwo_TT_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AC_fourier_mag_random_all, Xtwo_AG_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CA_fourier_mag_random_all, Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all, Xtwo_CT_fourier_mag_random_all,
                            Xtwo_GA_fourier_mag_random_all, Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all, Xtwo_GT_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TC_fourier_mag_random_all, Xtwo_TG_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AC_fourier_phase_random_all, Xtwo_AG_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CA_fourier_phase_random_all, Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all, Xtwo_CT_fourier_phase_random_all,
                            Xtwo_GA_fourier_phase_random_all, Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all, Xtwo_GT_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TC_fourier_phase_random_all, Xtwo_TG_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all)


temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.4593343
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.4465066
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4661394
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9802761
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7075332
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6402811
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()
dat_temp.feature_importance



# Nucleotides + Dinucleotides
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xtwo_AA_fourier_mag, Xtwo_AC_fourier_mag, Xtwo_AG_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CA_fourier_mag, Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag, Xtwo_CT_fourier_mag,
                 Xtwo_GA_fourier_mag, Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag, Xtwo_GT_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TC_fourier_mag, Xtwo_TG_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AC_fourier_phase, Xtwo_AG_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CA_fourier_phase, Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase, Xtwo_CT_fourier_phase,
                 Xtwo_GA_fourier_phase, Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase, Xtwo_GT_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TC_fourier_phase, Xtwo_TG_fourier_phase, Xtwo_TT_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xtwo_AA_fourier_mag_test, Xtwo_AC_fourier_mag_test, Xtwo_AG_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CA_fourier_mag_test, Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test, Xtwo_CT_fourier_mag_test,
                      Xtwo_GA_fourier_mag_test, Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test, Xtwo_GT_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TC_fourier_mag_test, Xtwo_TG_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AC_fourier_phase_test, Xtwo_AG_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CA_fourier_phase_test, Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test, Xtwo_CT_fourier_phase_test,
                      Xtwo_GA_fourier_phase_test, Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test, Xtwo_GT_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TC_fourier_phase_test, Xtwo_TG_fourier_phase_test, Xtwo_TT_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AC_fourier_mag_random_all, Xtwo_AG_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CA_fourier_mag_random_all, Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all, Xtwo_CT_fourier_mag_random_all,
                            Xtwo_GA_fourier_mag_random_all, Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all, Xtwo_GT_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TC_fourier_mag_random_all, Xtwo_TG_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AC_fourier_phase_random_all, Xtwo_AG_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CA_fourier_phase_random_all, Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all, Xtwo_CT_fourier_phase_random_all,
                            Xtwo_GA_fourier_phase_random_all, Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all, Xtwo_GT_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TC_fourier_phase_random_all, Xtwo_TG_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.49626
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.4835217
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.4813084
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9637695
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7170481
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6604852
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()
dat_temp.feature_importance


# Nucleotides + Pairs of Nucleotides + Dinucleotides
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorC_fourier_mag, Xone_AorG_fourier_mag, Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag, Xone_CorT_fourier_mag, Xone_GorT_fourier_mag,
                 Xone_AorC_fourier_phase, Xone_AorG_fourier_phase, Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, Xone_CorT_fourier_phase, Xone_GorT_fourier_phase,
                 Xtwo_AA_fourier_mag, Xtwo_AC_fourier_mag, Xtwo_AG_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CA_fourier_mag, Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag, Xtwo_CT_fourier_mag,
                 Xtwo_GA_fourier_mag, Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag, Xtwo_GT_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TC_fourier_mag, Xtwo_TG_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AC_fourier_phase, Xtwo_AG_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CA_fourier_phase, Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase, Xtwo_CT_fourier_phase,
                 Xtwo_GA_fourier_phase, Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase, Xtwo_GT_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TC_fourier_phase, Xtwo_TG_fourier_phase, Xtwo_TT_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorC_fourier_mag_test, Xone_AorG_fourier_mag_test, Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test, Xone_CorT_fourier_mag_test, Xone_GorT_fourier_mag_test,
                      Xone_AorC_fourier_phase_test, Xone_AorG_fourier_phase_test, Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, Xone_CorT_fourier_phase_test, Xone_GorT_fourier_phase_test,
                      Xtwo_AA_fourier_mag_test, Xtwo_AC_fourier_mag_test, Xtwo_AG_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CA_fourier_mag_test, Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test, Xtwo_CT_fourier_mag_test,
                      Xtwo_GA_fourier_mag_test, Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test, Xtwo_GT_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TC_fourier_mag_test, Xtwo_TG_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AC_fourier_phase_test, Xtwo_AG_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CA_fourier_phase_test, Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test, Xtwo_CT_fourier_phase_test,
                      Xtwo_GA_fourier_phase_test, Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test, Xtwo_GT_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TC_fourier_phase_test, Xtwo_TG_fourier_phase_test, Xtwo_TT_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorC_fourier_mag_random_all, Xone_AorG_fourier_mag_random_all, Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all, Xone_CorT_fourier_mag_random_all, Xone_GorT_fourier_mag_random_all,
                            Xone_AorC_fourier_phase_random_all, Xone_AorG_fourier_phase_random_all, Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, Xone_CorT_fourier_phase_random_all, Xone_GorT_fourier_phase_random_all,
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AC_fourier_mag_random_all, Xtwo_AG_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CA_fourier_mag_random_all, Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all, Xtwo_CT_fourier_mag_random_all,
                            Xtwo_GA_fourier_mag_random_all, Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all, Xtwo_GT_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TC_fourier_mag_random_all, Xtwo_TG_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AC_fourier_phase_random_all, Xtwo_AG_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CA_fourier_phase_random_all, Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all, Xtwo_CT_fourier_phase_random_all,
                            Xtwo_GA_fourier_phase_random_all, Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all, Xtwo_GT_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TC_fourier_phase_random_all, Xtwo_TG_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5795711
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.571136
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5411141
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9599999
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7339418
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6838103
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw()
dat_temp.feature_importance



# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides
dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase,
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all,
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5773291
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5697716
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5385766
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9706429
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7408798
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6898604
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides, LightGBM")
dat_temp.feature_importance




# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides

dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5775296
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5697376
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5385988
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9676123
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7507324
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.694905
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, Random, LightGBM")
abline(0,1)

high_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) > 1.5
plot(dat_temp.lgb_random_all_pred[high_resid_ids_random_all], y_random_all[high_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() +
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides, LightGBM")
dat_temp.feature_importance




# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + npolyA/T

dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase,
                 npolyA, npolyT, npolyAorT)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test,
                      npolyA_test, npolyT_test, npolyAorT_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all,
                            npolyA_random_all, npolyT_random_all, npolyAorT_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5942835
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5901965
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.546724
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9725579
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7535151
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6985298
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + npolyA/T, LightGBM")
dat_temp.feature_importance



# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + max polyA/T

dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase,
                 dat_max_poly)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test,
                      dat_max_poly_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all,
                            dat_max_poly_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5810199
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5756244
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.540403
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9599137
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7524906
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6932112
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + max polyA/T, LightGBM")
dat_temp.feature_importance



# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + all polyA/T

dat_temp = cbind(dat %>% 
                   select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                          all_of(dinucleotides), all_of(trinucleotides), C0_new),
                 # select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase,
                 dat_poly)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                               all_of(dinucleotides), all_of(trinucleotides), C0_new),
                      # select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test,
                      dat_poly_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                                     all_of(dinucleotides), all_of(trinucleotides), C0_new),
                            # select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all,
                            dat_poly_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.5947672
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.5909181
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5465432
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(3),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9517412
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7526449
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6951263
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T, LightGBM")
dat_temp.feature_importance

high_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) > 1.5
plot(dat_temp.lgb_random_all_pred[high_resid_ids_random_all], y_random_all[high_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")

temp_interpretation_high = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(high_resid_ids_random_all==1))
lgb.plot.interpretation(temp_interpretation_high[[8]])

low_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) < 0.0005
plot(dat_temp.lgb_random_all_pred[low_resid_ids_random_all], y_random_all[low_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")

temp_interpretation_low = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(low_resid_ids_random_all==1))
lgb.plot.interpretation(temp_interpretation_low[[7]])


# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides

dat_temp = cbind(dat %>% 
                   # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                   #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 # select(all_of(trinucleotides), C0),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase,
                 dat_poly, n_gap_di)

dat_temp_test = cbind(dat_test %>% 
                        # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                        #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      # select(all_of(trinucleotides), C0),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test,
                      dat_poly_test, n_gap_di_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                              #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            # select(all_of(trinucleotides), C0),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all,
                            dat_poly_random_all, n_gap_di_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.6146394
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.6183059
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5781801
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LM")
abline(0,1)

temp_gam_formula = as.formula(paste0("C0_new~", paste0(colnames(dat_temp%>%select(all_of(ps2),
                                                                                  all_of(trinucleotides),
                                                                                  all_of(colnames(dat_npoly)))), 
                                                       collapse="+"),
                                     "+s(", paste0(colnames(dat_temp%>%select(-C0_new, 
                                                                              -all_of(ps2),
                                                                              -all_of(trinucleotides),
                                                                              -all_of(colnames(dat_npoly)))), 
                                                   collapse=",k=6)+s("), ",k=6)"))
temp_gam = gam(temp_gam_formula, data=dat_temp)

cor(temp_gam$fitted.values, y)
# 0.6450205
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



set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9632829
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LightGBM")
title("Trained on Tiling, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7596548
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LightGBM")
title("Trained on Tiling, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6933248
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LightGBM")
title("Trained on Tiling, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  # ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, LightGBM")
  ggtitle("Trained on Tiling, LightGBM")
dat_temp.feature_importance

# high_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) > 1.5
# plot(dat_temp.lgb_random_all_pred[high_resid_ids_random_all], y_random_all[high_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_high = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(high_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_high[[8]])
# 
# low_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) < 0.0005
# plot(dat_temp.lgb_random_all_pred[low_resid_ids_random_all], y_random_all[low_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_low = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(low_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_low[[7]])


dat_temp_random_trval = dat_temp_random_all[random_train_indices,]
dat_temp_random_test = dat_temp_random_all[-random_train_indices,]

full_dat_temp_trval = rbind(dat_temp, dat_temp_random_trval)
full_dat_temp_test = rbind(dat_temp_test, dat_temp_random_test)

full_nrow = 85355
set.seed(50)
full_train_indices = sample(1:full_nrow, full_nrow*.9, replace=FALSE)

full_dat_temp_train = full_dat_temp_trval[full_train_indices,]
full_dat_temp_val = full_dat_temp_trval[-full_train_indices,]

set.seed(50)
full_dat_temp.lgb_train_set = full_dat_temp_train %>%
  select(-C0_new)
full_dat_temp.lgb_train_rules = lgb.convert_with_rules(data = full_dat_temp.lgb_train_set)
full_dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(full_dat_temp.lgb_train_rules$data),
                                           label = full_dat_temp_train$C0_new,
                                           categorical_feature = handle_categorical(colnames(full_dat_temp.lgb_train_rules$data)))
full_dat_temp.lgb_val_set = full_dat_temp_val %>%
  select(-C0_new)
full_dat_temp.lgb_val_matrix = as.matrix(lgb.convert_with_rules(data = full_dat_temp.lgb_val_set, rules = full_dat_temp.lgb_train_rules$rules)$data)
full_dat_temp.lgb_val_data = lgb.Dataset(data = full_dat_temp.lgb_val_matrix,
                                         label = full_dat_temp_val$C0_new,
                                         categorical_feature = handle_categorical(colnames(full_dat_temp.lgb_val_matrix)))
full_dat_temp.lgb_obj = lightgbm(data = full_dat_temp.lgb_train_data,
                                 params = list(learning_rate = c(0.05),
                                               objective = c("regression"),
                                               min_data_in_leaf = 1000,
                                               max_depth = c(5),
                                               num_leaves = c(1000),
                                               lambda_l2 = c(10),
                                               boosting = c("gbdt")),
                                 valids = list(valid = full_dat_temp.lgb_val_data),
                                 nrounds = 16000,
                                 early_stopping_rounds = 50)

full_dat_temp.lgb_train_pred = predict(full_dat_temp.lgb_obj, data = as.matrix(full_dat_temp.lgb_train_rules$data))
cor(full_dat_temp.lgb_train_pred, full_dat_temp_train$C0_new, method="pearson")
# 0.9704614
plot(full_dat_temp.lgb_train_pred, full_dat_temp_train$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LightGBM")
title("Trained on Tiling+Random, Train, LightGBM")
abline(0,1)

full_dat_temp.lgb_val_pred = predict(full_dat_temp.lgb_obj, data = full_dat_temp.lgb_val_matrix)

cor(full_dat_temp.lgb_val_pred, full_dat_temp_val$C0_new, method="pearson")
# 0.7512439
plot(full_dat_temp.lgb_val_pred, full_dat_temp_val$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LightGBM")
title("Trained on Tiling+Random, Validation, LightGBM")
abline(0,1)

full_dat_temp.lgb_test_set = full_dat_temp_test %>%
  select(-C0_new)
full_dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = full_dat_temp.lgb_test_set, rules = full_dat_temp.lgb_train_rules$rules)$data)
full_dat_temp.lgb_test_pred = predict(full_dat_temp.lgb_obj, data = full_dat_temp.lgb_test_matrix)

cor(full_dat_temp.lgb_test_pred, full_dat_temp_test$C0_new, method="pearson")
# 0.7547138
plot(full_dat_temp.lgb_test_pred, full_dat_temp_test$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LightGBM")
title("Trained on Tiling+Random, Test, LightGBM")
abline(0,1)

full_dat_temp.feature_importance = as.data.frame(lgb.importance(full_dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  # ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, LightGBM")
  ggtitle("Trained on Tiling+Random, LightGBM")
full_dat_temp.feature_importance

# high_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) > 1.5
# plot(dat_temp.lgb_random_all_pred[high_resid_ids_random_all], y_random_all[high_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_high = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(high_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_high[[8]])
# 
# low_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) < 0.0005
# plot(dat_temp.lgb_random_all_pred[low_resid_ids_random_all], y_random_all[low_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_low = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(low_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_low[[7]])




# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + Groups of Dinucleotides + all polyA/T + N-Gapped Dinucleotides

dat_temp = cbind(dat %>% 
                   # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                   #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
                 # select(all_of(trinucleotides), C0),
                 Xone_A_fourier_mag, Xone_C_fourier_mag, Xone_G_fourier_mag, Xone_T_fourier_mag,
                 Xone_A_fourier_phase, Xone_C_fourier_phase, Xone_G_fourier_phase, Xone_T_fourier_phase,
                 Xone_AorT_fourier_mag, 
                 Xone_CorG_fourier_mag,
                 Xone_AorT_fourier_phase, 
                 Xone_CorG_fourier_phase, 
                 Xtwo_AA_fourier_mag, Xtwo_AT_fourier_mag,
                 Xtwo_CC_fourier_mag, Xtwo_CG_fourier_mag,
                 Xtwo_GC_fourier_mag, Xtwo_GG_fourier_mag,
                 Xtwo_TA_fourier_mag, Xtwo_TT_fourier_mag,
                 Xtwo_AA_fourier_phase, Xtwo_AT_fourier_phase,
                 Xtwo_CC_fourier_phase, Xtwo_CG_fourier_phase,
                 Xtwo_GC_fourier_phase, Xtwo_GG_fourier_phase,
                 Xtwo_TA_fourier_phase, Xtwo_TT_fourier_phase,
                 Xtwo_AAorAT_fourier_mag, Xtwo_AAorTA_fourier_mag,
                 Xtwo_AAorTT_fourier_mag, Xtwo_ATorTA_fourier_mag,
                 Xtwo_ATorTT_fourier_mag, Xtwo_TAorTT_fourier_mag,
                 Xtwo_CCorCG_fourier_mag, Xtwo_CCorGC_fourier_mag,
                 Xtwo_CCorGG_fourier_mag, Xtwo_CGorGC_fourier_mag,
                 Xtwo_CGorGG_fourier_mag, Xtwo_GCorGG_fourier_mag,
                 Xtwo_AAorAT_fourier_phase, Xtwo_AAorTA_fourier_phase,
                 Xtwo_AAorTT_fourier_phase, Xtwo_ATorTA_fourier_phase,
                 Xtwo_ATorTT_fourier_phase, Xtwo_TAorTT_fourier_phase,
                 Xtwo_CCorCG_fourier_phase, Xtwo_CCorGC_fourier_phase,
                 Xtwo_CCorGG_fourier_phase, Xtwo_CGorGC_fourier_phase,
                 Xtwo_CGorGG_fourier_phase, Xtwo_GCorGG_fourier_phase,
                 Xtwo_AAorATorTAorTT_fourier_mag, Xtwo_CCorCGorGCorGG_fourier_mag,
                 Xtwo_AAorATorTAorTT_fourier_phase, Xtwo_CCorCGorGCorGG_fourier_phase,
                 dat_poly, n_gap_di)

dat_temp_test = cbind(dat_test %>% 
                        # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                        #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
                      # select(all_of(trinucleotides), C0),
                      Xone_A_fourier_mag_test, Xone_C_fourier_mag_test, Xone_G_fourier_mag_test, Xone_T_fourier_mag_test,
                      Xone_A_fourier_phase_test, Xone_C_fourier_phase_test, Xone_G_fourier_phase_test, Xone_T_fourier_phase_test,
                      Xone_AorT_fourier_mag_test, 
                      Xone_CorG_fourier_mag_test,
                      Xone_AorT_fourier_phase_test, 
                      Xone_CorG_fourier_phase_test, 
                      Xtwo_AA_fourier_mag_test, Xtwo_AT_fourier_mag_test,
                      Xtwo_CC_fourier_mag_test, Xtwo_CG_fourier_mag_test,
                      Xtwo_GC_fourier_mag_test, Xtwo_GG_fourier_mag_test,
                      Xtwo_TA_fourier_mag_test, Xtwo_TT_fourier_mag_test,
                      Xtwo_AA_fourier_phase_test, Xtwo_AT_fourier_phase_test,
                      Xtwo_CC_fourier_phase_test, Xtwo_CG_fourier_phase_test,
                      Xtwo_GC_fourier_phase_test, Xtwo_GG_fourier_phase_test,
                      Xtwo_TA_fourier_phase_test, Xtwo_TT_fourier_phase_test,
                      Xtwo_AAorAT_fourier_mag_test, Xtwo_AAorTA_fourier_mag_test,
                      Xtwo_AAorTT_fourier_mag_test, Xtwo_ATorTA_fourier_mag_test,
                      Xtwo_ATorTT_fourier_mag_test, Xtwo_TAorTT_fourier_mag_test,
                      Xtwo_CCorCG_fourier_mag_test, Xtwo_CCorGC_fourier_mag_test,
                      Xtwo_CCorGG_fourier_mag_test, Xtwo_CGorGC_fourier_mag_test,
                      Xtwo_CGorGG_fourier_mag_test, Xtwo_GCorGG_fourier_mag_test,
                      Xtwo_AAorAT_fourier_phase_test, Xtwo_AAorTA_fourier_phase_test,
                      Xtwo_AAorTT_fourier_phase_test, Xtwo_ATorTA_fourier_phase_test,
                      Xtwo_ATorTT_fourier_phase_test, Xtwo_TAorTT_fourier_phase_test,
                      Xtwo_CCorCG_fourier_phase_test, Xtwo_CCorGC_fourier_phase_test,
                      Xtwo_CCorGG_fourier_phase_test, Xtwo_CGorGC_fourier_phase_test,
                      Xtwo_CGorGG_fourier_phase_test, Xtwo_GCorGG_fourier_phase_test,
                      Xtwo_AAorATorTAorTT_fourier_mag_test, Xtwo_CCorCGorGCorGG_fourier_mag_test,
                      Xtwo_AAorATorTAorTT_fourier_phase_test, Xtwo_CCorCGorGCorGG_fourier_phase_test,
                      dat_poly_test, n_gap_di_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              # select(all_of(ps1), all_of(ps2), all_of(nucleotides), 
                              #        all_of(dinucleotides), all_of(trinucleotides), C0_new),
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
                            # select(all_of(trinucleotides), C0),
                            Xone_A_fourier_mag_random_all, Xone_C_fourier_mag_random_all, Xone_G_fourier_mag_random_all, Xone_T_fourier_mag_random_all,
                            Xone_A_fourier_phase_random_all, Xone_C_fourier_phase_random_all, Xone_G_fourier_phase_random_all, Xone_T_fourier_phase_random_all,
                            Xone_AorT_fourier_mag_random_all, 
                            Xone_CorG_fourier_mag_random_all,
                            Xone_AorT_fourier_phase_random_all, 
                            Xone_CorG_fourier_phase_random_all, 
                            Xtwo_AA_fourier_mag_random_all, Xtwo_AT_fourier_mag_random_all,
                            Xtwo_CC_fourier_mag_random_all, Xtwo_CG_fourier_mag_random_all,
                            Xtwo_GC_fourier_mag_random_all, Xtwo_GG_fourier_mag_random_all,
                            Xtwo_TA_fourier_mag_random_all, Xtwo_TT_fourier_mag_random_all,
                            Xtwo_AA_fourier_phase_random_all, Xtwo_AT_fourier_phase_random_all,
                            Xtwo_CC_fourier_phase_random_all, Xtwo_CG_fourier_phase_random_all,
                            Xtwo_GC_fourier_phase_random_all, Xtwo_GG_fourier_phase_random_all,
                            Xtwo_TA_fourier_phase_random_all, Xtwo_TT_fourier_phase_random_all,
                            Xtwo_AAorAT_fourier_mag_random_all, Xtwo_AAorTA_fourier_mag_random_all,
                            Xtwo_AAorTT_fourier_mag_random_all, Xtwo_ATorTA_fourier_mag_random_all,
                            Xtwo_ATorTT_fourier_mag_random_all, Xtwo_TAorTT_fourier_mag_random_all,
                            Xtwo_CCorCG_fourier_mag_random_all, Xtwo_CCorGC_fourier_mag_random_all,
                            Xtwo_CCorGG_fourier_mag_random_all, Xtwo_CGorGC_fourier_mag_random_all,
                            Xtwo_CGorGG_fourier_mag_random_all, Xtwo_GCorGG_fourier_mag_random_all,
                            Xtwo_AAorAT_fourier_phase_random_all, Xtwo_AAorTA_fourier_phase_random_all,
                            Xtwo_AAorTT_fourier_phase_random_all, Xtwo_ATorTA_fourier_phase_random_all,
                            Xtwo_ATorTT_fourier_phase_random_all, Xtwo_TAorTT_fourier_phase_random_all,
                            Xtwo_CCorCG_fourier_phase_random_all, Xtwo_CCorGC_fourier_phase_random_all,
                            Xtwo_CCorGG_fourier_phase_random_all, Xtwo_CGorGC_fourier_phase_random_all,
                            Xtwo_CGorGG_fourier_phase_random_all, Xtwo_GCorGG_fourier_phase_random_all,
                            Xtwo_AAorATorTAorTT_fourier_mag_random_all, Xtwo_CCorCGorGCorGG_fourier_mag_random_all,
                            Xtwo_AAorATorTAorTT_fourier_phase_random_all, Xtwo_CCorCGorGCorGG_fourier_phase_random_all,
                            dat_poly_random_all, n_gap_di_random_all)

temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.6425337
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.6472649
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.5917098
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LM")
abline(0,1)

temp_gam_formula = as.formula(paste0("C0_new~", paste0(colnames(dat_temp%>%select(all_of(ps2),
                                                                                  all_of(trinucleotides),
                                                                                  all_of(colnames(dat_npoly)))), 
                                                       collapse="+"),
                                     "+s(", paste0(colnames(dat_temp%>%select(-C0_new, 
                                                                              -all_of(ps2),
                                                                              -all_of(trinucleotides),
                                                                              -all_of(colnames(dat_npoly)))), 
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



set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9660618
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LightGBM")
title("Trained on Tiling, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7628257
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LightGBM")
title("Trained on Tiling, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6955792
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LightGBM")
title("Trained on Tiling, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  # ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, LightGBM")
  ggtitle("Trained on Tiling, LightGBM")
dat_temp.feature_importance

# high_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) > 1.5
# plot(dat_temp.lgb_random_all_pred[high_resid_ids_random_all], y_random_all[high_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_high = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(high_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_high[[8]])
# 
# low_resid_ids_random_all = abs(dat_temp.lgb_random_all_pred - y_random_all) < 0.0005
# plot(dat_temp.lgb_random_all_pred[low_resid_ids_random_all], y_random_all[low_resid_ids_random_all], xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col="red")
# 
# temp_interpretation_low = lgb.interprete(dat_temp.lgb_obj, dat_temp.lgb_random_all_matrix, which(low_resid_ids_random_all==1))
# lgb.plot.interpretation(temp_interpretation_low[[7]])


dat_temp_random_trval = dat_temp_random_all[random_train_indices,]
dat_temp_random_test = dat_temp_random_all[-random_train_indices,]

full_dat_temp_trval = rbind(dat_temp, dat_temp_random_trval)
full_dat_temp_test = rbind(dat_temp_test, dat_temp_random_test)

full_nrow = 85355
set.seed(50)
full_train_indices = sample(1:full_nrow, full_nrow*.9, replace=FALSE)

full_dat_temp_train = full_dat_temp_trval[full_train_indices,]
full_dat_temp_val = full_dat_temp_trval[-full_train_indices,]

set.seed(50)
full_dat_temp.lgb_train_set = full_dat_temp_train %>%
  select(-C0_new)
full_dat_temp.lgb_train_rules = lgb.convert_with_rules(data = full_dat_temp.lgb_train_set)
full_dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(full_dat_temp.lgb_train_rules$data),
                                           label = full_dat_temp_train$C0_new,
                                           categorical_feature = handle_categorical(colnames(full_dat_temp.lgb_train_rules$data)))
full_dat_temp.lgb_val_set = full_dat_temp_val %>%
  select(-C0_new)
full_dat_temp.lgb_val_matrix = as.matrix(lgb.convert_with_rules(data = full_dat_temp.lgb_val_set, rules = full_dat_temp.lgb_train_rules$rules)$data)
full_dat_temp.lgb_val_data = lgb.Dataset(data = full_dat_temp.lgb_val_matrix,
                                         label = full_dat_temp_val$C0_new,
                                         categorical_feature = handle_categorical(colnames(full_dat_temp.lgb_val_matrix)))
full_dat_temp.lgb_obj = lightgbm(data = full_dat_temp.lgb_train_data,
                                 params = list(learning_rate = c(0.05),
                                               objective = c("regression"),
                                               min_data_in_leaf = 1000,
                                               max_depth = c(5),
                                               num_leaves = c(1000),
                                               lambda_l2 = c(10),
                                               boosting = c("gbdt")),
                                 valids = list(valid = full_dat_temp.lgb_val_data),
                                 nrounds = 16000,
                                 early_stopping_rounds = 50)

full_dat_temp.lgb_train_pred = predict(full_dat_temp.lgb_obj, data = as.matrix(full_dat_temp.lgb_train_rules$data))
cor(full_dat_temp.lgb_train_pred, full_dat_temp_train$C0_new, method="pearson")
# 
plot(full_dat_temp.lgb_train_pred, full_dat_temp_train$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Train, LightGBM")
title("Trained on Tiling+Random, Train, LightGBM")
abline(0,1)

full_dat_temp.lgb_val_pred = predict(full_dat_temp.lgb_obj, data = full_dat_temp.lgb_val_matrix)

cor(full_dat_temp.lgb_val_pred, full_dat_temp_val$C0_new, method="pearson")
# 
plot(full_dat_temp.lgb_val_pred, full_dat_temp_val$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Test, LightGBM")
title("Trained on Tiling+Random, Validation, LightGBM")
abline(0,1)

full_dat_temp.lgb_test_set = full_dat_temp_test %>%
  select(-C0_new)
full_dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = full_dat_temp.lgb_test_set, rules = full_dat_temp.lgb_train_rules$rules)$data)
full_dat_temp.lgb_test_pred = predict(full_dat_temp.lgb_obj, data = full_dat_temp.lgb_test_matrix)

cor(full_dat_temp.lgb_test_pred, full_dat_temp_test$C0_new, method="pearson")
# 
plot(full_dat_temp.lgb_test_pred, full_dat_temp_test$C0_new, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1), xlab="Predicted", ylab="Observed")
# title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, Random, LightGBM")
title("Trained on Tiling+Random, Test, LightGBM")
abline(0,1)

full_dat_temp.feature_importance = as.data.frame(lgb.importance(full_dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() + 
  # ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + all polyA/T + N-Gapped Dinucleotides, LightGBM")
  ggtitle("Trained on Tiling+Random, LightGBM")
full_dat_temp.feature_importance


# Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + 
# (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides 
# and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides

dat_temp = cbind(dat %>% 
                   select(all_of(ps2), all_of(trinucleotides), C0_new),
                 # select(all_of(trinucleotides), C0_new),
                 # select(C0_new),
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
                 Xtwo_AAorAT_fourier_phase^2, Xtwo_AAorTA_fourier_phase^2,
                 Xtwo_AAorTT_fourier_phase^2, Xtwo_ATorTA_fourier_phase^2,
                 Xtwo_ATorTT_fourier_phase^2, Xtwo_TAorTT_fourier_phase^2,
                 Xtwo_CCorCG_fourier_phase^2, Xtwo_CCorGC_fourier_phase^2,
                 Xtwo_CCorGG_fourier_phase^2, Xtwo_CGorGC_fourier_phase^2,
                 Xtwo_CGorGG_fourier_phase^2, Xtwo_GCorGG_fourier_phase^2,
                 dist_freq_di_AorTdi_bid2, dist_freq_di_CorGdi_bid2)

dat_temp_test = cbind(dat_test %>% 
                        select(all_of(ps2), all_of(trinucleotides), C0_new),
                      # select(all_of(trinucleotides), C0_new),
                      # select(C0_new),
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
                      Xtwo_AAorAT_fourier_phase_test^2, Xtwo_AAorTA_fourier_phase_test^2,
                      Xtwo_AAorTT_fourier_phase_test^2, Xtwo_ATorTA_fourier_phase_test^2,
                      Xtwo_ATorTT_fourier_phase_test^2, Xtwo_TAorTT_fourier_phase_test^2,
                      Xtwo_CCorCG_fourier_phase_test^2, Xtwo_CCorGC_fourier_phase_test^2,
                      Xtwo_CCorGG_fourier_phase_test^2, Xtwo_CGorGC_fourier_phase_test^2,
                      Xtwo_CGorGG_fourier_phase_test^2, Xtwo_GCorGG_fourier_phase_test^2,
                      dist_freq_di_AorTdi_bid2_test, dist_freq_di_CorGdi_bid2_test)

dat_temp_random_all = cbind(dat_random_all %>% 
                              select(all_of(ps2), all_of(trinucleotides), C0_new),
                            # select(all_of(trinucleotides), C0_new),
                            # select(C0_new),
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
                            Xtwo_AAorAT_fourier_phase_random_all^2, Xtwo_AAorTA_fourier_phase_random_all^2,
                            Xtwo_AAorTT_fourier_phase_random_all^2, Xtwo_ATorTA_fourier_phase_random_all^2,
                            Xtwo_ATorTT_fourier_phase_random_all^2, Xtwo_TAorTT_fourier_phase_random_all^2,
                            Xtwo_CCorCG_fourier_phase_random_all^2, Xtwo_CCorGC_fourier_phase_random_all^2,
                            Xtwo_CCorGG_fourier_phase_random_all^2, Xtwo_CGorGC_fourier_phase_random_all^2,
                            Xtwo_CGorGG_fourier_phase_random_all^2, Xtwo_GCorGG_fourier_phase_random_all^2,
                            dist_freq_di_AorTdi_bid2_random_all, dist_freq_di_CorGdi_bid2_random_all)


temp_lm = lm(C0_new~., data=dat_temp)

cor(temp_lm$fitted.values, y)
# 0.7624294
plot(temp_lm$fitted.values, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Train, LM")
abline(0,1)

temp_lm_pred = predict(temp_lm, dat_temp_test)
cor(temp_lm_pred, y_test)
# 0.7558507
plot(temp_lm_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Test, LM")
abline(0,1)

temp_lm_random_all_pred = predict(temp_lm, dat_temp_random_all)
cor(temp_lm_random_all_pred, y_random_all)
# 0.6785853
plot(temp_lm_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Random, LM")
abline(0,1)


set.seed(50)
dat_temp.lgb_train_set = dat_temp %>%
  select(-C0_new)
dat_temp.lgb_train_rules = lgb.convert_with_rules(data = dat_temp.lgb_train_set)
dat_temp.lgb_train_data = lgb.Dataset(data = as.matrix(dat_temp.lgb_train_rules$data),
                                      label = dat_temp$C0_new,
                                      categorical_feature = handle_categorical(colnames(dat_temp.lgb_train_rules$data)))
dat_temp.lgb_test_set = dat_temp_test %>%
  select(-C0_new)
dat_temp.lgb_test_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_test_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_test_data = lgb.Dataset(data = dat_temp.lgb_test_matrix,
                                     label = dat_temp_test$C0_new,
                                     categorical_feature = handle_categorical(colnames(dat_temp.lgb_test_matrix)))
dat_temp.lgb_obj = lightgbm(data = dat_temp.lgb_train_data,
                            params = list(learning_rate = c(0.05),
                                          objective = c("regression"),
                                          min_data_in_leaf = 1000,
                                          max_depth = c(5),
                                          num_leaves = c(1000),
                                          lambda_l2 = c(10),
                                          boosting = c("gbdt")),
                            valids = list(valid = dat_temp.lgb_test_data),
                            nrounds = 16000,
                            early_stopping_rounds = 50)

dat_temp.lgb_train_pred = predict(dat_temp.lgb_obj, data = as.matrix(dat_temp.lgb_train_rules$data))
cor(dat_temp.lgb_train_pred, y, method="pearson")
# 0.9756495
plot(dat_temp.lgb_train_pred, y, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Train, LightGBM")
abline(0,1)

dat_temp.lgb_test_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_test_matrix)

cor(dat_temp.lgb_test_pred, y_test, method="pearson")
# 0.7830386
plot(dat_temp.lgb_test_pred, y_test, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Test, LightGBM")
abline(0,1)

dat_temp.lgb_random_all_set = dat_temp_random_all %>%
  select(-C0_new)
dat_temp.lgb_random_all_matrix = as.matrix(lgb.convert_with_rules(data = dat_temp.lgb_random_all_set, rules = dat_temp.lgb_train_rules$rules)$data)
dat_temp.lgb_random_all_pred = predict(dat_temp.lgb_obj, data = dat_temp.lgb_random_all_matrix)

cor(dat_temp.lgb_random_all_pred, y_random_all, method="pearson")
# 0.6865105
plot(dat_temp.lgb_random_all_pred, y_random_all, xlim=c(-3.5, 3.5), ylim=c(-3.5,3.5), col=alpha("blue", 0.1))
title("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, Random, LightGBM")
abline(0,1)

dat_temp.feature_importance = as.data.frame(lgb.importance(dat_temp.lgb_obj)) %>%
  mutate(Feature = fct_rev(fct_inorder(Feature))) %>%
  dplyr::slice(1:50) %>%
  ggplot(aes(x = Feature, y = Gain)) +
  geom_col(fill = "dodgerblue") +
  coord_flip() +
  theme_bw() +
  ggtitle("Nucleotides + (Selected) Pairs of Nucleotides + (Selected) Dinucleotides + (Selected) Pairs of Dinucleotides + All Distance Frequency Between Dinucleotides and AA/TT/AT/TA or CC/GG/CG/GC Dinucleotides, LightGBM")
dat_temp.feature_importance








