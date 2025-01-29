dat = readRDS("data/Created/processed_chrV_post_smooth_C0.rds")
dat_test = readRDS("data/Created/processed_chrV_post_smooth_C0_test.rds")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")
ps3 <- paste0("X", 1:48, "tri")


# Nucleotides:

find_fourier_ps1 = function(dat, nuc_pos, fft_indices, nuc_neg=NA) {
  Xone = dat %>% select(all_of(ps1))
  Xone_signal = matrix(0, nrow=nrow(Xone), ncol=ncol(Xone))
  colnames(Xone_signal) = colnames(Xone)
  for(nuc in nuc_pos) {
    Xone_signal[] = Xone_signal[] + ((Xone == nuc) %>% as.matrix() %>% as.numeric())
  }
  if(!is.na(nuc_neg)) {
    for(nuc in nuc_neg) {
      Xone_signal[] = Xone_signal[] - ((Xone ==nuc) %>% as.matrix() %>% as.numeric())
    }
  }
  Xone_fourier = t(apply(Xone_signal, 1, fft))[,fft_indices]
  
  Xone_fourier_mag = as.matrix(Mod(Xone_fourier))
  Xone_fourier_phase = as.matrix(Arg(Xone_fourier))
  
  nuc_pos_str = paste(nuc_pos, collapse="or")
  if(is.na(nuc_neg)) {
    colnames(Xone_fourier_mag) = paste0("fourier_mag_", nuc_pos_str, "_", fft_indices)
    colnames(Xone_fourier_phase) = paste0("fourier_phase_", nuc_pos_str, "_", fft_indices)
  }
  else {
    nuc_neg_str = paste(nuc_neg, collapse="or")
    colnames(Xone_fourier_mag) = paste0("fourier_mag_", nuc_pos_str, "_minus_", nuc_neg_str, "_", fft_indices)
    colnames(Xone_fourier_phase) = paste0("fourier_phase_", nuc_pos_str, "_minus_", nuc_neg_str, "_", fft_indices)
  }
  ret = list()
  ret$mag = Xone_fourier_mag
  ret$phase = Xone_fourier_phase
  return(ret)
}

fft_indices = c(1,2,3,4,5,6,11,21)

Xone_A_fourier = find_fourier_ps1(dat, c("A"), fft_indices)
Xone_C_fourier = find_fourier_ps1(dat, c("C"), fft_indices)
Xone_G_fourier = find_fourier_ps1(dat, c("G"), fft_indices)
Xone_T_fourier = find_fourier_ps1(dat, c("T"), fft_indices)

Xone_nuc_fourier_mag = cbind(Xone_A_fourier$mag, Xone_C_fourier$mag,
                             Xone_G_fourier$mag, Xone_T_fourier$mag)
saveRDS(Xone_nuc_fourier_mag, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_mag.rds")

Xone_nuc_fourier_phase = cbind(Xone_A_fourier$phase, Xone_C_fourier$phase,
                               Xone_G_fourier$phase, Xone_T_fourier$phase)
saveRDS(Xone_nuc_fourier_phase, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_phase.rds")


Xone_A_fourier_test = find_fourier_ps1(dat_test, c("A"), fft_indices)
Xone_C_fourier_test = find_fourier_ps1(dat_test, c("C"), fft_indices)
Xone_G_fourier_test = find_fourier_ps1(dat_test, c("G"), fft_indices)
Xone_T_fourier_test = find_fourier_ps1(dat_test, c("T"), fft_indices)

Xone_nuc_fourier_mag_test = cbind(Xone_A_fourier_test$mag, Xone_C_fourier_test$mag,
                                  Xone_G_fourier_test$mag, Xone_T_fourier_test$mag)
saveRDS(Xone_nuc_fourier_mag_test, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_mag_test.rds")

Xone_nuc_fourier_phase_test = cbind(Xone_A_fourier_test$phase, Xone_C_fourier_test$phase,
                                    Xone_G_fourier_test$phase, Xone_T_fourier_test$phase)
saveRDS(Xone_nuc_fourier_phase_test, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_phase_test.rds")



# Differences between single nucleotides:

Xone_AminusC_fourier = find_fourier_ps1(dat, c("A"), fft_indices, c("C"))
Xone_AminusG_fourier = find_fourier_ps1(dat, c("A"), fft_indices, c("G"))
Xone_AminusT_fourier = find_fourier_ps1(dat, c("A"), fft_indices, c("T"))
Xone_CminusG_fourier = find_fourier_ps1(dat, c("C"), fft_indices, c("G"))
Xone_CminusT_fourier = find_fourier_ps1(dat, c("C"), fft_indices, c("T"))
Xone_GminusT_fourier = find_fourier_ps1(dat, c("G"), fft_indices, c("T"))

Xone_nucminusnuc_fourier_mag = cbind(Xone_AminusC_fourier$mag, 
                                     Xone_AminusG_fourier$mag,
                                     Xone_AminusT_fourier$mag,
                                     Xone_CminusG_fourier$mag,
                                     Xone_CminusT_fourier$mag,
                                     Xone_GminusT_fourier$mag)
saveRDS(Xone_nucminusnuc_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_mag.rds")

Xone_nucminusnuc_fourier_phase = cbind(Xone_AminusC_fourier$phase, 
                                       Xone_AminusG_fourier$phase,
                                       Xone_AminusT_fourier$phase,
                                       Xone_CminusG_fourier$phase,
                                       Xone_CminusT_fourier$phase,
                                       Xone_GminusT_fourier$phase)
saveRDS(Xone_nucminusnuc_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_phase.rds")


Xone_AminusC_fourier_test = find_fourier_ps1(dat_test, c("A"), fft_indices, c("C"))
Xone_AminusG_fourier_test = find_fourier_ps1(dat_test, c("A"), fft_indices, c("G"))
Xone_AminusT_fourier_test = find_fourier_ps1(dat_test, c("A"), fft_indices, c("T"))
Xone_CminusG_fourier_test = find_fourier_ps1(dat_test, c("C"), fft_indices, c("G"))
Xone_CminusT_fourier_test = find_fourier_ps1(dat_test, c("C"), fft_indices, c("T"))
Xone_GminusT_fourier_test = find_fourier_ps1(dat_test, c("G"), fft_indices, c("T"))

Xone_nucminusnuc_fourier_mag_test = cbind(Xone_AminusC_fourier_test$mag, 
                                          Xone_AminusG_fourier_test$mag,
                                          Xone_AminusT_fourier_test$mag,
                                          Xone_CminusG_fourier_test$mag,
                                          Xone_CminusT_fourier_test$mag,
                                          Xone_GminusT_fourier_test$mag)
saveRDS(Xone_nucminusnuc_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_mag_test.rds")

Xone_nucminusnuc_fourier_phase_test = cbind(Xone_AminusC_fourier_test$phase, 
                                            Xone_AminusG_fourier_test$phase,
                                            Xone_AminusT_fourier_test$phase,
                                            Xone_CminusG_fourier_test$phase,
                                            Xone_CminusT_fourier_test$phase,
                                            Xone_GminusT_fourier_test$phase)
saveRDS(Xone_nucminusnuc_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_phase_test.rds")


# Pairs of nucleotides:

Xone_AorC_fourier = find_fourier_ps1(dat, c("A", "C"), fft_indices)
Xone_AorG_fourier = find_fourier_ps1(dat, c("A", "G"), fft_indices)
Xone_AorT_fourier = find_fourier_ps1(dat, c("A", "T"), fft_indices)
Xone_CorG_fourier = find_fourier_ps1(dat, c("C", "G"), fft_indices)
Xone_CorT_fourier = find_fourier_ps1(dat, c("C", "T"), fft_indices)
Xone_GorT_fourier = find_fourier_ps1(dat, c("G", "T"), fft_indices)

Xone_nucornuc_fourier_mag = cbind(Xone_AorC_fourier$mag, 
                                  Xone_AorG_fourier$mag,
                                  Xone_AorT_fourier$mag,
                                  Xone_CorG_fourier$mag,
                                  Xone_CorT_fourier$mag,
                                  Xone_GorT_fourier$mag)
saveRDS(Xone_nucornuc_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xone_nucornuc_fourier_mag.rds")

Xone_nucornuc_fourier_phase = cbind(Xone_AorC_fourier$phase, 
                                    Xone_AorG_fourier$phase,
                                    Xone_AorT_fourier$phase,
                                    Xone_CorG_fourier$phase,
                                    Xone_CorT_fourier$phase,
                                    Xone_GorT_fourier$phase)
saveRDS(Xone_nucornuc_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xone_nucornuc_fourier_phase.rds")



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

Xtwo_AA_fourier = t(apply(Xtwo_AA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AC_fourier = t(apply(Xtwo_AC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AG_fourier = t(apply(Xtwo_AG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AT_fourier = t(apply(Xtwo_AT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CA_fourier = t(apply(Xtwo_CA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CC_fourier = t(apply(Xtwo_CC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CG_fourier = t(apply(Xtwo_CG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CT_fourier = t(apply(Xtwo_CT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GA_fourier = t(apply(Xtwo_GA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GC_fourier = t(apply(Xtwo_GC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GG_fourier = t(apply(Xtwo_GG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GT_fourier = t(apply(Xtwo_GT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TA_fourier = t(apply(Xtwo_TA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TC_fourier = t(apply(Xtwo_TC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TG_fourier = t(apply(Xtwo_TG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TT_fourier = t(apply(Xtwo_TT, 1, fft))[,c(1,2,3,4,5,6,11,21)]

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

colnames(Xtwo_AA_fourier_mag) = paste0("fourier_mag_AA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AC_fourier_mag) = paste0("fourier_mag_AC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AG_fourier_mag) = paste0("fourier_mag_AG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AT_fourier_mag) = paste0("fourier_mag_AT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CA_fourier_mag) = paste0("fourier_mag_CA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CC_fourier_mag) = paste0("fourier_mag_CC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CG_fourier_mag) = paste0("fourier_mag_CG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CT_fourier_mag) = paste0("fourier_mag_CT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GA_fourier_mag) = paste0("fourier_mag_GA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GC_fourier_mag) = paste0("fourier_mag_GC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GG_fourier_mag) = paste0("fourier_mag_GG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GT_fourier_mag) = paste0("fourier_mag_GT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TA_fourier_mag) = paste0("fourier_mag_TA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TC_fourier_mag) = paste0("fourier_mag_TC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TG_fourier_mag) = paste0("fourier_mag_TG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TT_fourier_mag) = paste0("fourier_mag_TT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AA_fourier_phase) = paste0("fourier_phase_AA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AC_fourier_phase) = paste0("fourier_phase_AC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AG_fourier_phase) = paste0("fourier_phase_AG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AT_fourier_phase) = paste0("fourier_phase_AT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CA_fourier_phase) = paste0("fourier_phase_CA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CC_fourier_phase) = paste0("fourier_phase_CC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CG_fourier_phase) = paste0("fourier_phase_CG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CT_fourier_phase) = paste0("fourier_phase_CT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GA_fourier_phase) = paste0("fourier_phase_GA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GC_fourier_phase) = paste0("fourier_phase_GC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GG_fourier_phase) = paste0("fourier_phase_GG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GT_fourier_phase) = paste0("fourier_phase_GT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TA_fourier_phase) = paste0("fourier_phase_TA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TC_fourier_phase) = paste0("fourier_phase_TC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TG_fourier_phase) = paste0("fourier_phase_TG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TT_fourier_phase) = paste0("fourier_phase_TT_", c(1,2,3,4,5,6,11,21))

Xtwo_di_fourier_mag = cbind(Xtwo_AA_fourier_mag,
                            Xtwo_AC_fourier_mag,
                            Xtwo_AG_fourier_mag,
                            Xtwo_AT_fourier_mag,
                            Xtwo_CA_fourier_mag,
                            Xtwo_CC_fourier_mag,
                            Xtwo_CG_fourier_mag,
                            Xtwo_CT_fourier_mag,
                            Xtwo_GA_fourier_mag,
                            Xtwo_GC_fourier_mag,
                            Xtwo_GG_fourier_mag,
                            Xtwo_GT_fourier_mag,
                            Xtwo_TA_fourier_mag,
                            Xtwo_TC_fourier_mag,
                            Xtwo_TG_fourier_mag,
                            Xtwo_TT_fourier_mag)
saveRDS(Xtwo_di_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xtwo_di_fourier_mag.rds")

Xtwo_di_fourier_phase = cbind(Xtwo_AA_fourier_phase,
                              Xtwo_AC_fourier_phase,
                              Xtwo_AG_fourier_phase,
                              Xtwo_AT_fourier_phase,
                              Xtwo_CA_fourier_phase,
                              Xtwo_CC_fourier_phase,
                              Xtwo_CG_fourier_phase,
                              Xtwo_CT_fourier_phase,
                              Xtwo_GA_fourier_phase,
                              Xtwo_GC_fourier_phase,
                              Xtwo_GG_fourier_phase,
                              Xtwo_GT_fourier_phase,
                              Xtwo_TA_fourier_phase,
                              Xtwo_TC_fourier_phase,
                              Xtwo_TG_fourier_phase,
                              Xtwo_TT_fourier_phase)
saveRDS(Xtwo_di_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xtwo_di_fourier_phase.rds")


# Difference between single dinucleotides:

Xtwo_AAminusAC_fourier = Xtwo_AA_fourier - Xtwo_AC_fourier
Xtwo_AAminusAG_fourier = Xtwo_AA_fourier - Xtwo_AG_fourier
Xtwo_AAminusAT_fourier = Xtwo_AA_fourier - Xtwo_AT_fourier
Xtwo_AAminusCA_fourier = Xtwo_AA_fourier - Xtwo_CA_fourier
Xtwo_AAminusCC_fourier = Xtwo_AA_fourier - Xtwo_CC_fourier
Xtwo_AAminusCG_fourier = Xtwo_AA_fourier - Xtwo_CG_fourier
Xtwo_AAminusCT_fourier = Xtwo_AA_fourier - Xtwo_CT_fourier
Xtwo_AAminusGA_fourier = Xtwo_AA_fourier - Xtwo_GA_fourier
Xtwo_AAminusGC_fourier = Xtwo_AA_fourier - Xtwo_GC_fourier
Xtwo_AAminusGG_fourier = Xtwo_AA_fourier - Xtwo_GG_fourier
Xtwo_AAminusGT_fourier = Xtwo_AA_fourier - Xtwo_GT_fourier
Xtwo_AAminusTA_fourier = Xtwo_AA_fourier - Xtwo_TA_fourier
Xtwo_AAminusTC_fourier = Xtwo_AA_fourier - Xtwo_TC_fourier
Xtwo_AAminusTG_fourier = Xtwo_AA_fourier - Xtwo_TG_fourier
Xtwo_AAminusTT_fourier = Xtwo_AA_fourier - Xtwo_TT_fourier

Xtwo_ACminusAG_fourier = Xtwo_AC_fourier - Xtwo_AG_fourier
Xtwo_ACminusAT_fourier = Xtwo_AC_fourier - Xtwo_AT_fourier
Xtwo_ACminusCA_fourier = Xtwo_AC_fourier - Xtwo_CA_fourier
Xtwo_ACminusCC_fourier = Xtwo_AC_fourier - Xtwo_CC_fourier
Xtwo_ACminusCG_fourier = Xtwo_AC_fourier - Xtwo_CG_fourier
Xtwo_ACminusCT_fourier = Xtwo_AC_fourier - Xtwo_CT_fourier
Xtwo_ACminusGA_fourier = Xtwo_AC_fourier - Xtwo_GA_fourier
Xtwo_ACminusGC_fourier = Xtwo_AC_fourier - Xtwo_GC_fourier
Xtwo_ACminusGG_fourier = Xtwo_AC_fourier - Xtwo_GG_fourier
Xtwo_ACminusGT_fourier = Xtwo_AC_fourier - Xtwo_GT_fourier
Xtwo_ACminusTA_fourier = Xtwo_AC_fourier - Xtwo_TA_fourier
Xtwo_ACminusTC_fourier = Xtwo_AC_fourier - Xtwo_TC_fourier
Xtwo_ACminusTG_fourier = Xtwo_AC_fourier - Xtwo_TG_fourier
Xtwo_ACminusTT_fourier = Xtwo_AC_fourier - Xtwo_TT_fourier

Xtwo_AGminusAT_fourier = Xtwo_AG_fourier - Xtwo_AT_fourier
Xtwo_AGminusCA_fourier = Xtwo_AG_fourier - Xtwo_CA_fourier
Xtwo_AGminusCC_fourier = Xtwo_AG_fourier - Xtwo_CC_fourier
Xtwo_AGminusCG_fourier = Xtwo_AG_fourier - Xtwo_CG_fourier
Xtwo_AGminusCT_fourier = Xtwo_AG_fourier - Xtwo_CT_fourier
Xtwo_AGminusGA_fourier = Xtwo_AG_fourier - Xtwo_GA_fourier
Xtwo_AGminusGC_fourier = Xtwo_AG_fourier - Xtwo_GC_fourier
Xtwo_AGminusGG_fourier = Xtwo_AG_fourier - Xtwo_GG_fourier
Xtwo_AGminusGT_fourier = Xtwo_AG_fourier - Xtwo_GT_fourier
Xtwo_AGminusTA_fourier = Xtwo_AG_fourier - Xtwo_TA_fourier
Xtwo_AGminusTC_fourier = Xtwo_AG_fourier - Xtwo_TC_fourier
Xtwo_AGminusTG_fourier = Xtwo_AG_fourier - Xtwo_TG_fourier
Xtwo_AGminusTT_fourier = Xtwo_AG_fourier - Xtwo_TT_fourier

Xtwo_ATminusCA_fourier = Xtwo_AT_fourier - Xtwo_CA_fourier
Xtwo_ATminusCC_fourier = Xtwo_AT_fourier - Xtwo_CC_fourier
Xtwo_ATminusCG_fourier = Xtwo_AT_fourier - Xtwo_CG_fourier
Xtwo_ATminusCT_fourier = Xtwo_AT_fourier - Xtwo_CT_fourier
Xtwo_ATminusGA_fourier = Xtwo_AT_fourier - Xtwo_GA_fourier
Xtwo_ATminusGC_fourier = Xtwo_AT_fourier - Xtwo_GC_fourier
Xtwo_ATminusGG_fourier = Xtwo_AT_fourier - Xtwo_GG_fourier
Xtwo_ATminusGT_fourier = Xtwo_AT_fourier - Xtwo_GT_fourier
Xtwo_ATminusTA_fourier = Xtwo_AT_fourier - Xtwo_TA_fourier
Xtwo_ATminusTC_fourier = Xtwo_AT_fourier - Xtwo_TC_fourier
Xtwo_ATminusTG_fourier = Xtwo_AT_fourier - Xtwo_TG_fourier
Xtwo_ATminusTT_fourier = Xtwo_AT_fourier - Xtwo_TT_fourier

Xtwo_CAminusCC_fourier = Xtwo_CA_fourier - Xtwo_CC_fourier
Xtwo_CAminusCG_fourier = Xtwo_CA_fourier - Xtwo_CG_fourier
Xtwo_CAminusCT_fourier = Xtwo_CA_fourier - Xtwo_CT_fourier
Xtwo_CAminusGA_fourier = Xtwo_CA_fourier - Xtwo_GA_fourier
Xtwo_CAminusGC_fourier = Xtwo_CA_fourier - Xtwo_GC_fourier
Xtwo_CAminusGG_fourier = Xtwo_CA_fourier - Xtwo_GG_fourier
Xtwo_CAminusGT_fourier = Xtwo_CA_fourier - Xtwo_GT_fourier
Xtwo_CAminusTA_fourier = Xtwo_CA_fourier - Xtwo_TA_fourier
Xtwo_CAminusTC_fourier = Xtwo_CA_fourier - Xtwo_TC_fourier
Xtwo_CAminusTG_fourier = Xtwo_CA_fourier - Xtwo_TG_fourier
Xtwo_CAminusTT_fourier = Xtwo_CA_fourier - Xtwo_TT_fourier

Xtwo_CCminusCG_fourier = Xtwo_CC_fourier - Xtwo_CG_fourier
Xtwo_CCminusCT_fourier = Xtwo_CC_fourier - Xtwo_CT_fourier
Xtwo_CCminusGA_fourier = Xtwo_CC_fourier - Xtwo_GA_fourier
Xtwo_CCminusGC_fourier = Xtwo_CC_fourier - Xtwo_GC_fourier
Xtwo_CCminusGG_fourier = Xtwo_CC_fourier - Xtwo_GG_fourier
Xtwo_CCminusGT_fourier = Xtwo_CC_fourier - Xtwo_GT_fourier
Xtwo_CCminusTA_fourier = Xtwo_CC_fourier - Xtwo_TA_fourier
Xtwo_CCminusTC_fourier = Xtwo_CC_fourier - Xtwo_TC_fourier
Xtwo_CCminusTG_fourier = Xtwo_CC_fourier - Xtwo_TG_fourier
Xtwo_CCminusTT_fourier = Xtwo_CC_fourier - Xtwo_TT_fourier

Xtwo_CGminusCT_fourier = Xtwo_CG_fourier - Xtwo_CT_fourier
Xtwo_CGminusGA_fourier = Xtwo_CG_fourier - Xtwo_GA_fourier
Xtwo_CGminusGC_fourier = Xtwo_CG_fourier - Xtwo_GC_fourier
Xtwo_CGminusGG_fourier = Xtwo_CG_fourier - Xtwo_GG_fourier
Xtwo_CGminusGT_fourier = Xtwo_CG_fourier - Xtwo_GT_fourier
Xtwo_CGminusTA_fourier = Xtwo_CG_fourier - Xtwo_TA_fourier
Xtwo_CGminusTC_fourier = Xtwo_CG_fourier - Xtwo_TC_fourier
Xtwo_CGminusTG_fourier = Xtwo_CG_fourier - Xtwo_TG_fourier
Xtwo_CGminusTT_fourier = Xtwo_CG_fourier - Xtwo_TT_fourier

Xtwo_CTminusGA_fourier = Xtwo_CT_fourier - Xtwo_GA_fourier
Xtwo_CTminusGC_fourier = Xtwo_CT_fourier - Xtwo_GC_fourier
Xtwo_CTminusGG_fourier = Xtwo_CT_fourier - Xtwo_GG_fourier
Xtwo_CTminusGT_fourier = Xtwo_CT_fourier - Xtwo_GT_fourier
Xtwo_CTminusTA_fourier = Xtwo_CT_fourier - Xtwo_TA_fourier
Xtwo_CTminusTC_fourier = Xtwo_CT_fourier - Xtwo_TC_fourier
Xtwo_CTminusTG_fourier = Xtwo_CT_fourier - Xtwo_TG_fourier
Xtwo_CTminusTT_fourier = Xtwo_CT_fourier - Xtwo_TT_fourier

Xtwo_GAminusGC_fourier = Xtwo_GA_fourier - Xtwo_GC_fourier
Xtwo_GAminusGG_fourier = Xtwo_GA_fourier - Xtwo_GG_fourier
Xtwo_GAminusGT_fourier = Xtwo_GA_fourier - Xtwo_GT_fourier
Xtwo_GAminusTA_fourier = Xtwo_GA_fourier - Xtwo_TA_fourier
Xtwo_GAminusTC_fourier = Xtwo_GA_fourier - Xtwo_TC_fourier
Xtwo_GAminusTG_fourier = Xtwo_GA_fourier - Xtwo_TG_fourier
Xtwo_GAminusTT_fourier = Xtwo_GA_fourier - Xtwo_TT_fourier

Xtwo_GCminusGG_fourier = Xtwo_GC_fourier - Xtwo_GG_fourier
Xtwo_GCminusGT_fourier = Xtwo_GC_fourier - Xtwo_GT_fourier
Xtwo_GCminusTA_fourier = Xtwo_GC_fourier - Xtwo_TA_fourier
Xtwo_GCminusTC_fourier = Xtwo_GC_fourier - Xtwo_TC_fourier
Xtwo_GCminusTG_fourier = Xtwo_GC_fourier - Xtwo_TG_fourier
Xtwo_GCminusTT_fourier = Xtwo_GC_fourier - Xtwo_TT_fourier

Xtwo_GGminusGT_fourier = Xtwo_GG_fourier - Xtwo_GT_fourier
Xtwo_GGminusTA_fourier = Xtwo_GG_fourier - Xtwo_TA_fourier
Xtwo_GGminusTC_fourier = Xtwo_GG_fourier - Xtwo_TC_fourier
Xtwo_GGminusTG_fourier = Xtwo_GG_fourier - Xtwo_TG_fourier
Xtwo_GGminusTT_fourier = Xtwo_GG_fourier - Xtwo_TT_fourier

Xtwo_GTminusTA_fourier = Xtwo_GT_fourier - Xtwo_TA_fourier
Xtwo_GTminusTC_fourier = Xtwo_GT_fourier - Xtwo_TC_fourier
Xtwo_GTminusTG_fourier = Xtwo_GT_fourier - Xtwo_TG_fourier
Xtwo_GTminusTT_fourier = Xtwo_GT_fourier - Xtwo_TT_fourier

Xtwo_TAminusTC_fourier = Xtwo_TA_fourier - Xtwo_TC_fourier
Xtwo_TAminusTG_fourier = Xtwo_TA_fourier - Xtwo_TG_fourier
Xtwo_TAminusTT_fourier = Xtwo_TA_fourier - Xtwo_TT_fourier

Xtwo_TCminusTG_fourier = Xtwo_TC_fourier - Xtwo_TG_fourier
Xtwo_TCminusTT_fourier = Xtwo_TC_fourier - Xtwo_TT_fourier

Xtwo_TGminusTT_fourier = Xtwo_TG_fourier - Xtwo_TT_fourier

Xtwo_AAminusAC_fourier_mag = as.matrix(Mod(Xtwo_AAminusAC_fourier))
Xtwo_AAminusAC_fourier_mag = as.matrix(Mod(Xtwo_AAminusAC_fourier))
Xtwo_AAminusAG_fourier_mag = as.matrix(Mod(Xtwo_AAminusAG_fourier))
Xtwo_AAminusAT_fourier_mag = as.matrix(Mod(Xtwo_AAminusAT_fourier))
Xtwo_AAminusCA_fourier_mag = as.matrix(Mod(Xtwo_AAminusCA_fourier))
Xtwo_AAminusCC_fourier_mag = as.matrix(Mod(Xtwo_AAminusCC_fourier))
Xtwo_AAminusCG_fourier_mag = as.matrix(Mod(Xtwo_AAminusCG_fourier))
Xtwo_AAminusCT_fourier_mag = as.matrix(Mod(Xtwo_AAminusCT_fourier))
Xtwo_AAminusGA_fourier_mag = as.matrix(Mod(Xtwo_AAminusGA_fourier))
Xtwo_AAminusGC_fourier_mag = as.matrix(Mod(Xtwo_AAminusGC_fourier))
Xtwo_AAminusGG_fourier_mag = as.matrix(Mod(Xtwo_AAminusGG_fourier))
Xtwo_AAminusGT_fourier_mag = as.matrix(Mod(Xtwo_AAminusGT_fourier))
Xtwo_AAminusTA_fourier_mag = as.matrix(Mod(Xtwo_AAminusTA_fourier))
Xtwo_AAminusTC_fourier_mag = as.matrix(Mod(Xtwo_AAminusTC_fourier))
Xtwo_AAminusTG_fourier_mag = as.matrix(Mod(Xtwo_AAminusTG_fourier))
Xtwo_AAminusTT_fourier_mag = as.matrix(Mod(Xtwo_AAminusTT_fourier))

Xtwo_ACminusAG_fourier_mag = as.matrix(Mod(Xtwo_ACminusAG_fourier))
Xtwo_ACminusAT_fourier_mag = as.matrix(Mod(Xtwo_ACminusAT_fourier))
Xtwo_ACminusCA_fourier_mag = as.matrix(Mod(Xtwo_ACminusCA_fourier))
Xtwo_ACminusCC_fourier_mag = as.matrix(Mod(Xtwo_ACminusCC_fourier))
Xtwo_ACminusCG_fourier_mag = as.matrix(Mod(Xtwo_ACminusCG_fourier))
Xtwo_ACminusCT_fourier_mag = as.matrix(Mod(Xtwo_ACminusCT_fourier))
Xtwo_ACminusGA_fourier_mag = as.matrix(Mod(Xtwo_ACminusGA_fourier))
Xtwo_ACminusGC_fourier_mag = as.matrix(Mod(Xtwo_ACminusGC_fourier))
Xtwo_ACminusGG_fourier_mag = as.matrix(Mod(Xtwo_ACminusGG_fourier))
Xtwo_ACminusGT_fourier_mag = as.matrix(Mod(Xtwo_ACminusGT_fourier))
Xtwo_ACminusTA_fourier_mag = as.matrix(Mod(Xtwo_ACminusTA_fourier))
Xtwo_ACminusTC_fourier_mag = as.matrix(Mod(Xtwo_ACminusTC_fourier))
Xtwo_ACminusTG_fourier_mag = as.matrix(Mod(Xtwo_ACminusTG_fourier))
Xtwo_ACminusTT_fourier_mag = as.matrix(Mod(Xtwo_ACminusTT_fourier))

Xtwo_AGminusAT_fourier_mag = as.matrix(Mod(Xtwo_AGminusAT_fourier))
Xtwo_AGminusCA_fourier_mag = as.matrix(Mod(Xtwo_AGminusCA_fourier))
Xtwo_AGminusCC_fourier_mag = as.matrix(Mod(Xtwo_AGminusCC_fourier))
Xtwo_AGminusCG_fourier_mag = as.matrix(Mod(Xtwo_AGminusCG_fourier))
Xtwo_AGminusCT_fourier_mag = as.matrix(Mod(Xtwo_AGminusCT_fourier))
Xtwo_AGminusGA_fourier_mag = as.matrix(Mod(Xtwo_AGminusGA_fourier))
Xtwo_AGminusGC_fourier_mag = as.matrix(Mod(Xtwo_AGminusGC_fourier))
Xtwo_AGminusGG_fourier_mag = as.matrix(Mod(Xtwo_AGminusGG_fourier))
Xtwo_AGminusGT_fourier_mag = as.matrix(Mod(Xtwo_AGminusGT_fourier))
Xtwo_AGminusTA_fourier_mag = as.matrix(Mod(Xtwo_AGminusTA_fourier))
Xtwo_AGminusTC_fourier_mag = as.matrix(Mod(Xtwo_AGminusTC_fourier))
Xtwo_AGminusTG_fourier_mag = as.matrix(Mod(Xtwo_AGminusTG_fourier))
Xtwo_AGminusTT_fourier_mag = as.matrix(Mod(Xtwo_AGminusTT_fourier))

Xtwo_ATminusCA_fourier_mag = as.matrix(Mod(Xtwo_ATminusCA_fourier))
Xtwo_ATminusCC_fourier_mag = as.matrix(Mod(Xtwo_ATminusCC_fourier))
Xtwo_ATminusCG_fourier_mag = as.matrix(Mod(Xtwo_ATminusCG_fourier))
Xtwo_ATminusCT_fourier_mag = as.matrix(Mod(Xtwo_ATminusCT_fourier))
Xtwo_ATminusGA_fourier_mag = as.matrix(Mod(Xtwo_ATminusGA_fourier))
Xtwo_ATminusGC_fourier_mag = as.matrix(Mod(Xtwo_ATminusGC_fourier))
Xtwo_ATminusGG_fourier_mag = as.matrix(Mod(Xtwo_ATminusGG_fourier))
Xtwo_ATminusGT_fourier_mag = as.matrix(Mod(Xtwo_ATminusGT_fourier))
Xtwo_ATminusTA_fourier_mag = as.matrix(Mod(Xtwo_ATminusTA_fourier))
Xtwo_ATminusTC_fourier_mag = as.matrix(Mod(Xtwo_ATminusTC_fourier))
Xtwo_ATminusTG_fourier_mag = as.matrix(Mod(Xtwo_ATminusTG_fourier))
Xtwo_ATminusTT_fourier_mag = as.matrix(Mod(Xtwo_ATminusTT_fourier))

Xtwo_CAminusCC_fourier_mag = as.matrix(Mod(Xtwo_CAminusCC_fourier))
Xtwo_CAminusCG_fourier_mag = as.matrix(Mod(Xtwo_CAminusCG_fourier))
Xtwo_CAminusCT_fourier_mag = as.matrix(Mod(Xtwo_CAminusCT_fourier))
Xtwo_CAminusGA_fourier_mag = as.matrix(Mod(Xtwo_CAminusGA_fourier))
Xtwo_CAminusGC_fourier_mag = as.matrix(Mod(Xtwo_CAminusGC_fourier))
Xtwo_CAminusGG_fourier_mag = as.matrix(Mod(Xtwo_CAminusGG_fourier))
Xtwo_CAminusGT_fourier_mag = as.matrix(Mod(Xtwo_CAminusGT_fourier))
Xtwo_CAminusTA_fourier_mag = as.matrix(Mod(Xtwo_CAminusTA_fourier))
Xtwo_CAminusTC_fourier_mag = as.matrix(Mod(Xtwo_CAminusTC_fourier))
Xtwo_CAminusTG_fourier_mag = as.matrix(Mod(Xtwo_CAminusTG_fourier))
Xtwo_CAminusTT_fourier_mag = as.matrix(Mod(Xtwo_CAminusTT_fourier))

Xtwo_CCminusCG_fourier_mag = as.matrix(Mod(Xtwo_CCminusCG_fourier))
Xtwo_CCminusCT_fourier_mag = as.matrix(Mod(Xtwo_CCminusCT_fourier))
Xtwo_CCminusGA_fourier_mag = as.matrix(Mod(Xtwo_CCminusGA_fourier))
Xtwo_CCminusGC_fourier_mag = as.matrix(Mod(Xtwo_CCminusGC_fourier))
Xtwo_CCminusGG_fourier_mag = as.matrix(Mod(Xtwo_CCminusGG_fourier))
Xtwo_CCminusGT_fourier_mag = as.matrix(Mod(Xtwo_CCminusGT_fourier))
Xtwo_CCminusTA_fourier_mag = as.matrix(Mod(Xtwo_CCminusTA_fourier))
Xtwo_CCminusTC_fourier_mag = as.matrix(Mod(Xtwo_CCminusTC_fourier))
Xtwo_CCminusTG_fourier_mag = as.matrix(Mod(Xtwo_CCminusTG_fourier))
Xtwo_CCminusTT_fourier_mag = as.matrix(Mod(Xtwo_CCminusTT_fourier))

Xtwo_CGminusCT_fourier_mag = as.matrix(Mod(Xtwo_CGminusCT_fourier))
Xtwo_CGminusGA_fourier_mag = as.matrix(Mod(Xtwo_CGminusGA_fourier))
Xtwo_CGminusGC_fourier_mag = as.matrix(Mod(Xtwo_CGminusGC_fourier))
Xtwo_CGminusGG_fourier_mag = as.matrix(Mod(Xtwo_CGminusGG_fourier))
Xtwo_CGminusGT_fourier_mag = as.matrix(Mod(Xtwo_CGminusGT_fourier))
Xtwo_CGminusTA_fourier_mag = as.matrix(Mod(Xtwo_CGminusTA_fourier))
Xtwo_CGminusTC_fourier_mag = as.matrix(Mod(Xtwo_CGminusTC_fourier))
Xtwo_CGminusTG_fourier_mag = as.matrix(Mod(Xtwo_CGminusTG_fourier))
Xtwo_CGminusTT_fourier_mag = as.matrix(Mod(Xtwo_CGminusTT_fourier))

Xtwo_CTminusGA_fourier_mag = as.matrix(Mod(Xtwo_CTminusGA_fourier))
Xtwo_CTminusGC_fourier_mag = as.matrix(Mod(Xtwo_CTminusGC_fourier))
Xtwo_CTminusGG_fourier_mag = as.matrix(Mod(Xtwo_CTminusGG_fourier))
Xtwo_CTminusGT_fourier_mag = as.matrix(Mod(Xtwo_CTminusGT_fourier))
Xtwo_CTminusTA_fourier_mag = as.matrix(Mod(Xtwo_CTminusTA_fourier))
Xtwo_CTminusTC_fourier_mag = as.matrix(Mod(Xtwo_CTminusTC_fourier))
Xtwo_CTminusTG_fourier_mag = as.matrix(Mod(Xtwo_CTminusTG_fourier))
Xtwo_CTminusTT_fourier_mag = as.matrix(Mod(Xtwo_CTminusTT_fourier))

Xtwo_GAminusGC_fourier_mag = as.matrix(Mod(Xtwo_GAminusGC_fourier))
Xtwo_GAminusGG_fourier_mag = as.matrix(Mod(Xtwo_GAminusGG_fourier))
Xtwo_GAminusGT_fourier_mag = as.matrix(Mod(Xtwo_GAminusGT_fourier))
Xtwo_GAminusTA_fourier_mag = as.matrix(Mod(Xtwo_GAminusTA_fourier))
Xtwo_GAminusTC_fourier_mag = as.matrix(Mod(Xtwo_GAminusTC_fourier))
Xtwo_GAminusTG_fourier_mag = as.matrix(Mod(Xtwo_GAminusTG_fourier))
Xtwo_GAminusTT_fourier_mag = as.matrix(Mod(Xtwo_GAminusTT_fourier))

Xtwo_GCminusGG_fourier_mag = as.matrix(Mod(Xtwo_GCminusGG_fourier))
Xtwo_GCminusGT_fourier_mag = as.matrix(Mod(Xtwo_GCminusGT_fourier))
Xtwo_GCminusTA_fourier_mag = as.matrix(Mod(Xtwo_GCminusTA_fourier))
Xtwo_GCminusTC_fourier_mag = as.matrix(Mod(Xtwo_GCminusTC_fourier))
Xtwo_GCminusTG_fourier_mag = as.matrix(Mod(Xtwo_GCminusTG_fourier))
Xtwo_GCminusTT_fourier_mag = as.matrix(Mod(Xtwo_GCminusTT_fourier))

Xtwo_GGminusGT_fourier_mag = as.matrix(Mod(Xtwo_GGminusGT_fourier))
Xtwo_GGminusTA_fourier_mag = as.matrix(Mod(Xtwo_GGminusTA_fourier))
Xtwo_GGminusTC_fourier_mag = as.matrix(Mod(Xtwo_GGminusTC_fourier))
Xtwo_GGminusTG_fourier_mag = as.matrix(Mod(Xtwo_GGminusTG_fourier))
Xtwo_GGminusTT_fourier_mag = as.matrix(Mod(Xtwo_GGminusTT_fourier))

Xtwo_GTminusTA_fourier_mag = as.matrix(Mod(Xtwo_GTminusTA_fourier))
Xtwo_GTminusTC_fourier_mag = as.matrix(Mod(Xtwo_GTminusTC_fourier))
Xtwo_GTminusTG_fourier_mag = as.matrix(Mod(Xtwo_GTminusTG_fourier))
Xtwo_GTminusTT_fourier_mag = as.matrix(Mod(Xtwo_GTminusTT_fourier))

Xtwo_TAminusTC_fourier_mag = as.matrix(Mod(Xtwo_TAminusTC_fourier))
Xtwo_TAminusTG_fourier_mag = as.matrix(Mod(Xtwo_TAminusTG_fourier))
Xtwo_TAminusTT_fourier_mag = as.matrix(Mod(Xtwo_TAminusTT_fourier))

Xtwo_TCminusTG_fourier_mag = as.matrix(Mod(Xtwo_TCminusTG_fourier))
Xtwo_TCminusTT_fourier_mag = as.matrix(Mod(Xtwo_TCminusTT_fourier))

Xtwo_TGminusTT_fourier_mag = as.matrix(Mod(Xtwo_TGminusTT_fourier))

Xtwo_AAminusAC_fourier_phase = as.matrix(Arg(Xtwo_AAminusAC_fourier))
Xtwo_AAminusAC_fourier_phase = as.matrix(Arg(Xtwo_AAminusAC_fourier))
Xtwo_AAminusAG_fourier_phase = as.matrix(Arg(Xtwo_AAminusAG_fourier))
Xtwo_AAminusAT_fourier_phase = as.matrix(Arg(Xtwo_AAminusAT_fourier))
Xtwo_AAminusCA_fourier_phase = as.matrix(Arg(Xtwo_AAminusCA_fourier))
Xtwo_AAminusCC_fourier_phase = as.matrix(Arg(Xtwo_AAminusCC_fourier))
Xtwo_AAminusCG_fourier_phase = as.matrix(Arg(Xtwo_AAminusCG_fourier))
Xtwo_AAminusCT_fourier_phase = as.matrix(Arg(Xtwo_AAminusCT_fourier))
Xtwo_AAminusGA_fourier_phase = as.matrix(Arg(Xtwo_AAminusGA_fourier))
Xtwo_AAminusGC_fourier_phase = as.matrix(Arg(Xtwo_AAminusGC_fourier))
Xtwo_AAminusGG_fourier_phase = as.matrix(Arg(Xtwo_AAminusGG_fourier))
Xtwo_AAminusGT_fourier_phase = as.matrix(Arg(Xtwo_AAminusGT_fourier))
Xtwo_AAminusTA_fourier_phase = as.matrix(Arg(Xtwo_AAminusTA_fourier))
Xtwo_AAminusTC_fourier_phase = as.matrix(Arg(Xtwo_AAminusTC_fourier))
Xtwo_AAminusTG_fourier_phase = as.matrix(Arg(Xtwo_AAminusTG_fourier))
Xtwo_AAminusTT_fourier_phase = as.matrix(Arg(Xtwo_AAminusTT_fourier))

Xtwo_ACminusAG_fourier_phase = as.matrix(Arg(Xtwo_ACminusAG_fourier))
Xtwo_ACminusAT_fourier_phase = as.matrix(Arg(Xtwo_ACminusAT_fourier))
Xtwo_ACminusCA_fourier_phase = as.matrix(Arg(Xtwo_ACminusCA_fourier))
Xtwo_ACminusCC_fourier_phase = as.matrix(Arg(Xtwo_ACminusCC_fourier))
Xtwo_ACminusCG_fourier_phase = as.matrix(Arg(Xtwo_ACminusCG_fourier))
Xtwo_ACminusCT_fourier_phase = as.matrix(Arg(Xtwo_ACminusCT_fourier))
Xtwo_ACminusGA_fourier_phase = as.matrix(Arg(Xtwo_ACminusGA_fourier))
Xtwo_ACminusGC_fourier_phase = as.matrix(Arg(Xtwo_ACminusGC_fourier))
Xtwo_ACminusGG_fourier_phase = as.matrix(Arg(Xtwo_ACminusGG_fourier))
Xtwo_ACminusGT_fourier_phase = as.matrix(Arg(Xtwo_ACminusGT_fourier))
Xtwo_ACminusTA_fourier_phase = as.matrix(Arg(Xtwo_ACminusTA_fourier))
Xtwo_ACminusTC_fourier_phase = as.matrix(Arg(Xtwo_ACminusTC_fourier))
Xtwo_ACminusTG_fourier_phase = as.matrix(Arg(Xtwo_ACminusTG_fourier))
Xtwo_ACminusTT_fourier_phase = as.matrix(Arg(Xtwo_ACminusTT_fourier))

Xtwo_AGminusAT_fourier_phase = as.matrix(Arg(Xtwo_AGminusAT_fourier))
Xtwo_AGminusCA_fourier_phase = as.matrix(Arg(Xtwo_AGminusCA_fourier))
Xtwo_AGminusCC_fourier_phase = as.matrix(Arg(Xtwo_AGminusCC_fourier))
Xtwo_AGminusCG_fourier_phase = as.matrix(Arg(Xtwo_AGminusCG_fourier))
Xtwo_AGminusCT_fourier_phase = as.matrix(Arg(Xtwo_AGminusCT_fourier))
Xtwo_AGminusGA_fourier_phase = as.matrix(Arg(Xtwo_AGminusGA_fourier))
Xtwo_AGminusGC_fourier_phase = as.matrix(Arg(Xtwo_AGminusGC_fourier))
Xtwo_AGminusGG_fourier_phase = as.matrix(Arg(Xtwo_AGminusGG_fourier))
Xtwo_AGminusGT_fourier_phase = as.matrix(Arg(Xtwo_AGminusGT_fourier))
Xtwo_AGminusTA_fourier_phase = as.matrix(Arg(Xtwo_AGminusTA_fourier))
Xtwo_AGminusTC_fourier_phase = as.matrix(Arg(Xtwo_AGminusTC_fourier))
Xtwo_AGminusTG_fourier_phase = as.matrix(Arg(Xtwo_AGminusTG_fourier))
Xtwo_AGminusTT_fourier_phase = as.matrix(Arg(Xtwo_AGminusTT_fourier))

Xtwo_ATminusCA_fourier_phase = as.matrix(Arg(Xtwo_ATminusCA_fourier))
Xtwo_ATminusCC_fourier_phase = as.matrix(Arg(Xtwo_ATminusCC_fourier))
Xtwo_ATminusCG_fourier_phase = as.matrix(Arg(Xtwo_ATminusCG_fourier))
Xtwo_ATminusCT_fourier_phase = as.matrix(Arg(Xtwo_ATminusCT_fourier))
Xtwo_ATminusGA_fourier_phase = as.matrix(Arg(Xtwo_ATminusGA_fourier))
Xtwo_ATminusGC_fourier_phase = as.matrix(Arg(Xtwo_ATminusGC_fourier))
Xtwo_ATminusGG_fourier_phase = as.matrix(Arg(Xtwo_ATminusGG_fourier))
Xtwo_ATminusGT_fourier_phase = as.matrix(Arg(Xtwo_ATminusGT_fourier))
Xtwo_ATminusTA_fourier_phase = as.matrix(Arg(Xtwo_ATminusTA_fourier))
Xtwo_ATminusTC_fourier_phase = as.matrix(Arg(Xtwo_ATminusTC_fourier))
Xtwo_ATminusTG_fourier_phase = as.matrix(Arg(Xtwo_ATminusTG_fourier))
Xtwo_ATminusTT_fourier_phase = as.matrix(Arg(Xtwo_ATminusTT_fourier))

Xtwo_CAminusCC_fourier_phase = as.matrix(Arg(Xtwo_CAminusCC_fourier))
Xtwo_CAminusCG_fourier_phase = as.matrix(Arg(Xtwo_CAminusCG_fourier))
Xtwo_CAminusCT_fourier_phase = as.matrix(Arg(Xtwo_CAminusCT_fourier))
Xtwo_CAminusGA_fourier_phase = as.matrix(Arg(Xtwo_CAminusGA_fourier))
Xtwo_CAminusGC_fourier_phase = as.matrix(Arg(Xtwo_CAminusGC_fourier))
Xtwo_CAminusGG_fourier_phase = as.matrix(Arg(Xtwo_CAminusGG_fourier))
Xtwo_CAminusGT_fourier_phase = as.matrix(Arg(Xtwo_CAminusGT_fourier))
Xtwo_CAminusTA_fourier_phase = as.matrix(Arg(Xtwo_CAminusTA_fourier))
Xtwo_CAminusTC_fourier_phase = as.matrix(Arg(Xtwo_CAminusTC_fourier))
Xtwo_CAminusTG_fourier_phase = as.matrix(Arg(Xtwo_CAminusTG_fourier))
Xtwo_CAminusTT_fourier_phase = as.matrix(Arg(Xtwo_CAminusTT_fourier))

Xtwo_CCminusCG_fourier_phase = as.matrix(Arg(Xtwo_CCminusCG_fourier))
Xtwo_CCminusCT_fourier_phase = as.matrix(Arg(Xtwo_CCminusCT_fourier))
Xtwo_CCminusGA_fourier_phase = as.matrix(Arg(Xtwo_CCminusGA_fourier))
Xtwo_CCminusGC_fourier_phase = as.matrix(Arg(Xtwo_CCminusGC_fourier))
Xtwo_CCminusGG_fourier_phase = as.matrix(Arg(Xtwo_CCminusGG_fourier))
Xtwo_CCminusGT_fourier_phase = as.matrix(Arg(Xtwo_CCminusGT_fourier))
Xtwo_CCminusTA_fourier_phase = as.matrix(Arg(Xtwo_CCminusTA_fourier))
Xtwo_CCminusTC_fourier_phase = as.matrix(Arg(Xtwo_CCminusTC_fourier))
Xtwo_CCminusTG_fourier_phase = as.matrix(Arg(Xtwo_CCminusTG_fourier))
Xtwo_CCminusTT_fourier_phase = as.matrix(Arg(Xtwo_CCminusTT_fourier))

Xtwo_CGminusCT_fourier_phase = as.matrix(Arg(Xtwo_CGminusCT_fourier))
Xtwo_CGminusGA_fourier_phase = as.matrix(Arg(Xtwo_CGminusGA_fourier))
Xtwo_CGminusGC_fourier_phase = as.matrix(Arg(Xtwo_CGminusGC_fourier))
Xtwo_CGminusGG_fourier_phase = as.matrix(Arg(Xtwo_CGminusGG_fourier))
Xtwo_CGminusGT_fourier_phase = as.matrix(Arg(Xtwo_CGminusGT_fourier))
Xtwo_CGminusTA_fourier_phase = as.matrix(Arg(Xtwo_CGminusTA_fourier))
Xtwo_CGminusTC_fourier_phase = as.matrix(Arg(Xtwo_CGminusTC_fourier))
Xtwo_CGminusTG_fourier_phase = as.matrix(Arg(Xtwo_CGminusTG_fourier))
Xtwo_CGminusTT_fourier_phase = as.matrix(Arg(Xtwo_CGminusTT_fourier))

Xtwo_CTminusGA_fourier_phase = as.matrix(Arg(Xtwo_CTminusGA_fourier))
Xtwo_CTminusGC_fourier_phase = as.matrix(Arg(Xtwo_CTminusGC_fourier))
Xtwo_CTminusGG_fourier_phase = as.matrix(Arg(Xtwo_CTminusGG_fourier))
Xtwo_CTminusGT_fourier_phase = as.matrix(Arg(Xtwo_CTminusGT_fourier))
Xtwo_CTminusTA_fourier_phase = as.matrix(Arg(Xtwo_CTminusTA_fourier))
Xtwo_CTminusTC_fourier_phase = as.matrix(Arg(Xtwo_CTminusTC_fourier))
Xtwo_CTminusTG_fourier_phase = as.matrix(Arg(Xtwo_CTminusTG_fourier))
Xtwo_CTminusTT_fourier_phase = as.matrix(Arg(Xtwo_CTminusTT_fourier))

Xtwo_GAminusGC_fourier_phase = as.matrix(Arg(Xtwo_GAminusGC_fourier))
Xtwo_GAminusGG_fourier_phase = as.matrix(Arg(Xtwo_GAminusGG_fourier))
Xtwo_GAminusGT_fourier_phase = as.matrix(Arg(Xtwo_GAminusGT_fourier))
Xtwo_GAminusTA_fourier_phase = as.matrix(Arg(Xtwo_GAminusTA_fourier))
Xtwo_GAminusTC_fourier_phase = as.matrix(Arg(Xtwo_GAminusTC_fourier))
Xtwo_GAminusTG_fourier_phase = as.matrix(Arg(Xtwo_GAminusTG_fourier))
Xtwo_GAminusTT_fourier_phase = as.matrix(Arg(Xtwo_GAminusTT_fourier))

Xtwo_GCminusGG_fourier_phase = as.matrix(Arg(Xtwo_GCminusGG_fourier))
Xtwo_GCminusGT_fourier_phase = as.matrix(Arg(Xtwo_GCminusGT_fourier))
Xtwo_GCminusTA_fourier_phase = as.matrix(Arg(Xtwo_GCminusTA_fourier))
Xtwo_GCminusTC_fourier_phase = as.matrix(Arg(Xtwo_GCminusTC_fourier))
Xtwo_GCminusTG_fourier_phase = as.matrix(Arg(Xtwo_GCminusTG_fourier))
Xtwo_GCminusTT_fourier_phase = as.matrix(Arg(Xtwo_GCminusTT_fourier))

Xtwo_GGminusGT_fourier_phase = as.matrix(Arg(Xtwo_GGminusGT_fourier))
Xtwo_GGminusTA_fourier_phase = as.matrix(Arg(Xtwo_GGminusTA_fourier))
Xtwo_GGminusTC_fourier_phase = as.matrix(Arg(Xtwo_GGminusTC_fourier))
Xtwo_GGminusTG_fourier_phase = as.matrix(Arg(Xtwo_GGminusTG_fourier))
Xtwo_GGminusTT_fourier_phase = as.matrix(Arg(Xtwo_GGminusTT_fourier))

Xtwo_GTminusTA_fourier_phase = as.matrix(Arg(Xtwo_GTminusTA_fourier))
Xtwo_GTminusTC_fourier_phase = as.matrix(Arg(Xtwo_GTminusTC_fourier))
Xtwo_GTminusTG_fourier_phase = as.matrix(Arg(Xtwo_GTminusTG_fourier))
Xtwo_GTminusTT_fourier_phase = as.matrix(Arg(Xtwo_GTminusTT_fourier))

Xtwo_TAminusTC_fourier_phase = as.matrix(Arg(Xtwo_TAminusTC_fourier))
Xtwo_TAminusTG_fourier_phase = as.matrix(Arg(Xtwo_TAminusTG_fourier))
Xtwo_TAminusTT_fourier_phase = as.matrix(Arg(Xtwo_TAminusTT_fourier))

Xtwo_TCminusTG_fourier_phase = as.matrix(Arg(Xtwo_TCminusTG_fourier))
Xtwo_TCminusTT_fourier_phase = as.matrix(Arg(Xtwo_TCminusTT_fourier))

Xtwo_TGminusTT_fourier_phase = as.matrix(Arg(Xtwo_TGminusTT_fourier))

colnames(Xtwo_AAminusAC_fourier_mag) = paste0("fourier_mag_AAminusAC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAG_fourier_mag) = paste0("fourier_mag_AAminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAT_fourier_mag) = paste0("fourier_mag_AAminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCA_fourier_mag) = paste0("fourier_mag_AAminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCC_fourier_mag) = paste0("fourier_mag_AAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCG_fourier_mag) = paste0("fourier_mag_AAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCT_fourier_mag) = paste0("fourier_mag_AAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGA_fourier_mag) = paste0("fourier_mag_AAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGC_fourier_mag) = paste0("fourier_mag_AAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGG_fourier_mag) = paste0("fourier_mag_AAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGT_fourier_mag) = paste0("fourier_mag_AAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTA_fourier_mag) = paste0("fourier_mag_AAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTC_fourier_mag) = paste0("fourier_mag_AAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTG_fourier_mag) = paste0("fourier_mag_AAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTT_fourier_mag) = paste0("fourier_mag_AAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ACminusAG_fourier_mag) = paste0("fourier_mag_ACminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusAT_fourier_mag) = paste0("fourier_mag_ACminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCA_fourier_mag) = paste0("fourier_mag_ACminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCC_fourier_mag) = paste0("fourier_mag_ACminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCG_fourier_mag) = paste0("fourier_mag_ACminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCT_fourier_mag) = paste0("fourier_mag_ACminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGA_fourier_mag) = paste0("fourier_mag_ACminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGC_fourier_mag) = paste0("fourier_mag_ACminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGG_fourier_mag) = paste0("fourier_mag_ACminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGT_fourier_mag) = paste0("fourier_mag_ACminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTA_fourier_mag) = paste0("fourier_mag_ACminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTC_fourier_mag) = paste0("fourier_mag_ACminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTG_fourier_mag) = paste0("fourier_mag_ACminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTT_fourier_mag) = paste0("fourier_mag_ACminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AGminusAT_fourier_mag) = paste0("fourier_mag_AGminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCA_fourier_mag) = paste0("fourier_mag_AGminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCC_fourier_mag) = paste0("fourier_mag_AGminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCG_fourier_mag) = paste0("fourier_mag_AGminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCT_fourier_mag) = paste0("fourier_mag_AGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGA_fourier_mag) = paste0("fourier_mag_AGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGC_fourier_mag) = paste0("fourier_mag_AGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGG_fourier_mag) = paste0("fourier_mag_AGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGT_fourier_mag) = paste0("fourier_mag_AGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTA_fourier_mag) = paste0("fourier_mag_AGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTC_fourier_mag) = paste0("fourier_mag_AGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTG_fourier_mag) = paste0("fourier_mag_AGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTT_fourier_mag) = paste0("fourier_mag_AGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ATminusCA_fourier_mag) = paste0("fourier_mag_ATminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCC_fourier_mag) = paste0("fourier_mag_ATminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCG_fourier_mag) = paste0("fourier_mag_ATminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCT_fourier_mag) = paste0("fourier_mag_ATminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGA_fourier_mag) = paste0("fourier_mag_ATminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGC_fourier_mag) = paste0("fourier_mag_ATminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGG_fourier_mag) = paste0("fourier_mag_ATminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGT_fourier_mag) = paste0("fourier_mag_ATminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTA_fourier_mag) = paste0("fourier_mag_ATminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTC_fourier_mag) = paste0("fourier_mag_ATminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTG_fourier_mag) = paste0("fourier_mag_ATminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTT_fourier_mag) = paste0("fourier_mag_ATminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CAminusCC_fourier_mag) = paste0("fourier_mag_CAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCG_fourier_mag) = paste0("fourier_mag_CAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCT_fourier_mag) = paste0("fourier_mag_CAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGA_fourier_mag) = paste0("fourier_mag_CAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGC_fourier_mag) = paste0("fourier_mag_CAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGG_fourier_mag) = paste0("fourier_mag_CAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGT_fourier_mag) = paste0("fourier_mag_CAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTA_fourier_mag) = paste0("fourier_mag_CAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTC_fourier_mag) = paste0("fourier_mag_CAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTG_fourier_mag) = paste0("fourier_mag_CAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTT_fourier_mag) = paste0("fourier_mag_CAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CCminusCG_fourier_mag) = paste0("fourier_mag_CCminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusCT_fourier_mag) = paste0("fourier_mag_CCminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGA_fourier_mag) = paste0("fourier_mag_CCminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGC_fourier_mag) = paste0("fourier_mag_CCminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGG_fourier_mag) = paste0("fourier_mag_CCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGT_fourier_mag) = paste0("fourier_mag_CCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTA_fourier_mag) = paste0("fourier_mag_CCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTC_fourier_mag) = paste0("fourier_mag_CCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTG_fourier_mag) = paste0("fourier_mag_CCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTT_fourier_mag) = paste0("fourier_mag_CCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CGminusCT_fourier_mag) = paste0("fourier_mag_CGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGA_fourier_mag) = paste0("fourier_mag_CGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGC_fourier_mag) = paste0("fourier_mag_CGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGG_fourier_mag) = paste0("fourier_mag_CGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGT_fourier_mag) = paste0("fourier_mag_CGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTA_fourier_mag) = paste0("fourier_mag_CGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTC_fourier_mag) = paste0("fourier_mag_CGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTG_fourier_mag) = paste0("fourier_mag_CGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTT_fourier_mag) = paste0("fourier_mag_CGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CTminusGA_fourier_mag) = paste0("fourier_mag_CTminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGC_fourier_mag) = paste0("fourier_mag_CTminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGG_fourier_mag) = paste0("fourier_mag_CTminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGT_fourier_mag) = paste0("fourier_mag_CTminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTA_fourier_mag) = paste0("fourier_mag_CTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTC_fourier_mag) = paste0("fourier_mag_CTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTG_fourier_mag) = paste0("fourier_mag_CTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTT_fourier_mag) = paste0("fourier_mag_CTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GAminusGC_fourier_mag) = paste0("fourier_mag_GAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGG_fourier_mag) = paste0("fourier_mag_GAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGT_fourier_mag) = paste0("fourier_mag_GAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTA_fourier_mag) = paste0("fourier_mag_GAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTC_fourier_mag) = paste0("fourier_mag_GAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTG_fourier_mag) = paste0("fourier_mag_GAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTT_fourier_mag) = paste0("fourier_mag_GAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GCminusGG_fourier_mag) = paste0("fourier_mag_GCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusGT_fourier_mag) = paste0("fourier_mag_GCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTA_fourier_mag) = paste0("fourier_mag_GCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTC_fourier_mag) = paste0("fourier_mag_GCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTG_fourier_mag) = paste0("fourier_mag_GCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTT_fourier_mag) = paste0("fourier_mag_GCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GGminusGT_fourier_mag) = paste0("fourier_mag_GGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTA_fourier_mag) = paste0("fourier_mag_GGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTC_fourier_mag) = paste0("fourier_mag_GGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTG_fourier_mag) = paste0("fourier_mag_GGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTT_fourier_mag) = paste0("fourier_mag_GGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GTminusTA_fourier_mag) = paste0("fourier_mag_GTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTC_fourier_mag) = paste0("fourier_mag_GTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTG_fourier_mag) = paste0("fourier_mag_GTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTT_fourier_mag) = paste0("fourier_mag_GTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TAminusTC_fourier_mag) = paste0("fourier_mag_TAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTG_fourier_mag) = paste0("fourier_mag_TAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTT_fourier_mag) = paste0("fourier_mag_TAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TCminusTG_fourier_mag) = paste0("fourier_mag_TCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TCminusTT_fourier_mag) = paste0("fourier_mag_TCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TGminusTT_fourier_mag) = paste0("fourier_mag_TGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAminusAC_fourier_phase) = paste0("fourier_phase_AAminusAC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAG_fourier_phase) = paste0("fourier_phase_AAminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAT_fourier_phase) = paste0("fourier_phase_AAminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCA_fourier_phase) = paste0("fourier_phase_AAminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCC_fourier_phase) = paste0("fourier_phase_AAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCG_fourier_phase) = paste0("fourier_phase_AAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCT_fourier_phase) = paste0("fourier_phase_AAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGA_fourier_phase) = paste0("fourier_phase_AAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGC_fourier_phase) = paste0("fourier_phase_AAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGG_fourier_phase) = paste0("fourier_phase_AAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGT_fourier_phase) = paste0("fourier_phase_AAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTA_fourier_phase) = paste0("fourier_phase_AAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTC_fourier_phase) = paste0("fourier_phase_AAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTG_fourier_phase) = paste0("fourier_phase_AAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTT_fourier_phase) = paste0("fourier_phase_AAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ACminusAG_fourier_phase) = paste0("fourier_phase_ACminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusAT_fourier_phase) = paste0("fourier_phase_ACminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCA_fourier_phase) = paste0("fourier_phase_ACminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCC_fourier_phase) = paste0("fourier_phase_ACminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCG_fourier_phase) = paste0("fourier_phase_ACminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCT_fourier_phase) = paste0("fourier_phase_ACminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGA_fourier_phase) = paste0("fourier_phase_ACminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGC_fourier_phase) = paste0("fourier_phase_ACminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGG_fourier_phase) = paste0("fourier_phase_ACminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGT_fourier_phase) = paste0("fourier_phase_ACminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTA_fourier_phase) = paste0("fourier_phase_ACminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTC_fourier_phase) = paste0("fourier_phase_ACminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTG_fourier_phase) = paste0("fourier_phase_ACminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTT_fourier_phase) = paste0("fourier_phase_ACminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AGminusAT_fourier_phase) = paste0("fourier_phase_AGminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCA_fourier_phase) = paste0("fourier_phase_AGminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCC_fourier_phase) = paste0("fourier_phase_AGminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCG_fourier_phase) = paste0("fourier_phase_AGminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCT_fourier_phase) = paste0("fourier_phase_AGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGA_fourier_phase) = paste0("fourier_phase_AGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGC_fourier_phase) = paste0("fourier_phase_AGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGG_fourier_phase) = paste0("fourier_phase_AGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGT_fourier_phase) = paste0("fourier_phase_AGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTA_fourier_phase) = paste0("fourier_phase_AGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTC_fourier_phase) = paste0("fourier_phase_AGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTG_fourier_phase) = paste0("fourier_phase_AGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTT_fourier_phase) = paste0("fourier_phase_AGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ATminusCA_fourier_phase) = paste0("fourier_phase_ATminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCC_fourier_phase) = paste0("fourier_phase_ATminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCG_fourier_phase) = paste0("fourier_phase_ATminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCT_fourier_phase) = paste0("fourier_phase_ATminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGA_fourier_phase) = paste0("fourier_phase_ATminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGC_fourier_phase) = paste0("fourier_phase_ATminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGG_fourier_phase) = paste0("fourier_phase_ATminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGT_fourier_phase) = paste0("fourier_phase_ATminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTA_fourier_phase) = paste0("fourier_phase_ATminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTC_fourier_phase) = paste0("fourier_phase_ATminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTG_fourier_phase) = paste0("fourier_phase_ATminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTT_fourier_phase) = paste0("fourier_phase_ATminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CAminusCC_fourier_phase) = paste0("fourier_phase_CAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCG_fourier_phase) = paste0("fourier_phase_CAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCT_fourier_phase) = paste0("fourier_phase_CAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGA_fourier_phase) = paste0("fourier_phase_CAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGC_fourier_phase) = paste0("fourier_phase_CAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGG_fourier_phase) = paste0("fourier_phase_CAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGT_fourier_phase) = paste0("fourier_phase_CAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTA_fourier_phase) = paste0("fourier_phase_CAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTC_fourier_phase) = paste0("fourier_phase_CAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTG_fourier_phase) = paste0("fourier_phase_CAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTT_fourier_phase) = paste0("fourier_phase_CAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CCminusCG_fourier_phase) = paste0("fourier_phase_CCminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusCT_fourier_phase) = paste0("fourier_phase_CCminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGA_fourier_phase) = paste0("fourier_phase_CCminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGC_fourier_phase) = paste0("fourier_phase_CCminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGG_fourier_phase) = paste0("fourier_phase_CCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGT_fourier_phase) = paste0("fourier_phase_CCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTA_fourier_phase) = paste0("fourier_phase_CCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTC_fourier_phase) = paste0("fourier_phase_CCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTG_fourier_phase) = paste0("fourier_phase_CCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTT_fourier_phase) = paste0("fourier_phase_CCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CGminusCT_fourier_phase) = paste0("fourier_phase_CGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGA_fourier_phase) = paste0("fourier_phase_CGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGC_fourier_phase) = paste0("fourier_phase_CGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGG_fourier_phase) = paste0("fourier_phase_CGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGT_fourier_phase) = paste0("fourier_phase_CGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTA_fourier_phase) = paste0("fourier_phase_CGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTC_fourier_phase) = paste0("fourier_phase_CGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTG_fourier_phase) = paste0("fourier_phase_CGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTT_fourier_phase) = paste0("fourier_phase_CGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CTminusGA_fourier_phase) = paste0("fourier_phase_CTminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGC_fourier_phase) = paste0("fourier_phase_CTminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGG_fourier_phase) = paste0("fourier_phase_CTminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGT_fourier_phase) = paste0("fourier_phase_CTminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTA_fourier_phase) = paste0("fourier_phase_CTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTC_fourier_phase) = paste0("fourier_phase_CTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTG_fourier_phase) = paste0("fourier_phase_CTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTT_fourier_phase) = paste0("fourier_phase_CTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GAminusGC_fourier_phase) = paste0("fourier_phase_GAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGG_fourier_phase) = paste0("fourier_phase_GAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGT_fourier_phase) = paste0("fourier_phase_GAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTA_fourier_phase) = paste0("fourier_phase_GAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTC_fourier_phase) = paste0("fourier_phase_GAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTG_fourier_phase) = paste0("fourier_phase_GAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTT_fourier_phase) = paste0("fourier_phase_GAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GCminusGG_fourier_phase) = paste0("fourier_phase_GCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusGT_fourier_phase) = paste0("fourier_phase_GCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTA_fourier_phase) = paste0("fourier_phase_GCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTC_fourier_phase) = paste0("fourier_phase_GCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTG_fourier_phase) = paste0("fourier_phase_GCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTT_fourier_phase) = paste0("fourier_phase_GCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GGminusGT_fourier_phase) = paste0("fourier_phase_GGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTA_fourier_phase) = paste0("fourier_phase_GGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTC_fourier_phase) = paste0("fourier_phase_GGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTG_fourier_phase) = paste0("fourier_phase_GGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTT_fourier_phase) = paste0("fourier_phase_GGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GTminusTA_fourier_phase) = paste0("fourier_phase_GTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTC_fourier_phase) = paste0("fourier_phase_GTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTG_fourier_phase) = paste0("fourier_phase_GTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTT_fourier_phase) = paste0("fourier_phase_GTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TAminusTC_fourier_phase) = paste0("fourier_phase_TAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTG_fourier_phase) = paste0("fourier_phase_TAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTT_fourier_phase) = paste0("fourier_phase_TAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TCminusTG_fourier_phase) = paste0("fourier_phase_TCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TCminusTT_fourier_phase) = paste0("fourier_phase_TCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TGminusTT_fourier_phase) = paste0("fourier_phase_TGminusTT_", c(1,2,3,4,5,6,11,21))

Xtwo_diminusdi_fourier_mag = cbind(Xtwo_AAminusAC_fourier_mag,
                                   Xtwo_AAminusAG_fourier_mag,
                                   Xtwo_AAminusAT_fourier_mag,
                                   Xtwo_AAminusCA_fourier_mag,
                                   Xtwo_AAminusCC_fourier_mag,
                                   Xtwo_AAminusCG_fourier_mag,
                                   Xtwo_AAminusCT_fourier_mag,
                                   Xtwo_AAminusGA_fourier_mag,
                                   Xtwo_AAminusGC_fourier_mag,
                                   Xtwo_AAminusGG_fourier_mag,
                                   Xtwo_AAminusGT_fourier_mag,
                                   Xtwo_AAminusTA_fourier_mag,
                                   Xtwo_AAminusTC_fourier_mag,
                                   Xtwo_AAminusTG_fourier_mag,
                                   Xtwo_AAminusTT_fourier_mag,
                                   
                                   Xtwo_ACminusAG_fourier_mag,
                                   Xtwo_ACminusAT_fourier_mag,
                                   Xtwo_ACminusCA_fourier_mag,
                                   Xtwo_ACminusCC_fourier_mag,
                                   Xtwo_ACminusCG_fourier_mag,
                                   Xtwo_ACminusCT_fourier_mag,
                                   Xtwo_ACminusGA_fourier_mag,
                                   Xtwo_ACminusGC_fourier_mag,
                                   Xtwo_ACminusGG_fourier_mag,
                                   Xtwo_ACminusGT_fourier_mag,
                                   Xtwo_ACminusTA_fourier_mag,
                                   Xtwo_ACminusTC_fourier_mag,
                                   Xtwo_ACminusTG_fourier_mag,
                                   Xtwo_ACminusTT_fourier_mag,
                                   
                                   Xtwo_AGminusAT_fourier_mag,
                                   Xtwo_AGminusCA_fourier_mag,
                                   Xtwo_AGminusCC_fourier_mag,
                                   Xtwo_AGminusCG_fourier_mag,
                                   Xtwo_AGminusCT_fourier_mag,
                                   Xtwo_AGminusGA_fourier_mag,
                                   Xtwo_AGminusGC_fourier_mag,
                                   Xtwo_AGminusGG_fourier_mag,
                                   Xtwo_AGminusGT_fourier_mag,
                                   Xtwo_AGminusTA_fourier_mag,
                                   Xtwo_AGminusTC_fourier_mag,
                                   Xtwo_AGminusTG_fourier_mag,
                                   Xtwo_AGminusTT_fourier_mag,
                                   
                                   Xtwo_ATminusCA_fourier_mag,
                                   Xtwo_ATminusCC_fourier_mag,
                                   Xtwo_ATminusCG_fourier_mag,
                                   Xtwo_ATminusCT_fourier_mag,
                                   Xtwo_ATminusGA_fourier_mag,
                                   Xtwo_ATminusGC_fourier_mag,
                                   Xtwo_ATminusGG_fourier_mag,
                                   Xtwo_ATminusGT_fourier_mag,
                                   Xtwo_ATminusTA_fourier_mag,
                                   Xtwo_ATminusTC_fourier_mag,
                                   Xtwo_ATminusTG_fourier_mag,
                                   Xtwo_ATminusTT_fourier_mag,
                                   
                                   Xtwo_CAminusCC_fourier_mag,
                                   Xtwo_CAminusCG_fourier_mag,
                                   Xtwo_CAminusCT_fourier_mag,
                                   Xtwo_CAminusGA_fourier_mag,
                                   Xtwo_CAminusGC_fourier_mag,
                                   Xtwo_CAminusGG_fourier_mag,
                                   Xtwo_CAminusGT_fourier_mag,
                                   Xtwo_CAminusTA_fourier_mag,
                                   Xtwo_CAminusTC_fourier_mag,
                                   Xtwo_CAminusTG_fourier_mag,
                                   Xtwo_CAminusTT_fourier_mag,
                                   
                                   Xtwo_CCminusCG_fourier_mag,
                                   Xtwo_CCminusCT_fourier_mag,
                                   Xtwo_CCminusGA_fourier_mag,
                                   Xtwo_CCminusGC_fourier_mag,
                                   Xtwo_CCminusGG_fourier_mag,
                                   Xtwo_CCminusGT_fourier_mag,
                                   Xtwo_CCminusTA_fourier_mag,
                                   Xtwo_CCminusTC_fourier_mag,
                                   Xtwo_CCminusTG_fourier_mag,
                                   Xtwo_CCminusTT_fourier_mag,
                                   
                                   Xtwo_CGminusCT_fourier_mag,
                                   Xtwo_CGminusGA_fourier_mag,
                                   Xtwo_CGminusGC_fourier_mag,
                                   Xtwo_CGminusGG_fourier_mag,
                                   Xtwo_CGminusGT_fourier_mag,
                                   Xtwo_CGminusTA_fourier_mag,
                                   Xtwo_CGminusTC_fourier_mag,
                                   Xtwo_CGminusTG_fourier_mag,
                                   Xtwo_CGminusTT_fourier_mag,
                                   
                                   Xtwo_CTminusGA_fourier_mag,
                                   Xtwo_CTminusGC_fourier_mag,
                                   Xtwo_CTminusGG_fourier_mag,
                                   Xtwo_CTminusGT_fourier_mag,
                                   Xtwo_CTminusTA_fourier_mag,
                                   Xtwo_CTminusTC_fourier_mag,
                                   Xtwo_CTminusTG_fourier_mag,
                                   Xtwo_CTminusTT_fourier_mag,
                                   
                                   Xtwo_GAminusGC_fourier_mag,
                                   Xtwo_GAminusGG_fourier_mag,
                                   Xtwo_GAminusGT_fourier_mag,
                                   Xtwo_GAminusTA_fourier_mag,
                                   Xtwo_GAminusTC_fourier_mag,
                                   Xtwo_GAminusTG_fourier_mag,
                                   Xtwo_GAminusTT_fourier_mag,
                                   
                                   Xtwo_GCminusGG_fourier_mag,
                                   Xtwo_GCminusGT_fourier_mag,
                                   Xtwo_GCminusTA_fourier_mag,
                                   Xtwo_GCminusTC_fourier_mag,
                                   Xtwo_GCminusTG_fourier_mag,
                                   Xtwo_GCminusTT_fourier_mag,
                                   
                                   Xtwo_GGminusGT_fourier_mag,
                                   Xtwo_GGminusTA_fourier_mag,
                                   Xtwo_GGminusTC_fourier_mag,
                                   Xtwo_GGminusTG_fourier_mag,
                                   Xtwo_GGminusTT_fourier_mag,
                                   
                                   Xtwo_GTminusTA_fourier_mag,
                                   Xtwo_GTminusTC_fourier_mag,
                                   Xtwo_GTminusTG_fourier_mag,
                                   Xtwo_GTminusTT_fourier_mag,
                                   
                                   Xtwo_TAminusTC_fourier_mag,
                                   Xtwo_TAminusTG_fourier_mag,
                                   Xtwo_TAminusTT_fourier_mag,
                                   
                                   Xtwo_TCminusTG_fourier_mag,
                                   Xtwo_TCminusTT_fourier_mag,
                                   
                                   Xtwo_TGminusTT_fourier_mag)

Xtwo_diminusdi_fourier_phase = cbind(Xtwo_AAminusAC_fourier_phase,
                                     Xtwo_AAminusAG_fourier_phase,
                                     Xtwo_AAminusAT_fourier_phase,
                                     Xtwo_AAminusCA_fourier_phase,
                                     Xtwo_AAminusCC_fourier_phase,
                                     Xtwo_AAminusCG_fourier_phase,
                                     Xtwo_AAminusCT_fourier_phase,
                                     Xtwo_AAminusGA_fourier_phase,
                                     Xtwo_AAminusGC_fourier_phase,
                                     Xtwo_AAminusGG_fourier_phase,
                                     Xtwo_AAminusGT_fourier_phase,
                                     Xtwo_AAminusTA_fourier_phase,
                                     Xtwo_AAminusTC_fourier_phase,
                                     Xtwo_AAminusTG_fourier_phase,
                                     Xtwo_AAminusTT_fourier_phase,
                                     
                                     Xtwo_ACminusAG_fourier_phase,
                                     Xtwo_ACminusAT_fourier_phase,
                                     Xtwo_ACminusCA_fourier_phase,
                                     Xtwo_ACminusCC_fourier_phase,
                                     Xtwo_ACminusCG_fourier_phase,
                                     Xtwo_ACminusCT_fourier_phase,
                                     Xtwo_ACminusGA_fourier_phase,
                                     Xtwo_ACminusGC_fourier_phase,
                                     Xtwo_ACminusGG_fourier_phase,
                                     Xtwo_ACminusGT_fourier_phase,
                                     Xtwo_ACminusTA_fourier_phase,
                                     Xtwo_ACminusTC_fourier_phase,
                                     Xtwo_ACminusTG_fourier_phase,
                                     Xtwo_ACminusTT_fourier_phase,
                                     
                                     Xtwo_AGminusAT_fourier_phase,
                                     Xtwo_AGminusCA_fourier_phase,
                                     Xtwo_AGminusCC_fourier_phase,
                                     Xtwo_AGminusCG_fourier_phase,
                                     Xtwo_AGminusCT_fourier_phase,
                                     Xtwo_AGminusGA_fourier_phase,
                                     Xtwo_AGminusGC_fourier_phase,
                                     Xtwo_AGminusGG_fourier_phase,
                                     Xtwo_AGminusGT_fourier_phase,
                                     Xtwo_AGminusTA_fourier_phase,
                                     Xtwo_AGminusTC_fourier_phase,
                                     Xtwo_AGminusTG_fourier_phase,
                                     Xtwo_AGminusTT_fourier_phase,
                                     
                                     Xtwo_ATminusCA_fourier_phase,
                                     Xtwo_ATminusCC_fourier_phase,
                                     Xtwo_ATminusCG_fourier_phase,
                                     Xtwo_ATminusCT_fourier_phase,
                                     Xtwo_ATminusGA_fourier_phase,
                                     Xtwo_ATminusGC_fourier_phase,
                                     Xtwo_ATminusGG_fourier_phase,
                                     Xtwo_ATminusGT_fourier_phase,
                                     Xtwo_ATminusTA_fourier_phase,
                                     Xtwo_ATminusTC_fourier_phase,
                                     Xtwo_ATminusTG_fourier_phase,
                                     Xtwo_ATminusTT_fourier_phase,
                                     
                                     Xtwo_CAminusCC_fourier_phase,
                                     Xtwo_CAminusCG_fourier_phase,
                                     Xtwo_CAminusCT_fourier_phase,
                                     Xtwo_CAminusGA_fourier_phase,
                                     Xtwo_CAminusGC_fourier_phase,
                                     Xtwo_CAminusGG_fourier_phase,
                                     Xtwo_CAminusGT_fourier_phase,
                                     Xtwo_CAminusTA_fourier_phase,
                                     Xtwo_CAminusTC_fourier_phase,
                                     Xtwo_CAminusTG_fourier_phase,
                                     Xtwo_CAminusTT_fourier_phase,
                                     
                                     Xtwo_CCminusCG_fourier_phase,
                                     Xtwo_CCminusCT_fourier_phase,
                                     Xtwo_CCminusGA_fourier_phase,
                                     Xtwo_CCminusGC_fourier_phase,
                                     Xtwo_CCminusGG_fourier_phase,
                                     Xtwo_CCminusGT_fourier_phase,
                                     Xtwo_CCminusTA_fourier_phase,
                                     Xtwo_CCminusTC_fourier_phase,
                                     Xtwo_CCminusTG_fourier_phase,
                                     Xtwo_CCminusTT_fourier_phase,
                                     
                                     Xtwo_CGminusCT_fourier_phase,
                                     Xtwo_CGminusGA_fourier_phase,
                                     Xtwo_CGminusGC_fourier_phase,
                                     Xtwo_CGminusGG_fourier_phase,
                                     Xtwo_CGminusGT_fourier_phase,
                                     Xtwo_CGminusTA_fourier_phase,
                                     Xtwo_CGminusTC_fourier_phase,
                                     Xtwo_CGminusTG_fourier_phase,
                                     Xtwo_CGminusTT_fourier_phase,
                                     
                                     Xtwo_CTminusGA_fourier_phase,
                                     Xtwo_CTminusGC_fourier_phase,
                                     Xtwo_CTminusGG_fourier_phase,
                                     Xtwo_CTminusGT_fourier_phase,
                                     Xtwo_CTminusTA_fourier_phase,
                                     Xtwo_CTminusTC_fourier_phase,
                                     Xtwo_CTminusTG_fourier_phase,
                                     Xtwo_CTminusTT_fourier_phase,
                                     
                                     Xtwo_GAminusGC_fourier_phase,
                                     Xtwo_GAminusGG_fourier_phase,
                                     Xtwo_GAminusGT_fourier_phase,
                                     Xtwo_GAminusTA_fourier_phase,
                                     Xtwo_GAminusTC_fourier_phase,
                                     Xtwo_GAminusTG_fourier_phase,
                                     Xtwo_GAminusTT_fourier_phase,
                                     
                                     Xtwo_GCminusGG_fourier_phase,
                                     Xtwo_GCminusGT_fourier_phase,
                                     Xtwo_GCminusTA_fourier_phase,
                                     Xtwo_GCminusTC_fourier_phase,
                                     Xtwo_GCminusTG_fourier_phase,
                                     Xtwo_GCminusTT_fourier_phase,
                                     
                                     Xtwo_GGminusGT_fourier_phase,
                                     Xtwo_GGminusTA_fourier_phase,
                                     Xtwo_GGminusTC_fourier_phase,
                                     Xtwo_GGminusTG_fourier_phase,
                                     Xtwo_GGminusTT_fourier_phase,
                                     
                                     Xtwo_GTminusTA_fourier_phase,
                                     Xtwo_GTminusTC_fourier_phase,
                                     Xtwo_GTminusTG_fourier_phase,
                                     Xtwo_GTminusTT_fourier_phase,
                                     
                                     Xtwo_TAminusTC_fourier_phase,
                                     Xtwo_TAminusTG_fourier_phase,
                                     Xtwo_TAminusTT_fourier_phase,
                                     
                                     Xtwo_TCminusTG_fourier_phase,
                                     Xtwo_TCminusTT_fourier_phase,
                                     
                                     Xtwo_TGminusTT_fourier_phase)

saveRDS(Xtwo_diminusdi_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diminusdi_fourier_mag.rds")

saveRDS(Xtwo_diminusdi_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diminusdi_fourier_phase.rds")


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

Xtwo_AAorAT_fourier = t(apply(Xtwo_AAorAT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AAorTA_fourier = t(apply(Xtwo_AAorTA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AAorTT_fourier = t(apply(Xtwo_AAorTT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_ATorTA_fourier = t(apply(Xtwo_ATorTA, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_ATorTT_fourier = t(apply(Xtwo_ATorTT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TAorTT_fourier = t(apply(Xtwo_TAorTT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorCG_fourier = t(apply(Xtwo_CCorCG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorGC_fourier = t(apply(Xtwo_CCorGC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorGG_fourier = t(apply(Xtwo_CCorGG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CGorGC_fourier = t(apply(Xtwo_CGorGC, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CGorGG_fourier = t(apply(Xtwo_CGorGG, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GCorGG_fourier = t(apply(Xtwo_GCorGG, 1, fft))[,c(1,2,3,4,5,6,11,21)]

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


colnames(Xtwo_AAorAT_fourier_mag) = paste0("fourier_mag_AAorAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTA_fourier_mag) = paste0("fourier_mag_AAorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTT_fourier_mag) = paste0("fourier_mag_AAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTA_fourier_mag) = paste0("fourier_mag_ATorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTT_fourier_mag) = paste0("fourier_mag_ATorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAorTT_fourier_mag) = paste0("fourier_mag_TAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCG_fourier_mag) = paste0("fourier_mag_CCorCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGC_fourier_mag) = paste0("fourier_mag_CCorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGG_fourier_mag) = paste0("fourier_mag_CCorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGC_fourier_mag) = paste0("fourier_mag_CGorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGG_fourier_mag) = paste0("fourier_mag_CGorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCorGG_fourier_mag) = paste0("fourier_mag_GCorGG_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAorAT_fourier_phase) = paste0("fourier_phase_AAorAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTA_fourier_phase) = paste0("fourier_phase_AAorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTT_fourier_phase) = paste0("fourier_phase_AAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTA_fourier_phase) = paste0("fourier_phase_ATorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTT_fourier_phase) = paste0("fourier_phase_ATorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAorTT_fourier_phase) = paste0("fourier_phase_TAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCG_fourier_phase) = paste0("fourier_phase_CCorCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGC_fourier_phase) = paste0("fourier_phase_CCorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGG_fourier_phase) = paste0("fourier_phase_CCorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGC_fourier_phase) = paste0("fourier_phase_CGorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGG_fourier_phase) = paste0("fourier_phase_CGorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCorGG_fourier_phase) = paste0("fourier_phase_GCorGG_", c(1,2,3,4,5,6,11,21))

Xtwo_diordi_fourier_mag = cbind(Xtwo_AAorAT_fourier_mag,
                                Xtwo_AAorTA_fourier_mag,
                                Xtwo_AAorTT_fourier_mag,
                                Xtwo_ATorTA_fourier_mag,
                                Xtwo_ATorTT_fourier_mag,
                                Xtwo_TAorTT_fourier_mag,
                                Xtwo_CCorCG_fourier_mag,
                                Xtwo_CCorGC_fourier_mag,
                                Xtwo_CCorGG_fourier_mag,
                                Xtwo_CGorGC_fourier_mag,
                                Xtwo_CGorGG_fourier_mag,
                                Xtwo_GCorGG_fourier_mag)
saveRDS(Xtwo_diordi_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordi_fourier_mag.rds")

Xtwo_diordi_fourier_phase = cbind(Xtwo_AAorAT_fourier_phase,
                                  Xtwo_AAorTA_fourier_phase,
                                  Xtwo_AAorTT_fourier_phase,
                                  Xtwo_ATorTA_fourier_phase,
                                  Xtwo_ATorTT_fourier_phase,
                                  Xtwo_TAorTT_fourier_phase,
                                  Xtwo_CCorCG_fourier_phase,
                                  Xtwo_CCorGC_fourier_phase,
                                  Xtwo_CCorGG_fourier_phase,
                                  Xtwo_CGorGC_fourier_phase,
                                  Xtwo_CGorGG_fourier_phase,
                                  Xtwo_GCorGG_fourier_phase)
saveRDS(Xtwo_diordi_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordi_fourier_phase.rds")


# Groups of Dinucleotides:

Xtwo_AAorATorTAorTT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CCorCGorGCorGG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AAorATorTAorTT) = colnames(Xtwo)
colnames(Xtwo_CCorCGorGCorGG) = colnames(Xtwo)
Xtwo_AAorATorTAorTT[] = ((Xtwo == "AA") | (Xtwo == "AT") | 
                           (Xtwo == "TA") | (Xtwo == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG[] = ((Xtwo == "CC") | (Xtwo == "CG") |
                           (Xtwo == "GC") | (Xtwo == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_fourier = t(apply(Xtwo_AAorATorTAorTT, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorCGorGCorGG_fourier = t(apply(Xtwo_CCorCGorGCorGG, 1, fft))[,c(1,2,3,4,5,6,11,21)]

Xtwo_AAorATorTAorTT_fourier_mag = as.matrix(Mod(Xtwo_AAorATorTAorTT_fourier))
Xtwo_CCorCGorGCorGG_fourier_mag = as.matrix(Mod(Xtwo_CCorCGorGCorGG_fourier))

Xtwo_AAorATorTAorTT_fourier_phase = as.matrix(Arg(Xtwo_AAorATorTAorTT_fourier))
Xtwo_CCorCGorGCorGG_fourier_phase = as.matrix(Arg(Xtwo_CCorCGorGCorGG_fourier))


colnames(Xtwo_AAorATorTAorTT_fourier_mag) = paste0("fourier_mag_AAorATorTAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCGorGCorGG_fourier_mag) = paste0("fourier_mag_CCorCGorGCorGG_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAorATorTAorTT_fourier_phase) = paste0("fourier_phase_AAorATorTAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCGorGCorGG_fourier_phase) = paste0("fourier_phase_CCorCGorGCorGG_", c(1,2,3,4,5,6,11,21))

Xtwo_diordiordiordi_fourier_mag = cbind(Xtwo_AAorATorTAorTT_fourier_mag,
                                        Xtwo_CCorCGorGCorGG_fourier_mag)
saveRDS(Xtwo_diordiordiordi_fourier_mag, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordiordiordi_fourier_mag.rds")

Xtwo_diordiordiordi_fourier_phase = cbind(Xtwo_AAorATorTAorTT_fourier_phase,
                                          Xtwo_CCorCGorGCorGG_fourier_phase)
saveRDS(Xtwo_diordiordiordi_fourier_phase, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordiordiordi_fourier_phase.rds")




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

# Xone_A_fourier = t(apply(Xone_A, 1, fft))[,1:26]
# Xone_C_fourier = t(apply(Xone_C, 1, fft))[,1:26]
# Xone_G_fourier = t(apply(Xone_G, 1, fft))[,1:26]
# Xone_T_fourier = t(apply(Xone_T, 1, fft))[,1:26]

Xone_A_fourier_test = t(apply(Xone_A_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_C_fourier_test = t(apply(Xone_C_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_G_fourier_test = t(apply(Xone_G_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_T_fourier_test = t(apply(Xone_T_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]

Xone_A_fourier_mag_test = as.matrix(Mod(Xone_A_fourier_test))
Xone_C_fourier_mag_test = as.matrix(Mod(Xone_C_fourier_test))
Xone_G_fourier_mag_test = as.matrix(Mod(Xone_G_fourier_test))
Xone_T_fourier_mag_test = as.matrix(Mod(Xone_T_fourier_test))

Xone_A_fourier_phase_test = as.matrix(Arg(Xone_A_fourier_test))
Xone_C_fourier_phase_test = as.matrix(Arg(Xone_C_fourier_test))
Xone_G_fourier_phase_test = as.matrix(Arg(Xone_G_fourier_test))
Xone_T_fourier_phase_test = as.matrix(Arg(Xone_T_fourier_test))

# colnames(Xone_A_fourier_mag) = paste0("fourier_mag_A_", 1:26)
# colnames(Xone_C_fourier_mag) = paste0("fourier_mag_C_", 1:26)
# colnames(Xone_G_fourier_mag) = paste0("fourier_mag_G_", 1:26)
# colnames(Xone_T_fourier_mag) = paste0("fourier_mag_T_", 1:26)

# colnames(Xone_A_fourier_phase) = paste0("fourier_phase_A_", 1:26)
# colnames(Xone_C_fourier_phase) = paste0("fourier_phase_C_", 1:26)
# colnames(Xone_G_fourier_phase) = paste0("fourier_phase_G_", 1:26)
# colnames(Xone_T_fourier_phase) = paste0("fourier_phase_T_", 1:26)

colnames(Xone_A_fourier_mag_test) = paste0("fourier_mag_A_", c(1,2,3,4,5,6,11,21))
colnames(Xone_C_fourier_mag_test) = paste0("fourier_mag_C_", c(1,2,3,4,5,6,11,21))
colnames(Xone_G_fourier_mag_test) = paste0("fourier_mag_G_", c(1,2,3,4,5,6,11,21))
colnames(Xone_T_fourier_mag_test) = paste0("fourier_mag_T_", c(1,2,3,4,5,6,11,21))

colnames(Xone_A_fourier_phase_test) = paste0("fourier_phase_A_", c(1,2,3,4,5,6,11,21))
colnames(Xone_C_fourier_phase_test) = paste0("fourier_phase_C_", c(1,2,3,4,5,6,11,21))
colnames(Xone_G_fourier_phase_test) = paste0("fourier_phase_G_", c(1,2,3,4,5,6,11,21))
colnames(Xone_T_fourier_phase_test) = paste0("fourier_phase_T_", c(1,2,3,4,5,6,11,21))

Xone_nuc_fourier_mag_test = cbind(Xone_A_fourier_mag_test, Xone_C_fourier_mag_test,
                                  Xone_G_fourier_mag_test, Xone_T_fourier_mag_test)
saveRDS(Xone_nuc_fourier_mag_test, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_mag_test.rds")

Xone_nuc_fourier_phase_test = cbind(Xone_A_fourier_phase_test, Xone_C_fourier_phase_test,
                                    Xone_G_fourier_phase_test, Xone_T_fourier_phase_test)
saveRDS(Xone_nuc_fourier_phase_test, "data/Created/chrV_post_smooth_C0_Xone_nuc_fourier_phase_test.rds")



# Differences between single nucleotides:

Xone_AminusC_fourier_test = Xone_A_fourier_test - Xone_C_fourier_test
Xone_AminusG_fourier_test = Xone_A_fourier_test - Xone_G_fourier_test
Xone_AminusT_fourier_test = Xone_A_fourier_test - Xone_T_fourier_test
Xone_CminusG_fourier_test = Xone_C_fourier_test - Xone_G_fourier_test
Xone_CminusT_fourier_test = Xone_C_fourier_test - Xone_T_fourier_test
Xone_GminusT_fourier_test = Xone_G_fourier_test - Xone_T_fourier_test

Xone_AminusC_fourier_mag_test = as.matrix(Mod(Xone_AminusC_fourier_test))
Xone_AminusG_fourier_mag_test = as.matrix(Mod(Xone_AminusG_fourier_test))
Xone_AminusT_fourier_mag_test = as.matrix(Mod(Xone_AminusT_fourier_test))
Xone_CminusG_fourier_mag_test = as.matrix(Mod(Xone_CminusG_fourier_test))
Xone_CminusT_fourier_mag_test = as.matrix(Mod(Xone_CminusT_fourier_test))
Xone_GminusT_fourier_mag_test = as.matrix(Mod(Xone_GminusT_fourier_test))

Xone_AminusC_fourier_phase_test = as.matrix(Arg(Xone_AminusC_fourier_test))
Xone_AminusG_fourier_phase_test = as.matrix(Arg(Xone_AminusG_fourier_test))
Xone_AminusT_fourier_phase_test = as.matrix(Arg(Xone_AminusT_fourier_test))
Xone_CminusG_fourier_phase_test = as.matrix(Arg(Xone_CminusG_fourier_test))
Xone_CminusT_fourier_phase_test = as.matrix(Arg(Xone_CminusT_fourier_test))
Xone_GminusT_fourier_phase_test = as.matrix(Arg(Xone_GminusT_fourier_test))

colnames(Xone_AminusC_fourier_mag_test) = paste0("fourier_mag_AminusC_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AminusG_fourier_mag_test) = paste0("fourier_mag_AminusG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AminusT_fourier_mag_test) = paste0("fourier_mag_AminusT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CminusG_fourier_mag_test) = paste0("fourier_mag_CminusG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CminusT_fourier_mag_test) = paste0("fourier_mag_CminusT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_GminusT_fourier_mag_test) = paste0("fourier_mag_GminusT_", c(1,2,3,4,5,6,11,21))

colnames(Xone_AminusC_fourier_phase_test) = paste0("fourier_phase_AminusC_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AminusG_fourier_phase_test) = paste0("fourier_phase_AminusG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AminusT_fourier_phase_test) = paste0("fourier_phase_AminusT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CminusG_fourier_phase_test) = paste0("fourier_phase_CminusG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CminusT_fourier_phase_test) = paste0("fourier_phase_CminusT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_GminusT_fourier_phase_test) = paste0("fourier_phase_GminusT_", c(1,2,3,4,5,6,11,21))

Xone_nucminusnuc_fourier_mag_test = cbind(Xone_AminusC_fourier_mag_test, 
                                          Xone_AminusG_fourier_mag_test,
                                          Xone_AminusT_fourier_mag_test,
                                          Xone_CminusG_fourier_mag_test,
                                          Xone_CminusT_fourier_mag_test,
                                          Xone_GminusT_fourier_mag_test)
saveRDS(Xone_nucminusnuc_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_mag_test.rds")
Xone_nucminusnuc_fourier_phase_test = cbind(Xone_AminusC_fourier_phase_test, 
                                            Xone_AminusG_fourier_phase_test,
                                            Xone_AminusT_fourier_phase_test,
                                            Xone_CminusG_fourier_phase_test,
                                            Xone_CminusT_fourier_phase_test,
                                            Xone_GminusT_fourier_phase_test)
saveRDS(Xone_nucminusnuc_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucminusnuc_fourier_phase_test.rds")


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

# Xone_AorC_fourier = t(apply(Xone_AorC, 1, fft))[,1:26]
# Xone_AorG_fourier = t(apply(Xone_AorG, 1, fft))[,1:26]
# Xone_AorT_fourier = t(apply(Xone_AorT, 1, fft))[,1:26]
# Xone_CorG_fourier = t(apply(Xone_CorG, 1, fft))[,1:26]
# Xone_CorT_fourier = t(apply(Xone_CorT, 1, fft))[,1:26]
# Xone_GorT_fourier = t(apply(Xone_GorT, 1, fft))[,1:26]

Xone_AorC_fourier_test = t(apply(Xone_AorC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_AorG_fourier_test = t(apply(Xone_AorG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_AorT_fourier_test = t(apply(Xone_AorT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_CorG_fourier_test = t(apply(Xone_CorG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_CorT_fourier_test = t(apply(Xone_CorT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xone_GorT_fourier_test = t(apply(Xone_GorT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]

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

colnames(Xone_AorC_fourier_mag_test) = paste0("fourier_mag_AorC_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AorG_fourier_mag_test) = paste0("fourier_mag_AorG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AorT_fourier_mag_test) = paste0("fourier_mag_AorT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CorG_fourier_mag_test) = paste0("fourier_mag_CorG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CorT_fourier_mag_test) = paste0("fourier_mag_CorT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_GorT_fourier_mag_test) = paste0("fourier_mag_GorT_", c(1,2,3,4,5,6,11,21))

colnames(Xone_AorC_fourier_phase_test) = paste0("fourier_phase_AorC_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AorG_fourier_phase_test) = paste0("fourier_phase_AorG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_AorT_fourier_phase_test) = paste0("fourier_phase_AorT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CorG_fourier_phase_test) = paste0("fourier_phase_CorG_", c(1,2,3,4,5,6,11,21))
colnames(Xone_CorT_fourier_phase_test) = paste0("fourier_phase_CorT_", c(1,2,3,4,5,6,11,21))
colnames(Xone_GorT_fourier_phase_test) = paste0("fourier_phase_GorT_", c(1,2,3,4,5,6,11,21))

Xone_nucornuc_fourier_mag_test = cbind(Xone_AorC_fourier_mag_test, 
                                       Xone_AorG_fourier_mag_test,
                                       Xone_AorT_fourier_mag_test,
                                       Xone_CorG_fourier_mag_test,
                                       Xone_CorT_fourier_mag_test,
                                       Xone_GorT_fourier_mag_test)
saveRDS(Xone_nucornuc_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucornuc_fourier_mag_test.rds")

Xone_nucornuc_fourier_phase_test = cbind(Xone_AorC_fourier_phase_test, 
                                         Xone_AorG_fourier_phase_test,
                                         Xone_AorT_fourier_phase_test,
                                         Xone_CorG_fourier_phase_test,
                                         Xone_CorT_fourier_phase_test,
                                         Xone_GorT_fourier_phase_test)
saveRDS(Xone_nucornuc_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xone_nucornuc_fourier_phase_test.rds")



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

Xtwo_AA_fourier_test = t(apply(Xtwo_AA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AC_fourier_test = t(apply(Xtwo_AC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AG_fourier_test = t(apply(Xtwo_AG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AT_fourier_test = t(apply(Xtwo_AT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CA_fourier_test = t(apply(Xtwo_CA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CC_fourier_test = t(apply(Xtwo_CC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CG_fourier_test = t(apply(Xtwo_CG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CT_fourier_test = t(apply(Xtwo_CT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GA_fourier_test = t(apply(Xtwo_GA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GC_fourier_test = t(apply(Xtwo_GC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GG_fourier_test = t(apply(Xtwo_GG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GT_fourier_test = t(apply(Xtwo_GT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TA_fourier_test = t(apply(Xtwo_TA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TC_fourier_test = t(apply(Xtwo_TC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TG_fourier_test = t(apply(Xtwo_TG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TT_fourier_test = t(apply(Xtwo_TT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]

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

colnames(Xtwo_AA_fourier_mag_test) = paste0("fourier_mag_AA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AC_fourier_mag_test) = paste0("fourier_mag_AC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AG_fourier_mag_test) = paste0("fourier_mag_AG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AT_fourier_mag_test) = paste0("fourier_mag_AT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CA_fourier_mag_test) = paste0("fourier_mag_CA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CC_fourier_mag_test) = paste0("fourier_mag_CC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CG_fourier_mag_test) = paste0("fourier_mag_CG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CT_fourier_mag_test) = paste0("fourier_mag_CT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GA_fourier_mag_test) = paste0("fourier_mag_GA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GC_fourier_mag_test) = paste0("fourier_mag_GC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GG_fourier_mag_test) = paste0("fourier_mag_GG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GT_fourier_mag_test) = paste0("fourier_mag_GT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TA_fourier_mag_test) = paste0("fourier_mag_TA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TC_fourier_mag_test) = paste0("fourier_mag_TC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TG_fourier_mag_test) = paste0("fourier_mag_TG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TT_fourier_mag_test) = paste0("fourier_mag_TT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AA_fourier_phase_test) = paste0("fourier_phase_AA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AC_fourier_phase_test) = paste0("fourier_phase_AC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AG_fourier_phase_test) = paste0("fourier_phase_AG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AT_fourier_phase_test) = paste0("fourier_phase_AT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CA_fourier_phase_test) = paste0("fourier_phase_CA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CC_fourier_phase_test) = paste0("fourier_phase_CC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CG_fourier_phase_test) = paste0("fourier_phase_CG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CT_fourier_phase_test) = paste0("fourier_phase_CT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GA_fourier_phase_test) = paste0("fourier_phase_GA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GC_fourier_phase_test) = paste0("fourier_phase_GC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GG_fourier_phase_test) = paste0("fourier_phase_GG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GT_fourier_phase_test) = paste0("fourier_phase_GT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TA_fourier_phase_test) = paste0("fourier_phase_TA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TC_fourier_phase_test) = paste0("fourier_phase_TC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TG_fourier_phase_test) = paste0("fourier_phase_TG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TT_fourier_phase_test) = paste0("fourier_phase_TT_", c(1,2,3,4,5,6,11,21))

Xtwo_di_fourier_mag_test = cbind(Xtwo_AA_fourier_mag_test,
                                 Xtwo_AC_fourier_mag_test,
                                 Xtwo_AG_fourier_mag_test,
                                 Xtwo_AT_fourier_mag_test,
                                 Xtwo_CA_fourier_mag_test,
                                 Xtwo_CC_fourier_mag_test,
                                 Xtwo_CG_fourier_mag_test,
                                 Xtwo_CT_fourier_mag_test,
                                 Xtwo_GA_fourier_mag_test,
                                 Xtwo_GC_fourier_mag_test,
                                 Xtwo_GG_fourier_mag_test,
                                 Xtwo_GT_fourier_mag_test,
                                 Xtwo_TA_fourier_mag_test,
                                 Xtwo_TC_fourier_mag_test,
                                 Xtwo_TG_fourier_mag_test,
                                 Xtwo_TT_fourier_mag_test)
saveRDS(Xtwo_di_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_di_fourier_mag_test.rds")

Xtwo_di_fourier_phase_test = cbind(Xtwo_AA_fourier_phase_test,
                                   Xtwo_AC_fourier_phase_test,
                                   Xtwo_AG_fourier_phase_test,
                                   Xtwo_AT_fourier_phase_test,
                                   Xtwo_CA_fourier_phase_test,
                                   Xtwo_CC_fourier_phase_test,
                                   Xtwo_CG_fourier_phase_test,
                                   Xtwo_CT_fourier_phase_test,
                                   Xtwo_GA_fourier_phase_test,
                                   Xtwo_GC_fourier_phase_test,
                                   Xtwo_GG_fourier_phase_test,
                                   Xtwo_GT_fourier_phase_test,
                                   Xtwo_TA_fourier_phase_test,
                                   Xtwo_TC_fourier_phase_test,
                                   Xtwo_TG_fourier_phase_test,
                                   Xtwo_TT_fourier_phase_test)
saveRDS(Xtwo_di_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_di_fourier_phase_test.rds")


# Difference between single dinucleotides:

Xtwo_AAminusAC_fourier_test = Xtwo_AA_fourier_test - Xtwo_AC_fourier_test
Xtwo_AAminusAG_fourier_test = Xtwo_AA_fourier_test - Xtwo_AG_fourier_test
Xtwo_AAminusAT_fourier_test = Xtwo_AA_fourier_test - Xtwo_AT_fourier_test
Xtwo_AAminusCA_fourier_test = Xtwo_AA_fourier_test - Xtwo_CA_fourier_test
Xtwo_AAminusCC_fourier_test = Xtwo_AA_fourier_test - Xtwo_CC_fourier_test
Xtwo_AAminusCG_fourier_test = Xtwo_AA_fourier_test - Xtwo_CG_fourier_test
Xtwo_AAminusCT_fourier_test = Xtwo_AA_fourier_test - Xtwo_CT_fourier_test
Xtwo_AAminusGA_fourier_test = Xtwo_AA_fourier_test - Xtwo_GA_fourier_test
Xtwo_AAminusGC_fourier_test = Xtwo_AA_fourier_test - Xtwo_GC_fourier_test
Xtwo_AAminusGG_fourier_test = Xtwo_AA_fourier_test - Xtwo_GG_fourier_test
Xtwo_AAminusGT_fourier_test = Xtwo_AA_fourier_test - Xtwo_GT_fourier_test
Xtwo_AAminusTA_fourier_test = Xtwo_AA_fourier_test - Xtwo_TA_fourier_test
Xtwo_AAminusTC_fourier_test = Xtwo_AA_fourier_test - Xtwo_TC_fourier_test
Xtwo_AAminusTG_fourier_test = Xtwo_AA_fourier_test - Xtwo_TG_fourier_test
Xtwo_AAminusTT_fourier_test = Xtwo_AA_fourier_test - Xtwo_TT_fourier_test

Xtwo_ACminusAG_fourier_test = Xtwo_AC_fourier_test - Xtwo_AG_fourier_test
Xtwo_ACminusAT_fourier_test = Xtwo_AC_fourier_test - Xtwo_AT_fourier_test
Xtwo_ACminusCA_fourier_test = Xtwo_AC_fourier_test - Xtwo_CA_fourier_test
Xtwo_ACminusCC_fourier_test = Xtwo_AC_fourier_test - Xtwo_CC_fourier_test
Xtwo_ACminusCG_fourier_test = Xtwo_AC_fourier_test - Xtwo_CG_fourier_test
Xtwo_ACminusCT_fourier_test = Xtwo_AC_fourier_test - Xtwo_CT_fourier_test
Xtwo_ACminusGA_fourier_test = Xtwo_AC_fourier_test - Xtwo_GA_fourier_test
Xtwo_ACminusGC_fourier_test = Xtwo_AC_fourier_test - Xtwo_GC_fourier_test
Xtwo_ACminusGG_fourier_test = Xtwo_AC_fourier_test - Xtwo_GG_fourier_test
Xtwo_ACminusGT_fourier_test = Xtwo_AC_fourier_test - Xtwo_GT_fourier_test
Xtwo_ACminusTA_fourier_test = Xtwo_AC_fourier_test - Xtwo_TA_fourier_test
Xtwo_ACminusTC_fourier_test = Xtwo_AC_fourier_test - Xtwo_TC_fourier_test
Xtwo_ACminusTG_fourier_test = Xtwo_AC_fourier_test - Xtwo_TG_fourier_test
Xtwo_ACminusTT_fourier_test = Xtwo_AC_fourier_test - Xtwo_TT_fourier_test

Xtwo_AGminusAT_fourier_test = Xtwo_AG_fourier_test - Xtwo_AT_fourier_test
Xtwo_AGminusCA_fourier_test = Xtwo_AG_fourier_test - Xtwo_CA_fourier_test
Xtwo_AGminusCC_fourier_test = Xtwo_AG_fourier_test - Xtwo_CC_fourier_test
Xtwo_AGminusCG_fourier_test = Xtwo_AG_fourier_test - Xtwo_CG_fourier_test
Xtwo_AGminusCT_fourier_test = Xtwo_AG_fourier_test - Xtwo_CT_fourier_test
Xtwo_AGminusGA_fourier_test = Xtwo_AG_fourier_test - Xtwo_GA_fourier_test
Xtwo_AGminusGC_fourier_test = Xtwo_AG_fourier_test - Xtwo_GC_fourier_test
Xtwo_AGminusGG_fourier_test = Xtwo_AG_fourier_test - Xtwo_GG_fourier_test
Xtwo_AGminusGT_fourier_test = Xtwo_AG_fourier_test - Xtwo_GT_fourier_test
Xtwo_AGminusTA_fourier_test = Xtwo_AG_fourier_test - Xtwo_TA_fourier_test
Xtwo_AGminusTC_fourier_test = Xtwo_AG_fourier_test - Xtwo_TC_fourier_test
Xtwo_AGminusTG_fourier_test = Xtwo_AG_fourier_test - Xtwo_TG_fourier_test
Xtwo_AGminusTT_fourier_test = Xtwo_AG_fourier_test - Xtwo_TT_fourier_test

Xtwo_ATminusCA_fourier_test = Xtwo_AT_fourier_test - Xtwo_CA_fourier_test
Xtwo_ATminusCC_fourier_test = Xtwo_AT_fourier_test - Xtwo_CC_fourier_test
Xtwo_ATminusCG_fourier_test = Xtwo_AT_fourier_test - Xtwo_CG_fourier_test
Xtwo_ATminusCT_fourier_test = Xtwo_AT_fourier_test - Xtwo_CT_fourier_test
Xtwo_ATminusGA_fourier_test = Xtwo_AT_fourier_test - Xtwo_GA_fourier_test
Xtwo_ATminusGC_fourier_test = Xtwo_AT_fourier_test - Xtwo_GC_fourier_test
Xtwo_ATminusGG_fourier_test = Xtwo_AT_fourier_test - Xtwo_GG_fourier_test
Xtwo_ATminusGT_fourier_test = Xtwo_AT_fourier_test - Xtwo_GT_fourier_test
Xtwo_ATminusTA_fourier_test = Xtwo_AT_fourier_test - Xtwo_TA_fourier_test
Xtwo_ATminusTC_fourier_test = Xtwo_AT_fourier_test - Xtwo_TC_fourier_test
Xtwo_ATminusTG_fourier_test = Xtwo_AT_fourier_test - Xtwo_TG_fourier_test
Xtwo_ATminusTT_fourier_test = Xtwo_AT_fourier_test - Xtwo_TT_fourier_test

Xtwo_CAminusCC_fourier_test = Xtwo_CA_fourier_test - Xtwo_CC_fourier_test
Xtwo_CAminusCG_fourier_test = Xtwo_CA_fourier_test - Xtwo_CG_fourier_test
Xtwo_CAminusCT_fourier_test = Xtwo_CA_fourier_test - Xtwo_CT_fourier_test
Xtwo_CAminusGA_fourier_test = Xtwo_CA_fourier_test - Xtwo_GA_fourier_test
Xtwo_CAminusGC_fourier_test = Xtwo_CA_fourier_test - Xtwo_GC_fourier_test
Xtwo_CAminusGG_fourier_test = Xtwo_CA_fourier_test - Xtwo_GG_fourier_test
Xtwo_CAminusGT_fourier_test = Xtwo_CA_fourier_test - Xtwo_GT_fourier_test
Xtwo_CAminusTA_fourier_test = Xtwo_CA_fourier_test - Xtwo_TA_fourier_test
Xtwo_CAminusTC_fourier_test = Xtwo_CA_fourier_test - Xtwo_TC_fourier_test
Xtwo_CAminusTG_fourier_test = Xtwo_CA_fourier_test - Xtwo_TG_fourier_test
Xtwo_CAminusTT_fourier_test = Xtwo_CA_fourier_test - Xtwo_TT_fourier_test

Xtwo_CCminusCG_fourier_test = Xtwo_CC_fourier_test - Xtwo_CG_fourier_test
Xtwo_CCminusCT_fourier_test = Xtwo_CC_fourier_test - Xtwo_CT_fourier_test
Xtwo_CCminusGA_fourier_test = Xtwo_CC_fourier_test - Xtwo_GA_fourier_test
Xtwo_CCminusGC_fourier_test = Xtwo_CC_fourier_test - Xtwo_GC_fourier_test
Xtwo_CCminusGG_fourier_test = Xtwo_CC_fourier_test - Xtwo_GG_fourier_test
Xtwo_CCminusGT_fourier_test = Xtwo_CC_fourier_test - Xtwo_GT_fourier_test
Xtwo_CCminusTA_fourier_test = Xtwo_CC_fourier_test - Xtwo_TA_fourier_test
Xtwo_CCminusTC_fourier_test = Xtwo_CC_fourier_test - Xtwo_TC_fourier_test
Xtwo_CCminusTG_fourier_test = Xtwo_CC_fourier_test - Xtwo_TG_fourier_test
Xtwo_CCminusTT_fourier_test = Xtwo_CC_fourier_test - Xtwo_TT_fourier_test

Xtwo_CGminusCT_fourier_test = Xtwo_CG_fourier_test - Xtwo_CT_fourier_test
Xtwo_CGminusGA_fourier_test = Xtwo_CG_fourier_test - Xtwo_GA_fourier_test
Xtwo_CGminusGC_fourier_test = Xtwo_CG_fourier_test - Xtwo_GC_fourier_test
Xtwo_CGminusGG_fourier_test = Xtwo_CG_fourier_test - Xtwo_GG_fourier_test
Xtwo_CGminusGT_fourier_test = Xtwo_CG_fourier_test - Xtwo_GT_fourier_test
Xtwo_CGminusTA_fourier_test = Xtwo_CG_fourier_test - Xtwo_TA_fourier_test
Xtwo_CGminusTC_fourier_test = Xtwo_CG_fourier_test - Xtwo_TC_fourier_test
Xtwo_CGminusTG_fourier_test = Xtwo_CG_fourier_test - Xtwo_TG_fourier_test
Xtwo_CGminusTT_fourier_test = Xtwo_CG_fourier_test - Xtwo_TT_fourier_test

Xtwo_CTminusGA_fourier_test = Xtwo_CT_fourier_test - Xtwo_GA_fourier_test
Xtwo_CTminusGC_fourier_test = Xtwo_CT_fourier_test - Xtwo_GC_fourier_test
Xtwo_CTminusGG_fourier_test = Xtwo_CT_fourier_test - Xtwo_GG_fourier_test
Xtwo_CTminusGT_fourier_test = Xtwo_CT_fourier_test - Xtwo_GT_fourier_test
Xtwo_CTminusTA_fourier_test = Xtwo_CT_fourier_test - Xtwo_TA_fourier_test
Xtwo_CTminusTC_fourier_test = Xtwo_CT_fourier_test - Xtwo_TC_fourier_test
Xtwo_CTminusTG_fourier_test = Xtwo_CT_fourier_test - Xtwo_TG_fourier_test
Xtwo_CTminusTT_fourier_test = Xtwo_CT_fourier_test - Xtwo_TT_fourier_test

Xtwo_GAminusGC_fourier_test = Xtwo_GA_fourier_test - Xtwo_GC_fourier_test
Xtwo_GAminusGG_fourier_test = Xtwo_GA_fourier_test - Xtwo_GG_fourier_test
Xtwo_GAminusGT_fourier_test = Xtwo_GA_fourier_test - Xtwo_GT_fourier_test
Xtwo_GAminusTA_fourier_test = Xtwo_GA_fourier_test - Xtwo_TA_fourier_test
Xtwo_GAminusTC_fourier_test = Xtwo_GA_fourier_test - Xtwo_TC_fourier_test
Xtwo_GAminusTG_fourier_test = Xtwo_GA_fourier_test - Xtwo_TG_fourier_test
Xtwo_GAminusTT_fourier_test = Xtwo_GA_fourier_test - Xtwo_TT_fourier_test

Xtwo_GCminusGG_fourier_test = Xtwo_GC_fourier_test - Xtwo_GG_fourier_test
Xtwo_GCminusGT_fourier_test = Xtwo_GC_fourier_test - Xtwo_GT_fourier_test
Xtwo_GCminusTA_fourier_test = Xtwo_GC_fourier_test - Xtwo_TA_fourier_test
Xtwo_GCminusTC_fourier_test = Xtwo_GC_fourier_test - Xtwo_TC_fourier_test
Xtwo_GCminusTG_fourier_test = Xtwo_GC_fourier_test - Xtwo_TG_fourier_test
Xtwo_GCminusTT_fourier_test = Xtwo_GC_fourier_test - Xtwo_TT_fourier_test

Xtwo_GGminusGT_fourier_test = Xtwo_GG_fourier_test - Xtwo_GT_fourier_test
Xtwo_GGminusTA_fourier_test = Xtwo_GG_fourier_test - Xtwo_TA_fourier_test
Xtwo_GGminusTC_fourier_test = Xtwo_GG_fourier_test - Xtwo_TC_fourier_test
Xtwo_GGminusTG_fourier_test = Xtwo_GG_fourier_test - Xtwo_TG_fourier_test
Xtwo_GGminusTT_fourier_test = Xtwo_GG_fourier_test - Xtwo_TT_fourier_test

Xtwo_GTminusTA_fourier_test = Xtwo_GT_fourier_test - Xtwo_TA_fourier_test
Xtwo_GTminusTC_fourier_test = Xtwo_GT_fourier_test - Xtwo_TC_fourier_test
Xtwo_GTminusTG_fourier_test = Xtwo_GT_fourier_test - Xtwo_TG_fourier_test
Xtwo_GTminusTT_fourier_test = Xtwo_GT_fourier_test - Xtwo_TT_fourier_test

Xtwo_TAminusTC_fourier_test = Xtwo_TA_fourier_test - Xtwo_TC_fourier_test
Xtwo_TAminusTG_fourier_test = Xtwo_TA_fourier_test - Xtwo_TG_fourier_test
Xtwo_TAminusTT_fourier_test = Xtwo_TA_fourier_test - Xtwo_TT_fourier_test

Xtwo_TCminusTG_fourier_test = Xtwo_TC_fourier_test - Xtwo_TG_fourier_test
Xtwo_TCminusTT_fourier_test = Xtwo_TC_fourier_test - Xtwo_TT_fourier_test

Xtwo_TGminusTT_fourier_test = Xtwo_TG_fourier_test - Xtwo_TT_fourier_test

Xtwo_AAminusAC_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusAC_fourier_test))
Xtwo_AAminusAC_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusAC_fourier_test))
Xtwo_AAminusAG_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusAG_fourier_test))
Xtwo_AAminusAT_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusAT_fourier_test))
Xtwo_AAminusCA_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusCA_fourier_test))
Xtwo_AAminusCC_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusCC_fourier_test))
Xtwo_AAminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusCG_fourier_test))
Xtwo_AAminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusCT_fourier_test))
Xtwo_AAminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusGA_fourier_test))
Xtwo_AAminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusGC_fourier_test))
Xtwo_AAminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusGG_fourier_test))
Xtwo_AAminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusGT_fourier_test))
Xtwo_AAminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusTA_fourier_test))
Xtwo_AAminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusTC_fourier_test))
Xtwo_AAminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusTG_fourier_test))
Xtwo_AAminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_AAminusTT_fourier_test))

Xtwo_ACminusAG_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusAG_fourier_test))
Xtwo_ACminusAT_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusAT_fourier_test))
Xtwo_ACminusCA_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusCA_fourier_test))
Xtwo_ACminusCC_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusCC_fourier_test))
Xtwo_ACminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusCG_fourier_test))
Xtwo_ACminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusCT_fourier_test))
Xtwo_ACminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusGA_fourier_test))
Xtwo_ACminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusGC_fourier_test))
Xtwo_ACminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusGG_fourier_test))
Xtwo_ACminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusGT_fourier_test))
Xtwo_ACminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusTA_fourier_test))
Xtwo_ACminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusTC_fourier_test))
Xtwo_ACminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusTG_fourier_test))
Xtwo_ACminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_ACminusTT_fourier_test))

Xtwo_AGminusAT_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusAT_fourier_test))
Xtwo_AGminusCA_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusCA_fourier_test))
Xtwo_AGminusCC_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusCC_fourier_test))
Xtwo_AGminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusCG_fourier_test))
Xtwo_AGminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusCT_fourier_test))
Xtwo_AGminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusGA_fourier_test))
Xtwo_AGminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusGC_fourier_test))
Xtwo_AGminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusGG_fourier_test))
Xtwo_AGminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusGT_fourier_test))
Xtwo_AGminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusTA_fourier_test))
Xtwo_AGminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusTC_fourier_test))
Xtwo_AGminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusTG_fourier_test))
Xtwo_AGminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_AGminusTT_fourier_test))

Xtwo_ATminusCA_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusCA_fourier_test))
Xtwo_ATminusCC_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusCC_fourier_test))
Xtwo_ATminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusCG_fourier_test))
Xtwo_ATminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusCT_fourier_test))
Xtwo_ATminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusGA_fourier_test))
Xtwo_ATminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusGC_fourier_test))
Xtwo_ATminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusGG_fourier_test))
Xtwo_ATminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusGT_fourier_test))
Xtwo_ATminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusTA_fourier_test))
Xtwo_ATminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusTC_fourier_test))
Xtwo_ATminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusTG_fourier_test))
Xtwo_ATminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_ATminusTT_fourier_test))

Xtwo_CAminusCC_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusCC_fourier_test))
Xtwo_CAminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusCG_fourier_test))
Xtwo_CAminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusCT_fourier_test))
Xtwo_CAminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusGA_fourier_test))
Xtwo_CAminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusGC_fourier_test))
Xtwo_CAminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusGG_fourier_test))
Xtwo_CAminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusGT_fourier_test))
Xtwo_CAminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusTA_fourier_test))
Xtwo_CAminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusTC_fourier_test))
Xtwo_CAminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusTG_fourier_test))
Xtwo_CAminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_CAminusTT_fourier_test))

Xtwo_CCminusCG_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusCG_fourier_test))
Xtwo_CCminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusCT_fourier_test))
Xtwo_CCminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusGA_fourier_test))
Xtwo_CCminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusGC_fourier_test))
Xtwo_CCminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusGG_fourier_test))
Xtwo_CCminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusGT_fourier_test))
Xtwo_CCminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusTA_fourier_test))
Xtwo_CCminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusTC_fourier_test))
Xtwo_CCminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusTG_fourier_test))
Xtwo_CCminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_CCminusTT_fourier_test))

Xtwo_CGminusCT_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusCT_fourier_test))
Xtwo_CGminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusGA_fourier_test))
Xtwo_CGminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusGC_fourier_test))
Xtwo_CGminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusGG_fourier_test))
Xtwo_CGminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusGT_fourier_test))
Xtwo_CGminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusTA_fourier_test))
Xtwo_CGminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusTC_fourier_test))
Xtwo_CGminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusTG_fourier_test))
Xtwo_CGminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_CGminusTT_fourier_test))

Xtwo_CTminusGA_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusGA_fourier_test))
Xtwo_CTminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusGC_fourier_test))
Xtwo_CTminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusGG_fourier_test))
Xtwo_CTminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusGT_fourier_test))
Xtwo_CTminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusTA_fourier_test))
Xtwo_CTminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusTC_fourier_test))
Xtwo_CTminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusTG_fourier_test))
Xtwo_CTminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_CTminusTT_fourier_test))

Xtwo_GAminusGC_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusGC_fourier_test))
Xtwo_GAminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusGG_fourier_test))
Xtwo_GAminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusGT_fourier_test))
Xtwo_GAminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusTA_fourier_test))
Xtwo_GAminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusTC_fourier_test))
Xtwo_GAminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusTG_fourier_test))
Xtwo_GAminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_GAminusTT_fourier_test))

Xtwo_GCminusGG_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusGG_fourier_test))
Xtwo_GCminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusGT_fourier_test))
Xtwo_GCminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusTA_fourier_test))
Xtwo_GCminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusTC_fourier_test))
Xtwo_GCminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusTG_fourier_test))
Xtwo_GCminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_GCminusTT_fourier_test))

Xtwo_GGminusGT_fourier_mag_test = as.matrix(Mod(Xtwo_GGminusGT_fourier_test))
Xtwo_GGminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_GGminusTA_fourier_test))
Xtwo_GGminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_GGminusTC_fourier_test))
Xtwo_GGminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_GGminusTG_fourier_test))
Xtwo_GGminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_GGminusTT_fourier_test))

Xtwo_GTminusTA_fourier_mag_test = as.matrix(Mod(Xtwo_GTminusTA_fourier_test))
Xtwo_GTminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_GTminusTC_fourier_test))
Xtwo_GTminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_GTminusTG_fourier_test))
Xtwo_GTminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_GTminusTT_fourier_test))

Xtwo_TAminusTC_fourier_mag_test = as.matrix(Mod(Xtwo_TAminusTC_fourier_test))
Xtwo_TAminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_TAminusTG_fourier_test))
Xtwo_TAminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_TAminusTT_fourier_test))

Xtwo_TCminusTG_fourier_mag_test = as.matrix(Mod(Xtwo_TCminusTG_fourier_test))
Xtwo_TCminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_TCminusTT_fourier_test))

Xtwo_TGminusTT_fourier_mag_test = as.matrix(Mod(Xtwo_TGminusTT_fourier_test))

Xtwo_AAminusAC_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusAC_fourier_test))
Xtwo_AAminusAC_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusAC_fourier_test))
Xtwo_AAminusAG_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusAG_fourier_test))
Xtwo_AAminusAT_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusAT_fourier_test))
Xtwo_AAminusCA_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusCA_fourier_test))
Xtwo_AAminusCC_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusCC_fourier_test))
Xtwo_AAminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusCG_fourier_test))
Xtwo_AAminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusCT_fourier_test))
Xtwo_AAminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusGA_fourier_test))
Xtwo_AAminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusGC_fourier_test))
Xtwo_AAminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusGG_fourier_test))
Xtwo_AAminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusGT_fourier_test))
Xtwo_AAminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusTA_fourier_test))
Xtwo_AAminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusTC_fourier_test))
Xtwo_AAminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusTG_fourier_test))
Xtwo_AAminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_AAminusTT_fourier_test))

Xtwo_ACminusAG_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusAG_fourier_test))
Xtwo_ACminusAT_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusAT_fourier_test))
Xtwo_ACminusCA_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusCA_fourier_test))
Xtwo_ACminusCC_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusCC_fourier_test))
Xtwo_ACminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusCG_fourier_test))
Xtwo_ACminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusCT_fourier_test))
Xtwo_ACminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusGA_fourier_test))
Xtwo_ACminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusGC_fourier_test))
Xtwo_ACminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusGG_fourier_test))
Xtwo_ACminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusGT_fourier_test))
Xtwo_ACminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusTA_fourier_test))
Xtwo_ACminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusTC_fourier_test))
Xtwo_ACminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusTG_fourier_test))
Xtwo_ACminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_ACminusTT_fourier_test))

Xtwo_AGminusAT_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusAT_fourier_test))
Xtwo_AGminusCA_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusCA_fourier_test))
Xtwo_AGminusCC_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusCC_fourier_test))
Xtwo_AGminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusCG_fourier_test))
Xtwo_AGminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusCT_fourier_test))
Xtwo_AGminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusGA_fourier_test))
Xtwo_AGminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusGC_fourier_test))
Xtwo_AGminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusGG_fourier_test))
Xtwo_AGminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusGT_fourier_test))
Xtwo_AGminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusTA_fourier_test))
Xtwo_AGminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusTC_fourier_test))
Xtwo_AGminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusTG_fourier_test))
Xtwo_AGminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_AGminusTT_fourier_test))

Xtwo_ATminusCA_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusCA_fourier_test))
Xtwo_ATminusCC_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusCC_fourier_test))
Xtwo_ATminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusCG_fourier_test))
Xtwo_ATminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusCT_fourier_test))
Xtwo_ATminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusGA_fourier_test))
Xtwo_ATminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusGC_fourier_test))
Xtwo_ATminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusGG_fourier_test))
Xtwo_ATminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusGT_fourier_test))
Xtwo_ATminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusTA_fourier_test))
Xtwo_ATminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusTC_fourier_test))
Xtwo_ATminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusTG_fourier_test))
Xtwo_ATminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_ATminusTT_fourier_test))

Xtwo_CAminusCC_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusCC_fourier_test))
Xtwo_CAminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusCG_fourier_test))
Xtwo_CAminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusCT_fourier_test))
Xtwo_CAminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusGA_fourier_test))
Xtwo_CAminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusGC_fourier_test))
Xtwo_CAminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusGG_fourier_test))
Xtwo_CAminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusGT_fourier_test))
Xtwo_CAminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusTA_fourier_test))
Xtwo_CAminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusTC_fourier_test))
Xtwo_CAminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusTG_fourier_test))
Xtwo_CAminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_CAminusTT_fourier_test))

Xtwo_CCminusCG_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusCG_fourier_test))
Xtwo_CCminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusCT_fourier_test))
Xtwo_CCminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusGA_fourier_test))
Xtwo_CCminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusGC_fourier_test))
Xtwo_CCminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusGG_fourier_test))
Xtwo_CCminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusGT_fourier_test))
Xtwo_CCminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusTA_fourier_test))
Xtwo_CCminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusTC_fourier_test))
Xtwo_CCminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusTG_fourier_test))
Xtwo_CCminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_CCminusTT_fourier_test))

Xtwo_CGminusCT_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusCT_fourier_test))
Xtwo_CGminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusGA_fourier_test))
Xtwo_CGminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusGC_fourier_test))
Xtwo_CGminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusGG_fourier_test))
Xtwo_CGminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusGT_fourier_test))
Xtwo_CGminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusTA_fourier_test))
Xtwo_CGminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusTC_fourier_test))
Xtwo_CGminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusTG_fourier_test))
Xtwo_CGminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_CGminusTT_fourier_test))

Xtwo_CTminusGA_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusGA_fourier_test))
Xtwo_CTminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusGC_fourier_test))
Xtwo_CTminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusGG_fourier_test))
Xtwo_CTminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusGT_fourier_test))
Xtwo_CTminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusTA_fourier_test))
Xtwo_CTminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusTC_fourier_test))
Xtwo_CTminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusTG_fourier_test))
Xtwo_CTminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_CTminusTT_fourier_test))

Xtwo_GAminusGC_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusGC_fourier_test))
Xtwo_GAminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusGG_fourier_test))
Xtwo_GAminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusGT_fourier_test))
Xtwo_GAminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusTA_fourier_test))
Xtwo_GAminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusTC_fourier_test))
Xtwo_GAminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusTG_fourier_test))
Xtwo_GAminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_GAminusTT_fourier_test))

Xtwo_GCminusGG_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusGG_fourier_test))
Xtwo_GCminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusGT_fourier_test))
Xtwo_GCminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusTA_fourier_test))
Xtwo_GCminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusTC_fourier_test))
Xtwo_GCminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusTG_fourier_test))
Xtwo_GCminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_GCminusTT_fourier_test))

Xtwo_GGminusGT_fourier_phase_test = as.matrix(Arg(Xtwo_GGminusGT_fourier_test))
Xtwo_GGminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_GGminusTA_fourier_test))
Xtwo_GGminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_GGminusTC_fourier_test))
Xtwo_GGminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_GGminusTG_fourier_test))
Xtwo_GGminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_GGminusTT_fourier_test))

Xtwo_GTminusTA_fourier_phase_test = as.matrix(Arg(Xtwo_GTminusTA_fourier_test))
Xtwo_GTminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_GTminusTC_fourier_test))
Xtwo_GTminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_GTminusTG_fourier_test))
Xtwo_GTminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_GTminusTT_fourier_test))

Xtwo_TAminusTC_fourier_phase_test = as.matrix(Arg(Xtwo_TAminusTC_fourier_test))
Xtwo_TAminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_TAminusTG_fourier_test))
Xtwo_TAminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_TAminusTT_fourier_test))

Xtwo_TCminusTG_fourier_phase_test = as.matrix(Arg(Xtwo_TCminusTG_fourier_test))
Xtwo_TCminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_TCminusTT_fourier_test))

Xtwo_TGminusTT_fourier_phase_test = as.matrix(Arg(Xtwo_TGminusTT_fourier_test))

colnames(Xtwo_AAminusAC_fourier_mag_test) = paste0("fourier_mag_AAminusAC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAG_fourier_mag_test) = paste0("fourier_mag_AAminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAT_fourier_mag_test) = paste0("fourier_mag_AAminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCA_fourier_mag_test) = paste0("fourier_mag_AAminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCC_fourier_mag_test) = paste0("fourier_mag_AAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCG_fourier_mag_test) = paste0("fourier_mag_AAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCT_fourier_mag_test) = paste0("fourier_mag_AAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGA_fourier_mag_test) = paste0("fourier_mag_AAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGC_fourier_mag_test) = paste0("fourier_mag_AAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGG_fourier_mag_test) = paste0("fourier_mag_AAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGT_fourier_mag_test) = paste0("fourier_mag_AAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTA_fourier_mag_test) = paste0("fourier_mag_AAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTC_fourier_mag_test) = paste0("fourier_mag_AAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTG_fourier_mag_test) = paste0("fourier_mag_AAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTT_fourier_mag_test) = paste0("fourier_mag_AAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ACminusAG_fourier_mag_test) = paste0("fourier_mag_ACminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusAT_fourier_mag_test) = paste0("fourier_mag_ACminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCA_fourier_mag_test) = paste0("fourier_mag_ACminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCC_fourier_mag_test) = paste0("fourier_mag_ACminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCG_fourier_mag_test) = paste0("fourier_mag_ACminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCT_fourier_mag_test) = paste0("fourier_mag_ACminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGA_fourier_mag_test) = paste0("fourier_mag_ACminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGC_fourier_mag_test) = paste0("fourier_mag_ACminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGG_fourier_mag_test) = paste0("fourier_mag_ACminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGT_fourier_mag_test) = paste0("fourier_mag_ACminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTA_fourier_mag_test) = paste0("fourier_mag_ACminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTC_fourier_mag_test) = paste0("fourier_mag_ACminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTG_fourier_mag_test) = paste0("fourier_mag_ACminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTT_fourier_mag_test) = paste0("fourier_mag_ACminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AGminusAT_fourier_mag_test) = paste0("fourier_mag_AGminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCA_fourier_mag_test) = paste0("fourier_mag_AGminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCC_fourier_mag_test) = paste0("fourier_mag_AGminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCG_fourier_mag_test) = paste0("fourier_mag_AGminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCT_fourier_mag_test) = paste0("fourier_mag_AGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGA_fourier_mag_test) = paste0("fourier_mag_AGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGC_fourier_mag_test) = paste0("fourier_mag_AGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGG_fourier_mag_test) = paste0("fourier_mag_AGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGT_fourier_mag_test) = paste0("fourier_mag_AGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTA_fourier_mag_test) = paste0("fourier_mag_AGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTC_fourier_mag_test) = paste0("fourier_mag_AGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTG_fourier_mag_test) = paste0("fourier_mag_AGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTT_fourier_mag_test) = paste0("fourier_mag_AGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ATminusCA_fourier_mag_test) = paste0("fourier_mag_ATminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCC_fourier_mag_test) = paste0("fourier_mag_ATminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCG_fourier_mag_test) = paste0("fourier_mag_ATminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCT_fourier_mag_test) = paste0("fourier_mag_ATminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGA_fourier_mag_test) = paste0("fourier_mag_ATminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGC_fourier_mag_test) = paste0("fourier_mag_ATminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGG_fourier_mag_test) = paste0("fourier_mag_ATminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGT_fourier_mag_test) = paste0("fourier_mag_ATminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTA_fourier_mag_test) = paste0("fourier_mag_ATminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTC_fourier_mag_test) = paste0("fourier_mag_ATminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTG_fourier_mag_test) = paste0("fourier_mag_ATminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTT_fourier_mag_test) = paste0("fourier_mag_ATminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CAminusCC_fourier_mag_test) = paste0("fourier_mag_CAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCG_fourier_mag_test) = paste0("fourier_mag_CAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCT_fourier_mag_test) = paste0("fourier_mag_CAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGA_fourier_mag_test) = paste0("fourier_mag_CAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGC_fourier_mag_test) = paste0("fourier_mag_CAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGG_fourier_mag_test) = paste0("fourier_mag_CAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGT_fourier_mag_test) = paste0("fourier_mag_CAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTA_fourier_mag_test) = paste0("fourier_mag_CAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTC_fourier_mag_test) = paste0("fourier_mag_CAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTG_fourier_mag_test) = paste0("fourier_mag_CAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTT_fourier_mag_test) = paste0("fourier_mag_CAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CCminusCG_fourier_mag_test) = paste0("fourier_mag_CCminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusCT_fourier_mag_test) = paste0("fourier_mag_CCminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGA_fourier_mag_test) = paste0("fourier_mag_CCminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGC_fourier_mag_test) = paste0("fourier_mag_CCminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGG_fourier_mag_test) = paste0("fourier_mag_CCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGT_fourier_mag_test) = paste0("fourier_mag_CCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTA_fourier_mag_test) = paste0("fourier_mag_CCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTC_fourier_mag_test) = paste0("fourier_mag_CCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTG_fourier_mag_test) = paste0("fourier_mag_CCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTT_fourier_mag_test) = paste0("fourier_mag_CCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CGminusCT_fourier_mag_test) = paste0("fourier_mag_CGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGA_fourier_mag_test) = paste0("fourier_mag_CGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGC_fourier_mag_test) = paste0("fourier_mag_CGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGG_fourier_mag_test) = paste0("fourier_mag_CGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGT_fourier_mag_test) = paste0("fourier_mag_CGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTA_fourier_mag_test) = paste0("fourier_mag_CGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTC_fourier_mag_test) = paste0("fourier_mag_CGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTG_fourier_mag_test) = paste0("fourier_mag_CGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTT_fourier_mag_test) = paste0("fourier_mag_CGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CTminusGA_fourier_mag_test) = paste0("fourier_mag_CTminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGC_fourier_mag_test) = paste0("fourier_mag_CTminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGG_fourier_mag_test) = paste0("fourier_mag_CTminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGT_fourier_mag_test) = paste0("fourier_mag_CTminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTA_fourier_mag_test) = paste0("fourier_mag_CTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTC_fourier_mag_test) = paste0("fourier_mag_CTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTG_fourier_mag_test) = paste0("fourier_mag_CTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTT_fourier_mag_test) = paste0("fourier_mag_CTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GAminusGC_fourier_mag_test) = paste0("fourier_mag_GAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGG_fourier_mag_test) = paste0("fourier_mag_GAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGT_fourier_mag_test) = paste0("fourier_mag_GAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTA_fourier_mag_test) = paste0("fourier_mag_GAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTC_fourier_mag_test) = paste0("fourier_mag_GAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTG_fourier_mag_test) = paste0("fourier_mag_GAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTT_fourier_mag_test) = paste0("fourier_mag_GAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GCminusGG_fourier_mag_test) = paste0("fourier_mag_GCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusGT_fourier_mag_test) = paste0("fourier_mag_GCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTA_fourier_mag_test) = paste0("fourier_mag_GCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTC_fourier_mag_test) = paste0("fourier_mag_GCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTG_fourier_mag_test) = paste0("fourier_mag_GCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTT_fourier_mag_test) = paste0("fourier_mag_GCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GGminusGT_fourier_mag_test) = paste0("fourier_mag_GGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTA_fourier_mag_test) = paste0("fourier_mag_GGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTC_fourier_mag_test) = paste0("fourier_mag_GGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTG_fourier_mag_test) = paste0("fourier_mag_GGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTT_fourier_mag_test) = paste0("fourier_mag_GGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GTminusTA_fourier_mag_test) = paste0("fourier_mag_GTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTC_fourier_mag_test) = paste0("fourier_mag_GTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTG_fourier_mag_test) = paste0("fourier_mag_GTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTT_fourier_mag_test) = paste0("fourier_mag_GTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TAminusTC_fourier_mag_test) = paste0("fourier_mag_TAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTG_fourier_mag_test) = paste0("fourier_mag_TAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTT_fourier_mag_test) = paste0("fourier_mag_TAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TCminusTG_fourier_mag_test) = paste0("fourier_mag_TCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TCminusTT_fourier_mag_test) = paste0("fourier_mag_TCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TGminusTT_fourier_mag_test) = paste0("fourier_mag_TGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAminusAC_fourier_phase_test) = paste0("fourier_phase_AAminusAC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAG_fourier_phase_test) = paste0("fourier_phase_AAminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusAT_fourier_phase_test) = paste0("fourier_phase_AAminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCA_fourier_phase_test) = paste0("fourier_phase_AAminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCC_fourier_phase_test) = paste0("fourier_phase_AAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCG_fourier_phase_test) = paste0("fourier_phase_AAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusCT_fourier_phase_test) = paste0("fourier_phase_AAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGA_fourier_phase_test) = paste0("fourier_phase_AAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGC_fourier_phase_test) = paste0("fourier_phase_AAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGG_fourier_phase_test) = paste0("fourier_phase_AAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusGT_fourier_phase_test) = paste0("fourier_phase_AAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTA_fourier_phase_test) = paste0("fourier_phase_AAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTC_fourier_phase_test) = paste0("fourier_phase_AAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTG_fourier_phase_test) = paste0("fourier_phase_AAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAminusTT_fourier_phase_test) = paste0("fourier_phase_AAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ACminusAG_fourier_phase_test) = paste0("fourier_phase_ACminusAG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusAT_fourier_phase_test) = paste0("fourier_phase_ACminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCA_fourier_phase_test) = paste0("fourier_phase_ACminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCC_fourier_phase_test) = paste0("fourier_phase_ACminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCG_fourier_phase_test) = paste0("fourier_phase_ACminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusCT_fourier_phase_test) = paste0("fourier_phase_ACminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGA_fourier_phase_test) = paste0("fourier_phase_ACminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGC_fourier_phase_test) = paste0("fourier_phase_ACminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGG_fourier_phase_test) = paste0("fourier_phase_ACminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusGT_fourier_phase_test) = paste0("fourier_phase_ACminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTA_fourier_phase_test) = paste0("fourier_phase_ACminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTC_fourier_phase_test) = paste0("fourier_phase_ACminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTG_fourier_phase_test) = paste0("fourier_phase_ACminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ACminusTT_fourier_phase_test) = paste0("fourier_phase_ACminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AGminusAT_fourier_phase_test) = paste0("fourier_phase_AGminusAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCA_fourier_phase_test) = paste0("fourier_phase_AGminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCC_fourier_phase_test) = paste0("fourier_phase_AGminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCG_fourier_phase_test) = paste0("fourier_phase_AGminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusCT_fourier_phase_test) = paste0("fourier_phase_AGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGA_fourier_phase_test) = paste0("fourier_phase_AGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGC_fourier_phase_test) = paste0("fourier_phase_AGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGG_fourier_phase_test) = paste0("fourier_phase_AGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusGT_fourier_phase_test) = paste0("fourier_phase_AGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTA_fourier_phase_test) = paste0("fourier_phase_AGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTC_fourier_phase_test) = paste0("fourier_phase_AGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTG_fourier_phase_test) = paste0("fourier_phase_AGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AGminusTT_fourier_phase_test) = paste0("fourier_phase_AGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_ATminusCA_fourier_phase_test) = paste0("fourier_phase_ATminusCA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCC_fourier_phase_test) = paste0("fourier_phase_ATminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCG_fourier_phase_test) = paste0("fourier_phase_ATminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusCT_fourier_phase_test) = paste0("fourier_phase_ATminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGA_fourier_phase_test) = paste0("fourier_phase_ATminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGC_fourier_phase_test) = paste0("fourier_phase_ATminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGG_fourier_phase_test) = paste0("fourier_phase_ATminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusGT_fourier_phase_test) = paste0("fourier_phase_ATminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTA_fourier_phase_test) = paste0("fourier_phase_ATminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTC_fourier_phase_test) = paste0("fourier_phase_ATminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTG_fourier_phase_test) = paste0("fourier_phase_ATminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATminusTT_fourier_phase_test) = paste0("fourier_phase_ATminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CAminusCC_fourier_phase_test) = paste0("fourier_phase_CAminusCC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCG_fourier_phase_test) = paste0("fourier_phase_CAminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusCT_fourier_phase_test) = paste0("fourier_phase_CAminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGA_fourier_phase_test) = paste0("fourier_phase_CAminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGC_fourier_phase_test) = paste0("fourier_phase_CAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGG_fourier_phase_test) = paste0("fourier_phase_CAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusGT_fourier_phase_test) = paste0("fourier_phase_CAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTA_fourier_phase_test) = paste0("fourier_phase_CAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTC_fourier_phase_test) = paste0("fourier_phase_CAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTG_fourier_phase_test) = paste0("fourier_phase_CAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CAminusTT_fourier_phase_test) = paste0("fourier_phase_CAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CCminusCG_fourier_phase_test) = paste0("fourier_phase_CCminusCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusCT_fourier_phase_test) = paste0("fourier_phase_CCminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGA_fourier_phase_test) = paste0("fourier_phase_CCminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGC_fourier_phase_test) = paste0("fourier_phase_CCminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGG_fourier_phase_test) = paste0("fourier_phase_CCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusGT_fourier_phase_test) = paste0("fourier_phase_CCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTA_fourier_phase_test) = paste0("fourier_phase_CCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTC_fourier_phase_test) = paste0("fourier_phase_CCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTG_fourier_phase_test) = paste0("fourier_phase_CCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCminusTT_fourier_phase_test) = paste0("fourier_phase_CCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CGminusCT_fourier_phase_test) = paste0("fourier_phase_CGminusCT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGA_fourier_phase_test) = paste0("fourier_phase_CGminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGC_fourier_phase_test) = paste0("fourier_phase_CGminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGG_fourier_phase_test) = paste0("fourier_phase_CGminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusGT_fourier_phase_test) = paste0("fourier_phase_CGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTA_fourier_phase_test) = paste0("fourier_phase_CGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTC_fourier_phase_test) = paste0("fourier_phase_CGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTG_fourier_phase_test) = paste0("fourier_phase_CGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGminusTT_fourier_phase_test) = paste0("fourier_phase_CGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_CTminusGA_fourier_phase_test) = paste0("fourier_phase_CTminusGA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGC_fourier_phase_test) = paste0("fourier_phase_CTminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGG_fourier_phase_test) = paste0("fourier_phase_CTminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusGT_fourier_phase_test) = paste0("fourier_phase_CTminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTA_fourier_phase_test) = paste0("fourier_phase_CTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTC_fourier_phase_test) = paste0("fourier_phase_CTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTG_fourier_phase_test) = paste0("fourier_phase_CTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CTminusTT_fourier_phase_test) = paste0("fourier_phase_CTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GAminusGC_fourier_phase_test) = paste0("fourier_phase_GAminusGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGG_fourier_phase_test) = paste0("fourier_phase_GAminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusGT_fourier_phase_test) = paste0("fourier_phase_GAminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTA_fourier_phase_test) = paste0("fourier_phase_GAminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTC_fourier_phase_test) = paste0("fourier_phase_GAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTG_fourier_phase_test) = paste0("fourier_phase_GAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GAminusTT_fourier_phase_test) = paste0("fourier_phase_GAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GCminusGG_fourier_phase_test) = paste0("fourier_phase_GCminusGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusGT_fourier_phase_test) = paste0("fourier_phase_GCminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTA_fourier_phase_test) = paste0("fourier_phase_GCminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTC_fourier_phase_test) = paste0("fourier_phase_GCminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTG_fourier_phase_test) = paste0("fourier_phase_GCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCminusTT_fourier_phase_test) = paste0("fourier_phase_GCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GGminusGT_fourier_phase_test) = paste0("fourier_phase_GGminusGT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTA_fourier_phase_test) = paste0("fourier_phase_GGminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTC_fourier_phase_test) = paste0("fourier_phase_GGminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTG_fourier_phase_test) = paste0("fourier_phase_GGminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GGminusTT_fourier_phase_test) = paste0("fourier_phase_GGminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_GTminusTA_fourier_phase_test) = paste0("fourier_phase_GTminusTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTC_fourier_phase_test) = paste0("fourier_phase_GTminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTG_fourier_phase_test) = paste0("fourier_phase_GTminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GTminusTT_fourier_phase_test) = paste0("fourier_phase_GTminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TAminusTC_fourier_phase_test) = paste0("fourier_phase_TAminusTC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTG_fourier_phase_test) = paste0("fourier_phase_TAminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAminusTT_fourier_phase_test) = paste0("fourier_phase_TAminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TCminusTG_fourier_phase_test) = paste0("fourier_phase_TCminusTG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TCminusTT_fourier_phase_test) = paste0("fourier_phase_TCminusTT_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_TGminusTT_fourier_phase_test) = paste0("fourier_phase_TGminusTT_", c(1,2,3,4,5,6,11,21))

Xtwo_diminusdi_fourier_mag_test = cbind(Xtwo_AAminusAC_fourier_mag_test,
                                        Xtwo_AAminusAG_fourier_mag_test,
                                        Xtwo_AAminusAT_fourier_mag_test,
                                        Xtwo_AAminusCA_fourier_mag_test,
                                        Xtwo_AAminusCC_fourier_mag_test,
                                        Xtwo_AAminusCG_fourier_mag_test,
                                        Xtwo_AAminusCT_fourier_mag_test,
                                        Xtwo_AAminusGA_fourier_mag_test,
                                        Xtwo_AAminusGC_fourier_mag_test,
                                        Xtwo_AAminusGG_fourier_mag_test,
                                        Xtwo_AAminusGT_fourier_mag_test,
                                        Xtwo_AAminusTA_fourier_mag_test,
                                        Xtwo_AAminusTC_fourier_mag_test,
                                        Xtwo_AAminusTG_fourier_mag_test,
                                        Xtwo_AAminusTT_fourier_mag_test,
                                        
                                        Xtwo_ACminusAG_fourier_mag_test,
                                        Xtwo_ACminusAT_fourier_mag_test,
                                        Xtwo_ACminusCA_fourier_mag_test,
                                        Xtwo_ACminusCC_fourier_mag_test,
                                        Xtwo_ACminusCG_fourier_mag_test,
                                        Xtwo_ACminusCT_fourier_mag_test,
                                        Xtwo_ACminusGA_fourier_mag_test,
                                        Xtwo_ACminusGC_fourier_mag_test,
                                        Xtwo_ACminusGG_fourier_mag_test,
                                        Xtwo_ACminusGT_fourier_mag_test,
                                        Xtwo_ACminusTA_fourier_mag_test,
                                        Xtwo_ACminusTC_fourier_mag_test,
                                        Xtwo_ACminusTG_fourier_mag_test,
                                        Xtwo_ACminusTT_fourier_mag_test,
                                        
                                        Xtwo_AGminusAT_fourier_mag_test,
                                        Xtwo_AGminusCA_fourier_mag_test,
                                        Xtwo_AGminusCC_fourier_mag_test,
                                        Xtwo_AGminusCG_fourier_mag_test,
                                        Xtwo_AGminusCT_fourier_mag_test,
                                        Xtwo_AGminusGA_fourier_mag_test,
                                        Xtwo_AGminusGC_fourier_mag_test,
                                        Xtwo_AGminusGG_fourier_mag_test,
                                        Xtwo_AGminusGT_fourier_mag_test,
                                        Xtwo_AGminusTA_fourier_mag_test,
                                        Xtwo_AGminusTC_fourier_mag_test,
                                        Xtwo_AGminusTG_fourier_mag_test,
                                        Xtwo_AGminusTT_fourier_mag_test,
                                        
                                        Xtwo_ATminusCA_fourier_mag_test,
                                        Xtwo_ATminusCC_fourier_mag_test,
                                        Xtwo_ATminusCG_fourier_mag_test,
                                        Xtwo_ATminusCT_fourier_mag_test,
                                        Xtwo_ATminusGA_fourier_mag_test,
                                        Xtwo_ATminusGC_fourier_mag_test,
                                        Xtwo_ATminusGG_fourier_mag_test,
                                        Xtwo_ATminusGT_fourier_mag_test,
                                        Xtwo_ATminusTA_fourier_mag_test,
                                        Xtwo_ATminusTC_fourier_mag_test,
                                        Xtwo_ATminusTG_fourier_mag_test,
                                        Xtwo_ATminusTT_fourier_mag_test,
                                        
                                        Xtwo_CAminusCC_fourier_mag_test,
                                        Xtwo_CAminusCG_fourier_mag_test,
                                        Xtwo_CAminusCT_fourier_mag_test,
                                        Xtwo_CAminusGA_fourier_mag_test,
                                        Xtwo_CAminusGC_fourier_mag_test,
                                        Xtwo_CAminusGG_fourier_mag_test,
                                        Xtwo_CAminusGT_fourier_mag_test,
                                        Xtwo_CAminusTA_fourier_mag_test,
                                        Xtwo_CAminusTC_fourier_mag_test,
                                        Xtwo_CAminusTG_fourier_mag_test,
                                        Xtwo_CAminusTT_fourier_mag_test,
                                        
                                        Xtwo_CCminusCG_fourier_mag_test,
                                        Xtwo_CCminusCT_fourier_mag_test,
                                        Xtwo_CCminusGA_fourier_mag_test,
                                        Xtwo_CCminusGC_fourier_mag_test,
                                        Xtwo_CCminusGG_fourier_mag_test,
                                        Xtwo_CCminusGT_fourier_mag_test,
                                        Xtwo_CCminusTA_fourier_mag_test,
                                        Xtwo_CCminusTC_fourier_mag_test,
                                        Xtwo_CCminusTG_fourier_mag_test,
                                        Xtwo_CCminusTT_fourier_mag_test,
                                        
                                        Xtwo_CGminusCT_fourier_mag_test,
                                        Xtwo_CGminusGA_fourier_mag_test,
                                        Xtwo_CGminusGC_fourier_mag_test,
                                        Xtwo_CGminusGG_fourier_mag_test,
                                        Xtwo_CGminusGT_fourier_mag_test,
                                        Xtwo_CGminusTA_fourier_mag_test,
                                        Xtwo_CGminusTC_fourier_mag_test,
                                        Xtwo_CGminusTG_fourier_mag_test,
                                        Xtwo_CGminusTT_fourier_mag_test,
                                        
                                        Xtwo_CTminusGA_fourier_mag_test,
                                        Xtwo_CTminusGC_fourier_mag_test,
                                        Xtwo_CTminusGG_fourier_mag_test,
                                        Xtwo_CTminusGT_fourier_mag_test,
                                        Xtwo_CTminusTA_fourier_mag_test,
                                        Xtwo_CTminusTC_fourier_mag_test,
                                        Xtwo_CTminusTG_fourier_mag_test,
                                        Xtwo_CTminusTT_fourier_mag_test,
                                        
                                        Xtwo_GAminusGC_fourier_mag_test,
                                        Xtwo_GAminusGG_fourier_mag_test,
                                        Xtwo_GAminusGT_fourier_mag_test,
                                        Xtwo_GAminusTA_fourier_mag_test,
                                        Xtwo_GAminusTC_fourier_mag_test,
                                        Xtwo_GAminusTG_fourier_mag_test,
                                        Xtwo_GAminusTT_fourier_mag_test,
                                        
                                        Xtwo_GCminusGG_fourier_mag_test,
                                        Xtwo_GCminusGT_fourier_mag_test,
                                        Xtwo_GCminusTA_fourier_mag_test,
                                        Xtwo_GCminusTC_fourier_mag_test,
                                        Xtwo_GCminusTG_fourier_mag_test,
                                        Xtwo_GCminusTT_fourier_mag_test,
                                        
                                        Xtwo_GGminusGT_fourier_mag_test,
                                        Xtwo_GGminusTA_fourier_mag_test,
                                        Xtwo_GGminusTC_fourier_mag_test,
                                        Xtwo_GGminusTG_fourier_mag_test,
                                        Xtwo_GGminusTT_fourier_mag_test,
                                        
                                        Xtwo_GTminusTA_fourier_mag_test,
                                        Xtwo_GTminusTC_fourier_mag_test,
                                        Xtwo_GTminusTG_fourier_mag_test,
                                        Xtwo_GTminusTT_fourier_mag_test,
                                        
                                        Xtwo_TAminusTC_fourier_mag_test,
                                        Xtwo_TAminusTG_fourier_mag_test,
                                        Xtwo_TAminusTT_fourier_mag_test,
                                        
                                        Xtwo_TCminusTG_fourier_mag_test,
                                        Xtwo_TCminusTT_fourier_mag_test,
                                        
                                        Xtwo_TGminusTT_fourier_mag_test)

Xtwo_diminusdi_fourier_phase_test = cbind(Xtwo_AAminusAC_fourier_phase_test,
                                          Xtwo_AAminusAG_fourier_phase_test,
                                          Xtwo_AAminusAT_fourier_phase_test,
                                          Xtwo_AAminusCA_fourier_phase_test,
                                          Xtwo_AAminusCC_fourier_phase_test,
                                          Xtwo_AAminusCG_fourier_phase_test,
                                          Xtwo_AAminusCT_fourier_phase_test,
                                          Xtwo_AAminusGA_fourier_phase_test,
                                          Xtwo_AAminusGC_fourier_phase_test,
                                          Xtwo_AAminusGG_fourier_phase_test,
                                          Xtwo_AAminusGT_fourier_phase_test,
                                          Xtwo_AAminusTA_fourier_phase_test,
                                          Xtwo_AAminusTC_fourier_phase_test,
                                          Xtwo_AAminusTG_fourier_phase_test,
                                          Xtwo_AAminusTT_fourier_phase_test,
                                          
                                          Xtwo_ACminusAG_fourier_phase_test,
                                          Xtwo_ACminusAT_fourier_phase_test,
                                          Xtwo_ACminusCA_fourier_phase_test,
                                          Xtwo_ACminusCC_fourier_phase_test,
                                          Xtwo_ACminusCG_fourier_phase_test,
                                          Xtwo_ACminusCT_fourier_phase_test,
                                          Xtwo_ACminusGA_fourier_phase_test,
                                          Xtwo_ACminusGC_fourier_phase_test,
                                          Xtwo_ACminusGG_fourier_phase_test,
                                          Xtwo_ACminusGT_fourier_phase_test,
                                          Xtwo_ACminusTA_fourier_phase_test,
                                          Xtwo_ACminusTC_fourier_phase_test,
                                          Xtwo_ACminusTG_fourier_phase_test,
                                          Xtwo_ACminusTT_fourier_phase_test,
                                          
                                          Xtwo_AGminusAT_fourier_phase_test,
                                          Xtwo_AGminusCA_fourier_phase_test,
                                          Xtwo_AGminusCC_fourier_phase_test,
                                          Xtwo_AGminusCG_fourier_phase_test,
                                          Xtwo_AGminusCT_fourier_phase_test,
                                          Xtwo_AGminusGA_fourier_phase_test,
                                          Xtwo_AGminusGC_fourier_phase_test,
                                          Xtwo_AGminusGG_fourier_phase_test,
                                          Xtwo_AGminusGT_fourier_phase_test,
                                          Xtwo_AGminusTA_fourier_phase_test,
                                          Xtwo_AGminusTC_fourier_phase_test,
                                          Xtwo_AGminusTG_fourier_phase_test,
                                          Xtwo_AGminusTT_fourier_phase_test,
                                          
                                          Xtwo_ATminusCA_fourier_phase_test,
                                          Xtwo_ATminusCC_fourier_phase_test,
                                          Xtwo_ATminusCG_fourier_phase_test,
                                          Xtwo_ATminusCT_fourier_phase_test,
                                          Xtwo_ATminusGA_fourier_phase_test,
                                          Xtwo_ATminusGC_fourier_phase_test,
                                          Xtwo_ATminusGG_fourier_phase_test,
                                          Xtwo_ATminusGT_fourier_phase_test,
                                          Xtwo_ATminusTA_fourier_phase_test,
                                          Xtwo_ATminusTC_fourier_phase_test,
                                          Xtwo_ATminusTG_fourier_phase_test,
                                          Xtwo_ATminusTT_fourier_phase_test,
                                          
                                          Xtwo_CAminusCC_fourier_phase_test,
                                          Xtwo_CAminusCG_fourier_phase_test,
                                          Xtwo_CAminusCT_fourier_phase_test,
                                          Xtwo_CAminusGA_fourier_phase_test,
                                          Xtwo_CAminusGC_fourier_phase_test,
                                          Xtwo_CAminusGG_fourier_phase_test,
                                          Xtwo_CAminusGT_fourier_phase_test,
                                          Xtwo_CAminusTA_fourier_phase_test,
                                          Xtwo_CAminusTC_fourier_phase_test,
                                          Xtwo_CAminusTG_fourier_phase_test,
                                          Xtwo_CAminusTT_fourier_phase_test,
                                          
                                          Xtwo_CCminusCG_fourier_phase_test,
                                          Xtwo_CCminusCT_fourier_phase_test,
                                          Xtwo_CCminusGA_fourier_phase_test,
                                          Xtwo_CCminusGC_fourier_phase_test,
                                          Xtwo_CCminusGG_fourier_phase_test,
                                          Xtwo_CCminusGT_fourier_phase_test,
                                          Xtwo_CCminusTA_fourier_phase_test,
                                          Xtwo_CCminusTC_fourier_phase_test,
                                          Xtwo_CCminusTG_fourier_phase_test,
                                          Xtwo_CCminusTT_fourier_phase_test,
                                          
                                          Xtwo_CGminusCT_fourier_phase_test,
                                          Xtwo_CGminusGA_fourier_phase_test,
                                          Xtwo_CGminusGC_fourier_phase_test,
                                          Xtwo_CGminusGG_fourier_phase_test,
                                          Xtwo_CGminusGT_fourier_phase_test,
                                          Xtwo_CGminusTA_fourier_phase_test,
                                          Xtwo_CGminusTC_fourier_phase_test,
                                          Xtwo_CGminusTG_fourier_phase_test,
                                          Xtwo_CGminusTT_fourier_phase_test,
                                          
                                          Xtwo_CTminusGA_fourier_phase_test,
                                          Xtwo_CTminusGC_fourier_phase_test,
                                          Xtwo_CTminusGG_fourier_phase_test,
                                          Xtwo_CTminusGT_fourier_phase_test,
                                          Xtwo_CTminusTA_fourier_phase_test,
                                          Xtwo_CTminusTC_fourier_phase_test,
                                          Xtwo_CTminusTG_fourier_phase_test,
                                          Xtwo_CTminusTT_fourier_phase_test,
                                          
                                          Xtwo_GAminusGC_fourier_phase_test,
                                          Xtwo_GAminusGG_fourier_phase_test,
                                          Xtwo_GAminusGT_fourier_phase_test,
                                          Xtwo_GAminusTA_fourier_phase_test,
                                          Xtwo_GAminusTC_fourier_phase_test,
                                          Xtwo_GAminusTG_fourier_phase_test,
                                          Xtwo_GAminusTT_fourier_phase_test,
                                          
                                          Xtwo_GCminusGG_fourier_phase_test,
                                          Xtwo_GCminusGT_fourier_phase_test,
                                          Xtwo_GCminusTA_fourier_phase_test,
                                          Xtwo_GCminusTC_fourier_phase_test,
                                          Xtwo_GCminusTG_fourier_phase_test,
                                          Xtwo_GCminusTT_fourier_phase_test,
                                          
                                          Xtwo_GGminusGT_fourier_phase_test,
                                          Xtwo_GGminusTA_fourier_phase_test,
                                          Xtwo_GGminusTC_fourier_phase_test,
                                          Xtwo_GGminusTG_fourier_phase_test,
                                          Xtwo_GGminusTT_fourier_phase_test,
                                          
                                          Xtwo_GTminusTA_fourier_phase_test,
                                          Xtwo_GTminusTC_fourier_phase_test,
                                          Xtwo_GTminusTG_fourier_phase_test,
                                          Xtwo_GTminusTT_fourier_phase_test,
                                          
                                          Xtwo_TAminusTC_fourier_phase_test,
                                          Xtwo_TAminusTG_fourier_phase_test,
                                          Xtwo_TAminusTT_fourier_phase_test,
                                          
                                          Xtwo_TCminusTG_fourier_phase_test,
                                          Xtwo_TCminusTT_fourier_phase_test,
                                          
                                          Xtwo_TGminusTT_fourier_phase_test)
saveRDS(Xtwo_diminusdi_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diminusdi_fourier_mag_test.rds")

saveRDS(Xtwo_diminusdi_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diminusdi_fourier_phase_test.rds")


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

Xtwo_AAorAT_fourier_test = t(apply(Xtwo_AAorAT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AAorTA_fourier_test = t(apply(Xtwo_AAorTA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_AAorTT_fourier_test = t(apply(Xtwo_AAorTT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_ATorTA_fourier_test = t(apply(Xtwo_ATorTA_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_ATorTT_fourier_test = t(apply(Xtwo_ATorTT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_TAorTT_fourier_test = t(apply(Xtwo_TAorTT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorCG_fourier_test = t(apply(Xtwo_CCorCG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorGC_fourier_test = t(apply(Xtwo_CCorGC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorGG_fourier_test = t(apply(Xtwo_CCorGG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CGorGC_fourier_test = t(apply(Xtwo_CGorGC_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CGorGG_fourier_test = t(apply(Xtwo_CGorGG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_GCorGG_fourier_test = t(apply(Xtwo_GCorGG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]

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


colnames(Xtwo_AAorAT_fourier_mag_test) = paste0("fourier_mag_AAorAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTA_fourier_mag_test) = paste0("fourier_mag_AAorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTT_fourier_mag_test) = paste0("fourier_mag_AAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTA_fourier_mag_test) = paste0("fourier_mag_ATorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTT_fourier_mag_test) = paste0("fourier_mag_ATorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAorTT_fourier_mag_test) = paste0("fourier_mag_TAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCG_fourier_mag_test) = paste0("fourier_mag_CCorCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGC_fourier_mag_test) = paste0("fourier_mag_CCorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGG_fourier_mag_test) = paste0("fourier_mag_CCorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGC_fourier_mag_test) = paste0("fourier_mag_CGorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGG_fourier_mag_test) = paste0("fourier_mag_CGorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCorGG_fourier_mag_test) = paste0("fourier_mag_GCorGG_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAorAT_fourier_phase_test) = paste0("fourier_phase_AAorAT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTA_fourier_phase_test) = paste0("fourier_phase_AAorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_AAorTT_fourier_phase_test) = paste0("fourier_phase_AAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTA_fourier_phase_test) = paste0("fourier_phase_ATorTA_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_ATorTT_fourier_phase_test) = paste0("fourier_phase_ATorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_TAorTT_fourier_phase_test) = paste0("fourier_phase_TAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCG_fourier_phase_test) = paste0("fourier_phase_CCorCG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGC_fourier_phase_test) = paste0("fourier_phase_CCorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorGG_fourier_phase_test) = paste0("fourier_phase_CCorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGC_fourier_phase_test) = paste0("fourier_phase_CGorGC_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CGorGG_fourier_phase_test) = paste0("fourier_phase_CGorGG_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_GCorGG_fourier_phase_test) = paste0("fourier_phase_GCorGG_", c(1,2,3,4,5,6,11,21))

Xtwo_diordi_fourier_mag_test = cbind(Xtwo_AAorAT_fourier_mag_test,
                                     Xtwo_AAorTA_fourier_mag_test,
                                     Xtwo_AAorTT_fourier_mag_test,
                                     Xtwo_ATorTA_fourier_mag_test,
                                     Xtwo_ATorTT_fourier_mag_test,
                                     Xtwo_TAorTT_fourier_mag_test,
                                     Xtwo_CCorCG_fourier_mag_test,
                                     Xtwo_CCorGC_fourier_mag_test,
                                     Xtwo_CCorGG_fourier_mag_test,
                                     Xtwo_CGorGC_fourier_mag_test,
                                     Xtwo_CGorGG_fourier_mag_test,
                                     Xtwo_GCorGG_fourier_mag_test)
saveRDS(Xtwo_diordi_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordi_fourier_mag_test.rds")

Xtwo_diordi_fourier_phase_test = cbind(Xtwo_AAorAT_fourier_phase_test,
                                       Xtwo_AAorTA_fourier_phase_test,
                                       Xtwo_AAorTT_fourier_phase_test,
                                       Xtwo_ATorTA_fourier_phase_test,
                                       Xtwo_ATorTT_fourier_phase_test,
                                       Xtwo_TAorTT_fourier_phase_test,
                                       Xtwo_CCorCG_fourier_phase_test,
                                       Xtwo_CCorGC_fourier_phase_test,
                                       Xtwo_CCorGG_fourier_phase_test,
                                       Xtwo_CGorGC_fourier_phase_test,
                                       Xtwo_CGorGG_fourier_phase_test,
                                       Xtwo_GCorGG_fourier_phase_test)
saveRDS(Xtwo_diordi_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordi_fourier_phase_test.rds")


# Groups of Dinucleotides:

Xtwo_AAorATorTAorTT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CCorCGorGCorGG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AAorATorTAorTT_test) = colnames(Xtwo_test)
colnames(Xtwo_CCorCGorGCorGG_test) = colnames(Xtwo_test)
Xtwo_AAorATorTAorTT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "AT") | 
                                (Xtwo_test == "TA") | (Xtwo_test == "TT")) %>% as.matrix() %>% as.numeric()
Xtwo_CCorCGorGCorGG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "CG") |
                                (Xtwo_test == "GC") | (Xtwo_test == "GG")) %>% as.matrix() %>% as.numeric()

Xtwo_AAorATorTAorTT_fourier_test = t(apply(Xtwo_AAorATorTAorTT_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]
Xtwo_CCorCGorGCorGG_fourier_test = t(apply(Xtwo_CCorCGorGCorGG_test, 1, fft))[,c(1,2,3,4,5,6,11,21)]

Xtwo_AAorATorTAorTT_fourier_mag_test = as.matrix(Mod(Xtwo_AAorATorTAorTT_fourier_test))
Xtwo_CCorCGorGCorGG_fourier_mag_test = as.matrix(Mod(Xtwo_CCorCGorGCorGG_fourier_test))

Xtwo_AAorATorTAorTT_fourier_phase_test = as.matrix(Arg(Xtwo_AAorATorTAorTT_fourier_test))
Xtwo_CCorCGorGCorGG_fourier_phase_test = as.matrix(Arg(Xtwo_CCorCGorGCorGG_fourier_test))


colnames(Xtwo_AAorATorTAorTT_fourier_mag_test) = paste0("fourier_mag_AAorATorTAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCGorGCorGG_fourier_mag_test) = paste0("fourier_mag_CCorCGorGCorGG_", c(1,2,3,4,5,6,11,21))

colnames(Xtwo_AAorATorTAorTT_fourier_phase_test) = paste0("fourier_phase_AAorATorTAorTT_", c(1,2,3,4,5,6,11,21))
colnames(Xtwo_CCorCGorGCorGG_fourier_phase_test) = paste0("fourier_phase_CCorCGorGCorGG_", c(1,2,3,4,5,6,11,21))

Xtwo_diordiordiordi_fourier_mag_test = cbind(Xtwo_AAorATorTAorTT_fourier_mag_test,
                                             Xtwo_CCorCGorGCorGG_fourier_mag_test)
saveRDS(Xtwo_diordiordiordi_fourier_mag_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordiordiordi_fourier_mag_test.rds")

Xtwo_diordiordiordi_fourier_phase_test = cbind(Xtwo_AAorATorTAorTT_fourier_phase_test,
                                               Xtwo_CCorCGorGCorGG_fourier_phase_test)
saveRDS(Xtwo_diordiordiordi_fourier_phase_test, 
        "data/Created/chrV_post_smooth_C0_Xtwo_diordiordiordi_fourier_phase_test.rds")
