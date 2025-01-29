ps1_freq = apply(dat %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  C_freq = sum(col == "C")
  G_freq = sum(col == "G")
  T_freq = sum(col == "T")
  return(c(A_freq, C_freq, G_freq, T_freq))
})
rownames(ps1_freq) = nucleotides

ps2_freq = apply(dat %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  AC_freq = sum(col == "AC")
  AG_freq = sum(col == "AG")
  AT_freq = sum(col == "AT")
  CA_freq = sum(col == "CA")
  CC_freq = sum(col == "CC")
  CG_freq = sum(col == "CG")
  CT_freq = sum(col == "CT")
  GA_freq = sum(col == "GA")
  GC_freq = sum(col == "GC")
  GG_freq = sum(col == "GG")
  GT_freq = sum(col == "GT")
  TA_freq = sum(col == "TA")
  TC_freq = sum(col == "TC")
  TG_freq = sum(col == "TG")
  TT_freq = sum(col == "TT")
  return(c(AA_freq, AC_freq, AG_freq, AT_freq,
           CA_freq, CC_freq, CG_freq, CT_freq,
           GA_freq, GC_freq, GG_freq, GT_freq,
           TA_freq, TC_freq, TG_freq, TT_freq))
})
rownames(ps2_freq) = dinucleotides

exp_ps2_freq = matrix(nrow = 16, ncol = 49)
colnames(exp_ps2_freq) = paste0("X", (1:49), "di")
rownames(exp_ps2_freq) = dinucleotides
for(i in 1:49) {
  exp_ps2_mat = (ps1_freq[,i]/nrow(dat)) %*% (t(ps1_freq[,i+1])/nrow(dat))
  exp_ps2_freq[,i] = as.vector(t(exp_ps2_mat))*nrow(dat)
}

chisq = (ps2_freq - exp_ps2_freq)^2/(exp_ps2_freq)

apply(chisq, 2, sum)

1-dchisq(apply(chisq, 2, sum), df=11)









################################################################################
# Quartiles:
################################################################################

ps1_freq_q1 = apply(dat_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  C_freq = sum(col == "C")
  G_freq = sum(col == "G")
  T_freq = sum(col == "T")
  return(c(A_freq, C_freq, G_freq, T_freq))
})
rownames(ps1_freq_q1) = nucleotides

ps2_freq_q1 = apply(dat_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  AC_freq = sum(col == "AC")
  AG_freq = sum(col == "AG")
  AT_freq = sum(col == "AT")
  CA_freq = sum(col == "CA")
  CC_freq = sum(col == "CC")
  CG_freq = sum(col == "CG")
  CT_freq = sum(col == "CT")
  GA_freq = sum(col == "GA")
  GC_freq = sum(col == "GC")
  GG_freq = sum(col == "GG")
  GT_freq = sum(col == "GT")
  TA_freq = sum(col == "TA")
  TC_freq = sum(col == "TC")
  TG_freq = sum(col == "TG")
  TT_freq = sum(col == "TT")
  return(c(AA_freq, AC_freq, AG_freq, AT_freq,
           CA_freq, CC_freq, CG_freq, CT_freq,
           GA_freq, GC_freq, GG_freq, GT_freq,
           TA_freq, TC_freq, TG_freq, TT_freq))
})
rownames(ps2_freq_q1) = dinucleotides

exp_ps2_freq_q1 = matrix(nrow = 16, ncol = 49)
colnames(exp_ps2_freq_q1) = paste0("X", (1:49), "di")
rownames(exp_ps2_freq_q1) = dinucleotides
for(i in 1:49) {
  exp_ps2_mat_q1 = (ps1_freq_q1[,i]/nrow(dat_q1)) %*% (t(ps1_freq_q1[,i+1])/nrow(dat_q1))
  exp_ps2_freq_q1[,i] = as.vector(t(exp_ps2_mat_q1))*nrow(dat_q1)
}

chisq_q1 = (ps2_freq_q1 - exp_ps2_freq_q1)^2/(exp_ps2_freq_q1)

apply(chisq_q1, 2, sum)

1-pchisq(apply(chisq_q1, 2, sum), df=15)



ps1_freq_q2 = apply(dat_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  C_freq = sum(col == "C")
  G_freq = sum(col == "G")
  T_freq = sum(col == "T")
  return(c(A_freq, C_freq, G_freq, T_freq))
})
rownames(ps1_freq_q2) = nucleotides

ps2_freq_q2 = apply(dat_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  AC_freq = sum(col == "AC")
  AG_freq = sum(col == "AG")
  AT_freq = sum(col == "AT")
  CA_freq = sum(col == "CA")
  CC_freq = sum(col == "CC")
  CG_freq = sum(col == "CG")
  CT_freq = sum(col == "CT")
  GA_freq = sum(col == "GA")
  GC_freq = sum(col == "GC")
  GG_freq = sum(col == "GG")
  GT_freq = sum(col == "GT")
  TA_freq = sum(col == "TA")
  TC_freq = sum(col == "TC")
  TG_freq = sum(col == "TG")
  TT_freq = sum(col == "TT")
  return(c(AA_freq, AC_freq, AG_freq, AT_freq,
           CA_freq, CC_freq, CG_freq, CT_freq,
           GA_freq, GC_freq, GG_freq, GT_freq,
           TA_freq, TC_freq, TG_freq, TT_freq))
})
rownames(ps2_freq_q2) = dinucleotides

exp_ps2_freq_q2 = matrix(nrow = 16, ncol = 49)
colnames(exp_ps2_freq_q2) = paste0("X", (1:49), "di")
rownames(exp_ps2_freq_q2) = dinucleotides
for(i in 1:49) {
  exp_ps2_mat_q2 = (ps1_freq_q2[,i]/nrow(dat_q2)) %*% t((ps1_freq_q2[,i+1]/nrow(dat_q2)))
  exp_ps2_freq_q2[,i] = as.vector(t(exp_ps2_mat_q2))*nrow(dat_q2)
}

chisq_q2 = (ps2_freq_q2 - exp_ps2_freq_q2)^2/(exp_ps2_freq_q2)

apply(chisq_q2, 2, sum)

1-dchisq(sum(chisq_q2), df=11)




ps1_rel_freq_q3 = apply(dat_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  C_freq = sum(col == "C")
  G_freq = sum(col == "G")
  T_freq = sum(col == "T")
  return(c(A_freq, C_freq, G_freq, T_freq)/nrow(dat_q3))
})
rownames(ps1_rel_freq_q3) = nucleotides

ps2_rel_freq_q3 = apply(dat_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  AC_freq = sum(col == "AC")
  AG_freq = sum(col == "AG")
  AT_freq = sum(col == "AT")
  CA_freq = sum(col == "CA")
  CC_freq = sum(col == "CC")
  CG_freq = sum(col == "CG")
  CT_freq = sum(col == "CT")
  GA_freq = sum(col == "GA")
  GC_freq = sum(col == "GC")
  GG_freq = sum(col == "GG")
  GT_freq = sum(col == "GT")
  TA_freq = sum(col == "TA")
  TC_freq = sum(col == "TC")
  TG_freq = sum(col == "TG")
  TT_freq = sum(col == "TT")
  return(c(AA_freq, AC_freq, AG_freq, AT_freq,
           CA_freq, CC_freq, CG_freq, CT_freq,
           GA_freq, GC_freq, GG_freq, GT_freq,
           TA_freq, TC_freq, TG_freq, TT_freq)/nrow(dat_q3))
})
rownames(ps2_rel_freq_q3) = dinucleotides

exp_ps2_rel_freq_q3 = matrix(nrow = 16, ncol = 49)
colnames(exp_ps2_rel_freq_q3) = paste0("X", (1:49), "di")
rownames(exp_ps2_rel_freq_q3) = dinucleotides
for(i in 1:49) {
  exp_ps2_mat_q3 = ps1_rel_freq_q3[,i] %*% t(ps1_rel_freq_q3[,i+1])
  exp_ps2_rel_freq_q3[,i] = as.vector(t(exp_ps2_mat_q3))
}

chisq_q3 = (ps2_rel_freq_q3 - exp_ps2_rel_freq_q3)^2/(exp_ps2_rel_freq_q3)

sum(chisq_q3)

1-dchisq(sum(chisq_q3), df=11)



ps1_freq_q4 = apply(dat_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  C_freq = sum(col == "C")
  G_freq = sum(col == "G")
  T_freq = sum(col == "T")
  return(c(A_freq, C_freq, G_freq, T_freq))
})
rownames(ps1_rel_freq_q4) = nucleotides

ps2_freq_q4 = apply(dat_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  AC_freq = sum(col == "AC")
  AG_freq = sum(col == "AG")
  AT_freq = sum(col == "AT")
  CA_freq = sum(col == "CA")
  CC_freq = sum(col == "CC")
  CG_freq = sum(col == "CG")
  CT_freq = sum(col == "CT")
  GA_freq = sum(col == "GA")
  GC_freq = sum(col == "GC")
  GG_freq = sum(col == "GG")
  GT_freq = sum(col == "GT")
  TA_freq = sum(col == "TA")
  TC_freq = sum(col == "TC")
  TG_freq = sum(col == "TG")
  TT_freq = sum(col == "TT")
  return(c(AA_freq, AC_freq, AG_freq, AT_freq,
           CA_freq, CC_freq, CG_freq, CT_freq,
           GA_freq, GC_freq, GG_freq, GT_freq,
           TA_freq, TC_freq, TG_freq, TT_freq))
})
rownames(ps2_rel_freq_q4) = dinucleotides

exp_ps2_rel_freq_q4 = matrix(nrow = 16, ncol = 49)
colnames(exp_ps2_rel_freq_q4) = paste0("X", (1:49), "di")
rownames(exp_ps2_rel_freq_q4) = dinucleotides
for(i in 1:49) {
  exp_ps2_mat_q4 = ps1_rel_freq_q4[,i] %*% t(ps1_rel_freq_q4[,i+1])
  exp_ps2_rel_freq_q4[,i] = as.vector(t(exp_ps2_mat_q4))
}

chisq_q4 = (ps2_rel_freq_q4 - exp_ps2_rel_freq_q4)^2/(exp_ps2_rel_freq_q4)

apply(chisq_q4, 2, sum)

1-dchisq(sum(chisq_q4), df=11)
