########################################################################
# Tiling
########################################################################

# Train data
Xone = dat %>% select(all_of(ps1))
Xone_AT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
Xone_CG = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AT) = colnames(Xone)
colnames(Xone_CG) = colnames(Xone)
Xone_AT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
Xone_CG[] = ((Xone == "C") | (Xone == "G")) %>% as.matrix() %>% as.numeric()

Y = dat$C0_new

N = 50
# rows: n (frequency), cols: k (location)
# period: N/n
cos_matrix = matrix(nrow = N, ncol = N)
sin_matrix = matrix(nrow = N, ncol = N)
for (n in 1:N) {
  for (k in 1:N) {
    cos_matrix[n, k] = cos(2*pi*n*(k-1)/N)
    sin_matrix[n, k] = sin(2*pi*n*(k-1)/N)
  }
}
Xone_AT_fourier_cos = Xone_AT %*% t(cos_matrix)
Xone_AT_fourier_sin = Xone_AT %*% t(sin_matrix)

Xone_CG_fourier_cos = Xone_CG %*% t(cos_matrix)
Xone_CG_fourier_sin = Xone_CG %*% t(sin_matrix)

colnames(Xone_AT_fourier_cos) = paste0("fouriercosAT", 1:50)
colnames(Xone_AT_fourier_sin) = paste0("fouriersinAT", 1:50)

colnames(Xone_CG_fourier_cos) = paste0("fouriercosCG", 1:50)
colnames(Xone_CG_fourier_sin) = paste0("fouriersinCG", 1:50)

saveRDS(Xone_AT_fourier_cos, "data/Created/tiling_Xone_AT_fourier_cos.rds")
saveRDS(Xone_AT_fourier_sin, "data/Created/tiling_Xone_AT_fourier_sin.rds")

saveRDS(Xone_CG_fourier_cos, "data/Created/tiling_Xone_CG_fourier_cos.rds")
saveRDS(Xone_CG_fourier_sin, "data/Created/tiling_Xone_CG_fourier_sin.rds")


Xtwo = dat %>% select(all_of(ps2))
Xtwo_AT = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
Xtwo_CG = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
colnames(Xtwo_AT) = colnames(Xtwo)
colnames(Xtwo_CG) = colnames(Xtwo)
Xtwo_AT[] = ((Xtwo == "AA") | (Xtwo == "TT") | (Xtwo == "AT") | (Xtwo == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CG[] = ((Xtwo == "CC") | (Xtwo == "GG") | (Xtwo == "CG") | (Xtwo == "GC")) %>% 
  as.matrix() %>% as.numeric()

N2 = 49
# rows: n (frequency), cols: k (location)
# period: N/n
cos_matrix2 = matrix(nrow = N2, ncol = N2)
sin_matrix2 = matrix(nrow = N2, ncol = N2)
for (n in 1:N2) {
  for (k in 1:N2) {
    cos_matrix2[n, k] = cos(2*pi*n*(k-1)/N)
    sin_matrix2[n, k] = sin(2*pi*n*(k-1)/N)
  }
}
Xtwo_AT_fourier_cos = Xtwo_AT %*% t(cos_matrix2)
Xtwo_AT_fourier_sin = Xtwo_AT %*% t(sin_matrix2)

Xtwo_CG_fourier_cos = Xtwo_CG %*% t(cos_matrix2)
Xtwo_CG_fourier_sin = Xtwo_CG %*% t(sin_matrix2)

colnames(Xtwo_AT_fourier_cos) = paste0("fouriercosAT2_", 1:49)
colnames(Xtwo_AT_fourier_sin) = paste0("fouriersinAT2_", 1:49)

colnames(Xtwo_CG_fourier_cos) = paste0("fouriercosCG2_", 1:49)
colnames(Xtwo_CG_fourier_sin) = paste0("fouriersinCG2_", 1:49)

saveRDS(Xtwo_AT_fourier_cos, "data/Created/tiling_Xtwo_AT_fourier_cos.rds")
saveRDS(Xtwo_AT_fourier_sin, "data/Created/tiling_Xtwo_AT_fourier_sin.rds")

saveRDS(Xtwo_CG_fourier_cos, "data/Created/tiling_Xtwo_CG_fourier_cos.rds")
saveRDS(Xtwo_CG_fourier_sin, "data/Created/tiling_Xtwo_CG_fourier_sin.rds")


#FIXME
# Amplitude
Xtwo_AT_fourier_amp = sqrt(Xtwo_AT_fourier_cos^2 + Xtwo_AT_fourier_sin^2)
colnames(Xtwo_AT_fourier_amp) = paste0("fourierampAT2_", 1:49)

Xtwo_CG_fourier_amp = sqrt(Xtwo_CG_fourier_cos^2 + Xtwo_CG_fourier_sin^2)
colnames(Xtwo_CG_fourier_amp) = paste0("fourierampCG2_", 1:49)

saveRDS(Xtwo_AT_fourier_amp, "data/Created/tiling_Xtwo_AT_fourier_amp.rds")
saveRDS(Xtwo_CG_fourier_amp, "data/Created/tiling_Xtwo_CG_fourier_amp.rds")

#FIXME
# Anti-phase measure
Xtwo_fourier_antiphase = (Xtwo_AT_fourier_cos*Xtwo_CG_fourier_cos + 
                            Xtwo_AT_fourier_sin*Xtwo_CG_fourier_sin)/
  (Xtwo_AT_fourier_amp*Xtwo_CG_fourier_amp)
colnames(Xtwo_fourier_antiphase) = paste0("fourier_antiphase_", 1:49)

saveRDS(Xtwo_fourier_antiphase, "data/Created/tiling_Xtwo_fourier_antiphase.rds")


# Distance
Xtwo_fourier_distance = sqrt((Xtwo_AT_fourier_cos - Xtwo_CG_fourier_cos)^2 + 
                               (Xtwo_AT_fourier_sin - Xtwo_CG_fourier_sin)^2)
colnames(Xtwo_fourier_distance) = paste0("fourier2_distance_", 1:49)

saveRDS(Xtwo_fourier_distance, "data/Created/tiling_Xtwo_fourier_distance.rds")



Xthree = dat %>% select(all_of(ps3))
Xthree_AT = matrix(nrow=nrow(Xthree), ncol=ncol(Xthree))
Xthree_CG = matrix(nrow=nrow(Xthree), ncol=ncol(Xthree))
colnames(Xthree_AT) = colnames(Xthree)
colnames(Xthree_CG) = colnames(Xthree)
Xthree_AT[] = ((Xthree == "AAA") | (Xthree == "AAT") | (Xthree == "ATA") | (Xthree == "ATT") |
                 (Xthree == "TAA") | (Xthree == "TAT") | (Xthree == "TTA") | (Xthree == "TTT")) %>% 
  as.matrix() %>% as.numeric()
Xthree_CG[] = ((Xthree == "CCC") | (Xthree == "CCG") | (Xthree == "CGC") | (Xthree == "CGG") | 
                 (Xthree == "GCC") | (Xthree == "GCG") | (Xthree == "GGC") | (Xthree == "GGG")) %>% 
  as.matrix() %>% as.numeric()

N3 = 48
# rows: n (frequency), cols: k (location)
# period: N/n
cos_matrix3 = matrix(nrow = N3, ncol = N3)
sin_matrix3 = matrix(nrow = N3, ncol = N3)
for (n in 1:N3) {
  for (k in 1:N3) {
    cos_matrix3[n, k] = cos(2*pi*n*(k-1)/N)
    sin_matrix3[n, k] = sin(2*pi*n*(k-1)/N)
  }
}
Xthree_AT_fourier_cos = Xthree_AT %*% t(cos_matrix3)
Xthree_AT_fourier_sin = Xthree_AT %*% t(sin_matrix3)

Xthree_CG_fourier_cos = Xthree_CG %*% t(cos_matrix3)
Xthree_CG_fourier_sin = Xthree_CG %*% t(sin_matrix3)

colnames(Xthree_AT_fourier_cos) = paste0("fouriercosAT3_", 1:48)
colnames(Xthree_AT_fourier_sin) = paste0("fouriersinAT3_", 1:48)

colnames(Xthree_CG_fourier_cos) = paste0("fouriercosCG3_", 1:48)
colnames(Xthree_CG_fourier_sin) = paste0("fouriersinCG3_", 1:48)

saveRDS(Xthree_AT_fourier_cos, "data/Created/tiling_Xthree_AT_fourier_cos.rds")
saveRDS(Xthree_AT_fourier_sin, "data/Created/tiling_Xthree_AT_fourier_sin.rds")

saveRDS(Xthree_CG_fourier_cos, "data/Created/tiling_Xthree_CG_fourier_cos.rds")
saveRDS(Xthree_CG_fourier_sin, "data/Created/tiling_Xthree_CG_fourier_sin.rds")








# Test data:
Xone_test = dat_test %>% select(all_of(ps1))
Xone_AT_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
Xone_CG_test = matrix(nrow=nrow(Xone_test), ncol=ncol(Xone_test))
colnames(Xone_AT_test) = colnames(Xone_test)
colnames(Xone_CG_test) = colnames(Xone_test)
Xone_AT_test[] = ((Xone_test == "A") | (Xone_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CG_test[] = ((Xone_test == "C") | (Xone_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AT_fourier_cos_test = Xone_AT_test %*% t(cos_matrix)
Xone_AT_fourier_sin_test = Xone_AT_test %*% t(sin_matrix)

Xone_CG_fourier_cos_test = Xone_CG_test %*% t(cos_matrix)
Xone_CG_fourier_sin_test = Xone_CG_test %*% t(sin_matrix)

colnames(Xone_AT_fourier_cos_test) = paste0("fouriercosAT", 1:50)
colnames(Xone_AT_fourier_sin_test) = paste0("fouriersinAT", 1:50)

colnames(Xone_CG_fourier_cos_test) = paste0("fouriercosCG", 1:50)
colnames(Xone_CG_fourier_sin_test) = paste0("fouriersinCG", 1:50)

saveRDS(Xone_AT_fourier_cos_test, "data/Created/tiling_Xone_AT_fourier_cos_test.rds")
saveRDS(Xone_AT_fourier_sin_test, "data/Created/tiling_Xone_AT_fourier_sin_test.rds")

saveRDS(Xone_CG_fourier_cos_test, "data/Created/tiling_Xone_CG_fourier_cos_test.rds")
saveRDS(Xone_CG_fourier_sin_test, "data/Created/tiling_Xone_CG_fourier_sin_test.rds")



Xtwo_test = dat_test %>% select(all_of(ps2))
Xtwo_AT_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
Xtwo_CG_test = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
colnames(Xtwo_AT_test) = colnames(Xtwo_test)
colnames(Xtwo_CG_test) = colnames(Xtwo_test)
Xtwo_AT_test[] = ((Xtwo_test == "AA") | (Xtwo_test == "TT") | (Xtwo_test == "AT") | 
                    (Xtwo_test == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_CG_test[] = ((Xtwo_test == "CC") | (Xtwo_test == "GG") | (Xtwo_test == "CG") | 
                    (Xtwo_test == "GC")) %>% as.matrix() %>% as.numeric()

Xtwo_AT_fourier_cos_test = Xtwo_AT_test %*% t(cos_matrix2)
Xtwo_AT_fourier_sin_test = Xtwo_AT_test %*% t(sin_matrix2)

Xtwo_CG_fourier_cos_test = Xtwo_CG_test %*% t(cos_matrix2)
Xtwo_CG_fourier_sin_test = Xtwo_CG_test %*% t(sin_matrix2)

colnames(Xtwo_AT_fourier_cos_test) = paste0("fouriercosAT2_", 1:49)
colnames(Xtwo_AT_fourier_sin_test) = paste0("fouriersinAT2_", 1:49)

colnames(Xtwo_CG_fourier_cos_test) = paste0("fouriercosCG2_", 1:49)
colnames(Xtwo_CG_fourier_sin_test) = paste0("fouriersinCG2_", 1:49)

saveRDS(Xtwo_AT_fourier_cos_test, "data/Created/tiling_Xtwo_AT_fourier_cos_test.rds")
saveRDS(Xtwo_AT_fourier_sin_test, "data/Created/tiling_Xtwo_AT_fourier_sin_test.rds")

saveRDS(Xtwo_CG_fourier_cos_test, "data/Created/tiling_Xtwo_CG_fourier_cos_test.rds")
saveRDS(Xtwo_CG_fourier_sin_test, "data/Created/tiling_Xtwo_CG_fourier_sin_test.rds")



#FIXME
# Amplitude
Xtwo_AT_fourier_amp_test = sqrt(Xtwo_AT_fourier_cos_test^2 + Xtwo_AT_fourier_sin_test^2)
colnames(Xtwo_AT_fourier_amp_test) = paste0("fourierampAT2_", 1:49)

Xtwo_CG_fourier_amp_test = sqrt(Xtwo_CG_fourier_cos_test^2 + Xtwo_CG_fourier_sin_test^2)
colnames(Xtwo_CG_fourier_amp_test) = paste0("fourierampCG2_", 1:49)

saveRDS(Xtwo_AT_fourier_amp_test, "data/Created/tiling_Xtwo_AT_fourier_amp_test.rds")
saveRDS(Xtwo_CG_fourier_amp_test, "data/Created/tiling_Xtwo_CG_fourier_amp_test.rds")

#FIXME
# Anti-phase measure
Xtwo_fourier_antiphase_test = (Xtwo_AT_fourier_cos_test*Xtwo_CG_fourier_cos_test + 
                                 Xtwo_AT_fourier_sin_test*Xtwo_CG_fourier_sin_test)/
  (Xtwo_AT_fourier_amp_test*Xtwo_CG_fourier_amp_test)
colnames(Xtwo_fourier_antiphase_test) = paste0("fourier_antiphase_", 1:49)

saveRDS(Xtwo_fourier_antiphase_test, "data/Created/tiling_Xtwo_fourier_antiphase_test.rds")


# Distance
Xtwo_fourier_distance_test = sqrt((Xtwo_AT_fourier_cos_test - Xtwo_CG_fourier_cos_test)^2 + 
                                    (Xtwo_AT_fourier_sin_test - Xtwo_CG_fourier_sin_test)^2)
colnames(Xtwo_fourier_distance_test) = paste0("fourier2_distance_", 1:49)

saveRDS(Xtwo_fourier_distance_test, "data/Created/tiling_Xtwo_fourier_distance_test.rds")


Xthree_test = dat_test %>% select(all_of(ps3))
Xthree_AT_test = matrix(nrow=nrow(Xthree_test), ncol=ncol(Xthree_test))
Xthree_CG_test = matrix(nrow=nrow(Xthree_test), ncol=ncol(Xthree_test))
colnames(Xthree_AT_test) = colnames(Xthree_test)
colnames(Xthree_CG_test) = colnames(Xthree_test)

Xthree_AT_test[] = ((Xthree_test == "AAA") | (Xthree_test == "AAT") | (Xthree_test == "ATA") | (Xthree_test == "ATT") |
                      (Xthree_test == "TAA") | (Xthree_test == "TAT") | (Xthree_test == "TTA") | (Xthree_test == "TTT")) %>%
  as.matrix() %>% as.numeric()
Xthree_CG_test[] = ((Xthree_test == "CCC") | (Xthree_test == "CCG") | (Xthree_test == "CGC") | (Xthree_test == "CGG") | 
                      (Xthree_test == "GCC") | (Xthree_test == "GCG") | (Xthree_test == "GGC") | (Xthree_test == "GGG")) %>% 
  as.matrix() %>% as.numeric()

Xthree_AT_fourier_cos_test = Xthree_AT_test %*% t(cos_matrix3)
Xthree_AT_fourier_sin_test = Xthree_AT_test %*% t(sin_matrix3)

Xthree_CG_fourier_cos_test = Xthree_CG_test %*% t(cos_matrix3)
Xthree_CG_fourier_sin_test = Xthree_CG_test %*% t(sin_matrix3)

colnames(Xthree_AT_fourier_cos_test) = paste0("fouriercosAT3_", 1:48)
colnames(Xthree_AT_fourier_sin_test) = paste0("fouriersinAT3_", 1:48)

colnames(Xthree_CG_fourier_cos_test) = paste0("fouriercosCG3_", 1:48)
colnames(Xthree_CG_fourier_sin_test) = paste0("fouriersinCG3_", 1:48)

saveRDS(Xthree_AT_fourier_cos_test, "data/Created/tiling_Xthree_AT_fourier_cos_test.rds")
saveRDS(Xthree_AT_fourier_sin_test, "data/Created/tiling_Xthree_AT_fourier_sin_test.rds")

saveRDS(Xthree_CG_fourier_cos_test, "data/Created/tiling_Xthree_CG_fourier_cos_test.rds")
saveRDS(Xthree_CG_fourier_sin_test, "data/Created/tiling_Xthree_CG_fourier_sin_test.rds")





########################################################################
# Random
########################################################################

# Train data
Xone_random = dat_random %>% select(all_of(ps1))
Xone_AT_random = matrix(nrow=nrow(Xone_random), ncol=ncol(Xone_random))
Xone_CG_random = matrix(nrow=nrow(Xone_random), ncol=ncol(Xone_random))
colnames(Xone_AT_random) = colnames(Xone_random)
colnames(Xone_CG_random) = colnames(Xone_random)
Xone_AT_random[] = ((Xone_random == "A") | (Xone_random == "T")) %>% as.matrix() %>% as.numeric()
Xone_CG_random[] = ((Xone_random == "C") | (Xone_random == "G")) %>% as.matrix() %>% as.numeric()

Xone_AT_random_fourier_cos = Xone_AT_random %*% t(cos_matrix)
Xone_AT_random_fourier_sin = Xone_AT_random %*% t(sin_matrix)

Xone_CG_random_fourier_cos = Xone_CG_random %*% t(cos_matrix)
Xone_CG_random_fourier_sin = Xone_CG_random %*% t(sin_matrix)

colnames(Xone_AT_random_fourier_cos) = paste0("fouriercosAT", 1:50)
colnames(Xone_AT_random_fourier_sin) = paste0("fouriersinAT", 1:50)

colnames(Xone_CG_random_fourier_cos) = paste0("fouriercosCG", 1:50)
colnames(Xone_CG_random_fourier_sin) = paste0("fouriersinCG", 1:50)

saveRDS(Xone_AT_random_fourier_cos, "data/Created/random_Xone_AT_fourier_cos.rds")
saveRDS(Xone_AT_random_fourier_sin, "data/Created/random_Xone_AT_fourier_sin.rds")

saveRDS(Xone_CG_random_fourier_cos, "data/Created/random_Xone_CG_fourier_cos.rds")
saveRDS(Xone_CG_random_fourier_sin, "data/Created/random_Xone_CG_fourier_sin.rds")




Xtwo_random = dat_random %>% select(all_of(ps2))
Xtwo_AT_random = matrix(nrow=nrow(Xtwo_random), ncol=ncol(Xtwo_random))
Xtwo_CG_random = matrix(nrow=nrow(Xtwo_random), ncol=ncol(Xtwo_random))
colnames(Xtwo_AT_random) = colnames(Xtwo_random)
colnames(Xtwo_CG_random) = colnames(Xtwo_random)
Xtwo_AT_random[] = ((Xtwo_random == "AA") | (Xtwo_random == "TT") | (Xtwo_random == "AT") | (Xtwo_random == "TA")) %>% 
  as.matrix() %>% as.numeric()
Xtwo_CG_random[] = ((Xtwo_random == "CC") | (Xtwo_random == "GG") | (Xtwo_random == "CG") | (Xtwo_random == "GC")) %>% 
  as.matrix() %>% as.numeric()

Xtwo_AT_random_fourier_cos = Xtwo_AT_random %*% t(cos_matrix2)
Xtwo_AT_random_fourier_sin = Xtwo_AT_random %*% t(sin_matrix2)

Xtwo_CG_random_fourier_cos = Xtwo_CG_random %*% t(cos_matrix2)
Xtwo_CG_random_fourier_sin = Xtwo_CG_random %*% t(sin_matrix2)

colnames(Xtwo_AT_random_fourier_cos) = paste0("fouriercosAT2_", 1:49)
colnames(Xtwo_AT_random_fourier_sin) = paste0("fouriersinAT2_", 1:49)

colnames(Xtwo_CG_random_fourier_cos) = paste0("fouriercosCG2_", 1:49)
colnames(Xtwo_CG_random_fourier_sin) = paste0("fouriersinCG2_", 1:49)

saveRDS(Xtwo_AT_random_fourier_cos, "data/Created/random_Xtwo_AT_fourier_cos.rds")
saveRDS(Xtwo_AT_random_fourier_sin, "data/Created/random_Xtwo_AT_fourier_sin.rds")

saveRDS(Xtwo_CG_random_fourier_cos, "data/Created/random_Xtwo_CG_fourier_cos.rds")
saveRDS(Xtwo_CG_random_fourier_sin, "data/Created/random_Xtwo_CG_fourier_sin.rds")






# Test data:
Xone_random_test = dat_random_test %>% select(all_of(ps1))
Xone_AT_random_test = matrix(nrow=nrow(Xone_random_test), ncol=ncol(Xone_random_test))
Xone_CG_random_test = matrix(nrow=nrow(Xone_random_test), ncol=ncol(Xone_random_test))
colnames(Xone_AT_random_test) = colnames(Xone_random_test)
colnames(Xone_CG_random_test) = colnames(Xone_random_test)
Xone_AT_random_test[] = ((Xone_random_test == "A") | (Xone_random_test == "T")) %>% as.matrix() %>% as.numeric()
Xone_CG_random_test[] = ((Xone_random_test == "C") | (Xone_random_test == "G")) %>% as.matrix() %>% as.numeric()

Xone_AT_random_fourier_cos_test = Xone_AT_random_test %*% t(cos_matrix)
Xone_AT_random_fourier_sin_test = Xone_AT_random_test %*% t(sin_matrix)

Xone_CG_random_fourier_cos_test = Xone_CG_random_test %*% t(cos_matrix)
Xone_CG_random_fourier_sin_test = Xone_CG_random_test %*% t(sin_matrix)

colnames(Xone_AT_random_fourier_cos_test) = paste0("fouriercosAT", 1:50)
colnames(Xone_AT_random_fourier_sin_test) = paste0("fouriersinAT", 1:50)

colnames(Xone_CG_random_fourier_cos_test) = paste0("fouriercosCG", 1:50)
colnames(Xone_CG_random_fourier_sin_test) = paste0("fouriersinCG", 1:50)

saveRDS(Xone_AT_random_fourier_cos_test, "data/Created/random_Xone_AT_fourier_cos_test.rds")
saveRDS(Xone_AT_random_fourier_sin_test, "data/Created/random_Xone_AT_fourier_sin_test.rds")

saveRDS(Xone_CG_random_fourier_cos_test, "data/Created/random_Xone_CG_fourier_cos_test.rds")
saveRDS(Xone_CG_random_fourier_sin_test, "data/Created/random_Xone_CG_fourier_sin_test.rds")



Xtwo_random_test = dat_random_test %>% select(all_of(ps2))
Xtwo_AT_random_test = matrix(nrow=nrow(Xtwo_random_test), ncol=ncol(Xtwo_random_test))
Xtwo_CG_random_test = matrix(nrow=nrow(Xtwo_random_test), ncol=ncol(Xtwo_random_test))
colnames(Xtwo_AT_random_test) = colnames(Xtwo_random_test)
colnames(Xtwo_CG_random_test) = colnames(Xtwo_random_test)
Xtwo_AT_random_test[] = ((Xtwo_random_test == "AA") | (Xtwo_random_test == "TT") | (Xtwo_random_test == "AT") | 
                           (Xtwo_random_test == "TA")) %>% as.matrix() %>% as.numeric()
Xtwo_CG_random_test[] = ((Xtwo_random_test == "CC") | (Xtwo_random_test == "GG") | (Xtwo_random_test == "CG") | 
                           (Xtwo_random_test == "GC")) %>% as.matrix() %>% as.numeric()

Xtwo_AT_random_fourier_cos_test = Xtwo_AT_random_test %*% t(cos_matrix2)
Xtwo_AT_random_fourier_sin_test = Xtwo_AT_random_test %*% t(sin_matrix2)

Xtwo_CG_random_fourier_cos_test = Xtwo_CG_random_test %*% t(cos_matrix2)
Xtwo_CG_random_fourier_sin_test = Xtwo_CG_random_test %*% t(sin_matrix2)

colnames(Xtwo_AT_random_fourier_cos_test) = paste0("fouriercosAT2_", 1:49)
colnames(Xtwo_AT_random_fourier_sin_test) = paste0("fouriersinAT2_", 1:49)

colnames(Xtwo_CG_random_fourier_cos_test) = paste0("fouriercosCG2_", 1:49)
colnames(Xtwo_CG_random_fourier_sin_test) = paste0("fouriersinCG2_", 1:49)

saveRDS(Xtwo_AT_random_fourier_cos_test, "data/Created/random_Xtwo_AT_fourier_cos_test.rds")
saveRDS(Xtwo_AT_random_fourier_sin_test, "data/Created/random_Xtwo_AT_fourier_sin_test.rds")

saveRDS(Xtwo_CG_random_fourier_cos_test, "data/Created/random_Xtwo_CG_fourier_cos_test.rds")
saveRDS(Xtwo_CG_random_fourier_sin_test, "data/Created/random_Xtwo_CG_fourier_sin_test.rds")






########################################################################
# Tiling Distances
########################################################################

cos_vector2_5 = cos_matrix2[5,]
sin_vector2_5 = sin_matrix2[5,]

# Xtwo_fourier_distance_onehot = matrix(nrow = nrow(Xtwo))

for (i in 1:15) {
  di1 = dinucleotides[i]
  cur1 = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
  colnames(cur1) = colnames(Xtwo)
  cur1[] = (Xtwo == di1) %>% 
    as.matrix() %>% as.numeric()
  
  # rows: n (frequency), cols: k (location)
  # period: N/n
  cur1_fourier_cos = cur1 %*% cos_vector2_5
  cur1_fourier_sin = cur1 %*% sin_vector2_5
  
  for (j in i+1:16) {
    di2 = dinucleotides[j]
    cur2 = matrix(nrow=nrow(Xtwo), ncol=ncol(Xtwo))
    colnames(cur2) = colnames(Xtwo)
    cur2[] = (Xtwo == di2) %>%
      as.matrix() %>% as.numeric()
    
    cur2_fourier_cos = cur2 %*% cos_vector2_5
    cur2_fourier_sin = cur2 %*% sin_vector2_5
    
    dist = sqrt((cur1_fourier_cos - cur2_fourier_cos)^2 + 
                  (cur1_fourier_sin - cur2_fourier_sin)^2)
    
    colnames(dist) = paste0("fourier2_distance_", di1, "_", di2, "_", 5)

    assign(paste0("Xtwo_fourier_distance_", di1, "_", di2), dist)
    
    # Xtwo_fourier_distance_onehot = cbind(Xtwo_fourier_distance_onehot, dist)
  }
}

saveRDS(Xtwo_fourier_distance_onehot, "data/Created/tiling_Xtwo_fourier_distance_onehot.rds")


# Test data:

Xtwo_fourier_distance_onehot_test = matrix(nrow = nrow(Xtwo_test))

for (i in 1:15) {
  di1 = dinucleotides[i]
  cur1 = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
  colnames(cur1) = colnames(Xtwo_test)
  cur1[] = (Xtwo_test == di1) %>% 
    as.matrix() %>% as.numeric()
  
  # rows: n (frequency), cols: k (location)
  # period: N/n
  cur1_fourier_cos = cur1 %*% t(cos_matrix2)
  cur1_fourier_sin = cur1 %*% t(sin_matrix2)
  
  for (j in i+1:16) {
    di2 = dinucleotides[j]
    cur2 = matrix(nrow=nrow(Xtwo_test), ncol=ncol(Xtwo_test))
    colnames(cur2) = colnames(Xtwo_test)
    cur2[] = (Xtwo_test == di2) %>%
      as.matrix() %>% as.numeric()
    
    cur2_fourier_cos = cur2 %*% t(cos_matrix2)
    cur2_fourier_sin = cur2 %*% t(sin_matrix2)
    
    dist = sqrt((cur1_fourier_cos - cur2_fourier_cos)^2 + 
                  (cur1_fourier_sin - cur2_fourier_sin)^2)
    
    colnames(dist) = paste0("fourier2_distance_", di1, "_", di2, "_", 1:49)

    assign(paste0("Xtwo_fourier_distance_", di1, "_", di2, "_test"), dist)
    
    Xtwo_fourier_distance_onehot_test = cbind(Xtwo_fourier_distance_onehot_test, dist)
  }
}

saveRDS(Xtwo_fourier_distance_onehot_test, "data/Created/tiling_Xtwo_fourier_distance_onehot_test.rds")


