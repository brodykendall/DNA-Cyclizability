ps2 <- paste0("X", 1:49, "di")

# Fourier for each sequence individually:
find_all_ps2_fourier = function(data) {
  data = data %>% select(all_of(ps2))
  Xtwo_AorT_data = matrix(nrow=nrow(data), ncol=ncol(data))
  colnames(Xtwo_AorT_data) = colnames(data)
  Xtwo_AorT_data[] = ((data == "AA") | (data == "TT") | (data == "AT") | (data == "TA")) %>% 
    as.matrix() %>% as.numeric()
  
  temp_fourier = t(apply(Xtwo_AorT_data, 1, fft))[,6]
  return(temp_fourier)
}

all_ps2_fourier_random = find_all_ps2_fourier(random)
all_ps2_fourier_tiling = find_all_ps2_fourier(tiling)
all_ps2_fourier_chrv = find_all_ps2_fourier(chrv)

# Variance:
Cn_var_ps2_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})
Cn_var_ps2_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})
Cn_var_ps2_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})

cor(Mod(all_ps2_fourier_random), Cn_var_ps2_random)
# 0.2216326
plot(Mod(all_ps2_fourier_random), Cn_var_ps2_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AATT, random")

cor(Mod(all_ps2_fourier_tiling), Cn_var_ps2_tiling)
# 0.176731
plot(Mod(all_ps2_fourier_tiling), Cn_var_ps2_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AATT, tiling")
plot(Mod(all_ps2_fourier_tiling), Cn_var_ps2_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AATT, tiling")

cor(Mod(all_ps2_fourier_chrv), Cn_var_ps2_chrv)
# 0.06031859
plot(Mod(all_ps2_fourier_chrv), Cn_var_ps2_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AATT, chrv")
plot(Mod(all_ps2_fourier_chrv), Cn_var_ps2_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AATT, chrv")


# Range:
Cn_range_ps2_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})
Cn_range_ps2_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})
Cn_range_ps2_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})

cor(Mod(all_ps2_fourier_random), Cn_range_ps2_random)
# 0.2226965
plot(Mod(all_ps2_fourier_random), Cn_range_ps2_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AATT, random")

cor(Mod(all_ps2_fourier_tiling), Cn_range_ps2_tiling)
# 0.20332
plot(Mod(all_ps2_fourier_tiling), Cn_range_ps2_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AATT, tiling")
plot(Mod(all_ps2_fourier_tiling), Cn_range_ps2_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AATT, tiling")

cor(Mod(all_ps2_fourier_chrv), Cn_range_ps2_chrv)
# 0.08597398
plot(Mod(all_ps2_fourier_chrv), Cn_range_ps2_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AATT, chrv")
plot(Mod(all_ps2_fourier_chrv), Cn_range_ps2_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AATT, chrv")


# Mean:
# Cn_mean_ps2_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
#   return(mean(row))
# })
# Cn_mean_ps2_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
#   return(mean(row))
# })
# Cn_mean_ps2_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
#   return(mean(row))
# })

Cn_mean_ps2_random = apply(random %>% select(C26, C31), 1, function(row) {
  return(mean(row))
})
Cn_mean_ps2_tiling = apply(tiling %>% select(C26, C31), 1, function(row) {
  return(mean(row))
})
Cn_mean_ps2_chrv = apply(chrv %>% select(C26, C31), 1, function(row) {
  return(mean(row))
})

cor(Mod(all_ps2_fourier_random), Cn_mean_ps2_random)
# 0.2630076
plot(Mod(all_ps2_fourier_random), Cn_mean_ps2_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AATT, random")

cor(Mod(all_ps2_fourier_tiling), Cn_mean_ps2_tiling)
# 0.3262742
plot(Mod(all_ps2_fourier_tiling), Cn_mean_ps2_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AATT, tiling")
plot(Mod(all_ps2_fourier_tiling), Cn_mean_ps2_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AATT, tiling")

cor(Mod(all_ps2_fourier_chrv), Cn_mean_ps2_chrv)
# 0.2840405
plot(Mod(all_ps2_fourier_chrv), Cn_mean_ps2_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AATT, chrv")
plot(Mod(all_ps2_fourier_chrv), Cn_mean_ps2_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AATT, chrv")



# CV
plot(Mod(all_ps2_fourier_random), sqrt(Cn_var_ps2_random)/abs(Cn_mean_ps2_random), col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 CV", ylim=c(0, 100))
title("AATT, random")
plot(Mod(all_ps2_fourier_tiling), sqrt(Cn_var_ps2_tiling)/abs(Cn_mean_ps2_tiling), col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 CV", ylim=c(0, 100))
title("AATT, tiling")
plot(Mod(all_ps2_fourier_chrv), sqrt(Cn_var_ps2_chrv)/abs(Cn_mean_ps2_chrv), col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 CV", ylim=c(-100, 100))
title("AATT, chrv")





# Provided Amplitude
cor(Mod(all_ps2_fourier_random), abs(random$Amplitude))
# 0.2102925
plot(Mod(all_ps2_fourier_random), abs(random$Amplitude), col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AATT, random")

cor(Mod(all_ps2_fourier_tiling), abs(tiling$Amplitude))
# 0.2008195
plot(Mod(all_ps2_fourier_tiling), abs(tiling$Amplitude), col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AATT, tiling")
plot(Mod(all_ps2_fourier_tiling), abs(tiling$Amplitude), col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AATT, tiling")

cor(Mod(all_ps2_fourier_chrv), abs(chrv$Amplitude))
# 0.06706283
plot(Mod(all_ps2_fourier_chrv), abs(chrv$Amplitude), col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AATT, chrv")
plot(Mod(all_ps2_fourier_chrv), abs(chrv$Amplitude), col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AATT, chrv")





