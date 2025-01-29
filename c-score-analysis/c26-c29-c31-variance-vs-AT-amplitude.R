ps1 <- paste0("X", 1:50, "mono")

# Fourier for each sequence individually:
find_all_ps1_fourier = function(data) {
  data = data %>% select(all_of(ps1))
  Xone_AorT_data = matrix(nrow=nrow(data), ncol=ncol(data))
  colnames(Xone_AorT_data) = colnames(data)
  Xone_AorT_data[] = ((data == "A") | (data == "T")) %>% 
    as.matrix() %>% as.numeric()
  
  temp_fourier = t(apply(Xone_AorT_data, 1, fft))[,6]
  return(temp_fourier)
}

all_ps1_fourier_random = find_all_ps1_fourier(random)
all_ps1_fourier_tiling = find_all_ps1_fourier(tiling)
all_ps1_fourier_chrv = find_all_ps1_fourier(chrv)


# Variance:
Cn_var_ps1_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})
Cn_var_ps1_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})
Cn_var_ps1_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
  return(var(row))
})

cor(Mod(all_ps1_fourier_random), Cn_var_ps1_random)
# 0.2496933
plot(Mod(all_ps1_fourier_random), Cn_var_ps1_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AT, random")

cor(Mod(all_ps1_fourier_tiling), Cn_var_ps1_tiling)
# 0.1743162
plot(Mod(all_ps1_fourier_tiling), Cn_var_ps1_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AT, tiling")
plot(Mod(all_ps1_fourier_tiling), Cn_var_ps1_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AT, tiling")

cor(Mod(all_ps1_fourier_chrv), Cn_var_ps1_chrv)
# 0.09877792
plot(Mod(all_ps1_fourier_chrv), Cn_var_ps1_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AT, chrv")
plot(Mod(all_ps1_fourier_chrv), Cn_var_ps1_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Variance")
title("AT, chrv")


# Range:
Cn_range_ps1_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})
Cn_range_ps1_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})
Cn_range_ps1_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
  return(diff(range(row)))
})

cor(Mod(all_ps1_fourier_random), Cn_range_ps1_random)
# 0.2537833
plot(Mod(all_ps1_fourier_random), Cn_range_ps1_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AT, random")

cor(Mod(all_ps1_fourier_tiling), Cn_range_ps1_tiling)
# 0.2077106
plot(Mod(all_ps1_fourier_tiling), Cn_range_ps1_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AT, tiling")
plot(Mod(all_ps1_fourier_tiling), Cn_range_ps1_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AT, tiling")

cor(Mod(all_ps1_fourier_chrv), Cn_range_ps1_chrv)
# 0.1263523
plot(Mod(all_ps1_fourier_chrv), Cn_range_ps1_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AT, chrv")
plot(Mod(all_ps1_fourier_chrv), Cn_range_ps1_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Range")
title("AT, chrv")


# Mean:
Cn_mean_ps1_random = apply(random %>% select(C26, C29, C31), 1, function(row) {
  return(mean(row))
})
Cn_mean_ps1_tiling = apply(tiling %>% select(C26, C29, C31), 1, function(row) {
  return(mean(row))
})
Cn_mean_ps1_chrv = apply(chrv %>% select(C26, C29, C31), 1, function(row) {
  return(mean(row))
})

cor(Mod(all_ps1_fourier_random), Cn_mean_ps1_random)
# 0.3157814
plot(Mod(all_ps1_fourier_random), Cn_mean_ps1_random, col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AT, random")

cor(Mod(all_ps1_fourier_tiling), Cn_mean_ps1_tiling)
# 0.3870065
plot(Mod(all_ps1_fourier_tiling), Cn_mean_ps1_tiling, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AT, tiling")
plot(Mod(all_ps1_fourier_tiling), Cn_mean_ps1_tiling, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AT, tiling")

cor(Mod(all_ps1_fourier_chrv), Cn_mean_ps1_chrv)
# 0.3433733
plot(Mod(all_ps1_fourier_chrv), Cn_mean_ps1_chrv, col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AT, chrv")
plot(Mod(all_ps1_fourier_chrv), Cn_mean_ps1_chrv, col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="C26, C29, C31 Mean")
title("AT, chrv")





# Provided Amplitude
cor(Mod(all_ps1_fourier_random), abs(random$Amplitude))
# 0.2380436
plot(Mod(all_ps1_fourier_random), abs(random$Amplitude), col=alpha("black", 0.1),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AT, random")

cor(Mod(all_ps1_fourier_tiling), abs(tiling$Amplitude))
# 0.2009403
plot(Mod(all_ps1_fourier_tiling), abs(tiling$Amplitude), col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AT, tiling")
plot(Mod(all_ps1_fourier_tiling), abs(tiling$Amplitude), col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AT, tiling")

cor(Mod(all_ps1_fourier_chrv), abs(chrv$Amplitude))
# 0.1087601
plot(Mod(all_ps1_fourier_chrv), abs(chrv$Amplitude), col=alpha("black", 0.025),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AT, chrv")
plot(Mod(all_ps1_fourier_chrv), abs(chrv$Amplitude), col=alpha("black", 0.5),
     xlab="Fourier Amplitude", ylab="Original Phi (Amplitude)")
title("AT, chrv")



