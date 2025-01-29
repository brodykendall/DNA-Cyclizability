######## A/T:

#### Top 4000/25000:

# Random A/T A29/A26
Mod(find_ps1_fourier(C29, random %>% top_n(4000, C29)))/
  Mod(find_ps1_fourier(C26, random %>% top_n(4000, C26)))
# 1.216022
# Random A/T A31/A26
Mod(find_ps1_fourier(C31, random %>% top_n(4000, C31)))/
  Mod(find_ps1_fourier(C26, random %>% top_n(4000, C26)))
# 1.280732

# Tiling A/T A29/A26
Mod(find_ps1_fourier(C29, tiling %>% top_n(25000, C29)))/
  Mod(find_ps1_fourier(C26, tiling %>% top_n(25000, C26)))
# 1.249591
# Tiling A/T A31/A26
Mod(find_ps1_fourier(C31, tiling %>% top_n(25000, C31)))/
  Mod(find_ps1_fourier(C26, tiling %>% top_n(25000, C26)))
# 1.18748

# ChrV A/T A29/A26
Mod(find_ps1_fourier(C29, chrv %>% top_n(25000, C29)))/
  Mod(find_ps1_fourier(C26, chrv %>% top_n(25000, C26)))
# 1.312908
# ChrV A/T A31/A26
Mod(find_ps1_fourier(C31, chrv %>% top_n(25000, C31)))/
  Mod(find_ps1_fourier(C26, chrv %>% top_n(25000, C26)))
# 1.251674


#### Top 2000/15000:

# Random A/T A29/A26
Mod(find_ps1_fourier(C29, random %>% top_n(2000, C29)))/
  Mod(find_ps1_fourier(C26, random %>% top_n(2000, C26)))
# 1.166823
# Random A/T A31/A26
Mod(find_ps1_fourier(C31, random %>% top_n(2000, C31)))/
  Mod(find_ps1_fourier(C26, random %>% top_n(2000, C26)))
# 1.150479

# Tiling A/T A29/A26
Mod(find_ps1_fourier(C29, tiling %>% top_n(15000, C29)))/
  Mod(find_ps1_fourier(C26, tiling %>% top_n(15000, C26)))
# 1.21366
# Tiling A/T A31/A26
Mod(find_ps1_fourier(C31, tiling %>% top_n(15000, C31)))/
  Mod(find_ps1_fourier(C26, tiling %>% top_n(15000, C26)))
# 1.158194

# ChrV A/T A29/A26
Mod(find_ps1_fourier(C29, chrv %>% top_n(15000, C29)))/
  Mod(find_ps1_fourier(C26, chrv %>% top_n(15000, C26)))
# 1.241191
# ChrV A/T A31/A26
Mod(find_ps1_fourier(C31, chrv %>% top_n(15000, C31)))/
  Mod(find_ps1_fourier(C26, chrv %>% top_n(15000, C26)))
# 1.196313


######## AA/TT:

#### Top Quartiles:

# Random AA/TT A29/A26
Mod(fft(random_AA_TT_AT_TA_C29_q1)[6])/Mod(fft(random_AA_TT_AT_TA_C26_q1)[6])
# 1.421354
# Random AA/TT A31/A26
Mod(fft(random_AA_TT_AT_TA_C31_q1)[6])/Mod(fft(random_AA_TT_AT_TA_C26_q1)[6])
# 1.457906

# Tiling AA/TT A29/A26
Mod(fft(tiling_AA_TT_AT_TA_C29_q1)[6])/Mod(fft(tiling_AA_TT_AT_TA_C26_q1)[6])
# 1.354385
# Tiling AA/TT A31/A26
Mod(fft(tiling_AA_TT_AT_TA_C31_q1)[6])/Mod(fft(tiling_AA_TT_AT_TA_C26_q1)[6])
# 1.250266

# ChrV AA/TT A29/A26
Mod(fft(chrv_AA_TT_AT_TA_C29_q1)[6])/Mod(fft(chrv_AA_TT_AT_TA_C26_q1)[6])
# 1.41639
# ChrV AA/TT A31/A26
Mod(fft(chrv_AA_TT_AT_TA_C31_q1)[6])/Mod(fft(chrv_AA_TT_AT_TA_C26_q1)[6])
# 1.313074


#### Top 70%:

# Random AA/TT A29/A26
Mod(fft(random_AA_TT_AT_TA_C29_custom1)[6])/Mod(fft(random_AA_TT_AT_TA_C26_custom1)[6])
# 1.236593
# Random AA/TT A31/A26
Mod(fft(random_AA_TT_AT_TA_C31_custom1)[6])/Mod(fft(random_AA_TT_AT_TA_C26_custom1)[6])
# 1.201936

# Tiling AA/TT A29/A26
Mod(fft(tiling_AA_TT_AT_TA_C29_custom1)[6])/Mod(fft(tiling_AA_TT_AT_TA_C26_custom1)[6])
# 1.274177
# Tiling AA/TT A31/A26
Mod(fft(tiling_AA_TT_AT_TA_C31_custom1)[6])/Mod(fft(tiling_AA_TT_AT_TA_C26_custom1)[6])
# 1.185538

# ChrV AA/TT A29/A26
Mod(fft(chrv_AA_TT_AT_TA_C29_custom1)[6])/Mod(fft(chrv_AA_TT_AT_TA_C26_custom1)[6])
# 1.355702
# ChrV AA/TT A31/A26
Mod(fft(chrv_AA_TT_AT_TA_C31_custom1)[6])/Mod(fft(chrv_AA_TT_AT_TA_C26_custom1)[6])
# 1.26201


