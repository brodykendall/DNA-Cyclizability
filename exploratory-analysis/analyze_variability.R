# dat_nuc = read.csv("cycle1.txt")
# dat_random = read.csv("cycle3.txt")
# dat_tiling = read.csv("cycle5.txt")
# dat_chrv = read.csv("cycle6.txt")
# 
# colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
# colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
# colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
# colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")


dat_nuc = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_cycle1_best_fold_predictions.csv")
dat_random = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_cycle3_best_fold_predictions.csv")
dat_tiling = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_cycle5_best_fold_predictions.csv")
dat_chrv = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_cycle6_best_fold_predictions.csv")

colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "SmoothC0_pred")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "SmoothC0_pred")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "SmoothC0_pred")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "SmoothC0_pred")

dat_chrV = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11_best_fold_predictions.csv")
dat_chrV = dat_chrV %>%
  rename(x50mer = sequence,
         SmoothC0 = smoothC0,
         SmoothC0_predictions = smoothC0_predictions)
dat_chrV_cn_avg_pred = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_predictions.csv")
dat_chrV_cn_avg_pred = dat_chrV_cn_avg_pred %>%
  rename(C26 = n.26,
         C29 = n.29,
         C31 = n.31)

dat_chrV = cbind(dat_chrV, dat_chrV_cn_avg_pred %>% select(C26, C29, C31))

normal_mean = -0.011196041799376931
normal_std = 0.651684644408004

dat_nuc$SmoothC0_pred_unnorm = dat_nuc$SmoothC0_pred*normal_std + normal_mean
dat_random$SmoothC0_pred_unnorm = dat_random$SmoothC0_pred*normal_std + normal_mean
dat_tiling$SmoothC0_pred_unnorm = dat_tiling$SmoothC0_pred*normal_std + normal_mean
dat_chrv$SmoothC0_pred_unnorm = dat_chrv$SmoothC0_pred*normal_std + normal_mean

dat_chrV$SmoothC0_pred_unnorm = dat_chrV$SmoothC0_pred*normal_std + normal_mean

dat_nuc$C26_norm = (dat_nuc$C26 - mean(dat_nuc$C26))/sqrt(var(dat_nuc$C26))
dat_nuc$C29_norm = (dat_nuc$C29 - mean(dat_nuc$C29))/sqrt(var(dat_nuc$C29))
dat_nuc$C31_norm = (dat_nuc$C31 - mean(dat_nuc$C31))/sqrt(var(dat_nuc$C31))
dat_nuc$C0_norm = (dat_nuc$C0 - mean(dat_nuc$C0))/sqrt(var(dat_nuc$C0))

dat_random$C26_norm = (dat_random$C26 - mean(dat_random$C26))/sqrt(var(dat_random$C26))
dat_random$C29_norm = (dat_random$C29 - mean(dat_random$C29))/sqrt(var(dat_random$C29))
dat_random$C31_norm = (dat_random$C31 - mean(dat_random$C31))/sqrt(var(dat_random$C31))
dat_random$C0_norm = (dat_random$C0 - mean(dat_random$C0))/sqrt(var(dat_random$C0))

dat_tiling$C26_norm = (dat_tiling$C26 - mean(dat_tiling$C26))/sqrt(var(dat_tiling$C26))
dat_tiling$C29_norm = (dat_tiling$C29 - mean(dat_tiling$C29))/sqrt(var(dat_tiling$C29))
dat_tiling$C31_norm = (dat_tiling$C31 - mean(dat_tiling$C31))/sqrt(var(dat_tiling$C31))
dat_tiling$C0_norm = (dat_tiling$C0 - mean(dat_tiling$C0))/sqrt(var(dat_tiling$C0))

dat_chrv$C26_norm = (dat_chrv$C26 - mean(dat_chrv$C26))/sqrt(var(dat_chrv$C26))
dat_chrv$C29_norm = (dat_chrv$C29 - mean(dat_chrv$C29))/sqrt(var(dat_chrv$C29))
dat_chrv$C31_norm = (dat_chrv$C31 - mean(dat_chrv$C31))/sqrt(var(dat_chrv$C31))
dat_chrv$C0_norm = (dat_chrv$C0 - mean(dat_chrv$C0))/sqrt(var(dat_chrv$C0))

# # Nucleosome:
# dat_nuc = dat_nuc %>%
#   mutate(periodic_term_26 = Amplitude*sin(2*pi*26/10.4 + Phase - pi),
#          periodic_term_29 = Amplitude/0.82*sin(2*pi*29/10.4 + Phase - pi),
#          periodic_term_31 = Amplitude/0.7*sin(2*pi*31/10.4 + Phase - pi))
# 
# # Variance due to periodic term:
# var(dat_nuc$periodic_term_26)
# # 0.09489595
# var(dat_nuc$periodic_term_29)
# # 0.2598312
# var(dat_nuc$periodic_term_31)
# # 0.191203
# 
# # Total variance in each Cn:
# var(dat_nuc$C26)
# # 0.3693024
# var(dat_nuc$C29)
# # 0.4046268
# var(dat_nuc$C31)
# # 0.3281334
# 
# 
# # Percentage variance contribution:
# var(dat_nuc$periodic_term_26)/var(dat_nuc$C26)
# # 0.25696
# var(dat_nuc$periodic_term_29)/var(dat_nuc$C29)
# # 0.6421503
# var(dat_nuc$periodic_term_31)/var(dat_nuc$C31)
# # 0.5826988
# 
# 
# # Random:
# dat_random = dat_random %>%
#   mutate(periodic_term_26 = Amplitude*sin(2*pi*26/10.4 + Phase - pi),
#          periodic_term_29 = Amplitude/0.82*sin(2*pi*29/10.4 + Phase - pi),
#          periodic_term_31 = Amplitude/0.7*sin(2*pi*31/10.4 + Phase - pi))
# 
# # Variance due to periodic term:
# var(dat_random$periodic_term_26)
# # 0.07112159
# var(dat_random$periodic_term_29)
# # 0.2189252
# var(dat_random$periodic_term_31)
# # 0.1551205
# 
# # Total variance in each Cn:
# var(dat_random$C26)
# # 0.2251218
# var(dat_random$C29)
# # 0.3170235
# var(dat_random$C31)
# # 0.2757755
# 
# 
# # Percentage variance contribution:
# var(dat_random$periodic_term_26)/var(dat_random$C26)
# # 0.3159249
# var(dat_random$periodic_term_29)/var(dat_random$C29)
# # 0.6905645
# var(dat_random$periodic_term_31)/var(dat_random$C31)
# # 0.5624882
# 
# 
# 
# # Tiling
# dat_tiling = dat_tiling %>%
#   mutate(periodic_term_26 = Amplitude*sin(2*pi*26/10.4 + Phase - pi),
#          periodic_term_29 = Amplitude/0.82*sin(2*pi*29/10.4 + Phase - pi),
#          periodic_term_31 = Amplitude/0.7*sin(2*pi*31/10.4 + Phase - pi))
# 
# # Variance due to periodic term:
# var(dat_tiling$periodic_term_26)
# # 0.09509013
# var(dat_tiling$periodic_term_29)
# # 0.2061907
# var(dat_tiling$periodic_term_31)
# # 0.200384
# 
# # Total variance in each Cn:
# var(dat_tiling$C26)
# # 0.331956
# var(dat_tiling$C29)
# # 0.3468428
# var(dat_tiling$C31)
# # 0.4259002
# 
# 
# # Percentage variance contribution:
# var(dat_tiling$periodic_term_26)/var(dat_tiling$C26)
# # 0.2864541
# var(dat_tiling$periodic_term_29)/var(dat_tiling$C29)
# # 0.594479
# var(dat_tiling$periodic_term_31)/var(dat_tiling$C31)
# # 0.4704951
# 
# 
# 
# # ChrV:
# dat_chrv = dat_chrv %>%
#   mutate(periodic_term_26 = Amplitude*sin(2*pi*26/10.4 + Phase - pi),
#          periodic_term_29 = Amplitude/0.82*sin(2*pi*29/10.4 + Phase - pi),
#          periodic_term_31 = Amplitude/0.7*sin(2*pi*31/10.4 + Phase - pi))
# 
# # Variance due to periodic term:
# var(dat_chrv$periodic_term_26)
# # 0.1343606
# var(dat_chrv$periodic_term_29)
# # 0.4560969
# var(dat_chrv$periodic_term_31)
# # 0.2807901
# 
# # Total variance in each Cn:
# var(dat_chrv$C26)
# # 0.4674042
# var(dat_chrv$C29)
# # 0.5226066
# var(dat_chrv$C31)
# # 0.4896776
# 
# 
# # Percentage variance contribution:
# var(dat_chrv$periodic_term_26)/var(dat_chrv$C26)
# # 0.2874612
# var(dat_chrv$periodic_term_29)/var(dat_chrv$C29)
# # 0.8727347
# var(dat_chrv$periodic_term_31)/var(dat_chrv$C31)
# # 0.5734183
# 
# 
# 
# 
# 
# 
# 
# 

# ACTUALLY USE THIS:

# N-SPECIFIC ANALYSIS:

dat_nuc_Cbar = (mean(dat_nuc$C26) + mean(dat_nuc$C29) + mean(dat_nuc$C31))/3
dat_nuc_C26bar = mean(dat_nuc$C26)
dat_nuc_C29bar = mean(dat_nuc$C29)
dat_nuc_C31bar = mean(dat_nuc$C31)

nuc_V_C26_in = sum((dat_nuc$C26 - dat_nuc_C26bar)^2)
nuc_V_C26_periodic = sum((dat_nuc$C26 - dat_nuc$C0)^2)
nuc_V_C26_cross_terms = sum((dat_nuc$C26 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_C26bar))

2*nuc_V_C26_cross_terms/nuc_V_C26_in
# 0.1225938

nuc_V_C26_periodic/nuc_V_C26_in
# 0.2676762


nuc_V_C29_in = sum((dat_nuc$C29 - dat_nuc_C29bar)^2)
nuc_V_C29_periodic = sum((dat_nuc$C29 - dat_nuc$C0)^2)
nuc_V_C29_cross_terms = sum((dat_nuc$C29 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_C29bar))

2*nuc_V_C29_cross_terms/nuc_V_C29_in
# -0.1908883

nuc_V_C29_periodic/nuc_V_C29_in
# 0.6431597


nuc_V_C31_in = sum((dat_nuc$C31 - dat_nuc_C31bar)^2)
nuc_V_C31_periodic = sum((dat_nuc$C31 - dat_nuc$C0)^2)
nuc_V_C31_cross_terms = sum((dat_nuc$C31 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_C31bar))

2*nuc_V_C31_cross_terms/nuc_V_C31_in
# -0.3056724

nuc_V_C31_periodic/nuc_V_C31_in
# 0.6071014


# Predicted Smooth C0:
dat_nuc_Cbar_norm = (mean(dat_nuc$C26_norm) + mean(dat_nuc$C29_norm) + mean(dat_nuc$C31_norm))/3
dat_nuc_C26bar_norm = mean(dat_nuc$C26_norm)
dat_nuc_C29bar_norm = mean(dat_nuc$C29_norm)
dat_nuc_C31bar_norm = mean(dat_nuc$C31_norm)

nuc_V_C26_norm_in = sum((dat_nuc$C26_norm - dat_nuc_C26bar_norm)^2)
nuc_V_C26_norm_periodic = sum((dat_nuc$C26_norm - dat_nuc$SmoothC0_pred_unnorm)^2)
nuc_V_C26_norm_cross_terms = sum((dat_nuc$C26_norm - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_C26bar_norm))

2*nuc_V_C26_norm_cross_terms/nuc_V_C26_norm_in
# 0.2341831

nuc_V_C26_norm_periodic/nuc_V_C26_norm_in
# 0.3992378


nuc_V_C29_norm_in = sum((dat_nuc$C29_norm - dat_nuc_C29bar_norm)^2)
nuc_V_C29_norm_periodic = sum((dat_nuc$C29_norm - dat_nuc$SmoothC0_pred_unnorm)^2)
nuc_V_C29_norm_cross_terms = sum((dat_nuc$C29_norm - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_C29bar_norm))

2*nuc_V_C29_norm_cross_terms/nuc_V_C29_norm_in
# 0.02703025

nuc_V_C29_norm_periodic/nuc_V_C29_norm_in
# 0.6063906


nuc_V_C31_norm_in = sum((dat_nuc$C31_norm - dat_nuc_C31bar_norm)^2)
nuc_V_C31_norm_periodic = sum((dat_nuc$C31_norm - dat_nuc$SmoothC0_pred_unnorm)^2)
nuc_V_C31_norm_cross_terms = sum((dat_nuc$C31_norm - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_C31bar_norm))

2*nuc_V_C31_norm_cross_terms/nuc_V_C31_norm_in
# -0.09986392

nuc_V_C31_norm_periodic/nuc_V_C31_norm_in
# 0.7332848








dat_random_Cbar = (mean(dat_random$C26) + mean(dat_random$C29) + mean(dat_random$C31))/3
dat_random_C26bar = mean(dat_random$C26)
dat_random_C29bar = mean(dat_random$C29)
dat_random_C31bar = mean(dat_random$C31)

random_V_C26_in = sum((dat_random$C26 - dat_random_C26bar)^2)
random_V_C26_periodic = sum((dat_random$C26 - dat_random$C0)^2)
random_V_C26_cross_terms = sum((dat_random$C26 - dat_random$C0)*(dat_random$C0 - dat_random_C26bar))



random_V_C26_in = sum((dat_random$C26 - dat_random_C26bar)^2)
random_V_C26_periodic = sum((dat_random$C26 - dat_random$C0)^2)
random_V_C26_cross_terms = sum((dat_random$C26 - dat_random$C0)*(dat_random$C0 - dat_random_C26bar))

2*random_V_C26_cross_terms/random_V_C26_in
# 0.04834965

random_V_C26_periodic/random_V_C26_in
# 0.3168609


random_V_C29_in = sum((dat_random$C29 - dat_random_C29bar)^2)
random_V_C29_periodic = sum((dat_random$C29 - dat_random$C0)^2)
random_V_C29_cross_terms = sum((dat_random$C29 - dat_random$C0)*(dat_random$C0 - dat_random_C29bar))

2*random_V_C29_cross_terms/random_V_C29_in
# -0.1446832

random_V_C29_periodic/random_V_C29_in
# 0.6925708


random_V_C31_in = sum((dat_random$C31 - dat_random_C31bar)^2)
random_V_C31_periodic = sum((dat_random$C31 - dat_random$C0)^2)
random_V_C31_cross_terms = sum((dat_random$C31 - dat_random$C0)*(dat_random$C0 - dat_random_C31bar))

2*random_V_C31_cross_terms/random_V_C31_in
# -0.08396948

random_V_C31_periodic/random_V_C31_in
# 0.5645144


# Predicted Smooth C0:
dat_random_Cbar_norm = (mean(dat_random$C26_norm) + mean(dat_random$C29_norm) + mean(dat_random$C31_norm))/3
dat_random_C26bar_norm = mean(dat_random$C26_norm)
dat_random_C29bar_norm = mean(dat_random$C29_norm)
dat_random_C31bar_norm = mean(dat_random$C31_norm)

random_V_C26_norm_in = sum((dat_random$C26_norm - dat_random_C26bar_norm)^2)
random_V_C26_norm_periodic = sum((dat_random$C26_norm - dat_random$SmoothC0_pred_unnorm)^2)
random_V_C26_norm_cross_terms = sum((dat_random$C26_norm - dat_random$SmoothC0_pred_unnorm)*(dat_random$SmoothC0_pred_unnorm - dat_random_C26bar_norm))

2*random_V_C26_norm_cross_terms/random_V_C26_norm_in
# 0.2236042

random_V_C26_norm_periodic/random_V_C26_norm_in
# 0.5436577


random_V_C29_norm_in = sum((dat_random$C29_norm - dat_random_C29bar_norm)^2)
random_V_C29_norm_periodic = sum((dat_random$C29_norm - dat_random$SmoothC0_pred_unnorm)^2)
random_V_C29_norm_cross_terms = sum((dat_random$C29_norm - dat_random$SmoothC0_pred_unnorm)*(dat_random$SmoothC0_pred_unnorm - dat_random_C29bar_norm))

2*random_V_C29_norm_cross_terms/random_V_C29_norm_in
# 0.09448441

random_V_C29_norm_periodic/random_V_C29_norm_in
# 0.6727775


random_V_C31_norm_in = sum((dat_random$C31_norm - dat_random_C31bar_norm)^2)
random_V_C31_norm_periodic = sum((dat_random$C31_norm - dat_random$SmoothC0_pred_unnorm)^2)
random_V_C31_norm_cross_terms = sum((dat_random$C31_norm - dat_random$SmoothC0_pred_unnorm)*(dat_random$SmoothC0_pred_unnorm - dat_random_C31bar_norm))

2*random_V_C31_norm_cross_terms/random_V_C31_norm_in
# 0.1455555

random_V_C31_norm_periodic/random_V_C31_norm_in
# 0.6217064






dat_tiling_Cbar = (mean(dat_tiling$C26) + mean(dat_tiling$C29) + mean(dat_tiling$C31))/3
dat_tiling_C26bar = mean(dat_tiling$C26)
dat_tiling_C29bar = mean(dat_tiling$C29)
dat_tiling_C31bar = mean(dat_tiling$C31)

tiling_V_C26_in = sum((dat_tiling$C26 - dat_tiling_C26bar)^2)
tiling_V_C26_periodic = sum((dat_tiling$C26 - dat_tiling$C0)^2)
tiling_V_C26_cross_terms = sum((dat_tiling$C26 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_C26bar))



tiling_V_C26_in = sum((dat_tiling$C26 - dat_tiling_C26bar)^2)
tiling_V_C26_periodic = sum((dat_tiling$C26 - dat_tiling$C0)^2)
tiling_V_C26_cross_terms = sum((dat_tiling$C26 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_C26bar))

2*tiling_V_C26_cross_terms/tiling_V_C26_in
# -0.005018646

tiling_V_C26_periodic/tiling_V_C26_in
# 0.2871789


tiling_V_C29_in = sum((dat_tiling$C29 - dat_tiling_C29bar)^2)
tiling_V_C29_periodic = sum((dat_tiling$C29 - dat_tiling$C0)^2)
tiling_V_C29_cross_terms = sum((dat_tiling$C29 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_C29bar))

2*tiling_V_C29_cross_terms/tiling_V_C29_in
# -0.2817492

tiling_V_C29_periodic/tiling_V_C29_in
# 0.5949462


tiling_V_C31_in = sum((dat_tiling$C31 - dat_tiling_C31bar)^2)
tiling_V_C31_periodic = sum((dat_tiling$C31 - dat_tiling$C0)^2)
tiling_V_C31_cross_terms = sum((dat_tiling$C31 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_C31bar))

2*tiling_V_C31_cross_terms/tiling_V_C31_in
# -0.03121361

tiling_V_C31_periodic/tiling_V_C31_in
# 0.4713868


# Predicted Smooth C0:
dat_tiling_Cbar_norm = (mean(dat_tiling$C26_norm) + mean(dat_tiling$C29_norm) + mean(dat_tiling$C31_norm))/3
dat_tiling_C26bar_norm = mean(dat_tiling$C26_norm)
dat_tiling_C29bar_norm = mean(dat_tiling$C29_norm)
dat_tiling_C31bar_norm = mean(dat_tiling$C31_norm)

tiling_V_C26_norm_in = sum((dat_tiling$C26_norm - dat_tiling_C26bar_norm)^2)
tiling_V_C26_norm_periodic = sum((dat_tiling$C26_norm - dat_tiling$SmoothC0_pred_unnorm)^2)
tiling_V_C26_norm_cross_terms = sum((dat_tiling$C26_norm - dat_tiling$SmoothC0_pred_unnorm)*(dat_tiling$SmoothC0_pred_unnorm - dat_tiling_C26bar_norm))

2*tiling_V_C26_norm_cross_terms/tiling_V_C26_norm_in
# 0.1078958

tiling_V_C26_norm_periodic/tiling_V_C26_norm_in
# 0.4796978


tiling_V_C29_norm_in = sum((dat_tiling$C29_norm - dat_tiling_C29bar_norm)^2)
tiling_V_C29_norm_periodic = sum((dat_tiling$C29_norm - dat_tiling$SmoothC0_pred_unnorm)^2)
tiling_V_C29_norm_cross_terms = sum((dat_tiling$C29_norm - dat_tiling$SmoothC0_pred_unnorm)*(dat_tiling$SmoothC0_pred_unnorm - dat_tiling_C29bar_norm))

2*tiling_V_C29_norm_cross_terms/tiling_V_C29_norm_in
# 0.01956878

tiling_V_C29_norm_periodic/tiling_V_C29_norm_in
# 0.5680248


tiling_V_C31_norm_in = sum((dat_tiling$C31_norm - dat_tiling_C31bar_norm)^2)
tiling_V_C31_norm_periodic = sum((dat_tiling$C31_norm - dat_tiling$SmoothC0_pred_unnorm)^2)
tiling_V_C31_norm_cross_terms = sum((dat_tiling$C31_norm - dat_tiling$SmoothC0_pred_unnorm)*(dat_tiling$SmoothC0_pred_unnorm - dat_tiling_C31bar_norm))

2*tiling_V_C31_norm_cross_terms/tiling_V_C31_norm_in
# 0.05894481

tiling_V_C31_norm_periodic/tiling_V_C31_norm_in
# 0.5286487







dat_chrv_Cbar = (mean(dat_chrv$C26) + mean(dat_chrv$C29) + mean(dat_chrv$C31))/3
dat_chrv_C26bar = mean(dat_chrv$C26)
dat_chrv_C29bar = mean(dat_chrv$C29)
dat_chrv_C31bar = mean(dat_chrv$C31)

chrv_V_C26_in = sum((dat_chrv$C26 - dat_chrv_C26bar)^2)
chrv_V_C26_periodic = sum((dat_chrv$C26 - dat_chrv$C0)^2)
chrv_V_C26_cross_terms = sum((dat_chrv$C26 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_C26bar))



chrv_V_C26_in = sum((dat_chrv$C26 - dat_chrv_C26bar)^2)
chrv_V_C26_periodic = sum((dat_chrv$C26 - dat_chrv$C0)^2)
chrv_V_C26_cross_terms = sum((dat_chrv$C26 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_C26bar))

2*chrv_V_C26_cross_terms/chrv_V_C26_in
# 0.08165593

chrv_V_C26_periodic/chrv_V_C26_in
# 0.2875924


chrv_V_C29_in = sum((dat_chrv$C29 - dat_chrv_C29bar)^2)
chrv_V_C29_periodic = sum((dat_chrv$C29 - dat_chrv$C0)^2)
chrv_V_C29_cross_terms = sum((dat_chrv$C29 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_C29bar))

2*chrv_V_C29_cross_terms/chrv_V_C29_in
# -0.4381992

chrv_V_C29_periodic/chrv_V_C29_in
# 0.8734626


chrv_V_C31_in = sum((dat_chrv$C31 - dat_chrv_C31bar)^2)
chrv_V_C31_periodic = sum((dat_chrv$C31 - dat_chrv$C0)^2)
chrv_V_C31_cross_terms = sum((dat_chrv$C31 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_C31bar))

2*chrv_V_C31_cross_terms/chrv_V_C31_in
# -0.1761106

chrv_V_C31_periodic/chrv_V_C31_in
# 0.5737963


# Predicted Smooth C0:
dat_chrv_Cbar_norm = (mean(dat_chrv$C26_norm) + mean(dat_chrv$C29_norm) + mean(dat_chrv$C31_norm))/3
dat_chrv_C26bar_norm = mean(dat_chrv$C26_norm)
dat_chrv_C29bar_norm = mean(dat_chrv$C29_norm)
dat_chrv_C31bar_norm = mean(dat_chrv$C31_norm)

chrv_V_C26_norm_in = sum((dat_chrv$C26_norm - dat_chrv_C26bar_norm)^2)
chrv_V_C26_norm_periodic = sum((dat_chrv$C26_norm - dat_chrv$SmoothC0_pred_unnorm)^2)
chrv_V_C26_norm_cross_terms = sum((dat_chrv$C26_norm - dat_chrv$SmoothC0_pred_unnorm)*(dat_chrv$SmoothC0_pred_unnorm - dat_chrv_C26bar_norm))

2*chrv_V_C26_norm_cross_terms/chrv_V_C26_norm_in
# -0.02886606

chrv_V_C26_norm_periodic/chrv_V_C26_norm_in
# 0.6529388


chrv_V_C29_norm_in = sum((dat_chrv$C29_norm - dat_chrv_C29bar_norm)^2)
chrv_V_C29_norm_periodic = sum((dat_chrv$C29_norm - dat_chrv$SmoothC0_pred_unnorm)^2)
chrv_V_C29_norm_cross_terms = sum((dat_chrv$C29_norm - dat_chrv$SmoothC0_pred_unnorm)*(dat_chrv$SmoothC0_pred_unnorm - dat_chrv_C29bar_norm))

2*chrv_V_C29_norm_cross_terms/chrv_V_C29_norm_in
# -0.08145486

chrv_V_C29_norm_periodic/chrv_V_C29_norm_in
# 0.7055276


chrv_V_C31_norm_in = sum((dat_chrv$C31_norm - dat_chrv_C31bar_norm)^2)
chrv_V_C31_norm_periodic = sum((dat_chrv$C31_norm - dat_chrv$SmoothC0_pred_unnorm)^2)
chrv_V_C31_norm_cross_terms = sum((dat_chrv$C31_norm - dat_chrv$SmoothC0_pred_unnorm)*(dat_chrv$SmoothC0_pred_unnorm - dat_chrv_C31bar_norm))

2*chrv_V_C31_norm_cross_terms/chrv_V_C31_norm_in
# -0.03961947

chrv_V_C31_norm_periodic/chrv_V_C31_norm_in
# 0.6636922




# Predicted Smooth C0:
dat_chrV_Cbar_norm = (mean(dat_chrV$C26) + mean(dat_chrV$C29) + mean(dat_chrV$C31))/3
dat_chrV_C26bar_norm = mean(dat_chrV$C26)
dat_chrV_C29bar_norm = mean(dat_chrV$C29)
dat_chrV_C31bar_norm = mean(dat_chrV$C31)

chrV_V_C26_norm_in = sum((dat_chrV$C26 - dat_chrV_C26bar_norm)^2)
chrV_V_C26_norm_periodic = sum((dat_chrV$C26 - dat_chrV$SmoothC0_pred_unnorm)^2)
chrV_V_C26_norm_cross_terms = sum((dat_chrV$C26 - dat_chrV$SmoothC0_pred_unnorm)*(dat_chrV$SmoothC0_pred_unnorm - dat_chrV_C26bar_norm))

2*chrV_V_C26_norm_cross_terms/chrV_V_C26_norm_in
# 0.1258134

chrV_V_C26_norm_periodic/chrV_V_C26_norm_in
# 0.4042445


chrV_V_C29_norm_in = sum((dat_chrV$C29 - dat_chrV_C29bar_norm)^2)
chrV_V_C29_norm_periodic = sum((dat_chrV$C29 - dat_chrV$SmoothC0_pred_unnorm)^2)
chrV_V_C29_norm_cross_terms = sum((dat_chrV$C29 - dat_chrV$SmoothC0_pred_unnorm)*(dat_chrV$SmoothC0_pred_unnorm - dat_chrV_C29bar_norm))

2*chrV_V_C29_norm_cross_terms/chrV_V_C29_norm_in
# -0.007181394

chrV_V_C29_norm_periodic/chrV_V_C29_norm_in
# 0.5477508


chrV_V_C31_norm_in = sum((dat_chrV$C31 - dat_chrV_C31bar_norm)^2)
chrV_V_C31_norm_periodic = sum((dat_chrV$C31 - dat_chrV$SmoothC0_pred_unnorm)^2)
chrV_V_C31_norm_cross_terms = sum((dat_chrV$C31 - dat_chrV$SmoothC0_pred_unnorm)*(dat_chrV$SmoothC0_pred_unnorm - dat_chrV_C31bar_norm))

2*chrV_V_C31_norm_cross_terms/chrV_V_C31_norm_in
# 0.0290143

chrV_V_C31_norm_periodic/chrV_V_C31_norm_in
# 0.5213315



# Calculated Smooth C0:
dat_chrV_smoothC0 = dat_chrV[which(!is.na(dat_chrV$SmoothC0)),]
dat_chrV_smoothC0_Cbar_norm = (mean(dat_chrV_smoothC0$C26) + mean(dat_chrV_smoothC0$C29) + mean(dat_chrV_smoothC0$C31))/3
dat_chrV_smoothC0_C26bar_norm = mean(dat_chrV_smoothC0$C26)
dat_chrV_smoothC0_C29bar_norm = mean(dat_chrV_smoothC0$C29)
dat_chrV_smoothC0_C31bar_norm = mean(dat_chrV_smoothC0$C31)

chrV_smoothC0_V_C26_norm_in = sum((dat_chrV_smoothC0$C26 - dat_chrV_smoothC0_C26bar_norm)^2)
chrV_smoothC0_V_C26_norm_periodic = sum((dat_chrV_smoothC0$C26 - dat_chrV_smoothC0$SmoothC0)^2)
chrV_smoothC0_V_C26_norm_cross_terms = sum((dat_chrV_smoothC0$C26 - dat_chrV_smoothC0$SmoothC0)*(dat_chrV_smoothC0$SmoothC0 - dat_chrV_smoothC0_C26bar_norm))

2*chrV_smoothC0_V_C26_norm_cross_terms/chrV_smoothC0_V_C26_norm_in
# 0.1093096

chrV_smoothC0_V_C26_norm_periodic/chrV_smoothC0_V_C26_norm_in
# 0.4028625


chrV_smoothC0_V_C29_norm_in = sum((dat_chrV_smoothC0$C29 - dat_chrV_smoothC0_C29bar_norm)^2)
chrV_smoothC0_V_C29_norm_periodic = sum((dat_chrV_smoothC0$C29 - dat_chrV_smoothC0$SmoothC0_pred_unnorm)^2)
chrV_smoothC0_V_C29_norm_cross_terms = sum((dat_chrV_smoothC0$C29 - dat_chrV_smoothC0$SmoothC0_pred_unnorm)*(dat_chrV_smoothC0$SmoothC0_pred_unnorm - dat_chrV_smoothC0_C29bar_norm))

2*chrV_smoothC0_V_C29_norm_cross_terms/chrV_smoothC0_V_C29_norm_in
# -0.007192708

chrV_smoothC0_V_C29_norm_periodic/chrV_smoothC0_V_C29_norm_in
# 0.5477574


chrV_smoothC0_V_C31_norm_in = sum((dat_chrV_smoothC0$C31 - dat_chrV_smoothC0_C31bar_norm)^2)
chrV_smoothC0_V_C31_norm_periodic = sum((dat_chrV_smoothC0$C31 - dat_chrV_smoothC0$SmoothC0_pred_unnorm)^2)
chrV_smoothC0_V_C31_norm_cross_terms = sum((dat_chrV_smoothC0$C31 - dat_chrV_smoothC0$SmoothC0_pred_unnorm)*(dat_chrV_smoothC0$SmoothC0_pred_unnorm - dat_chrV_smoothC0_C31bar_norm))

2*chrV_smoothC0_V_C31_norm_cross_terms/chrV_smoothC0_V_C31_norm_in
# 0.02901329

chrV_smoothC0_V_C31_norm_periodic/chrV_smoothC0_V_C31_norm_in
# 0.5213319





# Calculated Smooth C0 For Example Region:
starting_point = 479209
sequence_length = 100
dat_chrV_paper_region = dat_chrV[(starting_point:(starting_point+sequence_length)),]
dat_chrV_paper_region_Cbar_norm = (mean(dat_chrV_paper_region$C26) + mean(dat_chrV_paper_region$C29) + mean(dat_chrV_paper_region$C31))/3
dat_chrV_paper_region_C26bar_norm = mean(dat_chrV_paper_region$C26)
dat_chrV_paper_region_C29bar_norm = mean(dat_chrV_paper_region$C29)
dat_chrV_paper_region_C31bar_norm = mean(dat_chrV_paper_region$C31)

chrV_paper_region_V_C26_norm_in = sum((dat_chrV_paper_region$C26 - dat_chrV_paper_region_C26bar_norm)^2)
chrV_paper_region_V_C26_norm_periodic = sum((dat_chrV_paper_region$C26 - dat_chrV_paper_region$SmoothC0)^2)
chrV_paper_region_V_C26_norm_cross_terms = sum((dat_chrV_paper_region$C26 - dat_chrV_paper_region$SmoothC0)*(dat_chrV_paper_region$SmoothC0 - dat_chrV_paper_region_C26bar_norm))

2*chrV_paper_region_V_C26_norm_cross_terms/chrV_paper_region_V_C26_norm_in
# 0.06525573

chrV_paper_region_V_C26_norm_periodic/chrV_paper_region_V_C26_norm_in
# 0.717425


chrV_paper_region_V_C29_norm_in = sum((dat_chrV_paper_region$C29 - dat_chrV_paper_region_C29bar_norm)^2)
chrV_paper_region_V_C29_norm_periodic = sum((dat_chrV_paper_region$C29 - dat_chrV_paper_region$SmoothC0_pred_unnorm)^2)
chrV_paper_region_V_C29_norm_cross_terms = sum((dat_chrV_paper_region$C29 - dat_chrV_paper_region$SmoothC0_pred_unnorm)*(dat_chrV_paper_region$SmoothC0_pred_unnorm - dat_chrV_paper_region_C29bar_norm))

2*chrV_paper_region_V_C29_norm_cross_terms/chrV_paper_region_V_C29_norm_in
# 0.01488491

chrV_paper_region_V_C29_norm_periodic/chrV_paper_region_V_C29_norm_in
# 0.8142971


chrV_paper_region_V_C31_norm_in = sum((dat_chrV_paper_region$C31 - dat_chrV_paper_region_C31bar_norm)^2)
chrV_paper_region_V_C31_norm_periodic = sum((dat_chrV_paper_region$C31 - dat_chrV_paper_region$SmoothC0_pred_unnorm)^2)
chrV_paper_region_V_C31_norm_cross_terms = sum((dat_chrV_paper_region$C31 - dat_chrV_paper_region$SmoothC0_pred_unnorm)*(dat_chrV_paper_region$SmoothC0_pred_unnorm - dat_chrV_paper_region_C31bar_norm))

2*chrV_paper_region_V_C31_norm_cross_terms/chrV_paper_region_V_C31_norm_in
# 0.06553433

chrV_paper_region_V_C31_norm_periodic/chrV_paper_region_V_C31_norm_in
# 0.7778999











# ACROSS ALL N:

nuc_V_C_in = sum((dat_nuc$C26 - dat_nuc_Cbar)^2) + 
  sum((dat_nuc$C29 - dat_nuc_Cbar)^2) + 
  sum((dat_nuc$C31 - dat_nuc_Cbar)^2)

nuc_V_C_in_norm = sum((dat_nuc$C26_norm - dat_nuc_Cbar_norm)^2) + 
  sum((dat_nuc$C29_norm - dat_nuc_Cbar_norm)^2) + 
  sum((dat_nuc$C31_norm - dat_nuc_Cbar_norm)^2)

nuc_V_C_periodic = sum((dat_nuc$C26 - dat_nuc$C0)^2) +
  sum((dat_nuc$C29 - dat_nuc$C0)^2) +
  sum((dat_nuc$C31 - dat_nuc$C0)^2)

nuc_V_C_periodic_smooth_norm = sum((dat_nuc$C26_norm - dat_nuc$SmoothC0_pred)^2) +
  sum((dat_nuc$C29_norm - dat_nuc$SmoothC0_pred)^2) +
  sum((dat_nuc$C31_norm - dat_nuc$SmoothC0_pred)^2)

nuc_V_C_periodic_smooth_unnorm = sum((dat_nuc$C26 - dat_nuc$SmoothC0_pred_unnorm)^2) +
  sum((dat_nuc$C29 - dat_nuc$SmoothC0_pred_unnorm)^2) +
  sum((dat_nuc$C31 - dat_nuc$SmoothC0_pred_unnorm)^2)

nuc_V_cross_terms = sum((dat_nuc$C26 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_Cbar)) + 
  sum((dat_nuc$C29 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_Cbar)) + 
  sum((dat_nuc$C31 - dat_nuc$C0)*(dat_nuc$C0 - dat_nuc_Cbar))

nuc_V_cross_terms_smooth_norm = sum((dat_nuc$C26_norm - dat_nuc$SmoothC0_pred)*(dat_nuc$SmoothC0_pred - dat_nuc_Cbar_norm)) + 
  sum((dat_nuc$C29_norm - dat_nuc$SmoothC0_pred)*(dat_nuc$SmoothC0_pred - dat_nuc_Cbar_norm)) + 
  sum((dat_nuc$C31_norm - dat_nuc$SmoothC0_pred)*(dat_nuc$SmoothC0_pred - dat_nuc_Cbar_norm))

nuc_V_cross_terms_smooth_unnorm = sum((dat_nuc$C26 - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_Cbar)) + 
  sum((dat_nuc$C29 - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_Cbar)) + 
  sum((dat_nuc$C31 - dat_nuc$SmoothC0_pred_unnorm)*(dat_nuc$SmoothC0_pred_unnorm - dat_nuc_Cbar))


2*nuc_V_cross_terms/nuc_V_C_in
# -0.09785206

nuc_V_C_periodic/nuc_V_C_in
# 0.5013022


# ORIGINAL: (Think I put smoothC0 on the scale of the original C0)
2*nuc_V_cross_terms_smooth/nuc_V_C_in
# -0.1821574

nuc_V_C_periodic_smooth/nuc_V_C_in
# 0.5856075



2*nuc_V_cross_terms_smooth_norm/nuc_V_C_in_norm
# -0.525

nuc_V_C_periodic_smooth_norm/nuc_V_C_in_norm
# 0.6587253


2*nuc_V_cross_terms_smooth_unnorm/nuc_V_C_in
# -0.8605847

nuc_V_C_periodic_smooth_unnorm/nuc_V_C_in
# 0.7864136






random_V_C_in = sum((dat_random$C26 - dat_random_Cbar)^2) + 
  sum((dat_random$C29 - dat_random_Cbar)^2) + 
  sum((dat_random$C31 - dat_random_Cbar)^2)

random_V_C_periodic = sum((dat_random$C26 - dat_random$C0)^2) +
  sum((dat_random$C29 - dat_random$C0)^2) +
  sum((dat_random$C31 - dat_random$C0)^2)

random_V_C_periodic_smooth = sum((dat_random$C26 - dat_random$SmoothC0_pred)^2) +
  sum((dat_random$C29 - dat_random$SmoothC0_pred)^2) +
  sum((dat_random$C31 - dat_random$SmoothC0_pred)^2)

random_V_cross_terms = sum((dat_random$C26 - dat_random$C0)*(dat_random$C0 - dat_random_Cbar)) + 
  sum((dat_random$C29 - dat_random$C0)*(dat_random$C0 - dat_random_Cbar)) + 
  sum((dat_random$C31 - dat_random$C0)*(dat_random$C0 - dat_random_Cbar))

random_V_cross_terms_smooth = sum((dat_random$C26 - dat_random$SmoothC0_pred)*(dat_random$SmoothC0_pred - dat_random_Cbar)) + 
  sum((dat_random$C29 - dat_random$SmoothC0_pred)*(dat_random$SmoothC0_pred - dat_random_Cbar)) + 
  sum((dat_random$C31 - dat_random$SmoothC0_pred)*(dat_random$SmoothC0_pred - dat_random_Cbar))


2*random_V_cross_terms/random_V_C_in
# -0.06852252

random_V_C_periodic/random_V_C_in
# 0.5453104


2*random_V_cross_terms_smooth/random_V_C_in
# -0.1095594

random_V_C_periodic_smooth/random_V_C_in
# 0.5863473








tiling_V_C_in = sum((dat_tiling$C26 - dat_tiling_Cbar)^2) + 
  sum((dat_tiling$C29 - dat_tiling_Cbar)^2) + 
  sum((dat_tiling$C31 - dat_tiling_Cbar)^2)

tiling_V_C_periodic = sum((dat_tiling$C26 - dat_tiling$C0)^2) +
  sum((dat_tiling$C29 - dat_tiling$C0)^2) +
  sum((dat_tiling$C31 - dat_tiling$C0)^2)

tiling_V_C_periodic_smooth = sum((dat_tiling$C26 - dat_tiling$SmoothC0_pred)^2) +
  sum((dat_tiling$C29 - dat_tiling$SmoothC0_pred)^2) +
  sum((dat_tiling$C31 - dat_tiling$SmoothC0_pred)^2)

tiling_V_cross_terms = sum((dat_tiling$C26 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_Cbar)) + 
  sum((dat_tiling$C29 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_Cbar)) + 
  sum((dat_tiling$C31 - dat_tiling$C0)*(dat_tiling$C0 - dat_tiling_Cbar))

tiling_V_cross_terms_smooth = sum((dat_tiling$C26 - dat_tiling$SmoothC0_pred)*(dat_tiling$SmoothC0_pred - dat_tiling_Cbar)) + 
  sum((dat_tiling$C29 - dat_tiling$SmoothC0_pred)*(dat_tiling$SmoothC0_pred - dat_tiling_Cbar)) + 
  sum((dat_tiling$C31 - dat_tiling$SmoothC0_pred)*(dat_tiling$SmoothC0_pred - dat_tiling_Cbar))


2*tiling_V_cross_terms/tiling_V_C_in
# -0.1005638

tiling_V_C_periodic/tiling_V_C_in
# 0.4545159


2*tiling_V_cross_terms_smooth/tiling_V_C_in
# -0.1784878

tiling_V_C_periodic_smooth/tiling_V_C_in
# 0.53244








chrv_V_C_in = sum((dat_chrv$C26 - dat_chrv_Cbar)^2) + 
  sum((dat_chrv$C29 - dat_chrv_Cbar)^2) + 
  sum((dat_chrv$C31 - dat_chrv_Cbar)^2)

chrv_V_C_periodic = sum((dat_chrv$C26 - dat_chrv$C0)^2) +
  sum((dat_chrv$C29 - dat_chrv$C0)^2) +
  sum((dat_chrv$C31 - dat_chrv$C0)^2)

chrv_V_C_periodic_smooth = sum((dat_chrv$C26 - dat_chrv$SmoothC0_pred)^2) +
  sum((dat_chrv$C29 - dat_chrv$SmoothC0_pred)^2) +
  sum((dat_chrv$C31 - dat_chrv$SmoothC0_pred)^2)

chrv_V_cross_terms = sum((dat_chrv$C26 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_Cbar)) + 
  sum((dat_chrv$C29 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_Cbar)) + 
  sum((dat_chrv$C31 - dat_chrv$C0)*(dat_chrv$C0 - dat_chrv_Cbar))

chrv_V_cross_terms_smooth = sum((dat_chrv$C26 - dat_chrv$SmoothC0_pred)*(dat_chrv$SmoothC0_pred - dat_chrv_Cbar)) + 
  sum((dat_chrv$C29 - dat_chrv$SmoothC0_pred)*(dat_chrv$SmoothC0_pred - dat_chrv_Cbar)) + 
  sum((dat_chrv$C31 - dat_chrv$SmoothC0_pred)*(dat_chrv$SmoothC0_pred - dat_chrv_Cbar))


2*chrv_V_cross_terms/chrv_V_C_in
# -0.1866422

chrv_V_C_periodic/chrv_V_C_in
# 0.5890637


2*chrv_V_cross_terms_smooth/chrv_V_C_in
# -0.3047252

chrv_V_C_periodic_smooth/chrv_V_C_in
# 0.7071467

