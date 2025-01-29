nuc <- read_csv("cycle1.txt")
random <- read_csv("cycle3.txt")
tiling <- read_csv("cycle5.txt")
chrv <- read_csv("cycle6.txt")

colnames(nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

source("c-score-analysis/new-c0-functions.R")

# n=26,29,31, Amplitude has no dependence on n
nuc_C0Aphi_new_v1 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1,1,1))
random_C0Aphi_new_v1 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1,1,1))
tiling_C0Aphi_new_v1 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1,1,1))
chrv_C0Aphi_new_v1 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1,1,1))
nuc$C0_new_v1 = nuc_C0Aphi_new_v1[,1]
random$C0_new_v1 = random_C0Aphi_new_v1[,1]
tiling$C0_new_v1 = tiling_C0Aphi_new_v1[,1]
chrv$C0_new_v1 = chrv_C0Aphi_new_v1[,1]

# n=26,29,31, Amplitude depends on n as claimed in supplementary material:
nuc_C0Aphi_new_v2 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1,0.82,0.7))
random_C0Aphi_new_v2 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1,0.82,0.7))
tiling_C0Aphi_new_v2 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1,0.82,0.7))
chrv_C0Aphi_new_v2 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1,0.82,0.7))
nuc$C0_new_v2 = nuc_C0Aphi_new_v2[,1]
random$C0_new_v2 = random_C0Aphi_new_v2[,1]
tiling$C0_new_v2 = tiling_C0Aphi_new_v2[,1]
chrv$C0_new_v2 = chrv_C0Aphi_new_v2[,1]

# n=26,29,31, Optimal amplitude ratios calculated from top 4000 of Random Library (A/T)
nuc_C0Aphi_new_v3 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1.00, 1.66, 1.40))
random_C0Aphi_new_v3 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1.00, 1.66, 1.40))
tiling_C0Aphi_new_v3 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1.00, 1.66, 1.40))
chrv_C0Aphi_new_v3 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1.00, 1.66, 1.40))
nuc$C0_new_v3 = nuc_C0Aphi_new_v3[,1]
random$C0_new_v3 = random_C0Aphi_new_v3[,1]
tiling$C0_new_v3 = tiling_C0Aphi_new_v3[,1]
chrv$C0_new_v3 = chrv_C0Aphi_new_v3[,1]

# n=26,29,31, Optimal amplitude ratios calculated from top 25000 of Tiling Library (A/T)
nuc_C0Aphi_new_v4 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1.00, 0.39, 1.18))
random_C0Aphi_new_v4 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1.00, 0.39, 1.18))
tiling_C0Aphi_new_v4 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1.00, 0.39, 1.18))
chrv_C0Aphi_new_v4 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1.00, 0.39, 1.18))
nuc$C0_new_v4 = nuc_C0Aphi_new_v4[,1]
random$C0_new_v4 = random_C0Aphi_new_v4[,1]
tiling$C0_new_v4 = tiling_C0Aphi_new_v4[,1]
chrv$C0_new_v4 = chrv_C0Aphi_new_v4[,1]

# n=26,29,31, Optimal amplitude ratios calculated from top 2000 of Random Library (A/T)
nuc_C0Aphi_new_v5 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1.00, 1.36, 1.34))
random_C0Aphi_new_v5 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1.00, 1.36, 1.34))
tiling_C0Aphi_new_v5 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1.00, 1.36, 1.34))
chrv_C0Aphi_new_v5 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1.00, 1.36, 1.34))
nuc$C0_new_v5 = nuc_C0Aphi_new_v5[,1]
random$C0_new_v5 = random_C0Aphi_new_v5[,1]
tiling$C0_new_v5 = tiling_C0Aphi_new_v5[,1]
chrv$C0_new_v5 = chrv_C0Aphi_new_v5[,1]

# n=26,29,31, Optimal amplitude ratios calculated from top 15000 of Tiling Library (A/T)
nuc_C0Aphi_new_v6 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1.00, 0.35, 1.11))
random_C0Aphi_new_v6 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1.00, 0.35, 1.11))
tiling_C0Aphi_new_v6 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1.00, 0.35, 1.11))
chrv_C0Aphi_new_v6 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1.00, 0.35, 1.11))
nuc$C0_new_v6 = nuc_C0Aphi_new_v6[,1]
random$C0_new_v6 = random_C0Aphi_new_v6[,1]
tiling$C0_new_v6 = tiling_C0Aphi_new_v6[,1]
chrv$C0_new_v6 = chrv_C0Aphi_new_v6[,1]

# n, Amplitude ratios calculated from fft of A/T on top 4000 of Random Library (A/T)
nuc_C0Aphi_new_v7 = find_c0Aphi(nuc[,2:4], n=c(9.72, 3.02, 4.55), aa=c(1, 1.22, 1.28))
random_C0Aphi_new_v7 = find_c0Aphi(random[,2:4], n=c(9.72, 3.02, 4.55), aa=c(1, 1.22, 1.28))
tiling_C0Aphi_new_v7 = find_c0Aphi(tiling[,2:4], n=c(9.72, 3.02, 4.55), aa=c(1, 1.22, 1.28))
chrv_C0Aphi_new_v7 = find_c0Aphi(chrv[,2:4], n=c(9.72, 3.02, 4.55), aa=c(1, 1.22, 1.28))
nuc$C0_new_v7 = nuc_C0Aphi_new_v7[,1]
random$C0_new_v7 = random_C0Aphi_new_v7[,1]
tiling$C0_new_v7 = tiling_C0Aphi_new_v7[,1]
chrv$C0_new_v7 = chrv_C0Aphi_new_v7[,1]

# n, Amplitude ratios calculated from fft of A/T on top 25000 of Tiling Library (A/T)
nuc_C0Aphi_new_v8 = find_c0Aphi(nuc[,2:4], n=c(9.89, 2.81, 4.32), aa=c(1, 1.25, 1.19))
random_C0Aphi_new_v8 = find_c0Aphi(random[,2:4], n=c(9.89, 2.81, 4.32), aa=c(1, 1.25, 1.19))
tiling_C0Aphi_new_v8 = find_c0Aphi(tiling[,2:4], n=c(9.89, 2.81, 4.32), aa=c(1, 1.25, 1.19))
chrv_C0Aphi_new_v8 = find_c0Aphi(chrv[,2:4], n=c(9.89, 2.81, 4.32), aa=c(1, 1.25, 1.19))
nuc$C0_new_v8 = nuc_C0Aphi_new_v8[,1]
random$C0_new_v8 = random_C0Aphi_new_v8[,1]
tiling$C0_new_v8 = tiling_C0Aphi_new_v8[,1]
chrv$C0_new_v8 = chrv_C0Aphi_new_v8[,1]

# n, Amplitude ratios calculated from fft of A/T on top 2000 of Random Library (A/T)
nuc_C0Aphi_new_v9 = find_c0Aphi(nuc[,2:4], n=c(9.71, 2.99, 4.48), aa=c(1, 1.17, 1.15))
random_C0Aphi_new_v9 = find_c0Aphi(random[,2:4], n=c(9.71, 2.99, 4.48), aa=c(1, 1.17, 1.15))
tiling_C0Aphi_new_v9 = find_c0Aphi(tiling[,2:4], n=c(9.71, 2.99, 4.48), aa=c(1, 1.17, 1.15))
chrv_C0Aphi_new_v9 = find_c0Aphi(chrv[,2:4], n=c(9.71, 2.99, 4.48), aa=c(1, 1.17, 1.15))
nuc$C0_new_v9 = nuc_C0Aphi_new_v9[,1]
random$C0_new_v9 = random_C0Aphi_new_v9[,1]
tiling$C0_new_v9 = tiling_C0Aphi_new_v9[,1]
chrv$C0_new_v9 = chrv_C0Aphi_new_v9[,1]

# n, Amplitude ratios calculated from fft of A/T on top 15000 of Tiling Library (A/T)
nuc_C0Aphi_new_v10 = find_c0Aphi(nuc[,2:4], n=c(9.87, 2.82, 4.27), aa=c(1, 1.21, 1.16))
random_C0Aphi_new_v10 = find_c0Aphi(random[,2:4], n=c(9.87, 2.82, 4.27), aa=c(1, 1.21, 1.16))
tiling_C0Aphi_new_v10 = find_c0Aphi(tiling[,2:4], n=c(9.87, 2.82, 4.27), aa=c(1, 1.21, 1.16))
chrv_C0Aphi_new_v10 = find_c0Aphi(chrv[,2:4], n=c(9.87, 2.82, 4.27), aa=c(1, 1.21, 1.16))
nuc$C0_new_v10 = nuc_C0Aphi_new_v10[,1]
random$C0_new_v10 = random_C0Aphi_new_v10[,1]
tiling$C0_new_v10 = tiling_C0Aphi_new_v10[,1]
chrv$C0_new_v10 = chrv_C0Aphi_new_v10[,1]

# n=26,29,31, Optimal amplitude ratios calculated from average of top 25000 of Tiling Library 
# and top 4000 of Random Library (A/T)
nuc_C0Aphi_new_v11 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1, 0.64, 1.29))
random_C0Aphi_new_v11 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1, 0.64, 1.29))
tiling_C0Aphi_new_v11 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1, 0.64, 1.29))
chrv_C0Aphi_new_v11 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1, 0.64, 1.29))
nuc$C0_new_v11 = nuc_C0Aphi_new_v11[,1]
random$C0_new_v11 = random_C0Aphi_new_v11[,1]
tiling$C0_new_v11 = tiling_C0Aphi_new_v11[,1]
chrv$C0_new_v11 = chrv_C0Aphi_new_v11[,1]

# # n=26,29,31, Optimal amplitude ratios calculated from average of top 15000 of Tiling Library 
# # and top 2000 of Random Library (A/T)
# nuc_C0Aphi_new_v12 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1, 0.68, 1.24))
# random_C0Aphi_new_v12 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1, 0.68, 1.24))
# tiling_C0Aphi_new_v12 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1, 0.68, 1.24))
# chrv_C0Aphi_new_v12 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1, 0.68, 1.24))
# nuc$C0_new_v12 = nuc_C0Aphi_new_v12[,1]
# random$C0_new_v12 = random_C0Aphi_new_v12[,1]
# tiling$C0_new_v12 = tiling_C0Aphi_new_v12[,1]
# chrv$C0_new_v12 = chrv_C0Aphi_new_v12[,1]

# n=26,29,31, Optimal amplitude ratios calculated from average of top 25000 of Tiling Library, 
# top 25000 of ChrV Library, and top 4000 of Random Library (A/T)
nuc_C0Aphi_new_v12 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1, 0.55, 1.22))
random_C0Aphi_new_v12 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1, 0.55, 1.22))
tiling_C0Aphi_new_v12 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1, 0.55, 1.22))
chrv_C0Aphi_new_v12 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1, 0.55, 1.22))
nuc$C0_new_v12 = nuc_C0Aphi_new_v12[,1]
random$C0_new_v12 = random_C0Aphi_new_v12[,1]
tiling$C0_new_v12 = tiling_C0Aphi_new_v12[,1]
chrv$C0_new_v12 = chrv_C0Aphi_new_v12[,1]

# n=26,29,31, Optimal amplitude ratios calculated from average of top 15000 of Tiling Library, 
# top 15000 of ChrV Library, and top 2000 of Random Library (A/T)
nuc_C0Aphi_new_v13 = find_c0Aphi(nuc[,2:4], n=c(26,29,31), aa=c(1, 0.59, 1.2))
random_C0Aphi_new_v13 = find_c0Aphi(random[,2:4], n=c(26,29,31), aa=c(1, 0.59, 1.2))
tiling_C0Aphi_new_v13 = find_c0Aphi(tiling[,2:4], n=c(26,29,31), aa=c(1, 0.59, 1.2))
chrv_C0Aphi_new_v13 = find_c0Aphi(chrv[,2:4], n=c(26,29,31), aa=c(1, 0.59, 1.2))
nuc$C0_new_v13 = nuc_C0Aphi_new_v13[,1]
random$C0_new_v13 = random_C0Aphi_new_v13[,1]
tiling$C0_new_v13 = tiling_C0Aphi_new_v13[,1]
chrv$C0_new_v13 = chrv_C0Aphi_new_v13[,1]

#FIXME:
# Best of for each:
# Nuc: FIXME (Same as Tiling)
# Random: FIXME (e.g. v7)
# Tiling: FIXME (e.g. v7)
# ChrV: FIXME (e.g. v7)
# nuc$C0_new_v13 = nuc$C0_new_v7
# random$C0_new_v13 = random$C0_new_v7
# tiling$C0_new_v13 = tiling$C0_new_v7
# chrv$C0_new_v3 = chrv$C0_new_v7

cor(nuc %>% select(C26, C29, C31, C0, C0_new_v12))
cor(random %>% select(C26, C29, C31, C0, C0_new_v12))
cor(tiling %>% select(C26, C29, C31, C0, C0_new_v12))
cor(chrv %>% select(C26, C29, C31, C0, C0_new_v12))

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

source("scripts/functions/sequence-functions.R")

sequence_1_df_nuc <- nuc %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_nuc <- nuc %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(nuc), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_nuc <- sequence_1_df_nuc %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_nuc <- sequence_2_df_nuc %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)


sequence_1_df_random <- random %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_random <- random %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(random), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_random <- sequence_1_df_random %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_random <- sequence_2_df_random %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)


sequence_1_df_tiling <- tiling %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_tiling <- tiling %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(tiling), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_tiling <- sequence_1_df_tiling %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_tiling <- sequence_2_df_tiling %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)


sequence_1_df_chrv <- chrv %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrv <- chrv %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(chrv), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_chrv <- sequence_1_df_chrv %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_chrv <- sequence_2_df_chrv %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)



nuc = cbind(nuc, sequence_1_factor_nuc, sequence_2_factor_nuc)
random = cbind(random, sequence_1_factor_random, sequence_2_factor_random)
tiling = cbind(tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
chrv = cbind(chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)


# Find cutoffs for first-fourth quartiles + custom of C26
cutoffs_C26_0_25_nuc = quantile(nuc$C26, c(0, 0.25)) 
cutoffs_C26_25_50_nuc = quantile(nuc$C26, c(0.25, 0.5)) 
cutoffs_C26_50_75_nuc = quantile(nuc$C26, c(0.5, 0.75)) 
cutoffs_C26_75_100_nuc = quantile(nuc$C26, c(0.75, 1)) 
cutoffs_C26_custom1_nuc = quantile(nuc$C26, c(0.7, 1)) 

cutoffs_C26_0_25_random = quantile(random$C26, c(0, 0.25)) 
cutoffs_C26_25_50_random = quantile(random$C26, c(0.25, 0.5)) 
cutoffs_C26_50_75_random = quantile(random$C26, c(0.5, 0.75)) 
cutoffs_C26_75_100_random = quantile(random$C26, c(0.75, 1)) 
cutoffs_C26_custom1_random = quantile(random$C26, c(0.7, 1)) 

cutoffs_C26_0_25_tiling = quantile(tiling$C26, c(0, 0.25)) 
cutoffs_C26_25_50_tiling = quantile(tiling$C26, c(0.25, 0.5)) 
cutoffs_C26_50_75_tiling = quantile(tiling$C26, c(0.5, 0.75)) 
cutoffs_C26_75_100_tiling = quantile(tiling$C26, c(0.75, 1)) 
cutoffs_C26_custom1_tiling = quantile(tiling$C26, c(0.7, 1)) 

cutoffs_C26_0_25_chrv = quantile(chrv$C26, c(0, 0.25)) 
cutoffs_C26_25_50_chrv = quantile(chrv$C26, c(0.25, 0.5)) 
cutoffs_C26_50_75_chrv = quantile(chrv$C26, c(0.5, 0.75)) 
cutoffs_C26_75_100_chrv = quantile(chrv$C26, c(0.75, 1)) 
cutoffs_C26_custom1_chrv = quantile(chrv$C26, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C29
cutoffs_C29_0_25_nuc = quantile(nuc$C29, c(0, 0.25)) 
cutoffs_C29_25_50_nuc = quantile(nuc$C29, c(0.25, 0.5)) 
cutoffs_C29_50_75_nuc = quantile(nuc$C29, c(0.5, 0.75)) 
cutoffs_C29_75_100_nuc = quantile(nuc$C29, c(0.75, 1)) 
cutoffs_C29_custom1_nuc = quantile(nuc$C29, c(0.7, 1)) 

cutoffs_C29_0_25_random = quantile(random$C29, c(0, 0.25)) 
cutoffs_C29_25_50_random = quantile(random$C29, c(0.25, 0.5)) 
cutoffs_C29_50_75_random = quantile(random$C29, c(0.5, 0.75)) 
cutoffs_C29_75_100_random = quantile(random$C29, c(0.75, 1)) 
cutoffs_C29_custom1_random = quantile(random$C29, c(0.7, 1)) 

cutoffs_C29_0_25_tiling = quantile(tiling$C29, c(0, 0.25)) 
cutoffs_C29_25_50_tiling = quantile(tiling$C29, c(0.25, 0.5)) 
cutoffs_C29_50_75_tiling = quantile(tiling$C29, c(0.5, 0.75)) 
cutoffs_C29_75_100_tiling = quantile(tiling$C29, c(0.75, 1)) 
cutoffs_C29_custom1_tiling = quantile(tiling$C29, c(0.7, 1)) 

cutoffs_C29_0_25_chrv = quantile(chrv$C29, c(0, 0.25)) 
cutoffs_C29_25_50_chrv = quantile(chrv$C29, c(0.25, 0.5)) 
cutoffs_C29_50_75_chrv = quantile(chrv$C29, c(0.5, 0.75)) 
cutoffs_C29_75_100_chrv = quantile(chrv$C29, c(0.75, 1)) 
cutoffs_C29_custom1_chrv = quantile(chrv$C29, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C31
cutoffs_C31_0_25_nuc = quantile(nuc$C31, c(0, 0.25)) 
cutoffs_C31_25_50_nuc = quantile(nuc$C31, c(0.25, 0.5)) 
cutoffs_C31_50_75_nuc = quantile(nuc$C31, c(0.5, 0.75)) 
cutoffs_C31_75_100_nuc = quantile(nuc$C31, c(0.75, 1)) 
cutoffs_C31_custom1_nuc = quantile(nuc$C31, c(0.7, 1)) 

cutoffs_C31_0_25_random = quantile(random$C31, c(0, 0.25)) 
cutoffs_C31_25_50_random = quantile(random$C31, c(0.25, 0.5)) 
cutoffs_C31_50_75_random = quantile(random$C31, c(0.5, 0.75)) 
cutoffs_C31_75_100_random = quantile(random$C31, c(0.75, 1)) 
cutoffs_C31_custom1_random = quantile(random$C31, c(0.7, 1)) 

cutoffs_C31_0_25_tiling = quantile(tiling$C31, c(0, 0.25)) 
cutoffs_C31_25_50_tiling = quantile(tiling$C31, c(0.25, 0.5)) 
cutoffs_C31_50_75_tiling = quantile(tiling$C31, c(0.5, 0.75)) 
cutoffs_C31_75_100_tiling = quantile(tiling$C31, c(0.75, 1)) 
cutoffs_C31_custom1_tiling = quantile(tiling$C31, c(0.7, 1)) 

cutoffs_C31_0_25_chrv = quantile(chrv$C31, c(0, 0.25)) 
cutoffs_C31_25_50_chrv = quantile(chrv$C31, c(0.25, 0.5)) 
cutoffs_C31_50_75_chrv = quantile(chrv$C31, c(0.5, 0.75)) 
cutoffs_C31_75_100_chrv = quantile(chrv$C31, c(0.75, 1)) 
cutoffs_C31_custom1_chrv = quantile(chrv$C31, c(0.7, 1)) 

# Find cutoffs for first-fourth quartiles + custom of C0
cutoffs_C0_0_25_nuc = quantile(nuc$C0, c(0, 0.25)) 
cutoffs_C0_25_50_nuc = quantile(nuc$C0, c(0.25, 0.5)) 
cutoffs_C0_50_75_nuc = quantile(nuc$C0, c(0.5, 0.75)) 
cutoffs_C0_75_100_nuc = quantile(nuc$C0, c(0.75, 1)) 
cutoffs_C0_custom1_nuc = quantile(nuc$C0, c(0.7, 1)) 

cutoffs_C0_0_25_random = quantile(random$C0, c(0, 0.25)) 
cutoffs_C0_25_50_random = quantile(random$C0, c(0.25, 0.5)) 
cutoffs_C0_50_75_random = quantile(random$C0, c(0.5, 0.75)) 
cutoffs_C0_75_100_random = quantile(random$C0, c(0.75, 1)) 
cutoffs_C0_custom1_random = quantile(random$C0, c(0.7, 1)) 

cutoffs_C0_0_25_tiling = quantile(tiling$C0, c(0, 0.25)) 
cutoffs_C0_25_50_tiling = quantile(tiling$C0, c(0.25, 0.5)) 
cutoffs_C0_50_75_tiling = quantile(tiling$C0, c(0.5, 0.75)) 
cutoffs_C0_75_100_tiling = quantile(tiling$C0, c(0.75, 1)) 
cutoffs_C0_custom1_tiling = quantile(tiling$C0, c(0.7, 1)) 

cutoffs_C0_0_25_chrv = quantile(chrv$C0, c(0, 0.25)) 
cutoffs_C0_25_50_chrv = quantile(chrv$C0, c(0.25, 0.5)) 
cutoffs_C0_50_75_chrv = quantile(chrv$C0, c(0.5, 0.75)) 
cutoffs_C0_75_100_chrv = quantile(chrv$C0, c(0.75, 1)) 
cutoffs_C0_custom1_chrv = quantile(chrv$C0, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v1
cutoffs_C0_new_v1_0_25_nuc = quantile(nuc$C0_new_v1, c(0, 0.25)) 
cutoffs_C0_new_v1_25_50_nuc = quantile(nuc$C0_new_v1, c(0.25, 0.5)) 
cutoffs_C0_new_v1_50_75_nuc = quantile(nuc$C0_new_v1, c(0.5, 0.75)) 
cutoffs_C0_new_v1_75_100_nuc = quantile(nuc$C0_new_v1, c(0.75, 1)) 
cutoffs_C0_new_v1_custom1_nuc = quantile(nuc$C0_new_v1, c(0.7, 1)) 

cutoffs_C0_new_v1_0_25_random = quantile(random$C0_new_v1, c(0, 0.25)) 
cutoffs_C0_new_v1_25_50_random = quantile(random$C0_new_v1, c(0.25, 0.5)) 
cutoffs_C0_new_v1_50_75_random = quantile(random$C0_new_v1, c(0.5, 0.75)) 
cutoffs_C0_new_v1_75_100_random = quantile(random$C0_new_v1, c(0.75, 1)) 
cutoffs_C0_new_v1_custom1_random = quantile(random$C0_new_v1, c(0.7, 1)) 

cutoffs_C0_new_v1_0_25_tiling = quantile(tiling$C0_new_v1, c(0, 0.25)) 
cutoffs_C0_new_v1_25_50_tiling = quantile(tiling$C0_new_v1, c(0.25, 0.5)) 
cutoffs_C0_new_v1_50_75_tiling = quantile(tiling$C0_new_v1, c(0.5, 0.75)) 
cutoffs_C0_new_v1_75_100_tiling = quantile(tiling$C0_new_v1, c(0.75, 1)) 
cutoffs_C0_new_v1_custom1_tiling = quantile(tiling$C0_new_v1, c(0.7, 1)) 

cutoffs_C0_new_v1_0_25_chrv = quantile(chrv$C0_new_v1, c(0, 0.25)) 
cutoffs_C0_new_v1_25_50_chrv = quantile(chrv$C0_new_v1, c(0.25, 0.5)) 
cutoffs_C0_new_v1_50_75_chrv = quantile(chrv$C0_new_v1, c(0.5, 0.75)) 
cutoffs_C0_new_v1_75_100_chrv = quantile(chrv$C0_new_v1, c(0.75, 1)) 
cutoffs_C0_new_v1_custom1_chrv = quantile(chrv$C0_new_v1, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v2
cutoffs_C0_new_v2_0_25_nuc = quantile(nuc$C0_new_v2, c(0, 0.25)) 
cutoffs_C0_new_v2_25_50_nuc = quantile(nuc$C0_new_v2, c(0.25, 0.5)) 
cutoffs_C0_new_v2_50_75_nuc = quantile(nuc$C0_new_v2, c(0.5, 0.75)) 
cutoffs_C0_new_v2_75_100_nuc = quantile(nuc$C0_new_v2, c(0.75, 1)) 
cutoffs_C0_new_v2_custom1_nuc = quantile(nuc$C0_new_v2, c(0.7, 1)) 

cutoffs_C0_new_v2_0_25_random = quantile(random$C0_new_v2, c(0, 0.25)) 
cutoffs_C0_new_v2_25_50_random = quantile(random$C0_new_v2, c(0.25, 0.5)) 
cutoffs_C0_new_v2_50_75_random = quantile(random$C0_new_v2, c(0.5, 0.75)) 
cutoffs_C0_new_v2_75_100_random = quantile(random$C0_new_v2, c(0.75, 1)) 
cutoffs_C0_new_v2_custom1_random = quantile(random$C0_new_v2, c(0.7, 1)) 

cutoffs_C0_new_v2_0_25_tiling = quantile(tiling$C0_new_v2, c(0, 0.25)) 
cutoffs_C0_new_v2_25_50_tiling = quantile(tiling$C0_new_v2, c(0.25, 0.5)) 
cutoffs_C0_new_v2_50_75_tiling = quantile(tiling$C0_new_v2, c(0.5, 0.75)) 
cutoffs_C0_new_v2_75_100_tiling = quantile(tiling$C0_new_v2, c(0.75, 1)) 
cutoffs_C0_new_v2_custom1_tiling = quantile(tiling$C0_new_v2, c(0.7, 1)) 

cutoffs_C0_new_v2_0_25_chrv = quantile(chrv$C0_new_v2, c(0, 0.25)) 
cutoffs_C0_new_v2_25_50_chrv = quantile(chrv$C0_new_v2, c(0.25, 0.5)) 
cutoffs_C0_new_v2_50_75_chrv = quantile(chrv$C0_new_v2, c(0.5, 0.75)) 
cutoffs_C0_new_v2_75_100_chrv = quantile(chrv$C0_new_v2, c(0.75, 1)) 
cutoffs_C0_new_v2_custom1_chrv = quantile(chrv$C0_new_v2, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v3
cutoffs_C0_new_v3_0_25_nuc = quantile(nuc$C0_new_v3, c(0, 0.25)) 
cutoffs_C0_new_v3_25_50_nuc = quantile(nuc$C0_new_v3, c(0.25, 0.5)) 
cutoffs_C0_new_v3_50_75_nuc = quantile(nuc$C0_new_v3, c(0.5, 0.75)) 
cutoffs_C0_new_v3_75_100_nuc = quantile(nuc$C0_new_v3, c(0.75, 1)) 
cutoffs_C0_new_v3_custom1_nuc = quantile(nuc$C0_new_v3, c(0.7, 1)) 

cutoffs_C0_new_v3_0_25_random = quantile(random$C0_new_v3, c(0, 0.25)) 
cutoffs_C0_new_v3_25_50_random = quantile(random$C0_new_v3, c(0.25, 0.5)) 
cutoffs_C0_new_v3_50_75_random = quantile(random$C0_new_v3, c(0.5, 0.75)) 
cutoffs_C0_new_v3_75_100_random = quantile(random$C0_new_v3, c(0.75, 1)) 
cutoffs_C0_new_v3_custom1_random = quantile(random$C0_new_v3, c(0.7, 1)) 

cutoffs_C0_new_v3_0_25_tiling = quantile(tiling$C0_new_v3, c(0, 0.25)) 
cutoffs_C0_new_v3_25_50_tiling = quantile(tiling$C0_new_v3, c(0.25, 0.5)) 
cutoffs_C0_new_v3_50_75_tiling = quantile(tiling$C0_new_v3, c(0.5, 0.75)) 
cutoffs_C0_new_v3_75_100_tiling = quantile(tiling$C0_new_v3, c(0.75, 1)) 
cutoffs_C0_new_v3_custom1_tiling = quantile(tiling$C0_new_v3, c(0.7, 1)) 

cutoffs_C0_new_v3_0_25_chrv = quantile(chrv$C0_new_v3, c(0, 0.25)) 
cutoffs_C0_new_v3_25_50_chrv = quantile(chrv$C0_new_v3, c(0.25, 0.5)) 
cutoffs_C0_new_v3_50_75_chrv = quantile(chrv$C0_new_v3, c(0.5, 0.75)) 
cutoffs_C0_new_v3_75_100_chrv = quantile(chrv$C0_new_v3, c(0.75, 1)) 
cutoffs_C0_new_v3_custom1_chrv = quantile(chrv$C0_new_v3, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v4
cutoffs_C0_new_v4_0_25_nuc = quantile(nuc$C0_new_v4, c(0, 0.25)) 
cutoffs_C0_new_v4_25_50_nuc = quantile(nuc$C0_new_v4, c(0.25, 0.5)) 
cutoffs_C0_new_v4_50_75_nuc = quantile(nuc$C0_new_v4, c(0.5, 0.75)) 
cutoffs_C0_new_v4_75_100_nuc = quantile(nuc$C0_new_v4, c(0.75, 1)) 
cutoffs_C0_new_v4_custom1_nuc = quantile(nuc$C0_new_v4, c(0.7, 1)) 

cutoffs_C0_new_v4_0_25_random = quantile(random$C0_new_v4, c(0, 0.25)) 
cutoffs_C0_new_v4_25_50_random = quantile(random$C0_new_v4, c(0.25, 0.5)) 
cutoffs_C0_new_v4_50_75_random = quantile(random$C0_new_v4, c(0.5, 0.75)) 
cutoffs_C0_new_v4_75_100_random = quantile(random$C0_new_v4, c(0.75, 1)) 
cutoffs_C0_new_v4_custom1_random = quantile(random$C0_new_v4, c(0.7, 1)) 

cutoffs_C0_new_v4_0_25_tiling = quantile(tiling$C0_new_v4, c(0, 0.25)) 
cutoffs_C0_new_v4_25_50_tiling = quantile(tiling$C0_new_v4, c(0.25, 0.5)) 
cutoffs_C0_new_v4_50_75_tiling = quantile(tiling$C0_new_v4, c(0.5, 0.75)) 
cutoffs_C0_new_v4_75_100_tiling = quantile(tiling$C0_new_v4, c(0.75, 1)) 
cutoffs_C0_new_v4_custom1_tiling = quantile(tiling$C0_new_v4, c(0.7, 1)) 

cutoffs_C0_new_v4_0_25_chrv = quantile(chrv$C0_new_v4, c(0, 0.25)) 
cutoffs_C0_new_v4_25_50_chrv = quantile(chrv$C0_new_v4, c(0.25, 0.5)) 
cutoffs_C0_new_v4_50_75_chrv = quantile(chrv$C0_new_v4, c(0.5, 0.75)) 
cutoffs_C0_new_v4_75_100_chrv = quantile(chrv$C0_new_v4, c(0.75, 1)) 
cutoffs_C0_new_v4_custom1_chrv = quantile(chrv$C0_new_v4, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v5
cutoffs_C0_new_v5_0_25_nuc = quantile(nuc$C0_new_v5, c(0, 0.25)) 
cutoffs_C0_new_v5_25_50_nuc = quantile(nuc$C0_new_v5, c(0.25, 0.5)) 
cutoffs_C0_new_v5_50_75_nuc = quantile(nuc$C0_new_v5, c(0.5, 0.75)) 
cutoffs_C0_new_v5_75_100_nuc = quantile(nuc$C0_new_v5, c(0.75, 1)) 
cutoffs_C0_new_v5_custom1_nuc = quantile(nuc$C0_new_v5, c(0.7, 1)) 

cutoffs_C0_new_v5_0_25_random = quantile(random$C0_new_v5, c(0, 0.25)) 
cutoffs_C0_new_v5_25_50_random = quantile(random$C0_new_v5, c(0.25, 0.5)) 
cutoffs_C0_new_v5_50_75_random = quantile(random$C0_new_v5, c(0.5, 0.75)) 
cutoffs_C0_new_v5_75_100_random = quantile(random$C0_new_v5, c(0.75, 1)) 
cutoffs_C0_new_v5_custom1_random = quantile(random$C0_new_v5, c(0.7, 1)) 

cutoffs_C0_new_v5_0_25_tiling = quantile(tiling$C0_new_v5, c(0, 0.25)) 
cutoffs_C0_new_v5_25_50_tiling = quantile(tiling$C0_new_v5, c(0.25, 0.5)) 
cutoffs_C0_new_v5_50_75_tiling = quantile(tiling$C0_new_v5, c(0.5, 0.75)) 
cutoffs_C0_new_v5_75_100_tiling = quantile(tiling$C0_new_v5, c(0.75, 1)) 
cutoffs_C0_new_v5_custom1_tiling = quantile(tiling$C0_new_v5, c(0.7, 1)) 

cutoffs_C0_new_v5_0_25_chrv = quantile(chrv$C0_new_v5, c(0, 0.25)) 
cutoffs_C0_new_v5_25_50_chrv = quantile(chrv$C0_new_v5, c(0.25, 0.5)) 
cutoffs_C0_new_v5_50_75_chrv = quantile(chrv$C0_new_v5, c(0.5, 0.75)) 
cutoffs_C0_new_v5_75_100_chrv = quantile(chrv$C0_new_v5, c(0.75, 1)) 
cutoffs_C0_new_v5_custom1_chrv = quantile(chrv$C0_new_v5, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v6
cutoffs_C0_new_v6_0_25_nuc = quantile(nuc$C0_new_v6, c(0, 0.25)) 
cutoffs_C0_new_v6_25_50_nuc = quantile(nuc$C0_new_v6, c(0.25, 0.5)) 
cutoffs_C0_new_v6_50_75_nuc = quantile(nuc$C0_new_v6, c(0.5, 0.75)) 
cutoffs_C0_new_v6_75_100_nuc = quantile(nuc$C0_new_v6, c(0.75, 1)) 
cutoffs_C0_new_v6_custom1_nuc = quantile(nuc$C0_new_v6, c(0.7, 1)) 

cutoffs_C0_new_v6_0_25_random = quantile(random$C0_new_v6, c(0, 0.25)) 
cutoffs_C0_new_v6_25_50_random = quantile(random$C0_new_v6, c(0.25, 0.5)) 
cutoffs_C0_new_v6_50_75_random = quantile(random$C0_new_v6, c(0.5, 0.75)) 
cutoffs_C0_new_v6_75_100_random = quantile(random$C0_new_v6, c(0.75, 1)) 
cutoffs_C0_new_v6_custom1_random = quantile(random$C0_new_v6, c(0.7, 1)) 

cutoffs_C0_new_v6_0_25_tiling = quantile(tiling$C0_new_v6, c(0, 0.25)) 
cutoffs_C0_new_v6_25_50_tiling = quantile(tiling$C0_new_v6, c(0.25, 0.5)) 
cutoffs_C0_new_v6_50_75_tiling = quantile(tiling$C0_new_v6, c(0.5, 0.75)) 
cutoffs_C0_new_v6_75_100_tiling = quantile(tiling$C0_new_v6, c(0.75, 1)) 
cutoffs_C0_new_v6_custom1_tiling = quantile(tiling$C0_new_v6, c(0.7, 1)) 

cutoffs_C0_new_v6_0_25_chrv = quantile(chrv$C0_new_v6, c(0, 0.25)) 
cutoffs_C0_new_v6_25_50_chrv = quantile(chrv$C0_new_v6, c(0.25, 0.5)) 
cutoffs_C0_new_v6_50_75_chrv = quantile(chrv$C0_new_v6, c(0.5, 0.75)) 
cutoffs_C0_new_v6_75_100_chrv = quantile(chrv$C0_new_v6, c(0.75, 1)) 
cutoffs_C0_new_v6_custom1_chrv = quantile(chrv$C0_new_v6, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v7
cutoffs_C0_new_v7_0_25_nuc = quantile(nuc$C0_new_v7, c(0, 0.25)) 
cutoffs_C0_new_v7_25_50_nuc = quantile(nuc$C0_new_v7, c(0.25, 0.5)) 
cutoffs_C0_new_v7_50_75_nuc = quantile(nuc$C0_new_v7, c(0.5, 0.75)) 
cutoffs_C0_new_v7_75_100_nuc = quantile(nuc$C0_new_v7, c(0.75, 1)) 
cutoffs_C0_new_v7_custom1_nuc = quantile(nuc$C0_new_v7, c(0.7, 1)) 

cutoffs_C0_new_v7_0_25_random = quantile(random$C0_new_v7, c(0, 0.25)) 
cutoffs_C0_new_v7_25_50_random = quantile(random$C0_new_v7, c(0.25, 0.5)) 
cutoffs_C0_new_v7_50_75_random = quantile(random$C0_new_v7, c(0.5, 0.75)) 
cutoffs_C0_new_v7_75_100_random = quantile(random$C0_new_v7, c(0.75, 1)) 
cutoffs_C0_new_v7_custom1_random = quantile(random$C0_new_v7, c(0.7, 1)) 

cutoffs_C0_new_v7_0_25_tiling = quantile(tiling$C0_new_v7, c(0, 0.25)) 
cutoffs_C0_new_v7_25_50_tiling = quantile(tiling$C0_new_v7, c(0.25, 0.5)) 
cutoffs_C0_new_v7_50_75_tiling = quantile(tiling$C0_new_v7, c(0.5, 0.75)) 
cutoffs_C0_new_v7_75_100_tiling = quantile(tiling$C0_new_v7, c(0.75, 1)) 
cutoffs_C0_new_v7_custom1_tiling = quantile(tiling$C0_new_v7, c(0.7, 1)) 

cutoffs_C0_new_v7_0_25_chrv = quantile(chrv$C0_new_v7, c(0, 0.25)) 
cutoffs_C0_new_v7_25_50_chrv = quantile(chrv$C0_new_v7, c(0.25, 0.5)) 
cutoffs_C0_new_v7_50_75_chrv = quantile(chrv$C0_new_v7, c(0.5, 0.75)) 
cutoffs_C0_new_v7_75_100_chrv = quantile(chrv$C0_new_v7, c(0.75, 1)) 
cutoffs_C0_new_v7_custom1_chrv = quantile(chrv$C0_new_v7, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v8
cutoffs_C0_new_v8_0_25_nuc = quantile(nuc$C0_new_v8, c(0, 0.25)) 
cutoffs_C0_new_v8_25_50_nuc = quantile(nuc$C0_new_v8, c(0.25, 0.5)) 
cutoffs_C0_new_v8_50_75_nuc = quantile(nuc$C0_new_v8, c(0.5, 0.75)) 
cutoffs_C0_new_v8_75_100_nuc = quantile(nuc$C0_new_v8, c(0.75, 1)) 
cutoffs_C0_new_v8_custom1_nuc = quantile(nuc$C0_new_v8, c(0.7, 1)) 

cutoffs_C0_new_v8_0_25_random = quantile(random$C0_new_v8, c(0, 0.25)) 
cutoffs_C0_new_v8_25_50_random = quantile(random$C0_new_v8, c(0.25, 0.5)) 
cutoffs_C0_new_v8_50_75_random = quantile(random$C0_new_v8, c(0.5, 0.75)) 
cutoffs_C0_new_v8_75_100_random = quantile(random$C0_new_v8, c(0.75, 1)) 
cutoffs_C0_new_v8_custom1_random = quantile(random$C0_new_v8, c(0.7, 1)) 

cutoffs_C0_new_v8_0_25_tiling = quantile(tiling$C0_new_v8, c(0, 0.25)) 
cutoffs_C0_new_v8_25_50_tiling = quantile(tiling$C0_new_v8, c(0.25, 0.5)) 
cutoffs_C0_new_v8_50_75_tiling = quantile(tiling$C0_new_v8, c(0.5, 0.75)) 
cutoffs_C0_new_v8_75_100_tiling = quantile(tiling$C0_new_v8, c(0.75, 1)) 
cutoffs_C0_new_v8_custom1_tiling = quantile(tiling$C0_new_v8, c(0.7, 1)) 

cutoffs_C0_new_v8_0_25_chrv = quantile(chrv$C0_new_v8, c(0, 0.25)) 
cutoffs_C0_new_v8_25_50_chrv = quantile(chrv$C0_new_v8, c(0.25, 0.5)) 
cutoffs_C0_new_v8_50_75_chrv = quantile(chrv$C0_new_v8, c(0.5, 0.75)) 
cutoffs_C0_new_v8_75_100_chrv = quantile(chrv$C0_new_v8, c(0.75, 1)) 
cutoffs_C0_new_v8_custom1_chrv = quantile(chrv$C0_new_v8, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v9
cutoffs_C0_new_v9_0_25_nuc = quantile(nuc$C0_new_v9, c(0, 0.25)) 
cutoffs_C0_new_v9_25_50_nuc = quantile(nuc$C0_new_v9, c(0.25, 0.5)) 
cutoffs_C0_new_v9_50_75_nuc = quantile(nuc$C0_new_v9, c(0.5, 0.75)) 
cutoffs_C0_new_v9_75_100_nuc = quantile(nuc$C0_new_v9, c(0.75, 1)) 
cutoffs_C0_new_v9_custom1_nuc = quantile(nuc$C0_new_v9, c(0.7, 1)) 

cutoffs_C0_new_v9_0_25_random = quantile(random$C0_new_v9, c(0, 0.25)) 
cutoffs_C0_new_v9_25_50_random = quantile(random$C0_new_v9, c(0.25, 0.5)) 
cutoffs_C0_new_v9_50_75_random = quantile(random$C0_new_v9, c(0.5, 0.75)) 
cutoffs_C0_new_v9_75_100_random = quantile(random$C0_new_v9, c(0.75, 1)) 
cutoffs_C0_new_v9_custom1_random = quantile(random$C0_new_v9, c(0.7, 1)) 

cutoffs_C0_new_v9_0_25_tiling = quantile(tiling$C0_new_v9, c(0, 0.25)) 
cutoffs_C0_new_v9_25_50_tiling = quantile(tiling$C0_new_v9, c(0.25, 0.5)) 
cutoffs_C0_new_v9_50_75_tiling = quantile(tiling$C0_new_v9, c(0.5, 0.75)) 
cutoffs_C0_new_v9_75_100_tiling = quantile(tiling$C0_new_v9, c(0.75, 1)) 
cutoffs_C0_new_v9_custom1_tiling = quantile(tiling$C0_new_v9, c(0.7, 1)) 

cutoffs_C0_new_v9_0_25_chrv = quantile(chrv$C0_new_v9, c(0, 0.25)) 
cutoffs_C0_new_v9_25_50_chrv = quantile(chrv$C0_new_v9, c(0.25, 0.5)) 
cutoffs_C0_new_v9_50_75_chrv = quantile(chrv$C0_new_v9, c(0.5, 0.75)) 
cutoffs_C0_new_v9_75_100_chrv = quantile(chrv$C0_new_v9, c(0.75, 1)) 
cutoffs_C0_new_v9_custom1_chrv = quantile(chrv$C0_new_v9, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v10
cutoffs_C0_new_v10_0_25_nuc = quantile(nuc$C0_new_v10, c(0, 0.25)) 
cutoffs_C0_new_v10_25_50_nuc = quantile(nuc$C0_new_v10, c(0.25, 0.5)) 
cutoffs_C0_new_v10_50_75_nuc = quantile(nuc$C0_new_v10, c(0.5, 0.75)) 
cutoffs_C0_new_v10_75_100_nuc = quantile(nuc$C0_new_v10, c(0.75, 1)) 
cutoffs_C0_new_v10_custom1_nuc = quantile(nuc$C0_new_v10, c(0.7, 1)) 

cutoffs_C0_new_v10_0_25_random = quantile(random$C0_new_v10, c(0, 0.25)) 
cutoffs_C0_new_v10_25_50_random = quantile(random$C0_new_v10, c(0.25, 0.5)) 
cutoffs_C0_new_v10_50_75_random = quantile(random$C0_new_v10, c(0.5, 0.75)) 
cutoffs_C0_new_v10_75_100_random = quantile(random$C0_new_v10, c(0.75, 1)) 
cutoffs_C0_new_v10_custom1_random = quantile(random$C0_new_v10, c(0.7, 1)) 

cutoffs_C0_new_v10_0_25_tiling = quantile(tiling$C0_new_v10, c(0, 0.25)) 
cutoffs_C0_new_v10_25_50_tiling = quantile(tiling$C0_new_v10, c(0.25, 0.5)) 
cutoffs_C0_new_v10_50_75_tiling = quantile(tiling$C0_new_v10, c(0.5, 0.75)) 
cutoffs_C0_new_v10_75_100_tiling = quantile(tiling$C0_new_v10, c(0.75, 1)) 
cutoffs_C0_new_v10_custom1_tiling = quantile(tiling$C0_new_v10, c(0.7, 1)) 

cutoffs_C0_new_v10_0_25_chrv = quantile(chrv$C0_new_v10, c(0, 0.25)) 
cutoffs_C0_new_v10_25_50_chrv = quantile(chrv$C0_new_v10, c(0.25, 0.5)) 
cutoffs_C0_new_v10_50_75_chrv = quantile(chrv$C0_new_v10, c(0.5, 0.75)) 
cutoffs_C0_new_v10_75_100_chrv = quantile(chrv$C0_new_v10, c(0.75, 1)) 
cutoffs_C0_new_v10_custom1_chrv = quantile(chrv$C0_new_v10, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v11
cutoffs_C0_new_v11_0_25_nuc = quantile(nuc$C0_new_v11, c(0, 0.25)) 
cutoffs_C0_new_v11_25_50_nuc = quantile(nuc$C0_new_v11, c(0.25, 0.5)) 
cutoffs_C0_new_v11_50_75_nuc = quantile(nuc$C0_new_v11, c(0.5, 0.75)) 
cutoffs_C0_new_v11_75_100_nuc = quantile(nuc$C0_new_v11, c(0.75, 1)) 
cutoffs_C0_new_v11_custom1_nuc = quantile(nuc$C0_new_v11, c(0.7, 1)) 

cutoffs_C0_new_v11_0_25_random = quantile(random$C0_new_v11, c(0, 0.25)) 
cutoffs_C0_new_v11_25_50_random = quantile(random$C0_new_v11, c(0.25, 0.5)) 
cutoffs_C0_new_v11_50_75_random = quantile(random$C0_new_v11, c(0.5, 0.75)) 
cutoffs_C0_new_v11_75_100_random = quantile(random$C0_new_v11, c(0.75, 1)) 
cutoffs_C0_new_v11_custom1_random = quantile(random$C0_new_v11, c(0.7, 1)) 

cutoffs_C0_new_v11_0_25_tiling = quantile(tiling$C0_new_v11, c(0, 0.25)) 
cutoffs_C0_new_v11_25_50_tiling = quantile(tiling$C0_new_v11, c(0.25, 0.5)) 
cutoffs_C0_new_v11_50_75_tiling = quantile(tiling$C0_new_v11, c(0.5, 0.75)) 
cutoffs_C0_new_v11_75_100_tiling = quantile(tiling$C0_new_v11, c(0.75, 1)) 
cutoffs_C0_new_v11_custom1_tiling = quantile(tiling$C0_new_v11, c(0.7, 1)) 

cutoffs_C0_new_v11_0_25_chrv = quantile(chrv$C0_new_v11, c(0, 0.25)) 
cutoffs_C0_new_v11_25_50_chrv = quantile(chrv$C0_new_v11, c(0.25, 0.5)) 
cutoffs_C0_new_v11_50_75_chrv = quantile(chrv$C0_new_v11, c(0.5, 0.75)) 
cutoffs_C0_new_v11_75_100_chrv = quantile(chrv$C0_new_v11, c(0.75, 1)) 
cutoffs_C0_new_v11_custom1_chrv = quantile(chrv$C0_new_v11, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v12
cutoffs_C0_new_v12_0_25_nuc = quantile(nuc$C0_new_v12, c(0, 0.25)) 
cutoffs_C0_new_v12_25_50_nuc = quantile(nuc$C0_new_v12, c(0.25, 0.5)) 
cutoffs_C0_new_v12_50_75_nuc = quantile(nuc$C0_new_v12, c(0.5, 0.75)) 
cutoffs_C0_new_v12_75_100_nuc = quantile(nuc$C0_new_v12, c(0.75, 1)) 
cutoffs_C0_new_v12_custom1_nuc = quantile(nuc$C0_new_v12, c(0.7, 1)) 

cutoffs_C0_new_v12_0_25_random = quantile(random$C0_new_v12, c(0, 0.25)) 
cutoffs_C0_new_v12_25_50_random = quantile(random$C0_new_v12, c(0.25, 0.5)) 
cutoffs_C0_new_v12_50_75_random = quantile(random$C0_new_v12, c(0.5, 0.75)) 
cutoffs_C0_new_v12_75_100_random = quantile(random$C0_new_v12, c(0.75, 1)) 
cutoffs_C0_new_v12_custom1_random = quantile(random$C0_new_v12, c(0.7, 1)) 

cutoffs_C0_new_v12_0_25_tiling = quantile(tiling$C0_new_v12, c(0, 0.25)) 
cutoffs_C0_new_v12_25_50_tiling = quantile(tiling$C0_new_v12, c(0.25, 0.5)) 
cutoffs_C0_new_v12_50_75_tiling = quantile(tiling$C0_new_v12, c(0.5, 0.75)) 
cutoffs_C0_new_v12_75_100_tiling = quantile(tiling$C0_new_v12, c(0.75, 1)) 
cutoffs_C0_new_v12_custom1_tiling = quantile(tiling$C0_new_v12, c(0.7, 1)) 

cutoffs_C0_new_v12_0_25_chrv = quantile(chrv$C0_new_v12, c(0, 0.25)) 
cutoffs_C0_new_v12_25_50_chrv = quantile(chrv$C0_new_v12, c(0.25, 0.5)) 
cutoffs_C0_new_v12_50_75_chrv = quantile(chrv$C0_new_v12, c(0.5, 0.75)) 
cutoffs_C0_new_v12_75_100_chrv = quantile(chrv$C0_new_v12, c(0.75, 1)) 
cutoffs_C0_new_v12_custom1_chrv = quantile(chrv$C0_new_v12, c(0.7, 1)) 


# Find cutoffs for first-fourth quartiles + custom of C0_new_v13
cutoffs_C0_new_v13_0_25_nuc = quantile(nuc$C0_new_v13, c(0, 0.25)) 
cutoffs_C0_new_v13_25_50_nuc = quantile(nuc$C0_new_v13, c(0.25, 0.5)) 
cutoffs_C0_new_v13_50_75_nuc = quantile(nuc$C0_new_v13, c(0.5, 0.75)) 
cutoffs_C0_new_v13_75_100_nuc = quantile(nuc$C0_new_v13, c(0.75, 1)) 
cutoffs_C0_new_v13_custom1_nuc = quantile(nuc$C0_new_v13, c(0.7, 1)) 

cutoffs_C0_new_v13_0_25_random = quantile(random$C0_new_v13, c(0, 0.25)) 
cutoffs_C0_new_v13_25_50_random = quantile(random$C0_new_v13, c(0.25, 0.5)) 
cutoffs_C0_new_v13_50_75_random = quantile(random$C0_new_v13, c(0.5, 0.75)) 
cutoffs_C0_new_v13_75_100_random = quantile(random$C0_new_v13, c(0.75, 1)) 
cutoffs_C0_new_v13_custom1_random = quantile(random$C0_new_v13, c(0.7, 1)) 

cutoffs_C0_new_v13_0_25_tiling = quantile(tiling$C0_new_v13, c(0, 0.25)) 
cutoffs_C0_new_v13_25_50_tiling = quantile(tiling$C0_new_v13, c(0.25, 0.5)) 
cutoffs_C0_new_v13_50_75_tiling = quantile(tiling$C0_new_v13, c(0.5, 0.75)) 
cutoffs_C0_new_v13_75_100_tiling = quantile(tiling$C0_new_v13, c(0.75, 1)) 
cutoffs_C0_new_v13_custom1_tiling = quantile(tiling$C0_new_v13, c(0.7, 1)) 

cutoffs_C0_new_v13_0_25_chrv = quantile(chrv$C0_new_v13, c(0, 0.25)) 
cutoffs_C0_new_v13_25_50_chrv = quantile(chrv$C0_new_v13, c(0.25, 0.5)) 
cutoffs_C0_new_v13_50_75_chrv = quantile(chrv$C0_new_v13, c(0.5, 0.75)) 
cutoffs_C0_new_v13_75_100_chrv = quantile(chrv$C0_new_v13, c(0.75, 1)) 
cutoffs_C0_new_v13_custom1_chrv = quantile(chrv$C0_new_v13, c(0.7, 1)) 



# Divide each library into quartiles + custom for C26
nuc_C26_q1 = nuc %>%
  filter(C26 >= cutoffs_C26_0_25_nuc[1] & C26 <= cutoffs_C26_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
nuc_C26_q2 = nuc %>%
  filter(C26 >= cutoffs_C26_25_50_nuc[1] & C26 <= cutoffs_C26_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
nuc_C26_q3 = nuc %>%
  filter(C26 >= cutoffs_C26_50_75_nuc[1] & C26 <= cutoffs_C26_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
nuc_C26_q4 = nuc %>%
  filter(C26 >= cutoffs_C26_75_100_nuc[1] & C26 <= cutoffs_C26_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
nuc_C26_custom1 = nuc %>%
  filter(C26 >= cutoffs_C26_custom1_nuc[1] & C26 <= cutoffs_C26_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)

random_C26_q1 = random %>%
  filter(C26 >= cutoffs_C26_0_25_random[1] & C26 <= cutoffs_C26_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
random_C26_q2 = random %>%
  filter(C26 >= cutoffs_C26_25_50_random[1] & C26 <= cutoffs_C26_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
random_C26_q3 = random %>%
  filter(C26 >= cutoffs_C26_50_75_random[1] & C26 <= cutoffs_C26_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
random_C26_q4 = random %>%
  filter(C26 >= cutoffs_C26_75_100_random[1] & C26 <= cutoffs_C26_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
random_C26_custom1 = random %>%
  filter(C26 >= cutoffs_C26_custom1_random[1] & C26 <= cutoffs_C26_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)

tiling_C26_q1 = tiling %>%
  filter(C26 >= cutoffs_C26_0_25_tiling[1] & C26 <= cutoffs_C26_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
tiling_C26_q2 = tiling %>%
  filter(C26 >= cutoffs_C26_25_50_tiling[1] & C26 <= cutoffs_C26_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
tiling_C26_q3 = tiling %>%
  filter(C26 >= cutoffs_C26_50_75_tiling[1] & C26 <= cutoffs_C26_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
tiling_C26_q4 = tiling %>%
  filter(C26 >= cutoffs_C26_75_100_tiling[1] & C26 <= cutoffs_C26_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
tiling_C26_custom1 = tiling %>%
  filter(C26 >= cutoffs_C26_custom1_tiling[1] & C26 <= cutoffs_C26_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)

chrv_C26_q1 = chrv %>%
  filter(C26 >= cutoffs_C26_0_25_chrv[1] & C26 <= cutoffs_C26_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
chrv_C26_q2 = chrv %>%
  filter(C26 >= cutoffs_C26_25_50_chrv[1] & C26 <= cutoffs_C26_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
chrv_C26_q3 = chrv %>%
  filter(C26 >= cutoffs_C26_50_75_chrv[1] & C26 <= cutoffs_C26_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
chrv_C26_q4 = chrv %>%
  filter(C26 >= cutoffs_C26_75_100_chrv[1] & C26 <= cutoffs_C26_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)
chrv_C26_custom1 = chrv %>%
  filter(C26 >= cutoffs_C26_custom1_chrv[1] & C26 <= cutoffs_C26_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C26)



# Divide each library into quartiles + custom for C29
nuc_C29_q1 = nuc %>%
  filter(C29 >= cutoffs_C29_0_25_nuc[1] & C29 <= cutoffs_C29_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
nuc_C29_q2 = nuc %>%
  filter(C29 >= cutoffs_C29_25_50_nuc[1] & C29 <= cutoffs_C29_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
nuc_C29_q3 = nuc %>%
  filter(C29 >= cutoffs_C29_50_75_nuc[1] & C29 <= cutoffs_C29_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
nuc_C29_q4 = nuc %>%
  filter(C29 >= cutoffs_C29_75_100_nuc[1] & C29 <= cutoffs_C29_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
nuc_C29_custom1 = nuc %>%
  filter(C29 >= cutoffs_C29_custom1_nuc[1] & C29 <= cutoffs_C29_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)

random_C29_q1 = random %>%
  filter(C29 >= cutoffs_C29_0_25_random[1] & C29 <= cutoffs_C29_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
random_C29_q2 = random %>%
  filter(C29 >= cutoffs_C29_25_50_random[1] & C29 <= cutoffs_C29_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
random_C29_q3 = random %>%
  filter(C29 >= cutoffs_C29_50_75_random[1] & C29 <= cutoffs_C29_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
random_C29_q4 = random %>%
  filter(C29 >= cutoffs_C29_75_100_random[1] & C29 <= cutoffs_C29_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
random_C29_custom1 = random %>%
  filter(C29 >= cutoffs_C29_custom1_random[1] & C29 <= cutoffs_C29_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)

tiling_C29_q1 = tiling %>%
  filter(C29 >= cutoffs_C29_0_25_tiling[1] & C29 <= cutoffs_C29_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
tiling_C29_q2 = tiling %>%
  filter(C29 >= cutoffs_C29_25_50_tiling[1] & C29 <= cutoffs_C29_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
tiling_C29_q3 = tiling %>%
  filter(C29 >= cutoffs_C29_50_75_tiling[1] & C29 <= cutoffs_C29_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
tiling_C29_q4 = tiling %>%
  filter(C29 >= cutoffs_C29_75_100_tiling[1] & C29 <= cutoffs_C29_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
tiling_C29_custom1 = tiling %>%
  filter(C29 >= cutoffs_C29_custom1_tiling[1] & C29 <= cutoffs_C29_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)

chrv_C29_q1 = chrv %>%
  filter(C29 >= cutoffs_C29_0_25_chrv[1] & C29 <= cutoffs_C29_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
chrv_C29_q2 = chrv %>%
  filter(C29 >= cutoffs_C29_25_50_chrv[1] & C29 <= cutoffs_C29_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
chrv_C29_q3 = chrv %>%
  filter(C29 >= cutoffs_C29_50_75_chrv[1] & C29 <= cutoffs_C29_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
chrv_C29_q4 = chrv %>%
  filter(C29 >= cutoffs_C29_75_100_chrv[1] & C29 <= cutoffs_C29_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)
chrv_C29_custom1 = chrv %>%
  filter(C29 >= cutoffs_C29_custom1_chrv[1] & C29 <= cutoffs_C29_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C29)



# Divide each library into quartiles + custom for C31
nuc_C31_q1 = nuc %>%
  filter(C31 >= cutoffs_C31_0_25_nuc[1] & C31 <= cutoffs_C31_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
nuc_C31_q2 = nuc %>%
  filter(C31 >= cutoffs_C31_25_50_nuc[1] & C31 <= cutoffs_C31_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
nuc_C31_q3 = nuc %>%
  filter(C31 >= cutoffs_C31_50_75_nuc[1] & C31 <= cutoffs_C31_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
nuc_C31_q4 = nuc %>%
  filter(C31 >= cutoffs_C31_75_100_nuc[1] & C31 <= cutoffs_C31_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
nuc_C31_custom1 = nuc %>%
  filter(C31 >= cutoffs_C31_custom1_nuc[1] & C31 <= cutoffs_C31_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)

random_C31_q1 = random %>%
  filter(C31 >= cutoffs_C31_0_25_random[1] & C31 <= cutoffs_C31_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
random_C31_q2 = random %>%
  filter(C31 >= cutoffs_C31_25_50_random[1] & C31 <= cutoffs_C31_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
random_C31_q3 = random %>%
  filter(C31 >= cutoffs_C31_50_75_random[1] & C31 <= cutoffs_C31_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
random_C31_q4 = random %>%
  filter(C31 >= cutoffs_C31_75_100_random[1] & C31 <= cutoffs_C31_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
random_C31_custom1 = random %>%
  filter(C31 >= cutoffs_C31_custom1_random[1] & C31 <= cutoffs_C31_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)

tiling_C31_q1 = tiling %>%
  filter(C31 >= cutoffs_C31_0_25_tiling[1] & C31 <= cutoffs_C31_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
tiling_C31_q2 = tiling %>%
  filter(C31 >= cutoffs_C31_25_50_tiling[1] & C31 <= cutoffs_C31_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
tiling_C31_q3 = tiling %>%
  filter(C31 >= cutoffs_C31_50_75_tiling[1] & C31 <= cutoffs_C31_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
tiling_C31_q4 = tiling %>%
  filter(C31 >= cutoffs_C31_75_100_tiling[1] & C31 <= cutoffs_C31_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
tiling_C31_custom1 = tiling %>%
  filter(C31 >= cutoffs_C31_custom1_tiling[1] & C31 <= cutoffs_C31_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)

chrv_C31_q1 = chrv %>%
  filter(C31 >= cutoffs_C31_0_25_chrv[1] & C31 <= cutoffs_C31_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
chrv_C31_q2 = chrv %>%
  filter(C31 >= cutoffs_C31_25_50_chrv[1] & C31 <= cutoffs_C31_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
chrv_C31_q3 = chrv %>%
  filter(C31 >= cutoffs_C31_50_75_chrv[1] & C31 <= cutoffs_C31_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
chrv_C31_q4 = chrv %>%
  filter(C31 >= cutoffs_C31_75_100_chrv[1] & C31 <= cutoffs_C31_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)
chrv_C31_custom1 = chrv %>%
  filter(C31 >= cutoffs_C31_custom1_chrv[1] & C31 <= cutoffs_C31_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C31)


# Divide each library into quartiles + custom for C0
nuc_C0_q1 = nuc %>%
  filter(C0 >= cutoffs_C0_0_25_nuc[1] & C0 <= cutoffs_C0_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
nuc_C0_q2 = nuc %>%
  filter(C0 >= cutoffs_C0_25_50_nuc[1] & C0 <= cutoffs_C0_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
nuc_C0_q3 = nuc %>%
  filter(C0 >= cutoffs_C0_50_75_nuc[1] & C0 <= cutoffs_C0_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
nuc_C0_q4 = nuc %>%
  filter(C0 >= cutoffs_C0_75_100_nuc[1] & C0 <= cutoffs_C0_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
nuc_C0_custom1 = nuc %>%
  filter(C0 >= cutoffs_C0_custom1_nuc[1] & C0 <= cutoffs_C0_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)

random_C0_q1 = random %>%
  filter(C0 >= cutoffs_C0_0_25_random[1] & C0 <= cutoffs_C0_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
random_C0_q2 = random %>%
  filter(C0 >= cutoffs_C0_25_50_random[1] & C0 <= cutoffs_C0_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
random_C0_q3 = random %>%
  filter(C0 >= cutoffs_C0_50_75_random[1] & C0 <= cutoffs_C0_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
random_C0_q4 = random %>%
  filter(C0 >= cutoffs_C0_75_100_random[1] & C0 <= cutoffs_C0_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
random_C0_custom1 = random %>%
  filter(C0 >= cutoffs_C0_custom1_random[1] & C0 <= cutoffs_C0_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)

tiling_C0_q1 = tiling %>%
  filter(C0 >= cutoffs_C0_0_25_tiling[1] & C0 <= cutoffs_C0_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
tiling_C0_q2 = tiling %>%
  filter(C0 >= cutoffs_C0_25_50_tiling[1] & C0 <= cutoffs_C0_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
tiling_C0_q3 = tiling %>%
  filter(C0 >= cutoffs_C0_50_75_tiling[1] & C0 <= cutoffs_C0_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
tiling_C0_q4 = tiling %>%
  filter(C0 >= cutoffs_C0_75_100_tiling[1] & C0 <= cutoffs_C0_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
tiling_C0_custom1 = tiling %>%
  filter(C0 >= cutoffs_C0_custom1_tiling[1] & C0 <= cutoffs_C0_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)

chrv_C0_q1 = chrv %>%
  filter(C0 >= cutoffs_C0_0_25_chrv[1] & C0 <= cutoffs_C0_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
chrv_C0_q2 = chrv %>%
  filter(C0 >= cutoffs_C0_25_50_chrv[1] & C0 <= cutoffs_C0_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
chrv_C0_q3 = chrv %>%
  filter(C0 >= cutoffs_C0_50_75_chrv[1] & C0 <= cutoffs_C0_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
chrv_C0_q4 = chrv %>%
  filter(C0 >= cutoffs_C0_75_100_chrv[1] & C0 <= cutoffs_C0_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
chrv_C0_custom1 = chrv %>%
  filter(C0 >= cutoffs_C0_custom1_chrv[1] & C0 <= cutoffs_C0_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)


# Divide each library into quartiles + custom for C0_new_v1
nuc_C0_new_v1_q1 = nuc %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_0_25_nuc[1] & C0_new_v1 <= cutoffs_C0_new_v1_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
nuc_C0_new_v1_q2 = nuc %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_25_50_nuc[1] & C0_new_v1 <= cutoffs_C0_new_v1_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
nuc_C0_new_v1_q3 = nuc %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_50_75_nuc[1] & C0_new_v1 <= cutoffs_C0_new_v1_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
nuc_C0_new_v1_q4 = nuc %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_75_100_nuc[1] & C0_new_v1 <= cutoffs_C0_new_v1_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
nuc_C0_new_v1_custom1 = nuc %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_custom1_nuc[1] & C0_new_v1 <= cutoffs_C0_new_v1_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)

random_C0_new_v1_q1 = random %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_0_25_random[1] & C0_new_v1 <= cutoffs_C0_new_v1_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
random_C0_new_v1_q2 = random %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_25_50_random[1] & C0_new_v1 <= cutoffs_C0_new_v1_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
random_C0_new_v1_q3 = random %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_50_75_random[1] & C0_new_v1 <= cutoffs_C0_new_v1_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
random_C0_new_v1_q4 = random %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_75_100_random[1] & C0_new_v1 <= cutoffs_C0_new_v1_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
random_C0_new_v1_custom1 = random %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_custom1_random[1] & C0_new_v1 <= cutoffs_C0_new_v1_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)

tiling_C0_new_v1_q1 = tiling %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_0_25_tiling[1] & C0_new_v1 <= cutoffs_C0_new_v1_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
tiling_C0_new_v1_q2 = tiling %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_25_50_tiling[1] & C0_new_v1 <= cutoffs_C0_new_v1_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
tiling_C0_new_v1_q3 = tiling %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_50_75_tiling[1] & C0_new_v1 <= cutoffs_C0_new_v1_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
tiling_C0_new_v1_q4 = tiling %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_75_100_tiling[1] & C0_new_v1 <= cutoffs_C0_new_v1_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
tiling_C0_new_v1_custom1 = tiling %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_custom1_tiling[1] & C0_new_v1 <= cutoffs_C0_new_v1_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)

chrv_C0_new_v1_q1 = chrv %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_0_25_chrv[1] & C0_new_v1 <= cutoffs_C0_new_v1_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
chrv_C0_new_v1_q2 = chrv %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_25_50_chrv[1] & C0_new_v1 <= cutoffs_C0_new_v1_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
chrv_C0_new_v1_q3 = chrv %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_50_75_chrv[1] & C0_new_v1 <= cutoffs_C0_new_v1_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
chrv_C0_new_v1_q4 = chrv %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_75_100_chrv[1] & C0_new_v1 <= cutoffs_C0_new_v1_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)
chrv_C0_new_v1_custom1 = chrv %>%
  filter(C0_new_v1 >= cutoffs_C0_new_v1_custom1_chrv[1] & C0_new_v1 <= cutoffs_C0_new_v1_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v1)


# Divide each library into quartiles + custom for C0_new_v2
nuc_C0_new_v2_q1 = nuc %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_0_25_nuc[1] & C0_new_v2 <= cutoffs_C0_new_v2_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
nuc_C0_new_v2_q2 = nuc %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_25_50_nuc[1] & C0_new_v2 <= cutoffs_C0_new_v2_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
nuc_C0_new_v2_q3 = nuc %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_50_75_nuc[1] & C0_new_v2 <= cutoffs_C0_new_v2_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
nuc_C0_new_v2_q4 = nuc %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_75_100_nuc[1] & C0_new_v2 <= cutoffs_C0_new_v2_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
nuc_C0_new_v2_custom1 = nuc %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_custom1_nuc[1] & C0_new_v2 <= cutoffs_C0_new_v2_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)

random_C0_new_v2_q1 = random %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_0_25_random[1] & C0_new_v2 <= cutoffs_C0_new_v2_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
random_C0_new_v2_q2 = random %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_25_50_random[1] & C0_new_v2 <= cutoffs_C0_new_v2_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
random_C0_new_v2_q3 = random %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_50_75_random[1] & C0_new_v2 <= cutoffs_C0_new_v2_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
random_C0_new_v2_q4 = random %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_75_100_random[1] & C0_new_v2 <= cutoffs_C0_new_v2_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
random_C0_new_v2_custom1 = random %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_custom1_random[1] & C0_new_v2 <= cutoffs_C0_new_v2_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)

tiling_C0_new_v2_q1 = tiling %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_0_25_tiling[1] & C0_new_v2 <= cutoffs_C0_new_v2_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
tiling_C0_new_v2_q2 = tiling %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_25_50_tiling[1] & C0_new_v2 <= cutoffs_C0_new_v2_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
tiling_C0_new_v2_q3 = tiling %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_50_75_tiling[1] & C0_new_v2 <= cutoffs_C0_new_v2_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
tiling_C0_new_v2_q4 = tiling %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_75_100_tiling[1] & C0_new_v2 <= cutoffs_C0_new_v2_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
tiling_C0_new_v2_custom1 = tiling %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_custom1_tiling[1] & C0_new_v2 <= cutoffs_C0_new_v2_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)

chrv_C0_new_v2_q1 = chrv %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_0_25_chrv[1] & C0_new_v2 <= cutoffs_C0_new_v2_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
chrv_C0_new_v2_q2 = chrv %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_25_50_chrv[1] & C0_new_v2 <= cutoffs_C0_new_v2_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
chrv_C0_new_v2_q3 = chrv %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_50_75_chrv[1] & C0_new_v2 <= cutoffs_C0_new_v2_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
chrv_C0_new_v2_q4 = chrv %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_75_100_chrv[1] & C0_new_v2 <= cutoffs_C0_new_v2_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)
chrv_C0_new_v2_custom1 = chrv %>%
  filter(C0_new_v2 >= cutoffs_C0_new_v2_custom1_chrv[1] & C0_new_v2 <= cutoffs_C0_new_v2_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v2)


# Divide each library into quartiles + custom for C0_new_v3
nuc_C0_new_v3_q1 = nuc %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_0_25_nuc[1] & C0_new_v3 <= cutoffs_C0_new_v3_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
nuc_C0_new_v3_q2 = nuc %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_25_50_nuc[1] & C0_new_v3 <= cutoffs_C0_new_v3_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
nuc_C0_new_v3_q3 = nuc %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_50_75_nuc[1] & C0_new_v3 <= cutoffs_C0_new_v3_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
nuc_C0_new_v3_q4 = nuc %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_75_100_nuc[1] & C0_new_v3 <= cutoffs_C0_new_v3_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
nuc_C0_new_v3_custom1 = nuc %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_custom1_nuc[1] & C0_new_v3 <= cutoffs_C0_new_v3_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)

random_C0_new_v3_q1 = random %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_0_25_random[1] & C0_new_v3 <= cutoffs_C0_new_v3_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
random_C0_new_v3_q2 = random %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_25_50_random[1] & C0_new_v3 <= cutoffs_C0_new_v3_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
random_C0_new_v3_q3 = random %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_50_75_random[1] & C0_new_v3 <= cutoffs_C0_new_v3_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
random_C0_new_v3_q4 = random %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_75_100_random[1] & C0_new_v3 <= cutoffs_C0_new_v3_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
random_C0_new_v3_custom1 = random %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_custom1_random[1] & C0_new_v3 <= cutoffs_C0_new_v3_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)

tiling_C0_new_v3_q1 = tiling %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_0_25_tiling[1] & C0_new_v3 <= cutoffs_C0_new_v3_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
tiling_C0_new_v3_q2 = tiling %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_25_50_tiling[1] & C0_new_v3 <= cutoffs_C0_new_v3_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
tiling_C0_new_v3_q3 = tiling %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_50_75_tiling[1] & C0_new_v3 <= cutoffs_C0_new_v3_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
tiling_C0_new_v3_q4 = tiling %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_75_100_tiling[1] & C0_new_v3 <= cutoffs_C0_new_v3_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
tiling_C0_new_v3_custom1 = tiling %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_custom1_tiling[1] & C0_new_v3 <= cutoffs_C0_new_v3_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)

chrv_C0_new_v3_q1 = chrv %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_0_25_chrv[1] & C0_new_v3 <= cutoffs_C0_new_v3_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
chrv_C0_new_v3_q2 = chrv %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_25_50_chrv[1] & C0_new_v3 <= cutoffs_C0_new_v3_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
chrv_C0_new_v3_q3 = chrv %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_50_75_chrv[1] & C0_new_v3 <= cutoffs_C0_new_v3_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
chrv_C0_new_v3_q4 = chrv %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_75_100_chrv[1] & C0_new_v3 <= cutoffs_C0_new_v3_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)
chrv_C0_new_v3_custom1 = chrv %>%
  filter(C0_new_v3 >= cutoffs_C0_new_v3_custom1_chrv[1] & C0_new_v3 <= cutoffs_C0_new_v3_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v3)


# Divide each library into quartiles + custom for C0_new_v4
nuc_C0_new_v4_q1 = nuc %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_0_25_nuc[1] & C0_new_v4 <= cutoffs_C0_new_v4_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
nuc_C0_new_v4_q2 = nuc %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_25_50_nuc[1] & C0_new_v4 <= cutoffs_C0_new_v4_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
nuc_C0_new_v4_q3 = nuc %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_50_75_nuc[1] & C0_new_v4 <= cutoffs_C0_new_v4_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
nuc_C0_new_v4_q4 = nuc %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_75_100_nuc[1] & C0_new_v4 <= cutoffs_C0_new_v4_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
nuc_C0_new_v4_custom1 = nuc %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_custom1_nuc[1] & C0_new_v4 <= cutoffs_C0_new_v4_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)

random_C0_new_v4_q1 = random %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_0_25_random[1] & C0_new_v4 <= cutoffs_C0_new_v4_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
random_C0_new_v4_q2 = random %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_25_50_random[1] & C0_new_v4 <= cutoffs_C0_new_v4_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
random_C0_new_v4_q3 = random %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_50_75_random[1] & C0_new_v4 <= cutoffs_C0_new_v4_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
random_C0_new_v4_q4 = random %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_75_100_random[1] & C0_new_v4 <= cutoffs_C0_new_v4_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
random_C0_new_v4_custom1 = random %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_custom1_random[1] & C0_new_v4 <= cutoffs_C0_new_v4_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)

tiling_C0_new_v4_q1 = tiling %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_0_25_tiling[1] & C0_new_v4 <= cutoffs_C0_new_v4_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
tiling_C0_new_v4_q2 = tiling %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_25_50_tiling[1] & C0_new_v4 <= cutoffs_C0_new_v4_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
tiling_C0_new_v4_q3 = tiling %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_50_75_tiling[1] & C0_new_v4 <= cutoffs_C0_new_v4_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
tiling_C0_new_v4_q4 = tiling %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_75_100_tiling[1] & C0_new_v4 <= cutoffs_C0_new_v4_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
tiling_C0_new_v4_custom1 = tiling %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_custom1_tiling[1] & C0_new_v4 <= cutoffs_C0_new_v4_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)

chrv_C0_new_v4_q1 = chrv %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_0_25_chrv[1] & C0_new_v4 <= cutoffs_C0_new_v4_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
chrv_C0_new_v4_q2 = chrv %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_25_50_chrv[1] & C0_new_v4 <= cutoffs_C0_new_v4_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
chrv_C0_new_v4_q3 = chrv %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_50_75_chrv[1] & C0_new_v4 <= cutoffs_C0_new_v4_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
chrv_C0_new_v4_q4 = chrv %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_75_100_chrv[1] & C0_new_v4 <= cutoffs_C0_new_v4_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)
chrv_C0_new_v4_custom1 = chrv %>%
  filter(C0_new_v4 >= cutoffs_C0_new_v4_custom1_chrv[1] & C0_new_v4 <= cutoffs_C0_new_v4_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v4)


# Divide each library into quartiles + custom for C0_new_v5
nuc_C0_new_v5_q1 = nuc %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_0_25_nuc[1] & C0_new_v5 <= cutoffs_C0_new_v5_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
nuc_C0_new_v5_q2 = nuc %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_25_50_nuc[1] & C0_new_v5 <= cutoffs_C0_new_v5_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
nuc_C0_new_v5_q3 = nuc %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_50_75_nuc[1] & C0_new_v5 <= cutoffs_C0_new_v5_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
nuc_C0_new_v5_q4 = nuc %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_75_100_nuc[1] & C0_new_v5 <= cutoffs_C0_new_v5_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
nuc_C0_new_v5_custom1 = nuc %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_custom1_nuc[1] & C0_new_v5 <= cutoffs_C0_new_v5_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)

random_C0_new_v5_q1 = random %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_0_25_random[1] & C0_new_v5 <= cutoffs_C0_new_v5_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
random_C0_new_v5_q2 = random %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_25_50_random[1] & C0_new_v5 <= cutoffs_C0_new_v5_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
random_C0_new_v5_q3 = random %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_50_75_random[1] & C0_new_v5 <= cutoffs_C0_new_v5_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
random_C0_new_v5_q4 = random %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_75_100_random[1] & C0_new_v5 <= cutoffs_C0_new_v5_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
random_C0_new_v5_custom1 = random %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_custom1_random[1] & C0_new_v5 <= cutoffs_C0_new_v5_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)

tiling_C0_new_v5_q1 = tiling %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_0_25_tiling[1] & C0_new_v5 <= cutoffs_C0_new_v5_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
tiling_C0_new_v5_q2 = tiling %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_25_50_tiling[1] & C0_new_v5 <= cutoffs_C0_new_v5_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
tiling_C0_new_v5_q3 = tiling %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_50_75_tiling[1] & C0_new_v5 <= cutoffs_C0_new_v5_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
tiling_C0_new_v5_q4 = tiling %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_75_100_tiling[1] & C0_new_v5 <= cutoffs_C0_new_v5_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
tiling_C0_new_v5_custom1 = tiling %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_custom1_tiling[1] & C0_new_v5 <= cutoffs_C0_new_v5_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)

chrv_C0_new_v5_q1 = chrv %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_0_25_chrv[1] & C0_new_v5 <= cutoffs_C0_new_v5_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
chrv_C0_new_v5_q2 = chrv %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_25_50_chrv[1] & C0_new_v5 <= cutoffs_C0_new_v5_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
chrv_C0_new_v5_q3 = chrv %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_50_75_chrv[1] & C0_new_v5 <= cutoffs_C0_new_v5_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
chrv_C0_new_v5_q4 = chrv %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_75_100_chrv[1] & C0_new_v5 <= cutoffs_C0_new_v5_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)
chrv_C0_new_v5_custom1 = chrv %>%
  filter(C0_new_v5 >= cutoffs_C0_new_v5_custom1_chrv[1] & C0_new_v5 <= cutoffs_C0_new_v5_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v5)


# Divide each library into quartiles + custom for C0_new_v6
nuc_C0_new_v6_q1 = nuc %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_0_25_nuc[1] & C0_new_v6 <= cutoffs_C0_new_v6_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
nuc_C0_new_v6_q2 = nuc %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_25_50_nuc[1] & C0_new_v6 <= cutoffs_C0_new_v6_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
nuc_C0_new_v6_q3 = nuc %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_50_75_nuc[1] & C0_new_v6 <= cutoffs_C0_new_v6_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
nuc_C0_new_v6_q4 = nuc %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_75_100_nuc[1] & C0_new_v6 <= cutoffs_C0_new_v6_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
nuc_C0_new_v6_custom1 = nuc %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_custom1_nuc[1] & C0_new_v6 <= cutoffs_C0_new_v6_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)

random_C0_new_v6_q1 = random %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_0_25_random[1] & C0_new_v6 <= cutoffs_C0_new_v6_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
random_C0_new_v6_q2 = random %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_25_50_random[1] & C0_new_v6 <= cutoffs_C0_new_v6_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
random_C0_new_v6_q3 = random %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_50_75_random[1] & C0_new_v6 <= cutoffs_C0_new_v6_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
random_C0_new_v6_q4 = random %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_75_100_random[1] & C0_new_v6 <= cutoffs_C0_new_v6_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
random_C0_new_v6_custom1 = random %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_custom1_random[1] & C0_new_v6 <= cutoffs_C0_new_v6_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)

tiling_C0_new_v6_q1 = tiling %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_0_25_tiling[1] & C0_new_v6 <= cutoffs_C0_new_v6_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
tiling_C0_new_v6_q2 = tiling %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_25_50_tiling[1] & C0_new_v6 <= cutoffs_C0_new_v6_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
tiling_C0_new_v6_q3 = tiling %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_50_75_tiling[1] & C0_new_v6 <= cutoffs_C0_new_v6_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
tiling_C0_new_v6_q4 = tiling %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_75_100_tiling[1] & C0_new_v6 <= cutoffs_C0_new_v6_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
tiling_C0_new_v6_custom1 = tiling %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_custom1_tiling[1] & C0_new_v6 <= cutoffs_C0_new_v6_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)

chrv_C0_new_v6_q1 = chrv %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_0_25_chrv[1] & C0_new_v6 <= cutoffs_C0_new_v6_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
chrv_C0_new_v6_q2 = chrv %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_25_50_chrv[1] & C0_new_v6 <= cutoffs_C0_new_v6_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
chrv_C0_new_v6_q3 = chrv %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_50_75_chrv[1] & C0_new_v6 <= cutoffs_C0_new_v6_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
chrv_C0_new_v6_q4 = chrv %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_75_100_chrv[1] & C0_new_v6 <= cutoffs_C0_new_v6_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)
chrv_C0_new_v6_custom1 = chrv %>%
  filter(C0_new_v6 >= cutoffs_C0_new_v6_custom1_chrv[1] & C0_new_v6 <= cutoffs_C0_new_v6_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v6)


# Divide each library into quartiles + custom for C0_new_v7
nuc_C0_new_v7_q1 = nuc %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_0_25_nuc[1] & C0_new_v7 <= cutoffs_C0_new_v7_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
nuc_C0_new_v7_q2 = nuc %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_25_50_nuc[1] & C0_new_v7 <= cutoffs_C0_new_v7_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
nuc_C0_new_v7_q3 = nuc %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_50_75_nuc[1] & C0_new_v7 <= cutoffs_C0_new_v7_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
nuc_C0_new_v7_q4 = nuc %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_75_100_nuc[1] & C0_new_v7 <= cutoffs_C0_new_v7_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
nuc_C0_new_v7_custom1 = nuc %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_custom1_nuc[1] & C0_new_v7 <= cutoffs_C0_new_v7_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)

random_C0_new_v7_q1 = random %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_0_25_random[1] & C0_new_v7 <= cutoffs_C0_new_v7_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
random_C0_new_v7_q2 = random %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_25_50_random[1] & C0_new_v7 <= cutoffs_C0_new_v7_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
random_C0_new_v7_q3 = random %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_50_75_random[1] & C0_new_v7 <= cutoffs_C0_new_v7_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
random_C0_new_v7_q4 = random %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_75_100_random[1] & C0_new_v7 <= cutoffs_C0_new_v7_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
random_C0_new_v7_custom1 = random %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_custom1_random[1] & C0_new_v7 <= cutoffs_C0_new_v7_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)

tiling_C0_new_v7_q1 = tiling %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_0_25_tiling[1] & C0_new_v7 <= cutoffs_C0_new_v7_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
tiling_C0_new_v7_q2 = tiling %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_25_50_tiling[1] & C0_new_v7 <= cutoffs_C0_new_v7_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
tiling_C0_new_v7_q3 = tiling %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_50_75_tiling[1] & C0_new_v7 <= cutoffs_C0_new_v7_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
tiling_C0_new_v7_q4 = tiling %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_75_100_tiling[1] & C0_new_v7 <= cutoffs_C0_new_v7_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
tiling_C0_new_v7_custom1 = tiling %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_custom1_tiling[1] & C0_new_v7 <= cutoffs_C0_new_v7_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)

chrv_C0_new_v7_q1 = chrv %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_0_25_chrv[1] & C0_new_v7 <= cutoffs_C0_new_v7_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
chrv_C0_new_v7_q2 = chrv %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_25_50_chrv[1] & C0_new_v7 <= cutoffs_C0_new_v7_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
chrv_C0_new_v7_q3 = chrv %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_50_75_chrv[1] & C0_new_v7 <= cutoffs_C0_new_v7_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
chrv_C0_new_v7_q4 = chrv %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_75_100_chrv[1] & C0_new_v7 <= cutoffs_C0_new_v7_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)
chrv_C0_new_v7_custom1 = chrv %>%
  filter(C0_new_v7 >= cutoffs_C0_new_v7_custom1_chrv[1] & C0_new_v7 <= cutoffs_C0_new_v7_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v7)


# Divide each library into quartiles + custom for C0_new_v8
nuc_C0_new_v8_q1 = nuc %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_0_25_nuc[1] & C0_new_v8 <= cutoffs_C0_new_v8_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
nuc_C0_new_v8_q2 = nuc %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_25_50_nuc[1] & C0_new_v8 <= cutoffs_C0_new_v8_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
nuc_C0_new_v8_q3 = nuc %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_50_75_nuc[1] & C0_new_v8 <= cutoffs_C0_new_v8_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
nuc_C0_new_v8_q4 = nuc %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_75_100_nuc[1] & C0_new_v8 <= cutoffs_C0_new_v8_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
nuc_C0_new_v8_custom1 = nuc %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_custom1_nuc[1] & C0_new_v8 <= cutoffs_C0_new_v8_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)

random_C0_new_v8_q1 = random %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_0_25_random[1] & C0_new_v8 <= cutoffs_C0_new_v8_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
random_C0_new_v8_q2 = random %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_25_50_random[1] & C0_new_v8 <= cutoffs_C0_new_v8_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
random_C0_new_v8_q3 = random %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_50_75_random[1] & C0_new_v8 <= cutoffs_C0_new_v8_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
random_C0_new_v8_q4 = random %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_75_100_random[1] & C0_new_v8 <= cutoffs_C0_new_v8_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
random_C0_new_v8_custom1 = random %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_custom1_random[1] & C0_new_v8 <= cutoffs_C0_new_v8_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)

tiling_C0_new_v8_q1 = tiling %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_0_25_tiling[1] & C0_new_v8 <= cutoffs_C0_new_v8_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
tiling_C0_new_v8_q2 = tiling %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_25_50_tiling[1] & C0_new_v8 <= cutoffs_C0_new_v8_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
tiling_C0_new_v8_q3 = tiling %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_50_75_tiling[1] & C0_new_v8 <= cutoffs_C0_new_v8_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
tiling_C0_new_v8_q4 = tiling %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_75_100_tiling[1] & C0_new_v8 <= cutoffs_C0_new_v8_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
tiling_C0_new_v8_custom1 = tiling %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_custom1_tiling[1] & C0_new_v8 <= cutoffs_C0_new_v8_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)

chrv_C0_new_v8_q1 = chrv %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_0_25_chrv[1] & C0_new_v8 <= cutoffs_C0_new_v8_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
chrv_C0_new_v8_q2 = chrv %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_25_50_chrv[1] & C0_new_v8 <= cutoffs_C0_new_v8_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
chrv_C0_new_v8_q3 = chrv %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_50_75_chrv[1] & C0_new_v8 <= cutoffs_C0_new_v8_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
chrv_C0_new_v8_q4 = chrv %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_75_100_chrv[1] & C0_new_v8 <= cutoffs_C0_new_v8_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)
chrv_C0_new_v8_custom1 = chrv %>%
  filter(C0_new_v8 >= cutoffs_C0_new_v8_custom1_chrv[1] & C0_new_v8 <= cutoffs_C0_new_v8_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v8)


# Divide each library into quartiles + custom for C0_new_v9
nuc_C0_new_v9_q1 = nuc %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_0_25_nuc[1] & C0_new_v9 <= cutoffs_C0_new_v9_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
nuc_C0_new_v9_q2 = nuc %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_25_50_nuc[1] & C0_new_v9 <= cutoffs_C0_new_v9_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
nuc_C0_new_v9_q3 = nuc %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_50_75_nuc[1] & C0_new_v9 <= cutoffs_C0_new_v9_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
nuc_C0_new_v9_q4 = nuc %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_75_100_nuc[1] & C0_new_v9 <= cutoffs_C0_new_v9_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
nuc_C0_new_v9_custom1 = nuc %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_custom1_nuc[1] & C0_new_v9 <= cutoffs_C0_new_v9_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)

random_C0_new_v9_q1 = random %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_0_25_random[1] & C0_new_v9 <= cutoffs_C0_new_v9_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
random_C0_new_v9_q2 = random %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_25_50_random[1] & C0_new_v9 <= cutoffs_C0_new_v9_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
random_C0_new_v9_q3 = random %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_50_75_random[1] & C0_new_v9 <= cutoffs_C0_new_v9_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
random_C0_new_v9_q4 = random %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_75_100_random[1] & C0_new_v9 <= cutoffs_C0_new_v9_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
random_C0_new_v9_custom1 = random %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_custom1_random[1] & C0_new_v9 <= cutoffs_C0_new_v9_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)

tiling_C0_new_v9_q1 = tiling %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_0_25_tiling[1] & C0_new_v9 <= cutoffs_C0_new_v9_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
tiling_C0_new_v9_q2 = tiling %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_25_50_tiling[1] & C0_new_v9 <= cutoffs_C0_new_v9_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
tiling_C0_new_v9_q3 = tiling %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_50_75_tiling[1] & C0_new_v9 <= cutoffs_C0_new_v9_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
tiling_C0_new_v9_q4 = tiling %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_75_100_tiling[1] & C0_new_v9 <= cutoffs_C0_new_v9_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
tiling_C0_new_v9_custom1 = tiling %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_custom1_tiling[1] & C0_new_v9 <= cutoffs_C0_new_v9_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)

chrv_C0_new_v9_q1 = chrv %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_0_25_chrv[1] & C0_new_v9 <= cutoffs_C0_new_v9_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
chrv_C0_new_v9_q2 = chrv %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_25_50_chrv[1] & C0_new_v9 <= cutoffs_C0_new_v9_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
chrv_C0_new_v9_q3 = chrv %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_50_75_chrv[1] & C0_new_v9 <= cutoffs_C0_new_v9_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
chrv_C0_new_v9_q4 = chrv %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_75_100_chrv[1] & C0_new_v9 <= cutoffs_C0_new_v9_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)
chrv_C0_new_v9_custom1 = chrv %>%
  filter(C0_new_v9 >= cutoffs_C0_new_v9_custom1_chrv[1] & C0_new_v9 <= cutoffs_C0_new_v9_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v9)


# Divide each library into quartiles + custom for C0_new_v10
nuc_C0_new_v10_q1 = nuc %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_0_25_nuc[1] & C0_new_v10 <= cutoffs_C0_new_v10_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
nuc_C0_new_v10_q2 = nuc %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_25_50_nuc[1] & C0_new_v10 <= cutoffs_C0_new_v10_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
nuc_C0_new_v10_q3 = nuc %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_50_75_nuc[1] & C0_new_v10 <= cutoffs_C0_new_v10_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
nuc_C0_new_v10_q4 = nuc %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_75_100_nuc[1] & C0_new_v10 <= cutoffs_C0_new_v10_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
nuc_C0_new_v10_custom1 = nuc %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_custom1_nuc[1] & C0_new_v10 <= cutoffs_C0_new_v10_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)

random_C0_new_v10_q1 = random %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_0_25_random[1] & C0_new_v10 <= cutoffs_C0_new_v10_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
random_C0_new_v10_q2 = random %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_25_50_random[1] & C0_new_v10 <= cutoffs_C0_new_v10_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
random_C0_new_v10_q3 = random %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_50_75_random[1] & C0_new_v10 <= cutoffs_C0_new_v10_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
random_C0_new_v10_q4 = random %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_75_100_random[1] & C0_new_v10 <= cutoffs_C0_new_v10_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
random_C0_new_v10_custom1 = random %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_custom1_random[1] & C0_new_v10 <= cutoffs_C0_new_v10_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)

tiling_C0_new_v10_q1 = tiling %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_0_25_tiling[1] & C0_new_v10 <= cutoffs_C0_new_v10_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
tiling_C0_new_v10_q2 = tiling %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_25_50_tiling[1] & C0_new_v10 <= cutoffs_C0_new_v10_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
tiling_C0_new_v10_q3 = tiling %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_50_75_tiling[1] & C0_new_v10 <= cutoffs_C0_new_v10_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
tiling_C0_new_v10_q4 = tiling %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_75_100_tiling[1] & C0_new_v10 <= cutoffs_C0_new_v10_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
tiling_C0_new_v10_custom1 = tiling %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_custom1_tiling[1] & C0_new_v10 <= cutoffs_C0_new_v10_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)

chrv_C0_new_v10_q1 = chrv %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_0_25_chrv[1] & C0_new_v10 <= cutoffs_C0_new_v10_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
chrv_C0_new_v10_q2 = chrv %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_25_50_chrv[1] & C0_new_v10 <= cutoffs_C0_new_v10_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
chrv_C0_new_v10_q3 = chrv %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_50_75_chrv[1] & C0_new_v10 <= cutoffs_C0_new_v10_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
chrv_C0_new_v10_q4 = chrv %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_75_100_chrv[1] & C0_new_v10 <= cutoffs_C0_new_v10_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)
chrv_C0_new_v10_custom1 = chrv %>%
  filter(C0_new_v10 >= cutoffs_C0_new_v10_custom1_chrv[1] & C0_new_v10 <= cutoffs_C0_new_v10_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v10)


# Divide each library into quartiles + custom for C0_new_v11
nuc_C0_new_v11_q1 = nuc %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_0_25_nuc[1] & C0_new_v11 <= cutoffs_C0_new_v11_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
nuc_C0_new_v11_q2 = nuc %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_25_50_nuc[1] & C0_new_v11 <= cutoffs_C0_new_v11_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
nuc_C0_new_v11_q3 = nuc %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_50_75_nuc[1] & C0_new_v11 <= cutoffs_C0_new_v11_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
nuc_C0_new_v11_q4 = nuc %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_75_100_nuc[1] & C0_new_v11 <= cutoffs_C0_new_v11_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
nuc_C0_new_v11_custom1 = nuc %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_custom1_nuc[1] & C0_new_v11 <= cutoffs_C0_new_v11_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)

random_C0_new_v11_q1 = random %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_0_25_random[1] & C0_new_v11 <= cutoffs_C0_new_v11_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
random_C0_new_v11_q2 = random %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_25_50_random[1] & C0_new_v11 <= cutoffs_C0_new_v11_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
random_C0_new_v11_q3 = random %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_50_75_random[1] & C0_new_v11 <= cutoffs_C0_new_v11_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
random_C0_new_v11_q4 = random %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_75_100_random[1] & C0_new_v11 <= cutoffs_C0_new_v11_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
random_C0_new_v11_custom1 = random %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_custom1_random[1] & C0_new_v11 <= cutoffs_C0_new_v11_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)

tiling_C0_new_v11_q1 = tiling %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_0_25_tiling[1] & C0_new_v11 <= cutoffs_C0_new_v11_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
tiling_C0_new_v11_q2 = tiling %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_25_50_tiling[1] & C0_new_v11 <= cutoffs_C0_new_v11_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
tiling_C0_new_v11_q3 = tiling %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_50_75_tiling[1] & C0_new_v11 <= cutoffs_C0_new_v11_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
tiling_C0_new_v11_q4 = tiling %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_75_100_tiling[1] & C0_new_v11 <= cutoffs_C0_new_v11_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
tiling_C0_new_v11_custom1 = tiling %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_custom1_tiling[1] & C0_new_v11 <= cutoffs_C0_new_v11_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)

chrv_C0_new_v11_q1 = chrv %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_0_25_chrv[1] & C0_new_v11 <= cutoffs_C0_new_v11_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
chrv_C0_new_v11_q2 = chrv %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_25_50_chrv[1] & C0_new_v11 <= cutoffs_C0_new_v11_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
chrv_C0_new_v11_q3 = chrv %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_50_75_chrv[1] & C0_new_v11 <= cutoffs_C0_new_v11_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
chrv_C0_new_v11_q4 = chrv %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_75_100_chrv[1] & C0_new_v11 <= cutoffs_C0_new_v11_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)
chrv_C0_new_v11_custom1 = chrv %>%
  filter(C0_new_v11 >= cutoffs_C0_new_v11_custom1_chrv[1] & C0_new_v11 <= cutoffs_C0_new_v11_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v11)


# Divide each library into quartiles + custom for C0_new_v12
nuc_C0_new_v12_q1 = nuc %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_0_25_nuc[1] & C0_new_v12 <= cutoffs_C0_new_v12_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
nuc_C0_new_v12_q2 = nuc %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_25_50_nuc[1] & C0_new_v12 <= cutoffs_C0_new_v12_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
nuc_C0_new_v12_q3 = nuc %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_50_75_nuc[1] & C0_new_v12 <= cutoffs_C0_new_v12_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
nuc_C0_new_v12_q4 = nuc %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_75_100_nuc[1] & C0_new_v12 <= cutoffs_C0_new_v12_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
nuc_C0_new_v12_custom1 = nuc %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_custom1_nuc[1] & C0_new_v12 <= cutoffs_C0_new_v12_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)

random_C0_new_v12_q1 = random %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_0_25_random[1] & C0_new_v12 <= cutoffs_C0_new_v12_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
random_C0_new_v12_q2 = random %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_25_50_random[1] & C0_new_v12 <= cutoffs_C0_new_v12_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
random_C0_new_v12_q3 = random %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_50_75_random[1] & C0_new_v12 <= cutoffs_C0_new_v12_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
random_C0_new_v12_q4 = random %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_75_100_random[1] & C0_new_v12 <= cutoffs_C0_new_v12_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
random_C0_new_v12_custom1 = random %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_custom1_random[1] & C0_new_v12 <= cutoffs_C0_new_v12_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)

tiling_C0_new_v12_q1 = tiling %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_0_25_tiling[1] & C0_new_v12 <= cutoffs_C0_new_v12_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
tiling_C0_new_v12_q2 = tiling %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_25_50_tiling[1] & C0_new_v12 <= cutoffs_C0_new_v12_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
tiling_C0_new_v12_q3 = tiling %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_50_75_tiling[1] & C0_new_v12 <= cutoffs_C0_new_v12_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
tiling_C0_new_v12_q4 = tiling %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_75_100_tiling[1] & C0_new_v12 <= cutoffs_C0_new_v12_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
tiling_C0_new_v12_custom1 = tiling %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_custom1_tiling[1] & C0_new_v12 <= cutoffs_C0_new_v12_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)

chrv_C0_new_v12_q1 = chrv %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_0_25_chrv[1] & C0_new_v12 <= cutoffs_C0_new_v12_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
chrv_C0_new_v12_q2 = chrv %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_25_50_chrv[1] & C0_new_v12 <= cutoffs_C0_new_v12_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
chrv_C0_new_v12_q3 = chrv %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_50_75_chrv[1] & C0_new_v12 <= cutoffs_C0_new_v12_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
chrv_C0_new_v12_q4 = chrv %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_75_100_chrv[1] & C0_new_v12 <= cutoffs_C0_new_v12_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)
chrv_C0_new_v12_custom1 = chrv %>%
  filter(C0_new_v12 >= cutoffs_C0_new_v12_custom1_chrv[1] & C0_new_v12 <= cutoffs_C0_new_v12_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v12)


# Divide each library into quartiles + custom for C0_new_v13
nuc_C0_new_v13_q1 = nuc %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_0_25_nuc[1] & C0_new_v13 <= cutoffs_C0_new_v13_0_25_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
nuc_C0_new_v13_q2 = nuc %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_25_50_nuc[1] & C0_new_v13 <= cutoffs_C0_new_v13_25_50_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
nuc_C0_new_v13_q3 = nuc %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_50_75_nuc[1] & C0_new_v13 <= cutoffs_C0_new_v13_50_75_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
nuc_C0_new_v13_q4 = nuc %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_75_100_nuc[1] & C0_new_v13 <= cutoffs_C0_new_v13_75_100_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
nuc_C0_new_v13_custom1 = nuc %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_custom1_nuc[1] & C0_new_v13 <= cutoffs_C0_new_v13_custom1_nuc[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)

random_C0_new_v13_q1 = random %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_0_25_random[1] & C0_new_v13 <= cutoffs_C0_new_v13_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
random_C0_new_v13_q2 = random %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_25_50_random[1] & C0_new_v13 <= cutoffs_C0_new_v13_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
random_C0_new_v13_q3 = random %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_50_75_random[1] & C0_new_v13 <= cutoffs_C0_new_v13_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
random_C0_new_v13_q4 = random %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_75_100_random[1] & C0_new_v13 <= cutoffs_C0_new_v13_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
random_C0_new_v13_custom1 = random %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_custom1_random[1] & C0_new_v13 <= cutoffs_C0_new_v13_custom1_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)

tiling_C0_new_v13_q1 = tiling %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_0_25_tiling[1] & C0_new_v13 <= cutoffs_C0_new_v13_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
tiling_C0_new_v13_q2 = tiling %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_25_50_tiling[1] & C0_new_v13 <= cutoffs_C0_new_v13_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
tiling_C0_new_v13_q3 = tiling %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_50_75_tiling[1] & C0_new_v13 <= cutoffs_C0_new_v13_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
tiling_C0_new_v13_q4 = tiling %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_75_100_tiling[1] & C0_new_v13 <= cutoffs_C0_new_v13_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
tiling_C0_new_v13_custom1 = tiling %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_custom1_tiling[1] & C0_new_v13 <= cutoffs_C0_new_v13_custom1_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)

chrv_C0_new_v13_q1 = chrv %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_0_25_chrv[1] & C0_new_v13 <= cutoffs_C0_new_v13_0_25_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
chrv_C0_new_v13_q2 = chrv %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_25_50_chrv[1] & C0_new_v13 <= cutoffs_C0_new_v13_25_50_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
chrv_C0_new_v13_q3 = chrv %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_50_75_chrv[1] & C0_new_v13 <= cutoffs_C0_new_v13_50_75_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
chrv_C0_new_v13_q4 = chrv %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_75_100_chrv[1] & C0_new_v13 <= cutoffs_C0_new_v13_75_100_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)
chrv_C0_new_v13_custom1 = chrv %>%
  filter(C0_new_v13 >= cutoffs_C0_new_v13_custom1_chrv[1] & C0_new_v13 <= cutoffs_C0_new_v13_custom1_chrv[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new_v13)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C26
nuc_A_T_C26_q1 = apply(nuc_C26_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C26_q1))
})
nuc_A_T_C26_q2 = apply(nuc_C26_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C26_q2))
})
nuc_A_T_C26_q3 = apply(nuc_C26_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C26_q3))
})
nuc_A_T_C26_q4 = apply(nuc_C26_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C26_q4))
})
nuc_A_T_C26_custom1 = apply(nuc_C26_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C26_custom1))
})

random_A_T_C26_q1 = apply(random_C26_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C26_q1))
})
random_A_T_C26_q2 = apply(random_C26_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C26_q2))
})
random_A_T_C26_q3 = apply(random_C26_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C26_q3))
})
random_A_T_C26_q4 = apply(random_C26_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C26_q4))
})
random_A_T_C26_custom1 = apply(random_C26_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C26_custom1))
})

tiling_A_T_C26_q1 = apply(tiling_C26_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C26_q1))
})
tiling_A_T_C26_q2 = apply(tiling_C26_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C26_q2))
})
tiling_A_T_C26_q3 = apply(tiling_C26_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C26_q3))
})
tiling_A_T_C26_q4 = apply(tiling_C26_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C26_q4))
})
tiling_A_T_C26_custom1 = apply(tiling_C26_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C26_custom1))
})

chrv_A_T_C26_q1 = apply(chrv_C26_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C26_q1))
})
chrv_A_T_C26_q2 = apply(chrv_C26_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C26_q2))
})
chrv_A_T_C26_q3 = apply(chrv_C26_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C26_q3))
})
chrv_A_T_C26_q4 = apply(chrv_C26_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C26_q4))
})
chrv_A_T_C26_custom1 = apply(chrv_C26_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C26_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C26_q = data.frame(q1 = nuc_A_T_C26_q1,
                           q2 = nuc_A_T_C26_q2,
                           q3 = nuc_A_T_C26_q3,
                           q4 = nuc_A_T_C26_q4, 
                           custom1 = nuc_A_T_C26_custom1)

random_A_T_C26_q = data.frame(q1 = random_A_T_C26_q1,
                              q2 = random_A_T_C26_q2,
                              q3 = random_A_T_C26_q3,
                              q4 = random_A_T_C26_q4, 
                              custom1 = random_A_T_C26_custom1)

tiling_A_T_C26_q = data.frame(q1 = tiling_A_T_C26_q1,
                              q2 = tiling_A_T_C26_q2,
                              q3 = tiling_A_T_C26_q3,
                              q4 = tiling_A_T_C26_q4, 
                              custom1 = tiling_A_T_C26_custom1)

chrv_A_T_C26_q = data.frame(q1 = chrv_A_T_C26_q1,
                            q2 = chrv_A_T_C26_q2,
                            q3 = chrv_A_T_C26_q3,
                            q4 = chrv_A_T_C26_q4, 
                            custom1 = chrv_A_T_C26_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C29
nuc_A_T_C29_q1 = apply(nuc_C29_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C29_q1))
})
nuc_A_T_C29_q2 = apply(nuc_C29_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C29_q2))
})
nuc_A_T_C29_q3 = apply(nuc_C29_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C29_q3))
})
nuc_A_T_C29_q4 = apply(nuc_C29_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C29_q4))
})
nuc_A_T_C29_custom1 = apply(nuc_C29_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C29_custom1))
})

random_A_T_C29_q1 = apply(random_C29_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C29_q1))
})
random_A_T_C29_q2 = apply(random_C29_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C29_q2))
})
random_A_T_C29_q3 = apply(random_C29_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C29_q3))
})
random_A_T_C29_q4 = apply(random_C29_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C29_q4))
})
random_A_T_C29_custom1 = apply(random_C29_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C29_custom1))
})

tiling_A_T_C29_q1 = apply(tiling_C29_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C29_q1))
})
tiling_A_T_C29_q2 = apply(tiling_C29_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C29_q2))
})
tiling_A_T_C29_q3 = apply(tiling_C29_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C29_q3))
})
tiling_A_T_C29_q4 = apply(tiling_C29_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C29_q4))
})
tiling_A_T_C29_custom1 = apply(tiling_C29_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C29_custom1))
})

chrv_A_T_C29_q1 = apply(chrv_C29_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C29_q1))
})
chrv_A_T_C29_q2 = apply(chrv_C29_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C29_q2))
})
chrv_A_T_C29_q3 = apply(chrv_C29_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C29_q3))
})
chrv_A_T_C29_q4 = apply(chrv_C29_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C29_q4))
})
chrv_A_T_C29_custom1 = apply(chrv_C29_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C29_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C29_q = data.frame(q1 = nuc_A_T_C29_q1,
                           q2 = nuc_A_T_C29_q2,
                           q3 = nuc_A_T_C29_q3,
                           q4 = nuc_A_T_C29_q4, 
                           custom1 = nuc_A_T_C29_custom1)

random_A_T_C29_q = data.frame(q1 = random_A_T_C29_q1,
                              q2 = random_A_T_C29_q2,
                              q3 = random_A_T_C29_q3,
                              q4 = random_A_T_C29_q4, 
                              custom1 = random_A_T_C29_custom1)

tiling_A_T_C29_q = data.frame(q1 = tiling_A_T_C29_q1,
                              q2 = tiling_A_T_C29_q2,
                              q3 = tiling_A_T_C29_q3,
                              q4 = tiling_A_T_C29_q4, 
                              custom1 = tiling_A_T_C29_custom1)

chrv_A_T_C29_q = data.frame(q1 = chrv_A_T_C29_q1,
                            q2 = chrv_A_T_C29_q2,
                            q3 = chrv_A_T_C29_q3,
                            q4 = chrv_A_T_C29_q4, 
                            custom1 = chrv_A_T_C29_custom1)




# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C31
nuc_A_T_C31_q1 = apply(nuc_C31_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C31_q1))
})
nuc_A_T_C31_q2 = apply(nuc_C31_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C31_q2))
})
nuc_A_T_C31_q3 = apply(nuc_C31_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C31_q3))
})
nuc_A_T_C31_q4 = apply(nuc_C31_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C31_q4))
})
nuc_A_T_C31_custom1 = apply(nuc_C31_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C31_custom1))
})

random_A_T_C31_q1 = apply(random_C31_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C31_q1))
})
random_A_T_C31_q2 = apply(random_C31_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C31_q2))
})
random_A_T_C31_q3 = apply(random_C31_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C31_q3))
})
random_A_T_C31_q4 = apply(random_C31_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C31_q4))
})
random_A_T_C31_custom1 = apply(random_C31_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C31_custom1))
})

tiling_A_T_C31_q1 = apply(tiling_C31_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C31_q1))
})
tiling_A_T_C31_q2 = apply(tiling_C31_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C31_q2))
})
tiling_A_T_C31_q3 = apply(tiling_C31_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C31_q3))
})
tiling_A_T_C31_q4 = apply(tiling_C31_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C31_q4))
})
tiling_A_T_C31_custom1 = apply(tiling_C31_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C31_custom1))
})

chrv_A_T_C31_q1 = apply(chrv_C31_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C31_q1))
})
chrv_A_T_C31_q2 = apply(chrv_C31_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C31_q2))
})
chrv_A_T_C31_q3 = apply(chrv_C31_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C31_q3))
})
chrv_A_T_C31_q4 = apply(chrv_C31_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C31_q4))
})
chrv_A_T_C31_custom1 = apply(chrv_C31_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C31_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C31_q = data.frame(q1 = nuc_A_T_C31_q1,
                           q2 = nuc_A_T_C31_q2,
                           q3 = nuc_A_T_C31_q3,
                           q4 = nuc_A_T_C31_q4, 
                           custom1 = nuc_A_T_C31_custom1)

random_A_T_C31_q = data.frame(q1 = random_A_T_C31_q1,
                              q2 = random_A_T_C31_q2,
                              q3 = random_A_T_C31_q3,
                              q4 = random_A_T_C31_q4, 
                              custom1 = random_A_T_C31_custom1)

tiling_A_T_C31_q = data.frame(q1 = tiling_A_T_C31_q1,
                              q2 = tiling_A_T_C31_q2,
                              q3 = tiling_A_T_C31_q3,
                              q4 = tiling_A_T_C31_q4, 
                              custom1 = tiling_A_T_C31_custom1)

chrv_A_T_C31_q = data.frame(q1 = chrv_A_T_C31_q1,
                            q2 = chrv_A_T_C31_q2,
                            q3 = chrv_A_T_C31_q3,
                            q4 = chrv_A_T_C31_q4, 
                            custom1 = chrv_A_T_C31_custom1)




# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0
nuc_A_T_C0_q1 = apply(nuc_C0_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_q1))
})
nuc_A_T_C0_q2 = apply(nuc_C0_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_q2))
})
nuc_A_T_C0_q3 = apply(nuc_C0_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_q3))
})
nuc_A_T_C0_q4 = apply(nuc_C0_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_q4))
})
nuc_A_T_C0_custom1 = apply(nuc_C0_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_custom1))
})

random_A_T_C0_q1 = apply(random_C0_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_q1))
})
random_A_T_C0_q2 = apply(random_C0_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_q2))
})
random_A_T_C0_q3 = apply(random_C0_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_q3))
})
random_A_T_C0_q4 = apply(random_C0_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_q4))
})
random_A_T_C0_custom1 = apply(random_C0_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_custom1))
})

tiling_A_T_C0_q1 = apply(tiling_C0_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_q1))
})
tiling_A_T_C0_q2 = apply(tiling_C0_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_q2))
})
tiling_A_T_C0_q3 = apply(tiling_C0_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_q3))
})
tiling_A_T_C0_q4 = apply(tiling_C0_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_q4))
})
tiling_A_T_C0_custom1 = apply(tiling_C0_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_custom1))
})

chrv_A_T_C0_q1 = apply(chrv_C0_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_q1))
})
chrv_A_T_C0_q2 = apply(chrv_C0_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_q2))
})
chrv_A_T_C0_q3 = apply(chrv_C0_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_q3))
})
chrv_A_T_C0_q4 = apply(chrv_C0_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_q4))
})
chrv_A_T_C0_custom1 = apply(chrv_C0_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_q = data.frame(q1 = nuc_A_T_C0_q1,
                          q2 = nuc_A_T_C0_q2,
                          q3 = nuc_A_T_C0_q3,
                          q4 = nuc_A_T_C0_q4, 
                          custom1 = nuc_A_T_C0_custom1)

random_A_T_C0_q = data.frame(q1 = random_A_T_C0_q1,
                             q2 = random_A_T_C0_q2,
                             q3 = random_A_T_C0_q3,
                             q4 = random_A_T_C0_q4, 
                             custom1 = random_A_T_C0_custom1)

tiling_A_T_C0_q = data.frame(q1 = tiling_A_T_C0_q1,
                             q2 = tiling_A_T_C0_q2,
                             q3 = tiling_A_T_C0_q3,
                             q4 = tiling_A_T_C0_q4, 
                             custom1 = tiling_A_T_C0_custom1)

chrv_A_T_C0_q = data.frame(q1 = chrv_A_T_C0_q1,
                           q2 = chrv_A_T_C0_q2,
                           q3 = chrv_A_T_C0_q3,
                           q4 = chrv_A_T_C0_q4, 
                           custom1 = chrv_A_T_C0_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v1
nuc_A_T_C0_new_v1_q1 = apply(nuc_C0_new_v1_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v1_q1))
})
nuc_A_T_C0_new_v1_q2 = apply(nuc_C0_new_v1_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v1_q2))
})
nuc_A_T_C0_new_v1_q3 = apply(nuc_C0_new_v1_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v1_q3))
})
nuc_A_T_C0_new_v1_q4 = apply(nuc_C0_new_v1_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v1_q4))
})
nuc_A_T_C0_new_v1_custom1 = apply(nuc_C0_new_v1_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v1_custom1))
})

random_A_T_C0_new_v1_q1 = apply(random_C0_new_v1_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v1_q1))
})
random_A_T_C0_new_v1_q2 = apply(random_C0_new_v1_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v1_q2))
})
random_A_T_C0_new_v1_q3 = apply(random_C0_new_v1_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v1_q3))
})
random_A_T_C0_new_v1_q4 = apply(random_C0_new_v1_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v1_q4))
})
random_A_T_C0_new_v1_custom1 = apply(random_C0_new_v1_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v1_custom1))
})

tiling_A_T_C0_new_v1_q1 = apply(tiling_C0_new_v1_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v1_q1))
})
tiling_A_T_C0_new_v1_q2 = apply(tiling_C0_new_v1_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v1_q2))
})
tiling_A_T_C0_new_v1_q3 = apply(tiling_C0_new_v1_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v1_q3))
})
tiling_A_T_C0_new_v1_q4 = apply(tiling_C0_new_v1_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v1_q4))
})
tiling_A_T_C0_new_v1_custom1 = apply(tiling_C0_new_v1_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v1_custom1))
})

chrv_A_T_C0_new_v1_q1 = apply(chrv_C0_new_v1_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v1_q1))
})
chrv_A_T_C0_new_v1_q2 = apply(chrv_C0_new_v1_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v1_q2))
})
chrv_A_T_C0_new_v1_q3 = apply(chrv_C0_new_v1_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v1_q3))
})
chrv_A_T_C0_new_v1_q4 = apply(chrv_C0_new_v1_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v1_q4))
})
chrv_A_T_C0_new_v1_custom1 = apply(chrv_C0_new_v1_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v1_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v1_q = data.frame(q1 = nuc_A_T_C0_new_v1_q1,
                                 q2 = nuc_A_T_C0_new_v1_q2,
                                 q3 = nuc_A_T_C0_new_v1_q3,
                                 q4 = nuc_A_T_C0_new_v1_q4, 
                                 custom1 = nuc_A_T_C0_new_v1_custom1)

random_A_T_C0_new_v1_q = data.frame(q1 = random_A_T_C0_new_v1_q1,
                                    q2 = random_A_T_C0_new_v1_q2,
                                    q3 = random_A_T_C0_new_v1_q3,
                                    q4 = random_A_T_C0_new_v1_q4, 
                                    custom1 = random_A_T_C0_new_v1_custom1)

tiling_A_T_C0_new_v1_q = data.frame(q1 = tiling_A_T_C0_new_v1_q1,
                                    q2 = tiling_A_T_C0_new_v1_q2,
                                    q3 = tiling_A_T_C0_new_v1_q3,
                                    q4 = tiling_A_T_C0_new_v1_q4, 
                                    custom1 = tiling_A_T_C0_new_v1_custom1)

chrv_A_T_C0_new_v1_q = data.frame(q1 = chrv_A_T_C0_new_v1_q1,
                                  q2 = chrv_A_T_C0_new_v1_q2,
                                  q3 = chrv_A_T_C0_new_v1_q3,
                                  q4 = chrv_A_T_C0_new_v1_q4, 
                                  custom1 = chrv_A_T_C0_new_v1_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v2
nuc_A_T_C0_new_v2_q1 = apply(nuc_C0_new_v2_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v2_q1))
})
nuc_A_T_C0_new_v2_q2 = apply(nuc_C0_new_v2_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v2_q2))
})
nuc_A_T_C0_new_v2_q3 = apply(nuc_C0_new_v2_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v2_q3))
})
nuc_A_T_C0_new_v2_q4 = apply(nuc_C0_new_v2_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v2_q4))
})
nuc_A_T_C0_new_v2_custom1 = apply(nuc_C0_new_v2_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v2_custom1))
})

random_A_T_C0_new_v2_q1 = apply(random_C0_new_v2_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v2_q1))
})
random_A_T_C0_new_v2_q2 = apply(random_C0_new_v2_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v2_q2))
})
random_A_T_C0_new_v2_q3 = apply(random_C0_new_v2_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v2_q3))
})
random_A_T_C0_new_v2_q4 = apply(random_C0_new_v2_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v2_q4))
})
random_A_T_C0_new_v2_custom1 = apply(random_C0_new_v2_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v2_custom1))
})

tiling_A_T_C0_new_v2_q1 = apply(tiling_C0_new_v2_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v2_q1))
})
tiling_A_T_C0_new_v2_q2 = apply(tiling_C0_new_v2_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v2_q2))
})
tiling_A_T_C0_new_v2_q3 = apply(tiling_C0_new_v2_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v2_q3))
})
tiling_A_T_C0_new_v2_q4 = apply(tiling_C0_new_v2_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v2_q4))
})
tiling_A_T_C0_new_v2_custom1 = apply(tiling_C0_new_v2_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v2_custom1))
})

chrv_A_T_C0_new_v2_q1 = apply(chrv_C0_new_v2_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v2_q1))
})
chrv_A_T_C0_new_v2_q2 = apply(chrv_C0_new_v2_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v2_q2))
})
chrv_A_T_C0_new_v2_q3 = apply(chrv_C0_new_v2_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v2_q3))
})
chrv_A_T_C0_new_v2_q4 = apply(chrv_C0_new_v2_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v2_q4))
})
chrv_A_T_C0_new_v2_custom1 = apply(chrv_C0_new_v2_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v2_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v2_q = data.frame(q1 = nuc_A_T_C0_new_v2_q1,
                                 q2 = nuc_A_T_C0_new_v2_q2,
                                 q3 = nuc_A_T_C0_new_v2_q3,
                                 q4 = nuc_A_T_C0_new_v2_q4, 
                                 custom1 = nuc_A_T_C0_new_v2_custom1)

random_A_T_C0_new_v2_q = data.frame(q1 = random_A_T_C0_new_v2_q1,
                                    q2 = random_A_T_C0_new_v2_q2,
                                    q3 = random_A_T_C0_new_v2_q3,
                                    q4 = random_A_T_C0_new_v2_q4, 
                                    custom1 = random_A_T_C0_new_v2_custom1)

tiling_A_T_C0_new_v2_q = data.frame(q1 = tiling_A_T_C0_new_v2_q1,
                                    q2 = tiling_A_T_C0_new_v2_q2,
                                    q3 = tiling_A_T_C0_new_v2_q3,
                                    q4 = tiling_A_T_C0_new_v2_q4, 
                                    custom1 = tiling_A_T_C0_new_v2_custom1)

chrv_A_T_C0_new_v2_q = data.frame(q1 = chrv_A_T_C0_new_v2_q1,
                                  q2 = chrv_A_T_C0_new_v2_q2,
                                  q3 = chrv_A_T_C0_new_v2_q3,
                                  q4 = chrv_A_T_C0_new_v2_q4, 
                                  custom1 = chrv_A_T_C0_new_v2_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v3
nuc_A_T_C0_new_v3_q1 = apply(nuc_C0_new_v3_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v3_q1))
})
nuc_A_T_C0_new_v3_q2 = apply(nuc_C0_new_v3_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v3_q2))
})
nuc_A_T_C0_new_v3_q3 = apply(nuc_C0_new_v3_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v3_q3))
})
nuc_A_T_C0_new_v3_q4 = apply(nuc_C0_new_v3_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v3_q4))
})
nuc_A_T_C0_new_v3_custom1 = apply(nuc_C0_new_v3_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v3_custom1))
})

random_A_T_C0_new_v3_q1 = apply(random_C0_new_v3_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v3_q1))
})
random_A_T_C0_new_v3_q2 = apply(random_C0_new_v3_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v3_q2))
})
random_A_T_C0_new_v3_q3 = apply(random_C0_new_v3_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v3_q3))
})
random_A_T_C0_new_v3_q4 = apply(random_C0_new_v3_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v3_q4))
})
random_A_T_C0_new_v3_custom1 = apply(random_C0_new_v3_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v3_custom1))
})

tiling_A_T_C0_new_v3_q1 = apply(tiling_C0_new_v3_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v3_q1))
})
tiling_A_T_C0_new_v3_q2 = apply(tiling_C0_new_v3_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v3_q2))
})
tiling_A_T_C0_new_v3_q3 = apply(tiling_C0_new_v3_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v3_q3))
})
tiling_A_T_C0_new_v3_q4 = apply(tiling_C0_new_v3_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v3_q4))
})
tiling_A_T_C0_new_v3_custom1 = apply(tiling_C0_new_v3_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v3_custom1))
})

chrv_A_T_C0_new_v3_q1 = apply(chrv_C0_new_v3_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v3_q1))
})
chrv_A_T_C0_new_v3_q2 = apply(chrv_C0_new_v3_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v3_q2))
})
chrv_A_T_C0_new_v3_q3 = apply(chrv_C0_new_v3_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v3_q3))
})
chrv_A_T_C0_new_v3_q4 = apply(chrv_C0_new_v3_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v3_q4))
})
chrv_A_T_C0_new_v3_custom1 = apply(chrv_C0_new_v3_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v3_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v3_q = data.frame(q1 = nuc_A_T_C0_new_v3_q1,
                                 q2 = nuc_A_T_C0_new_v3_q2,
                                 q3 = nuc_A_T_C0_new_v3_q3,
                                 q4 = nuc_A_T_C0_new_v3_q4, 
                                 custom1 = nuc_A_T_C0_new_v3_custom1)

random_A_T_C0_new_v3_q = data.frame(q1 = random_A_T_C0_new_v3_q1,
                                    q2 = random_A_T_C0_new_v3_q2,
                                    q3 = random_A_T_C0_new_v3_q3,
                                    q4 = random_A_T_C0_new_v3_q4, 
                                    custom1 = random_A_T_C0_new_v3_custom1)

tiling_A_T_C0_new_v3_q = data.frame(q1 = tiling_A_T_C0_new_v3_q1,
                                    q2 = tiling_A_T_C0_new_v3_q2,
                                    q3 = tiling_A_T_C0_new_v3_q3,
                                    q4 = tiling_A_T_C0_new_v3_q4, 
                                    custom1 = tiling_A_T_C0_new_v3_custom1)

chrv_A_T_C0_new_v3_q = data.frame(q1 = chrv_A_T_C0_new_v3_q1,
                                  q2 = chrv_A_T_C0_new_v3_q2,
                                  q3 = chrv_A_T_C0_new_v3_q3,
                                  q4 = chrv_A_T_C0_new_v3_q4, 
                                  custom1 = chrv_A_T_C0_new_v3_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v4
nuc_A_T_C0_new_v4_q1 = apply(nuc_C0_new_v4_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v4_q1))
})
nuc_A_T_C0_new_v4_q2 = apply(nuc_C0_new_v4_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v4_q2))
})
nuc_A_T_C0_new_v4_q3 = apply(nuc_C0_new_v4_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v4_q3))
})
nuc_A_T_C0_new_v4_q4 = apply(nuc_C0_new_v4_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v4_q4))
})
nuc_A_T_C0_new_v4_custom1 = apply(nuc_C0_new_v4_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v4_custom1))
})

random_A_T_C0_new_v4_q1 = apply(random_C0_new_v4_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v4_q1))
})
random_A_T_C0_new_v4_q2 = apply(random_C0_new_v4_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v4_q2))
})
random_A_T_C0_new_v4_q3 = apply(random_C0_new_v4_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v4_q3))
})
random_A_T_C0_new_v4_q4 = apply(random_C0_new_v4_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v4_q4))
})
random_A_T_C0_new_v4_custom1 = apply(random_C0_new_v4_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v4_custom1))
})

tiling_A_T_C0_new_v4_q1 = apply(tiling_C0_new_v4_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v4_q1))
})
tiling_A_T_C0_new_v4_q2 = apply(tiling_C0_new_v4_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v4_q2))
})
tiling_A_T_C0_new_v4_q3 = apply(tiling_C0_new_v4_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v4_q3))
})
tiling_A_T_C0_new_v4_q4 = apply(tiling_C0_new_v4_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v4_q4))
})
tiling_A_T_C0_new_v4_custom1 = apply(tiling_C0_new_v4_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v4_custom1))
})

chrv_A_T_C0_new_v4_q1 = apply(chrv_C0_new_v4_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v4_q1))
})
chrv_A_T_C0_new_v4_q2 = apply(chrv_C0_new_v4_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v4_q2))
})
chrv_A_T_C0_new_v4_q3 = apply(chrv_C0_new_v4_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v4_q3))
})
chrv_A_T_C0_new_v4_q4 = apply(chrv_C0_new_v4_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v4_q4))
})
chrv_A_T_C0_new_v4_custom1 = apply(chrv_C0_new_v4_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v4_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v4_q = data.frame(q1 = nuc_A_T_C0_new_v4_q1,
                                 q2 = nuc_A_T_C0_new_v4_q2,
                                 q3 = nuc_A_T_C0_new_v4_q3,
                                 q4 = nuc_A_T_C0_new_v4_q4, 
                                 custom1 = nuc_A_T_C0_new_v4_custom1)

random_A_T_C0_new_v4_q = data.frame(q1 = random_A_T_C0_new_v4_q1,
                                    q2 = random_A_T_C0_new_v4_q2,
                                    q3 = random_A_T_C0_new_v4_q3,
                                    q4 = random_A_T_C0_new_v4_q4, 
                                    custom1 = random_A_T_C0_new_v4_custom1)

tiling_A_T_C0_new_v4_q = data.frame(q1 = tiling_A_T_C0_new_v4_q1,
                                    q2 = tiling_A_T_C0_new_v4_q2,
                                    q3 = tiling_A_T_C0_new_v4_q3,
                                    q4 = tiling_A_T_C0_new_v4_q4, 
                                    custom1 = tiling_A_T_C0_new_v4_custom1)

chrv_A_T_C0_new_v4_q = data.frame(q1 = chrv_A_T_C0_new_v4_q1,
                                  q2 = chrv_A_T_C0_new_v4_q2,
                                  q3 = chrv_A_T_C0_new_v4_q3,
                                  q4 = chrv_A_T_C0_new_v4_q4, 
                                  custom1 = chrv_A_T_C0_new_v4_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v5
nuc_A_T_C0_new_v5_q1 = apply(nuc_C0_new_v5_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v5_q1))
})
nuc_A_T_C0_new_v5_q2 = apply(nuc_C0_new_v5_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v5_q2))
})
nuc_A_T_C0_new_v5_q3 = apply(nuc_C0_new_v5_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v5_q3))
})
nuc_A_T_C0_new_v5_q4 = apply(nuc_C0_new_v5_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v5_q4))
})
nuc_A_T_C0_new_v5_custom1 = apply(nuc_C0_new_v5_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v5_custom1))
})

random_A_T_C0_new_v5_q1 = apply(random_C0_new_v5_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v5_q1))
})
random_A_T_C0_new_v5_q2 = apply(random_C0_new_v5_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v5_q2))
})
random_A_T_C0_new_v5_q3 = apply(random_C0_new_v5_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v5_q3))
})
random_A_T_C0_new_v5_q4 = apply(random_C0_new_v5_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v5_q4))
})
random_A_T_C0_new_v5_custom1 = apply(random_C0_new_v5_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v5_custom1))
})

tiling_A_T_C0_new_v5_q1 = apply(tiling_C0_new_v5_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v5_q1))
})
tiling_A_T_C0_new_v5_q2 = apply(tiling_C0_new_v5_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v5_q2))
})
tiling_A_T_C0_new_v5_q3 = apply(tiling_C0_new_v5_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v5_q3))
})
tiling_A_T_C0_new_v5_q4 = apply(tiling_C0_new_v5_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v5_q4))
})
tiling_A_T_C0_new_v5_custom1 = apply(tiling_C0_new_v5_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v5_custom1))
})

chrv_A_T_C0_new_v5_q1 = apply(chrv_C0_new_v5_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v5_q1))
})
chrv_A_T_C0_new_v5_q2 = apply(chrv_C0_new_v5_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v5_q2))
})
chrv_A_T_C0_new_v5_q3 = apply(chrv_C0_new_v5_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v5_q3))
})
chrv_A_T_C0_new_v5_q4 = apply(chrv_C0_new_v5_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v5_q4))
})
chrv_A_T_C0_new_v5_custom1 = apply(chrv_C0_new_v5_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v5_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v5_q = data.frame(q1 = nuc_A_T_C0_new_v5_q1,
                                 q2 = nuc_A_T_C0_new_v5_q2,
                                 q3 = nuc_A_T_C0_new_v5_q3,
                                 q4 = nuc_A_T_C0_new_v5_q4, 
                                 custom1 = nuc_A_T_C0_new_v5_custom1)

random_A_T_C0_new_v5_q = data.frame(q1 = random_A_T_C0_new_v5_q1,
                                    q2 = random_A_T_C0_new_v5_q2,
                                    q3 = random_A_T_C0_new_v5_q3,
                                    q4 = random_A_T_C0_new_v5_q4, 
                                    custom1 = random_A_T_C0_new_v5_custom1)

tiling_A_T_C0_new_v5_q = data.frame(q1 = tiling_A_T_C0_new_v5_q1,
                                    q2 = tiling_A_T_C0_new_v5_q2,
                                    q3 = tiling_A_T_C0_new_v5_q3,
                                    q4 = tiling_A_T_C0_new_v5_q4, 
                                    custom1 = tiling_A_T_C0_new_v5_custom1)

chrv_A_T_C0_new_v5_q = data.frame(q1 = chrv_A_T_C0_new_v5_q1,
                                  q2 = chrv_A_T_C0_new_v5_q2,
                                  q3 = chrv_A_T_C0_new_v5_q3,
                                  q4 = chrv_A_T_C0_new_v5_q4, 
                                  custom1 = chrv_A_T_C0_new_v5_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v6
nuc_A_T_C0_new_v6_q1 = apply(nuc_C0_new_v6_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v6_q1))
})
nuc_A_T_C0_new_v6_q2 = apply(nuc_C0_new_v6_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v6_q2))
})
nuc_A_T_C0_new_v6_q3 = apply(nuc_C0_new_v6_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v6_q3))
})
nuc_A_T_C0_new_v6_q4 = apply(nuc_C0_new_v6_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v6_q4))
})
nuc_A_T_C0_new_v6_custom1 = apply(nuc_C0_new_v6_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v6_custom1))
})

random_A_T_C0_new_v6_q1 = apply(random_C0_new_v6_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v6_q1))
})
random_A_T_C0_new_v6_q2 = apply(random_C0_new_v6_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v6_q2))
})
random_A_T_C0_new_v6_q3 = apply(random_C0_new_v6_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v6_q3))
})
random_A_T_C0_new_v6_q4 = apply(random_C0_new_v6_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v6_q4))
})
random_A_T_C0_new_v6_custom1 = apply(random_C0_new_v6_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v6_custom1))
})

tiling_A_T_C0_new_v6_q1 = apply(tiling_C0_new_v6_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v6_q1))
})
tiling_A_T_C0_new_v6_q2 = apply(tiling_C0_new_v6_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v6_q2))
})
tiling_A_T_C0_new_v6_q3 = apply(tiling_C0_new_v6_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v6_q3))
})
tiling_A_T_C0_new_v6_q4 = apply(tiling_C0_new_v6_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v6_q4))
})
tiling_A_T_C0_new_v6_custom1 = apply(tiling_C0_new_v6_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v6_custom1))
})

chrv_A_T_C0_new_v6_q1 = apply(chrv_C0_new_v6_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v6_q1))
})
chrv_A_T_C0_new_v6_q2 = apply(chrv_C0_new_v6_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v6_q2))
})
chrv_A_T_C0_new_v6_q3 = apply(chrv_C0_new_v6_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v6_q3))
})
chrv_A_T_C0_new_v6_q4 = apply(chrv_C0_new_v6_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v6_q4))
})
chrv_A_T_C0_new_v6_custom1 = apply(chrv_C0_new_v6_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v6_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v6_q = data.frame(q1 = nuc_A_T_C0_new_v6_q1,
                                 q2 = nuc_A_T_C0_new_v6_q2,
                                 q3 = nuc_A_T_C0_new_v6_q3,
                                 q4 = nuc_A_T_C0_new_v6_q4, 
                                 custom1 = nuc_A_T_C0_new_v6_custom1)

random_A_T_C0_new_v6_q = data.frame(q1 = random_A_T_C0_new_v6_q1,
                                    q2 = random_A_T_C0_new_v6_q2,
                                    q3 = random_A_T_C0_new_v6_q3,
                                    q4 = random_A_T_C0_new_v6_q4, 
                                    custom1 = random_A_T_C0_new_v6_custom1)

tiling_A_T_C0_new_v6_q = data.frame(q1 = tiling_A_T_C0_new_v6_q1,
                                    q2 = tiling_A_T_C0_new_v6_q2,
                                    q3 = tiling_A_T_C0_new_v6_q3,
                                    q4 = tiling_A_T_C0_new_v6_q4, 
                                    custom1 = tiling_A_T_C0_new_v6_custom1)

chrv_A_T_C0_new_v6_q = data.frame(q1 = chrv_A_T_C0_new_v6_q1,
                                  q2 = chrv_A_T_C0_new_v6_q2,
                                  q3 = chrv_A_T_C0_new_v6_q3,
                                  q4 = chrv_A_T_C0_new_v6_q4, 
                                  custom1 = chrv_A_T_C0_new_v6_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v7
nuc_A_T_C0_new_v7_q1 = apply(nuc_C0_new_v7_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v7_q1))
})
nuc_A_T_C0_new_v7_q2 = apply(nuc_C0_new_v7_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v7_q2))
})
nuc_A_T_C0_new_v7_q3 = apply(nuc_C0_new_v7_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v7_q3))
})
nuc_A_T_C0_new_v7_q4 = apply(nuc_C0_new_v7_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v7_q4))
})
nuc_A_T_C0_new_v7_custom1 = apply(nuc_C0_new_v7_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v7_custom1))
})

random_A_T_C0_new_v7_q1 = apply(random_C0_new_v7_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v7_q1))
})
random_A_T_C0_new_v7_q2 = apply(random_C0_new_v7_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v7_q2))
})
random_A_T_C0_new_v7_q3 = apply(random_C0_new_v7_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v7_q3))
})
random_A_T_C0_new_v7_q4 = apply(random_C0_new_v7_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v7_q4))
})
random_A_T_C0_new_v7_custom1 = apply(random_C0_new_v7_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v7_custom1))
})

tiling_A_T_C0_new_v7_q1 = apply(tiling_C0_new_v7_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v7_q1))
})
tiling_A_T_C0_new_v7_q2 = apply(tiling_C0_new_v7_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v7_q2))
})
tiling_A_T_C0_new_v7_q3 = apply(tiling_C0_new_v7_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v7_q3))
})
tiling_A_T_C0_new_v7_q4 = apply(tiling_C0_new_v7_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v7_q4))
})
tiling_A_T_C0_new_v7_custom1 = apply(tiling_C0_new_v7_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v7_custom1))
})

chrv_A_T_C0_new_v7_q1 = apply(chrv_C0_new_v7_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v7_q1))
})
chrv_A_T_C0_new_v7_q2 = apply(chrv_C0_new_v7_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v7_q2))
})
chrv_A_T_C0_new_v7_q3 = apply(chrv_C0_new_v7_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v7_q3))
})
chrv_A_T_C0_new_v7_q4 = apply(chrv_C0_new_v7_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v7_q4))
})
chrv_A_T_C0_new_v7_custom1 = apply(chrv_C0_new_v7_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v7_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v7_q = data.frame(q1 = nuc_A_T_C0_new_v7_q1,
                                 q2 = nuc_A_T_C0_new_v7_q2,
                                 q3 = nuc_A_T_C0_new_v7_q3,
                                 q4 = nuc_A_T_C0_new_v7_q4, 
                                 custom1 = nuc_A_T_C0_new_v7_custom1)

random_A_T_C0_new_v7_q = data.frame(q1 = random_A_T_C0_new_v7_q1,
                                    q2 = random_A_T_C0_new_v7_q2,
                                    q3 = random_A_T_C0_new_v7_q3,
                                    q4 = random_A_T_C0_new_v7_q4, 
                                    custom1 = random_A_T_C0_new_v7_custom1)

tiling_A_T_C0_new_v7_q = data.frame(q1 = tiling_A_T_C0_new_v7_q1,
                                    q2 = tiling_A_T_C0_new_v7_q2,
                                    q3 = tiling_A_T_C0_new_v7_q3,
                                    q4 = tiling_A_T_C0_new_v7_q4, 
                                    custom1 = tiling_A_T_C0_new_v7_custom1)

chrv_A_T_C0_new_v7_q = data.frame(q1 = chrv_A_T_C0_new_v7_q1,
                                  q2 = chrv_A_T_C0_new_v7_q2,
                                  q3 = chrv_A_T_C0_new_v7_q3,
                                  q4 = chrv_A_T_C0_new_v7_q4, 
                                  custom1 = chrv_A_T_C0_new_v7_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v8
nuc_A_T_C0_new_v8_q1 = apply(nuc_C0_new_v8_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v8_q1))
})
nuc_A_T_C0_new_v8_q2 = apply(nuc_C0_new_v8_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v8_q2))
})
nuc_A_T_C0_new_v8_q3 = apply(nuc_C0_new_v8_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v8_q3))
})
nuc_A_T_C0_new_v8_q4 = apply(nuc_C0_new_v8_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v8_q4))
})
nuc_A_T_C0_new_v8_custom1 = apply(nuc_C0_new_v8_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v8_custom1))
})

random_A_T_C0_new_v8_q1 = apply(random_C0_new_v8_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v8_q1))
})
random_A_T_C0_new_v8_q2 = apply(random_C0_new_v8_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v8_q2))
})
random_A_T_C0_new_v8_q3 = apply(random_C0_new_v8_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v8_q3))
})
random_A_T_C0_new_v8_q4 = apply(random_C0_new_v8_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v8_q4))
})
random_A_T_C0_new_v8_custom1 = apply(random_C0_new_v8_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v8_custom1))
})

tiling_A_T_C0_new_v8_q1 = apply(tiling_C0_new_v8_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v8_q1))
})
tiling_A_T_C0_new_v8_q2 = apply(tiling_C0_new_v8_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v8_q2))
})
tiling_A_T_C0_new_v8_q3 = apply(tiling_C0_new_v8_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v8_q3))
})
tiling_A_T_C0_new_v8_q4 = apply(tiling_C0_new_v8_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v8_q4))
})
tiling_A_T_C0_new_v8_custom1 = apply(tiling_C0_new_v8_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v8_custom1))
})

chrv_A_T_C0_new_v8_q1 = apply(chrv_C0_new_v8_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v8_q1))
})
chrv_A_T_C0_new_v8_q2 = apply(chrv_C0_new_v8_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v8_q2))
})
chrv_A_T_C0_new_v8_q3 = apply(chrv_C0_new_v8_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v8_q3))
})
chrv_A_T_C0_new_v8_q4 = apply(chrv_C0_new_v8_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v8_q4))
})
chrv_A_T_C0_new_v8_custom1 = apply(chrv_C0_new_v8_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v8_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v8_q = data.frame(q1 = nuc_A_T_C0_new_v8_q1,
                                 q2 = nuc_A_T_C0_new_v8_q2,
                                 q3 = nuc_A_T_C0_new_v8_q3,
                                 q4 = nuc_A_T_C0_new_v8_q4, 
                                 custom1 = nuc_A_T_C0_new_v8_custom1)

random_A_T_C0_new_v8_q = data.frame(q1 = random_A_T_C0_new_v8_q1,
                                    q2 = random_A_T_C0_new_v8_q2,
                                    q3 = random_A_T_C0_new_v8_q3,
                                    q4 = random_A_T_C0_new_v8_q4, 
                                    custom1 = random_A_T_C0_new_v8_custom1)

tiling_A_T_C0_new_v8_q = data.frame(q1 = tiling_A_T_C0_new_v8_q1,
                                    q2 = tiling_A_T_C0_new_v8_q2,
                                    q3 = tiling_A_T_C0_new_v8_q3,
                                    q4 = tiling_A_T_C0_new_v8_q4, 
                                    custom1 = tiling_A_T_C0_new_v8_custom1)

chrv_A_T_C0_new_v8_q = data.frame(q1 = chrv_A_T_C0_new_v8_q1,
                                  q2 = chrv_A_T_C0_new_v8_q2,
                                  q3 = chrv_A_T_C0_new_v8_q3,
                                  q4 = chrv_A_T_C0_new_v8_q4, 
                                  custom1 = chrv_A_T_C0_new_v8_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v9
nuc_A_T_C0_new_v9_q1 = apply(nuc_C0_new_v9_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v9_q1))
})
nuc_A_T_C0_new_v9_q2 = apply(nuc_C0_new_v9_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v9_q2))
})
nuc_A_T_C0_new_v9_q3 = apply(nuc_C0_new_v9_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v9_q3))
})
nuc_A_T_C0_new_v9_q4 = apply(nuc_C0_new_v9_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v9_q4))
})
nuc_A_T_C0_new_v9_custom1 = apply(nuc_C0_new_v9_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v9_custom1))
})

random_A_T_C0_new_v9_q1 = apply(random_C0_new_v9_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v9_q1))
})
random_A_T_C0_new_v9_q2 = apply(random_C0_new_v9_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v9_q2))
})
random_A_T_C0_new_v9_q3 = apply(random_C0_new_v9_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v9_q3))
})
random_A_T_C0_new_v9_q4 = apply(random_C0_new_v9_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v9_q4))
})
random_A_T_C0_new_v9_custom1 = apply(random_C0_new_v9_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v9_custom1))
})

tiling_A_T_C0_new_v9_q1 = apply(tiling_C0_new_v9_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v9_q1))
})
tiling_A_T_C0_new_v9_q2 = apply(tiling_C0_new_v9_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v9_q2))
})
tiling_A_T_C0_new_v9_q3 = apply(tiling_C0_new_v9_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v9_q3))
})
tiling_A_T_C0_new_v9_q4 = apply(tiling_C0_new_v9_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v9_q4))
})
tiling_A_T_C0_new_v9_custom1 = apply(tiling_C0_new_v9_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v9_custom1))
})

chrv_A_T_C0_new_v9_q1 = apply(chrv_C0_new_v9_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v9_q1))
})
chrv_A_T_C0_new_v9_q2 = apply(chrv_C0_new_v9_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v9_q2))
})
chrv_A_T_C0_new_v9_q3 = apply(chrv_C0_new_v9_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v9_q3))
})
chrv_A_T_C0_new_v9_q4 = apply(chrv_C0_new_v9_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v9_q4))
})
chrv_A_T_C0_new_v9_custom1 = apply(chrv_C0_new_v9_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v9_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v9_q = data.frame(q1 = nuc_A_T_C0_new_v9_q1,
                                 q2 = nuc_A_T_C0_new_v9_q2,
                                 q3 = nuc_A_T_C0_new_v9_q3,
                                 q4 = nuc_A_T_C0_new_v9_q4, 
                                 custom1 = nuc_A_T_C0_new_v9_custom1)

random_A_T_C0_new_v9_q = data.frame(q1 = random_A_T_C0_new_v9_q1,
                                    q2 = random_A_T_C0_new_v9_q2,
                                    q3 = random_A_T_C0_new_v9_q3,
                                    q4 = random_A_T_C0_new_v9_q4, 
                                    custom1 = random_A_T_C0_new_v9_custom1)

tiling_A_T_C0_new_v9_q = data.frame(q1 = tiling_A_T_C0_new_v9_q1,
                                    q2 = tiling_A_T_C0_new_v9_q2,
                                    q3 = tiling_A_T_C0_new_v9_q3,
                                    q4 = tiling_A_T_C0_new_v9_q4, 
                                    custom1 = tiling_A_T_C0_new_v9_custom1)

chrv_A_T_C0_new_v9_q = data.frame(q1 = chrv_A_T_C0_new_v9_q1,
                                  q2 = chrv_A_T_C0_new_v9_q2,
                                  q3 = chrv_A_T_C0_new_v9_q3,
                                  q4 = chrv_A_T_C0_new_v9_q4, 
                                  custom1 = chrv_A_T_C0_new_v9_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v10
nuc_A_T_C0_new_v10_q1 = apply(nuc_C0_new_v10_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v10_q1))
})
nuc_A_T_C0_new_v10_q2 = apply(nuc_C0_new_v10_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v10_q2))
})
nuc_A_T_C0_new_v10_q3 = apply(nuc_C0_new_v10_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v10_q3))
})
nuc_A_T_C0_new_v10_q4 = apply(nuc_C0_new_v10_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v10_q4))
})
nuc_A_T_C0_new_v10_custom1 = apply(nuc_C0_new_v10_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v10_custom1))
})

random_A_T_C0_new_v10_q1 = apply(random_C0_new_v10_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v10_q1))
})
random_A_T_C0_new_v10_q2 = apply(random_C0_new_v10_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v10_q2))
})
random_A_T_C0_new_v10_q3 = apply(random_C0_new_v10_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v10_q3))
})
random_A_T_C0_new_v10_q4 = apply(random_C0_new_v10_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v10_q4))
})
random_A_T_C0_new_v10_custom1 = apply(random_C0_new_v10_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v10_custom1))
})

tiling_A_T_C0_new_v10_q1 = apply(tiling_C0_new_v10_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v10_q1))
})
tiling_A_T_C0_new_v10_q2 = apply(tiling_C0_new_v10_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v10_q2))
})
tiling_A_T_C0_new_v10_q3 = apply(tiling_C0_new_v10_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v10_q3))
})
tiling_A_T_C0_new_v10_q4 = apply(tiling_C0_new_v10_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v10_q4))
})
tiling_A_T_C0_new_v10_custom1 = apply(tiling_C0_new_v10_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v10_custom1))
})

chrv_A_T_C0_new_v10_q1 = apply(chrv_C0_new_v10_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v10_q1))
})
chrv_A_T_C0_new_v10_q2 = apply(chrv_C0_new_v10_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v10_q2))
})
chrv_A_T_C0_new_v10_q3 = apply(chrv_C0_new_v10_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v10_q3))
})
chrv_A_T_C0_new_v10_q4 = apply(chrv_C0_new_v10_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v10_q4))
})
chrv_A_T_C0_new_v10_custom1 = apply(chrv_C0_new_v10_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v10_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v10_q = data.frame(q1 = nuc_A_T_C0_new_v10_q1,
                                 q2 = nuc_A_T_C0_new_v10_q2,
                                 q3 = nuc_A_T_C0_new_v10_q3,
                                 q4 = nuc_A_T_C0_new_v10_q4, 
                                 custom1 = nuc_A_T_C0_new_v10_custom1)

random_A_T_C0_new_v10_q = data.frame(q1 = random_A_T_C0_new_v10_q1,
                                    q2 = random_A_T_C0_new_v10_q2,
                                    q3 = random_A_T_C0_new_v10_q3,
                                    q4 = random_A_T_C0_new_v10_q4, 
                                    custom1 = random_A_T_C0_new_v10_custom1)

tiling_A_T_C0_new_v10_q = data.frame(q1 = tiling_A_T_C0_new_v10_q1,
                                    q2 = tiling_A_T_C0_new_v10_q2,
                                    q3 = tiling_A_T_C0_new_v10_q3,
                                    q4 = tiling_A_T_C0_new_v10_q4, 
                                    custom1 = tiling_A_T_C0_new_v10_custom1)

chrv_A_T_C0_new_v10_q = data.frame(q1 = chrv_A_T_C0_new_v10_q1,
                                  q2 = chrv_A_T_C0_new_v10_q2,
                                  q3 = chrv_A_T_C0_new_v10_q3,
                                  q4 = chrv_A_T_C0_new_v10_q4, 
                                  custom1 = chrv_A_T_C0_new_v10_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v11
nuc_A_T_C0_new_v11_q1 = apply(nuc_C0_new_v11_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v11_q1))
})
nuc_A_T_C0_new_v11_q2 = apply(nuc_C0_new_v11_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v11_q2))
})
nuc_A_T_C0_new_v11_q3 = apply(nuc_C0_new_v11_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v11_q3))
})
nuc_A_T_C0_new_v11_q4 = apply(nuc_C0_new_v11_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v11_q4))
})
nuc_A_T_C0_new_v11_custom1 = apply(nuc_C0_new_v11_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v11_custom1))
})

random_A_T_C0_new_v11_q1 = apply(random_C0_new_v11_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v11_q1))
})
random_A_T_C0_new_v11_q2 = apply(random_C0_new_v11_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v11_q2))
})
random_A_T_C0_new_v11_q3 = apply(random_C0_new_v11_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v11_q3))
})
random_A_T_C0_new_v11_q4 = apply(random_C0_new_v11_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v11_q4))
})
random_A_T_C0_new_v11_custom1 = apply(random_C0_new_v11_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v11_custom1))
})

tiling_A_T_C0_new_v11_q1 = apply(tiling_C0_new_v11_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v11_q1))
})
tiling_A_T_C0_new_v11_q2 = apply(tiling_C0_new_v11_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v11_q2))
})
tiling_A_T_C0_new_v11_q3 = apply(tiling_C0_new_v11_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v11_q3))
})
tiling_A_T_C0_new_v11_q4 = apply(tiling_C0_new_v11_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v11_q4))
})
tiling_A_T_C0_new_v11_custom1 = apply(tiling_C0_new_v11_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v11_custom1))
})

chrv_A_T_C0_new_v11_q1 = apply(chrv_C0_new_v11_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v11_q1))
})
chrv_A_T_C0_new_v11_q2 = apply(chrv_C0_new_v11_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v11_q2))
})
chrv_A_T_C0_new_v11_q3 = apply(chrv_C0_new_v11_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v11_q3))
})
chrv_A_T_C0_new_v11_q4 = apply(chrv_C0_new_v11_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v11_q4))
})
chrv_A_T_C0_new_v11_custom1 = apply(chrv_C0_new_v11_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v11_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v11_q = data.frame(q1 = nuc_A_T_C0_new_v11_q1,
                                 q2 = nuc_A_T_C0_new_v11_q2,
                                 q3 = nuc_A_T_C0_new_v11_q3,
                                 q4 = nuc_A_T_C0_new_v11_q4, 
                                 custom1 = nuc_A_T_C0_new_v11_custom1)

random_A_T_C0_new_v11_q = data.frame(q1 = random_A_T_C0_new_v11_q1,
                                    q2 = random_A_T_C0_new_v11_q2,
                                    q3 = random_A_T_C0_new_v11_q3,
                                    q4 = random_A_T_C0_new_v11_q4, 
                                    custom1 = random_A_T_C0_new_v11_custom1)

tiling_A_T_C0_new_v11_q = data.frame(q1 = tiling_A_T_C0_new_v11_q1,
                                    q2 = tiling_A_T_C0_new_v11_q2,
                                    q3 = tiling_A_T_C0_new_v11_q3,
                                    q4 = tiling_A_T_C0_new_v11_q4, 
                                    custom1 = tiling_A_T_C0_new_v11_custom1)

chrv_A_T_C0_new_v11_q = data.frame(q1 = chrv_A_T_C0_new_v11_q1,
                                  q2 = chrv_A_T_C0_new_v11_q2,
                                  q3 = chrv_A_T_C0_new_v11_q3,
                                  q4 = chrv_A_T_C0_new_v11_q4, 
                                  custom1 = chrv_A_T_C0_new_v11_custom1)


# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v12
nuc_A_T_C0_new_v12_q1 = apply(nuc_C0_new_v12_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v12_q1))
})
nuc_A_T_C0_new_v12_q2 = apply(nuc_C0_new_v12_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v12_q2))
})
nuc_A_T_C0_new_v12_q3 = apply(nuc_C0_new_v12_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v12_q3))
})
nuc_A_T_C0_new_v12_q4 = apply(nuc_C0_new_v12_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v12_q4))
})
nuc_A_T_C0_new_v12_custom1 = apply(nuc_C0_new_v12_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v12_custom1))
})

random_A_T_C0_new_v12_q1 = apply(random_C0_new_v12_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v12_q1))
})
random_A_T_C0_new_v12_q2 = apply(random_C0_new_v12_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v12_q2))
})
random_A_T_C0_new_v12_q3 = apply(random_C0_new_v12_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v12_q3))
})
random_A_T_C0_new_v12_q4 = apply(random_C0_new_v12_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v12_q4))
})
random_A_T_C0_new_v12_custom1 = apply(random_C0_new_v12_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v12_custom1))
})

tiling_A_T_C0_new_v12_q1 = apply(tiling_C0_new_v12_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v12_q1))
})
tiling_A_T_C0_new_v12_q2 = apply(tiling_C0_new_v12_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v12_q2))
})
tiling_A_T_C0_new_v12_q3 = apply(tiling_C0_new_v12_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v12_q3))
})
tiling_A_T_C0_new_v12_q4 = apply(tiling_C0_new_v12_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v12_q4))
})
tiling_A_T_C0_new_v12_custom1 = apply(tiling_C0_new_v12_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v12_custom1))
})

chrv_A_T_C0_new_v12_q1 = apply(chrv_C0_new_v12_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v12_q1))
})
chrv_A_T_C0_new_v12_q2 = apply(chrv_C0_new_v12_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v12_q2))
})
chrv_A_T_C0_new_v12_q3 = apply(chrv_C0_new_v12_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v12_q3))
})
chrv_A_T_C0_new_v12_q4 = apply(chrv_C0_new_v12_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v12_q4))
})
chrv_A_T_C0_new_v12_custom1 = apply(chrv_C0_new_v12_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v12_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v12_q = data.frame(q1 = nuc_A_T_C0_new_v12_q1,
                                 q2 = nuc_A_T_C0_new_v12_q2,
                                 q3 = nuc_A_T_C0_new_v12_q3,
                                 q4 = nuc_A_T_C0_new_v12_q4, 
                                 custom1 = nuc_A_T_C0_new_v12_custom1)

random_A_T_C0_new_v12_q = data.frame(q1 = random_A_T_C0_new_v12_q1,
                                    q2 = random_A_T_C0_new_v12_q2,
                                    q3 = random_A_T_C0_new_v12_q3,
                                    q4 = random_A_T_C0_new_v12_q4, 
                                    custom1 = random_A_T_C0_new_v12_custom1)

tiling_A_T_C0_new_v12_q = data.frame(q1 = tiling_A_T_C0_new_v12_q1,
                                    q2 = tiling_A_T_C0_new_v12_q2,
                                    q3 = tiling_A_T_C0_new_v12_q3,
                                    q4 = tiling_A_T_C0_new_v12_q4, 
                                    custom1 = tiling_A_T_C0_new_v12_custom1)

chrv_A_T_C0_new_v12_q = data.frame(q1 = chrv_A_T_C0_new_v12_q1,
                                  q2 = chrv_A_T_C0_new_v12_q2,
                                  q3 = chrv_A_T_C0_new_v12_q3,
                                  q4 = chrv_A_T_C0_new_v12_q4, 
                                  custom1 = chrv_A_T_C0_new_v12_custom1)



# Find the relative frequencies of A and T at each position (1-50) for each quartile + custom/library for C0_new_v13
nuc_A_T_C0_new_v13_q1 = apply(nuc_C0_new_v13_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v13_q1))
})
nuc_A_T_C0_new_v13_q2 = apply(nuc_C0_new_v13_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v13_q2))
})
nuc_A_T_C0_new_v13_q3 = apply(nuc_C0_new_v13_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v13_q3))
})
nuc_A_T_C0_new_v13_q4 = apply(nuc_C0_new_v13_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v13_q4))
})
nuc_A_T_C0_new_v13_custom1 = apply(nuc_C0_new_v13_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(nuc_C0_new_v13_custom1))
})

random_A_T_C0_new_v13_q1 = apply(random_C0_new_v13_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v13_q1))
})
random_A_T_C0_new_v13_q2 = apply(random_C0_new_v13_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v13_q2))
})
random_A_T_C0_new_v13_q3 = apply(random_C0_new_v13_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v13_q3))
})
random_A_T_C0_new_v13_q4 = apply(random_C0_new_v13_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v13_q4))
})
random_A_T_C0_new_v13_custom1 = apply(random_C0_new_v13_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(random_C0_new_v13_custom1))
})

tiling_A_T_C0_new_v13_q1 = apply(tiling_C0_new_v13_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v13_q1))
})
tiling_A_T_C0_new_v13_q2 = apply(tiling_C0_new_v13_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v13_q2))
})
tiling_A_T_C0_new_v13_q3 = apply(tiling_C0_new_v13_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v13_q3))
})
tiling_A_T_C0_new_v13_q4 = apply(tiling_C0_new_v13_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v13_q4))
})
tiling_A_T_C0_new_v13_custom1 = apply(tiling_C0_new_v13_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(tiling_C0_new_v13_custom1))
})

chrv_A_T_C0_new_v13_q1 = apply(chrv_C0_new_v13_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v13_q1))
})
chrv_A_T_C0_new_v13_q2 = apply(chrv_C0_new_v13_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v13_q2))
})
chrv_A_T_C0_new_v13_q3 = apply(chrv_C0_new_v13_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v13_q3))
})
chrv_A_T_C0_new_v13_q4 = apply(chrv_C0_new_v13_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v13_q4))
})
chrv_A_T_C0_new_v13_custom1 = apply(chrv_C0_new_v13_custom1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(chrv_C0_new_v13_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
nuc_A_T_C0_new_v13_q = data.frame(q1 = nuc_A_T_C0_new_v13_q1,
                                  q2 = nuc_A_T_C0_new_v13_q2,
                                  q3 = nuc_A_T_C0_new_v13_q3,
                                  q4 = nuc_A_T_C0_new_v13_q4, 
                                  custom1 = nuc_A_T_C0_new_v13_custom1)

random_A_T_C0_new_v13_q = data.frame(q1 = random_A_T_C0_new_v13_q1,
                                     q2 = random_A_T_C0_new_v13_q2,
                                     q3 = random_A_T_C0_new_v13_q3,
                                     q4 = random_A_T_C0_new_v13_q4, 
                                     custom1 = random_A_T_C0_new_v13_custom1)

tiling_A_T_C0_new_v13_q = data.frame(q1 = tiling_A_T_C0_new_v13_q1,
                                     q2 = tiling_A_T_C0_new_v13_q2,
                                     q3 = tiling_A_T_C0_new_v13_q3,
                                     q4 = tiling_A_T_C0_new_v13_q4, 
                                     custom1 = tiling_A_T_C0_new_v13_custom1)

chrv_A_T_C0_new_v13_q = data.frame(q1 = chrv_A_T_C0_new_v13_q1,
                                   q2 = chrv_A_T_C0_new_v13_q2,
                                   q3 = chrv_A_T_C0_new_v13_q3,
                                   q4 = chrv_A_T_C0_new_v13_q4, 
                                   custom1 = chrv_A_T_C0_new_v13_custom1)




# Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Nucleosome Library

nuc_AT_all_C_max = max(nuc_A_T_C26_q, nuc_A_T_C29_q, nuc_A_T_C31_q, nuc_A_T_C0_q,
                       nuc_A_T_C0_new_v1_q, nuc_A_T_C0_new_v2_q, nuc_A_T_C0_new_v3_q, nuc_A_T_C0_new_v4_q,
                       nuc_A_T_C0_new_v5_q, nuc_A_T_C0_new_v6_q, nuc_A_T_C0_new_v7_q, nuc_A_T_C0_new_v8_q,
                       nuc_A_T_C0_new_v9_q, nuc_A_T_C0_new_v10_q, nuc_A_T_C0_new_v11_q, nuc_A_T_C0_new_v12_q,
                       nuc_A_T_C0_new_v13_q)

nuc_AT_all_C_min = min(nuc_A_T_C26_q, nuc_A_T_C29_q, nuc_A_T_C31_q, nuc_A_T_C0_q,
                       nuc_A_T_C0_new_v1_q, nuc_A_T_C0_new_v2_q, nuc_A_T_C0_new_v3_q, nuc_A_T_C0_new_v4_q,
                       nuc_A_T_C0_new_v5_q, nuc_A_T_C0_new_v6_q, nuc_A_T_C0_new_v7_q, nuc_A_T_C0_new_v8_q,
                       nuc_A_T_C0_new_v9_q, nuc_A_T_C0_new_v10_q, nuc_A_T_C0_new_v11_q, nuc_A_T_C0_new_v12_q,
                       nuc_A_T_C0_new_v13_q)

nuc_AT_q1vq4_C26 = ggplot(data = nuc_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C26
ggsave("figures/periodicity/nuc_AT_q1vq4_C26.png", plot = nuc_AT_q1vq4_C26)

nuc_AT_q1vq4_C29 = ggplot(data = nuc_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C29
ggsave("figures/periodicity/nuc_AT_q1vq4_C29.png", plot = nuc_AT_q1vq4_C29)

nuc_AT_q1vq4_C31 = ggplot(data = nuc_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C31
ggsave("figures/periodicity/nuc_AT_q1vq4_C31.png", plot = nuc_AT_q1vq4_C31)

nuc_AT_q1vq4_C0 = ggplot(data = nuc_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0
ggsave("figures/periodicity/nuc_AT_q1vq4_C0.png", plot = nuc_AT_q1vq4_C0)

nuc_AT_q1vq4_C0_new_v1 = ggplot(data = nuc_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v1
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v1.png", plot = nuc_AT_q1vq4_C0_new_v1)

nuc_AT_q1vq4_C0_new_v2 = ggplot(data = nuc_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v2
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v2.png", plot = nuc_AT_q1vq4_C0_new_v2)

nuc_AT_q1vq4_C0_new_v3 = ggplot(data = nuc_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v3
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v3.png", plot = nuc_AT_q1vq4_C0_new_v3)

nuc_AT_q1vq4_C0_new_v4 = ggplot(data = nuc_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v4
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v4.png", plot = nuc_AT_q1vq4_C0_new_v4)

nuc_AT_q1vq4_C0_new_v5 = ggplot(data = nuc_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v5
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v5.png", plot = nuc_AT_q1vq4_C0_new_v5)

nuc_AT_q1vq4_C0_new_v6 = ggplot(data = nuc_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v6
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v6.png", plot = nuc_AT_q1vq4_C0_new_v6)

nuc_AT_q1vq4_C0_new_v7 = ggplot(data = nuc_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v7
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v7.png", plot = nuc_AT_q1vq4_C0_new_v7)

nuc_AT_q1vq4_C0_new_v8 = ggplot(data = nuc_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v8
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v8.png", plot = nuc_AT_q1vq4_C0_new_v8)

nuc_AT_q1vq4_C0_new_v9 = ggplot(data = nuc_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v9
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v9.png", plot = nuc_AT_q1vq4_C0_new_v9)

nuc_AT_q1vq4_C0_new_v10 = ggplot(data = nuc_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v10
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v10.png", plot = nuc_AT_q1vq4_C0_new_v10)

nuc_AT_q1vq4_C0_new_v11 = ggplot(data = nuc_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v11
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v11.png", plot = nuc_AT_q1vq4_C0_new_v11)

nuc_AT_q1vq4_C0_new_v12 = ggplot(data = nuc_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v12
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v12.png", plot = nuc_AT_q1vq4_C0_new_v12)

nuc_AT_q1vq4_C0_new_v13 = ggplot(data = nuc_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q1vq4_C0_new_v13
ggsave("figures/periodicity/nuc_AT_q1vq4_C0_new_v13.png", plot = nuc_AT_q1vq4_C0_new_v13)


nuc_AT_q2vq3_C26 = ggplot(data = nuc_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C26
ggsave("figures/periodicity/nuc_AT_q2vq3_C26.png", plot = nuc_AT_q2vq3_C26)

nuc_AT_q2vq3_C29 = ggplot(data = nuc_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C29
ggsave("figures/periodicity/nuc_AT_q2vq3_C29.png", plot = nuc_AT_q2vq3_C29)

nuc_AT_q2vq3_C31 = ggplot(data = nuc_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C31
ggsave("figures/periodicity/nuc_AT_q2vq3_C31.png", plot = nuc_AT_q2vq3_C31)

nuc_AT_q2vq3_C0 = ggplot(data = nuc_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0
ggsave("figures/periodicity/nuc_AT_q2vq3_C0.png", plot = nuc_AT_q2vq3_C0)

nuc_AT_q2vq3_C0_new_v1 = ggplot(data = nuc_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v1
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v1.png", plot = nuc_AT_q2vq3_C0_new_v1)

nuc_AT_q2vq3_C0_new_v2 = ggplot(data = nuc_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v2
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v2.png", plot = nuc_AT_q2vq3_C0_new_v2)

nuc_AT_q2vq3_C0_new_v3 = ggplot(data = nuc_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v3
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v3.png", plot = nuc_AT_q2vq3_C0_new_v3)

nuc_AT_q2vq3_C0_new_v4 = ggplot(data = nuc_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v4
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v4.png", plot = nuc_AT_q2vq3_C0_new_v4)

nuc_AT_q2vq3_C0_new_v5 = ggplot(data = nuc_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v5
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v5.png", plot = nuc_AT_q2vq3_C0_new_v5)

nuc_AT_q2vq3_C0_new_v6 = ggplot(data = nuc_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v6
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v6.png", plot = nuc_AT_q2vq3_C0_new_v6)

nuc_AT_q2vq3_C0_new_v7 = ggplot(data = nuc_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v7
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v7.png", plot = nuc_AT_q2vq3_C0_new_v7)

nuc_AT_q2vq3_C0_new_v8 = ggplot(data = nuc_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v8
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v8.png", plot = nuc_AT_q2vq3_C0_new_v8)

nuc_AT_q2vq3_C0_new_v9 = ggplot(data = nuc_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v9
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v9.png", plot = nuc_AT_q2vq3_C0_new_v9)

nuc_AT_q2vq3_C0_new_v10 = ggplot(data = nuc_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v10
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v10.png", plot = nuc_AT_q2vq3_C0_new_v10)

nuc_AT_q2vq3_C0_new_v11 = ggplot(data = nuc_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v11
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v11.png", plot = nuc_AT_q2vq3_C0_new_v11)

nuc_AT_q2vq3_C0_new_v12 = ggplot(data = nuc_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v12
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v12.png", plot = nuc_AT_q2vq3_C0_new_v12)

nuc_AT_q2vq3_C0_new_v13 = ggplot(data = nuc_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_q2vq3_C0_new_v13
ggsave("figures/periodicity/nuc_AT_q2vq3_C0_new_v13.png", plot = nuc_AT_q2vq3_C0_new_v13)


nuc_AT_custom1_C26 = ggplot(data = nuc_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C26
ggsave("figures/periodicity/nuc_AT_custom1_C26.png", plot = nuc_AT_custom1_C26)

nuc_AT_custom1_C29 = ggplot(data = nuc_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C29
ggsave("figures/periodicity/nuc_AT_custom1_C29.png", plot = nuc_AT_custom1_C29)

nuc_AT_custom1_C31 = ggplot(data = nuc_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C31
ggsave("figures/periodicity/nuc_AT_custom1_C31.png", plot = nuc_AT_custom1_C31)

nuc_AT_custom1_C0 = ggplot(data = nuc_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0
ggsave("figures/periodicity/nuc_AT_custom1_C0.png", plot = nuc_AT_custom1_C0)

nuc_AT_custom1_C0_new_v1 = ggplot(data = nuc_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v1
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v1.png", plot = nuc_AT_custom1_C0_new_v1)

nuc_AT_custom1_C0_new_v2 = ggplot(data = nuc_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v2
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v2.png", plot = nuc_AT_custom1_C0_new_v2)

nuc_AT_custom1_C0_new_v3 = ggplot(data = nuc_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v3
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v3.png", plot = nuc_AT_custom1_C0_new_v3)

nuc_AT_custom1_C0_new_v4 = ggplot(data = nuc_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v4
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v4.png", plot = nuc_AT_custom1_C0_new_v4)

nuc_AT_custom1_C0_new_v5 = ggplot(data = nuc_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v5
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v5.png", plot = nuc_AT_custom1_C0_new_v5)

nuc_AT_custom1_C0_new_v6 = ggplot(data = nuc_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v6
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v6.png", plot = nuc_AT_custom1_C0_new_v6)

nuc_AT_custom1_C0_new_v7 = ggplot(data = nuc_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v7
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v7.png", plot = nuc_AT_custom1_C0_new_v7)

nuc_AT_custom1_C0_new_v8 = ggplot(data = nuc_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v8
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v8.png", plot = nuc_AT_custom1_C0_new_v8)

nuc_AT_custom1_C0_new_v9 = ggplot(data = nuc_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v9
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v9.png", plot = nuc_AT_custom1_C0_new_v9)

nuc_AT_custom1_C0_new_v10 = ggplot(data = nuc_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v10
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v10.png", plot = nuc_AT_custom1_C0_new_v10)

nuc_AT_custom1_C0_new_v11 = ggplot(data = nuc_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v11
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v11.png", plot = nuc_AT_custom1_C0_new_v11)

nuc_AT_custom1_C0_new_v12 = ggplot(data = nuc_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v12
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v12.png", plot = nuc_AT_custom1_C0_new_v12)

nuc_AT_custom1_C0_new_v13 = ggplot(data = nuc_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AT_all_C_min, nuc_AT_all_C_max)
nuc_AT_custom1_C0_new_v13
ggsave("figures/periodicity/nuc_AT_custom1_C0_new_v13.png", plot = nuc_AT_custom1_C0_new_v13)




# Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Random Library

random_AT_all_C_max = max(random_A_T_C26_q, random_A_T_C29_q, random_A_T_C31_q, random_A_T_C0_q,
                          random_A_T_C0_new_v1_q, random_A_T_C0_new_v2_q, random_A_T_C0_new_v3_q, random_A_T_C0_new_v4_q,
                          random_A_T_C0_new_v5_q, random_A_T_C0_new_v6_q, random_A_T_C0_new_v7_q, random_A_T_C0_new_v8_q,
                          random_A_T_C0_new_v9_q, random_A_T_C0_new_v10_q, random_A_T_C0_new_v11_q, random_A_T_C0_new_v12_q,
                          random_A_T_C0_new_v13_q)
random_AT_all_C_min = min(random_A_T_C26_q, random_A_T_C29_q, random_A_T_C31_q, random_A_T_C0_q,
                          random_A_T_C0_new_v1_q, random_A_T_C0_new_v2_q, random_A_T_C0_new_v3_q, random_A_T_C0_new_v4_q,
                          random_A_T_C0_new_v5_q, random_A_T_C0_new_v6_q, random_A_T_C0_new_v7_q, random_A_T_C0_new_v8_q,
                          random_A_T_C0_new_v9_q, random_A_T_C0_new_v10_q, random_A_T_C0_new_v11_q, random_A_T_C0_new_v12_q,
                          random_A_T_C0_new_v13_q)

random_AT_q1vq4_C26 = ggplot(data = random_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C26") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C26
ggsave("figures/periodicity/random_AT_q1vq4_C26.png", plot = random_AT_q1vq4_C26)

random_AT_q1vq4_C29 = ggplot(data = random_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C29") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C29
ggsave("figures/periodicity/random_AT_q1vq4_C29.png", plot = random_AT_q1vq4_C29)

random_AT_q1vq4_C31 = ggplot(data = random_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C31") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C31
ggsave("figures/periodicity/random_AT_q1vq4_C31.png", plot = random_AT_q1vq4_C31)

random_AT_q1vq4_C0 = ggplot(data = random_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0
ggsave("figures/periodicity/random_AT_q1vq4_C0.png", plot = random_AT_q1vq4_C0)

random_AT_q1vq4_C0_new_v1 = ggplot(data = random_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v1") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v1
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v1.png", plot = random_AT_q1vq4_C0_new_v1)

random_AT_q1vq4_C0_new_v2 = ggplot(data = random_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v2") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v2
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v2.png", plot = random_AT_q1vq4_C0_new_v2)

random_AT_q1vq4_C0_new_v3 = ggplot(data = random_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v3") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v3
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v3.png", plot = random_AT_q1vq4_C0_new_v3)

random_AT_q1vq4_C0_new_v4 = ggplot(data = random_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v4") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v4
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v4.png", plot = random_AT_q1vq4_C0_new_v4)

random_AT_q1vq4_C0_new_v5 = ggplot(data = random_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v5") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v5
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v5.png", plot = random_AT_q1vq4_C0_new_v5)

random_AT_q1vq4_C0_new_v6 = ggplot(data = random_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v6") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v6
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v6.png", plot = random_AT_q1vq4_C0_new_v6)

random_AT_q1vq4_C0_new_v7 = ggplot(data = random_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v7") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v7
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v7.png", plot = random_AT_q1vq4_C0_new_v7)

random_AT_q1vq4_C0_new_v8 = ggplot(data = random_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v8") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v8
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v8.png", plot = random_AT_q1vq4_C0_new_v8)

random_AT_q1vq4_C0_new_v9 = ggplot(data = random_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v9") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v9
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v9.png", plot = random_AT_q1vq4_C0_new_v9)

random_AT_q1vq4_C0_new_v10 = ggplot(data = random_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v10") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v10
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v10.png", plot = random_AT_q1vq4_C0_new_v10)

random_AT_q1vq4_C0_new_v11 = ggplot(data = random_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v11") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v11
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v11.png", plot = random_AT_q1vq4_C0_new_v11)

random_AT_q1vq4_C0_new_v12 = ggplot(data = random_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v12") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v12
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v12.png", plot = random_AT_q1vq4_C0_new_v12)

random_AT_q1vq4_C0_new_v13 = ggplot(data = random_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v13") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q1vq4_C0_new_v13
ggsave("figures/periodicity/random_AT_q1vq4_C0_new_v13.png", plot = random_AT_q1vq4_C0_new_v13)


random_AT_q2vq3_C26 = ggplot(data = random_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C26") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C26
ggsave("figures/periodicity/random_AT_q2vq3_C26.png", plot = random_AT_q2vq3_C26)

random_AT_q2vq3_C29 = ggplot(data = random_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C29") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C29
ggsave("figures/periodicity/random_AT_q2vq3_C29.png", plot = random_AT_q2vq3_C29)

random_AT_q2vq3_C31 = ggplot(data = random_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C31") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C31
ggsave("figures/periodicity/random_AT_q2vq3_C31.png", plot = random_AT_q2vq3_C31)

random_AT_q2vq3_C0 = ggplot(data = random_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0
ggsave("figures/periodicity/random_AT_q2vq3_C0.png", plot = random_AT_q2vq3_C0)

random_AT_q2vq3_C0_new_v1 = ggplot(data = random_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v1") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v1
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v1.png", plot = random_AT_q2vq3_C0_new_v1)

random_AT_q2vq3_C0_new_v2 = ggplot(data = random_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v2") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v2
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v2.png", plot = random_AT_q2vq3_C0_new_v2)

random_AT_q2vq3_C0_new_v3 = ggplot(data = random_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v3") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v3
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v3.png", plot = random_AT_q2vq3_C0_new_v3)

random_AT_q2vq3_C0_new_v4 = ggplot(data = random_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v4") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v4
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v4.png", plot = random_AT_q2vq3_C0_new_v4)

random_AT_q2vq3_C0_new_v5 = ggplot(data = random_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v5") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v5
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v5.png", plot = random_AT_q2vq3_C0_new_v5)

random_AT_q2vq3_C0_new_v6 = ggplot(data = random_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v6") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v6
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v6.png", plot = random_AT_q2vq3_C0_new_v6)

random_AT_q2vq3_C0_new_v7 = ggplot(data = random_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v7") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v7
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v7.png", plot = random_AT_q2vq3_C0_new_v7)

random_AT_q2vq3_C0_new_v8 = ggplot(data = random_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v8") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v8
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v8.png", plot = random_AT_q2vq3_C0_new_v8)

random_AT_q2vq3_C0_new_v9 = ggplot(data = random_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v9") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v9
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v9.png", plot = random_AT_q2vq3_C0_new_v9)

random_AT_q2vq3_C0_new_v10 = ggplot(data = random_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v10") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v10
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v10.png", plot = random_AT_q2vq3_C0_new_v10)

random_AT_q2vq3_C0_new_v11 = ggplot(data = random_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v11") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v11
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v11.png", plot = random_AT_q2vq3_C0_new_v11)

random_AT_q2vq3_C0_new_v12 = ggplot(data = random_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v12") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v12
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v12.png", plot = random_AT_q2vq3_C0_new_v12)

random_AT_q2vq3_C0_new_v13 = ggplot(data = random_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v13") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_q2vq3_C0_new_v13
ggsave("figures/periodicity/random_AT_q2vq3_C0_new_v13.png", plot = random_AT_q2vq3_C0_new_v13)

random_AT_custom1_C26 = ggplot(data = random_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C26") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C26
ggsave("figures/periodicity/random_AT_custom1_C26.png", plot = random_AT_custom1_C26)

random_AT_custom1_C29 = ggplot(data = random_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C29") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C29
ggsave("figures/periodicity/random_AT_custom1_C29.png", plot = random_AT_custom1_C29)

random_AT_custom1_C31 = ggplot(data = random_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C31") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C31
ggsave("figures/periodicity/random_AT_custom1_C31.png", plot = random_AT_custom1_C31)

random_AT_custom1_C0 = ggplot(data = random_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0
ggsave("figures/periodicity/random_AT_custom1_C0.png", plot = random_AT_custom1_C0)

random_AT_custom1_C0_new_v1 = ggplot(data = random_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v1") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v1
ggsave("figures/periodicity/random_AT_custom1_C0_new_v1.png", plot = random_AT_custom1_C0_new_v1)

random_AT_custom1_C0_new_v2 = ggplot(data = random_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v2") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v2
ggsave("figures/periodicity/random_AT_custom1_C0_new_v2.png", plot = random_AT_custom1_C0_new_v2)

random_AT_custom1_C0_new_v3 = ggplot(data = random_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v3") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v3
ggsave("figures/periodicity/random_AT_custom1_C0_new_v3.png", plot = random_AT_custom1_C0_new_v3)

random_AT_custom1_C0_new_v4 = ggplot(data = random_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v4") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v4
ggsave("figures/periodicity/random_AT_custom1_C0_new_v4.png", plot = random_AT_custom1_C0_new_v4)

random_AT_custom1_C0_new_v5 = ggplot(data = random_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v5") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v5
ggsave("figures/periodicity/random_AT_custom1_C0_new_v5.png", plot = random_AT_custom1_C0_new_v5)

random_AT_custom1_C0_new_v6 = ggplot(data = random_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v6") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v6
ggsave("figures/periodicity/random_AT_custom1_C0_new_v6.png", plot = random_AT_custom1_C0_new_v6)

random_AT_custom1_C0_new_v7 = ggplot(data = random_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v7") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v7
ggsave("figures/periodicity/random_AT_custom1_C0_new_v7.png", plot = random_AT_custom1_C0_new_v7)

random_AT_custom1_C0_new_v8 = ggplot(data = random_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v8") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v8
ggsave("figures/periodicity/random_AT_custom1_C0_new_v8.png", plot = random_AT_custom1_C0_new_v8)

random_AT_custom1_C0_new_v9 = ggplot(data = random_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v9") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v9
ggsave("figures/periodicity/random_AT_custom1_C0_new_v9.png", plot = random_AT_custom1_C0_new_v9)

random_AT_custom1_C0_new_v10 = ggplot(data = random_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v10") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v10
ggsave("figures/periodicity/random_AT_custom1_C0_new_v10.png", plot = random_AT_custom1_C0_new_v10)

random_AT_custom1_C0_new_v11 = ggplot(data = random_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v11") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v11
ggsave("figures/periodicity/random_AT_custom1_C0_new_v11.png", plot = random_AT_custom1_C0_new_v11)

random_AT_custom1_C0_new_v12 = ggplot(data = random_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v12") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v12
ggsave("figures/periodicity/random_AT_custom1_C0_new_v12.png", plot = random_AT_custom1_C0_new_v12)

random_AT_custom1_C0_new_v13 = ggplot(data = random_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, C0_new_v13") +
  ylim(random_AT_all_C_min, random_AT_all_C_max)
random_AT_custom1_C0_new_v13
ggsave("figures/periodicity/random_AT_custom1_C0_new_v13.png", plot = random_AT_custom1_C0_new_v13)





# Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Tiling Library

tiling_AT_all_C_max = max(tiling_A_T_C26_q, tiling_A_T_C29_q, tiling_A_T_C31_q, tiling_A_T_C0_q,
                          tiling_A_T_C0_new_v1_q, tiling_A_T_C0_new_v2_q, tiling_A_T_C0_new_v3_q, tiling_A_T_C0_new_v4_q,
                          tiling_A_T_C0_new_v5_q, tiling_A_T_C0_new_v6_q, tiling_A_T_C0_new_v7_q, tiling_A_T_C0_new_v8_q,
                          tiling_A_T_C0_new_v9_q, tiling_A_T_C0_new_v10_q, tiling_A_T_C0_new_v11_q, tiling_A_T_C0_new_v12_q,
                          tiling_A_T_C0_new_v13_q)
tiling_AT_all_C_min = min(tiling_A_T_C26_q, tiling_A_T_C29_q, tiling_A_T_C31_q, tiling_A_T_C0_q,
                          tiling_A_T_C0_new_v1_q, tiling_A_T_C0_new_v2_q, tiling_A_T_C0_new_v3_q, tiling_A_T_C0_new_v4_q,
                          tiling_A_T_C0_new_v5_q, tiling_A_T_C0_new_v6_q, tiling_A_T_C0_new_v7_q, tiling_A_T_C0_new_v8_q,
                          tiling_A_T_C0_new_v9_q, tiling_A_T_C0_new_v10_q, tiling_A_T_C0_new_v11_q, tiling_A_T_C0_new_v12_q,
                          tiling_A_T_C0_new_v13_q)

tiling_AT_q1vq4_C26 = ggplot(data = tiling_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C26") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C26
ggsave("figures/periodicity/tiling_AT_q1vq4_C26.png", plot = tiling_AT_q1vq4_C26)

tiling_AT_q1vq4_C29 = ggplot(data = tiling_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C29") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C29
ggsave("figures/periodicity/tiling_AT_q1vq4_C29.png", plot = tiling_AT_q1vq4_C29)

tiling_AT_q1vq4_C31 = ggplot(data = tiling_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C31") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C31
ggsave("figures/periodicity/tiling_AT_q1vq4_C31.png", plot = tiling_AT_q1vq4_C31)

tiling_AT_q1vq4_C0 = ggplot(data = tiling_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0
ggsave("figures/periodicity/tiling_AT_q1vq4_C0.png", plot = tiling_AT_q1vq4_C0)

tiling_AT_q1vq4_C0_new_v1 = ggplot(data = tiling_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v1
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v1.png", plot = tiling_AT_q1vq4_C0_new_v1)

tiling_AT_q1vq4_C0_new_v2 = ggplot(data = tiling_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v2
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v2.png", plot = tiling_AT_q1vq4_C0_new_v2)

tiling_AT_q1vq4_C0_new_v3 = ggplot(data = tiling_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v3
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v3.png", plot = tiling_AT_q1vq4_C0_new_v3)

tiling_AT_q1vq4_C0_new_v4 = ggplot(data = tiling_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v4
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v4.png", plot = tiling_AT_q1vq4_C0_new_v4)

tiling_AT_q1vq4_C0_new_v5 = ggplot(data = tiling_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v5
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v5.png", plot = tiling_AT_q1vq4_C0_new_v5)

tiling_AT_q1vq4_C0_new_v6 = ggplot(data = tiling_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v6
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v6.png", plot = tiling_AT_q1vq4_C0_new_v6)

tiling_AT_q1vq4_C0_new_v7 = ggplot(data = tiling_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v7
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v7.png", plot = tiling_AT_q1vq4_C0_new_v7)

tiling_AT_q1vq4_C0_new_v8 = ggplot(data = tiling_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v8
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v8.png", plot = tiling_AT_q1vq4_C0_new_v8)

tiling_AT_q1vq4_C0_new_v9 = ggplot(data = tiling_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v9
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v9.png", plot = tiling_AT_q1vq4_C0_new_v9)

tiling_AT_q1vq4_C0_new_v10 = ggplot(data = tiling_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v10
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v10.png", plot = tiling_AT_q1vq4_C0_new_v10)

tiling_AT_q1vq4_C0_new_v11 = ggplot(data = tiling_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v11
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v11.png", plot = tiling_AT_q1vq4_C0_new_v11)

tiling_AT_q1vq4_C0_new_v12 = ggplot(data = tiling_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v12
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v12.png", plot = tiling_AT_q1vq4_C0_new_v12)

tiling_AT_q1vq4_C0_new_v13 = ggplot(data = tiling_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q1vq4_C0_new_v13
ggsave("figures/periodicity/tiling_AT_q1vq4_C0_new_v13.png", plot = tiling_AT_q1vq4_C0_new_v13)


tiling_AT_q2vq3_C26 = ggplot(data = tiling_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C26") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C26
ggsave("figures/periodicity/tiling_AT_q2vq3_C26.png", plot = tiling_AT_q2vq3_C26)

tiling_AT_q2vq3_C29 = ggplot(data = tiling_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C29") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C29
ggsave("figures/periodicity/tiling_AT_q2vq3_C29.png", plot = tiling_AT_q2vq3_C29)

tiling_AT_q2vq3_C31 = ggplot(data = tiling_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C31") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C31
ggsave("figures/periodicity/tiling_AT_q2vq3_C31.png", plot = tiling_AT_q2vq3_C31)

tiling_AT_q2vq3_C0 = ggplot(data = tiling_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0
ggsave("figures/periodicity/tiling_AT_q2vq3_C0.png", plot = tiling_AT_q2vq3_C0)

tiling_AT_q2vq3_C0_new_v1 = ggplot(data = tiling_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v1
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v1.png", plot = tiling_AT_q2vq3_C0_new_v1)

tiling_AT_q2vq3_C0_new_v2 = ggplot(data = tiling_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v2
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v2.png", plot = tiling_AT_q2vq3_C0_new_v2)

tiling_AT_q2vq3_C0_new_v3 = ggplot(data = tiling_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v3
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v3.png", plot = tiling_AT_q2vq3_C0_new_v3)

tiling_AT_q2vq3_C0_new_v4 = ggplot(data = tiling_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v4
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v4.png", plot = tiling_AT_q2vq3_C0_new_v4)

tiling_AT_q2vq3_C0_new_v5 = ggplot(data = tiling_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v5
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v5.png", plot = tiling_AT_q2vq3_C0_new_v5)

tiling_AT_q2vq3_C0_new_v6 = ggplot(data = tiling_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v6
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v6.png", plot = tiling_AT_q2vq3_C0_new_v6)

tiling_AT_q2vq3_C0_new_v7 = ggplot(data = tiling_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v7
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v7.png", plot = tiling_AT_q2vq3_C0_new_v7)

tiling_AT_q2vq3_C0_new_v8 = ggplot(data = tiling_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v8
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v8.png", plot = tiling_AT_q2vq3_C0_new_v8)

tiling_AT_q2vq3_C0_new_v9 = ggplot(data = tiling_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v9
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v9.png", plot = tiling_AT_q2vq3_C0_new_v9)

tiling_AT_q2vq3_C0_new_v10 = ggplot(data = tiling_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v10
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v10.png", plot = tiling_AT_q2vq3_C0_new_v10)

tiling_AT_q2vq3_C0_new_v11 = ggplot(data = tiling_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v11
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v11.png", plot = tiling_AT_q2vq3_C0_new_v11)

tiling_AT_q2vq3_C0_new_v12 = ggplot(data = tiling_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v12
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v12.png", plot = tiling_AT_q2vq3_C0_new_v12)

tiling_AT_q2vq3_C0_new_v13 = ggplot(data = tiling_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_q2vq3_C0_new_v13
ggsave("figures/periodicity/tiling_AT_q2vq3_C0_new_v13.png", plot = tiling_AT_q2vq3_C0_new_v13)



tiling_AT_custom1_C26 = ggplot(data = tiling_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C26") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C26
ggsave("figures/periodicity/tiling_AT_custom1_C26.png", plot = tiling_AT_custom1_C26)

tiling_AT_custom1_C29 = ggplot(data = tiling_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C29") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C29
ggsave("figures/periodicity/tiling_AT_custom1_C29.png", plot = tiling_AT_custom1_C29)

tiling_AT_custom1_C31 = ggplot(data = tiling_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C31") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C31
ggsave("figures/periodicity/tiling_AT_custom1_C31.png", plot = tiling_AT_custom1_C31)

tiling_AT_custom1_C0 = ggplot(data = tiling_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0
ggsave("figures/periodicity/tiling_AT_custom1_C0.png", plot = tiling_AT_custom1_C0)

tiling_AT_custom1_C0_new_v1 = ggplot(data = tiling_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v1
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v1.png", plot = tiling_AT_custom1_C0_new_v1)

tiling_AT_custom1_C0_new_v2 = ggplot(data = tiling_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v2
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v2.png", plot = tiling_AT_custom1_C0_new_v2)

tiling_AT_custom1_C0_new_v3 = ggplot(data = tiling_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v3
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v3.png", plot = tiling_AT_custom1_C0_new_v3)

tiling_AT_custom1_C0_new_v4 = ggplot(data = tiling_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v4
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v4.png", plot = tiling_AT_custom1_C0_new_v4)

tiling_AT_custom1_C0_new_v5 = ggplot(data = tiling_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v5
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v5.png", plot = tiling_AT_custom1_C0_new_v5)

tiling_AT_custom1_C0_new_v6 = ggplot(data = tiling_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v6
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v6.png", plot = tiling_AT_custom1_C0_new_v6)

tiling_AT_custom1_C0_new_v7 = ggplot(data = tiling_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v7
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v7.png", plot = tiling_AT_custom1_C0_new_v7)

tiling_AT_custom1_C0_new_v8 = ggplot(data = tiling_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v8
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v8.png", plot = tiling_AT_custom1_C0_new_v8)

tiling_AT_custom1_C0_new_v9 = ggplot(data = tiling_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v9
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v9.png", plot = tiling_AT_custom1_C0_new_v9)

tiling_AT_custom1_C0_new_v10 = ggplot(data = tiling_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v10
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v10.png", plot = tiling_AT_custom1_C0_new_v10)

tiling_AT_custom1_C0_new_v11 = ggplot(data = tiling_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v11
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v11.png", plot = tiling_AT_custom1_C0_new_v11)

tiling_AT_custom1_C0_new_v12 = ggplot(data = tiling_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v12
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v12.png", plot = tiling_AT_custom1_C0_new_v12)

tiling_AT_custom1_C0_new_v13 = ggplot(data = tiling_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AT_all_C_min, tiling_AT_all_C_max)
tiling_AT_custom1_C0_new_v13
ggsave("figures/periodicity/tiling_AT_custom1_C0_new_v13.png", plot = tiling_AT_custom1_C0_new_v13)





# Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for ChrV Library

chrv_AT_all_C_max = max(chrv_A_T_C26_q, chrv_A_T_C29_q, chrv_A_T_C31_q, chrv_A_T_C0_q,
                        chrv_A_T_C0_new_v1_q, chrv_A_T_C0_new_v2_q, chrv_A_T_C0_new_v3_q, chrv_A_T_C0_new_v4_q,
                        chrv_A_T_C0_new_v5_q, chrv_A_T_C0_new_v6_q, chrv_A_T_C0_new_v7_q, chrv_A_T_C0_new_v8_q,
                        chrv_A_T_C0_new_v9_q, chrv_A_T_C0_new_v10_q, chrv_A_T_C0_new_v11_q, chrv_A_T_C0_new_v12_q,
                        chrv_A_T_C0_new_v13_q)
chrv_AT_all_C_min = min(chrv_A_T_C26_q, chrv_A_T_C29_q, chrv_A_T_C31_q, chrv_A_T_C0_q,
                        chrv_A_T_C0_new_v1_q, chrv_A_T_C0_new_v2_q, chrv_A_T_C0_new_v3_q, chrv_A_T_C0_new_v4_q,
                        chrv_A_T_C0_new_v5_q, chrv_A_T_C0_new_v6_q, chrv_A_T_C0_new_v7_q, chrv_A_T_C0_new_v8_q,
                        chrv_A_T_C0_new_v9_q, chrv_A_T_C0_new_v10_q, chrv_A_T_C0_new_v11_q, chrv_A_T_C0_new_v12_q,
                        chrv_A_T_C0_new_v13_q)

chrv_AT_q1vq4_C26 = ggplot(data = chrv_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C26") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C26
ggsave("figures/periodicity/chrv_AT_q1vq4_C26.png", plot = chrv_AT_q1vq4_C26)

chrv_AT_q1vq4_C29 = ggplot(data = chrv_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C29") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C29
ggsave("figures/periodicity/chrv_AT_q1vq4_C29.png", plot = chrv_AT_q1vq4_C29)

chrv_AT_q1vq4_C31 = ggplot(data = chrv_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C31") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C31
ggsave("figures/periodicity/chrv_AT_q1vq4_C31.png", plot = chrv_AT_q1vq4_C31)

chrv_AT_q1vq4_C0 = ggplot(data = chrv_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0
ggsave("figures/periodicity/chrv_AT_q1vq4_C0.png", plot = chrv_AT_q1vq4_C0)

chrv_AT_q1vq4_C0_new_v1 = ggplot(data = chrv_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v1
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v1.png", plot = chrv_AT_q1vq4_C0_new_v1)

chrv_AT_q1vq4_C0_new_v2 = ggplot(data = chrv_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v2
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v2.png", plot = chrv_AT_q1vq4_C0_new_v2)

chrv_AT_q1vq4_C0_new_v3 = ggplot(data = chrv_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v3
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v3.png", plot = chrv_AT_q1vq4_C0_new_v3)

chrv_AT_q1vq4_C0_new_v4 = ggplot(data = chrv_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v4
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v4.png", plot = chrv_AT_q1vq4_C0_new_v4)

chrv_AT_q1vq4_C0_new_v5 = ggplot(data = chrv_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v5
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v5.png", plot = chrv_AT_q1vq4_C0_new_v5)

chrv_AT_q1vq4_C0_new_v6 = ggplot(data = chrv_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v6
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v6.png", plot = chrv_AT_q1vq4_C0_new_v6)

chrv_AT_q1vq4_C0_new_v7 = ggplot(data = chrv_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v7
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v7.png", plot = chrv_AT_q1vq4_C0_new_v7)

chrv_AT_q1vq4_C0_new_v8 = ggplot(data = chrv_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v8
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v8.png", plot = chrv_AT_q1vq4_C0_new_v8)

chrv_AT_q1vq4_C0_new_v9 = ggplot(data = chrv_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v9
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v9.png", plot = chrv_AT_q1vq4_C0_new_v9)

chrv_AT_q1vq4_C0_new_v10 = ggplot(data = chrv_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v10
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v10.png", plot = chrv_AT_q1vq4_C0_new_v10)

chrv_AT_q1vq4_C0_new_v11 = ggplot(data = chrv_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v11
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v11.png", plot = chrv_AT_q1vq4_C0_new_v11)

chrv_AT_q1vq4_C0_new_v12 = ggplot(data = chrv_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v12
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v12.png", plot = chrv_AT_q1vq4_C0_new_v12)

chrv_AT_q1vq4_C0_new_v13 = ggplot(data = chrv_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q1vq4_C0_new_v13
ggsave("figures/periodicity/chrv_AT_q1vq4_C0_new_v13.png", plot = chrv_AT_q1vq4_C0_new_v13)


chrv_AT_q2vq3_C26 = ggplot(data = chrv_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C26") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C26
ggsave("figures/periodicity/chrv_AT_q2vq3_C26.png", plot = chrv_AT_q2vq3_C26)

chrv_AT_q2vq3_C29 = ggplot(data = chrv_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C29") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C29
ggsave("figures/periodicity/chrv_AT_q2vq3_C29.png", plot = chrv_AT_q2vq3_C29)

chrv_AT_q2vq3_C31 = ggplot(data = chrv_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C31") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C31
ggsave("figures/periodicity/chrv_AT_q2vq3_C31.png", plot = chrv_AT_q2vq3_C31)

chrv_AT_q2vq3_C0 = ggplot(data = chrv_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0
ggsave("figures/periodicity/chrv_AT_q2vq3_C0.png", plot = chrv_AT_q2vq3_C0)

chrv_AT_q2vq3_C0_new_v1 = ggplot(data = chrv_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v1
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v1.png", plot = chrv_AT_q2vq3_C0_new_v1)

chrv_AT_q2vq3_C0_new_v2 = ggplot(data = chrv_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v2
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v2.png", plot = chrv_AT_q2vq3_C0_new_v2)

chrv_AT_q2vq3_C0_new_v3 = ggplot(data = chrv_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v3
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v3.png", plot = chrv_AT_q2vq3_C0_new_v3)

chrv_AT_q2vq3_C0_new_v4 = ggplot(data = chrv_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v4
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v4.png", plot = chrv_AT_q2vq3_C0_new_v4)

chrv_AT_q2vq3_C0_new_v5 = ggplot(data = chrv_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v5
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v5.png", plot = chrv_AT_q2vq3_C0_new_v5)

chrv_AT_q2vq3_C0_new_v6 = ggplot(data = chrv_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v6
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v6.png", plot = chrv_AT_q2vq3_C0_new_v6)

chrv_AT_q2vq3_C0_new_v7 = ggplot(data = chrv_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v7
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v7.png", plot = chrv_AT_q2vq3_C0_new_v7)

chrv_AT_q2vq3_C0_new_v8 = ggplot(data = chrv_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v8
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v8.png", plot = chrv_AT_q2vq3_C0_new_v8)

chrv_AT_q2vq3_C0_new_v9 = ggplot(data = chrv_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v9
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v9.png", plot = chrv_AT_q2vq3_C0_new_v9)

chrv_AT_q2vq3_C0_new_v10 = ggplot(data = chrv_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v10
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v10.png", plot = chrv_AT_q2vq3_C0_new_v10)

chrv_AT_q2vq3_C0_new_v11 = ggplot(data = chrv_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v11
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v11.png", plot = chrv_AT_q2vq3_C0_new_v11)

chrv_AT_q2vq3_C0_new_v12 = ggplot(data = chrv_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v12
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v12.png", plot = chrv_AT_q2vq3_C0_new_v12)

chrv_AT_q2vq3_C0_new_v13 = ggplot(data = chrv_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_q2vq3_C0_new_v13
ggsave("figures/periodicity/chrv_AT_q2vq3_C0_new_v13.png", plot = chrv_AT_q2vq3_C0_new_v13)


chrv_AT_custom1_C26 = ggplot(data = chrv_A_T_C26_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C26") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C26
ggsave("figures/periodicity/chrv_AT_custom1_C26.png", plot = chrv_AT_custom1_C26)

chrv_AT_custom1_C29 = ggplot(data = chrv_A_T_C29_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C29") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C29
ggsave("figures/periodicity/chrv_AT_custom1_C29.png", plot = chrv_AT_custom1_C29)

chrv_AT_custom1_C31 = ggplot(data = chrv_A_T_C31_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C31") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C31
ggsave("figures/periodicity/chrv_AT_custom1_C31.png", plot = chrv_AT_custom1_C31)

chrv_AT_custom1_C0 = ggplot(data = chrv_A_T_C0_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0
ggsave("figures/periodicity/chrv_AT_custom1_C0.png", plot = chrv_AT_custom1_C0)

chrv_AT_custom1_C0_new_v1 = ggplot(data = chrv_A_T_C0_new_v1_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v1
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v1.png", plot = chrv_AT_custom1_C0_new_v1)

chrv_AT_custom1_C0_new_v2 = ggplot(data = chrv_A_T_C0_new_v2_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v2
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v2.png", plot = chrv_AT_custom1_C0_new_v2)

chrv_AT_custom1_C0_new_v3 = ggplot(data = chrv_A_T_C0_new_v3_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v3
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v3.png", plot = chrv_AT_custom1_C0_new_v3)

chrv_AT_custom1_C0_new_v4 = ggplot(data = chrv_A_T_C0_new_v4_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v4
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v4.png", plot = chrv_AT_custom1_C0_new_v4)

chrv_AT_custom1_C0_new_v5 = ggplot(data = chrv_A_T_C0_new_v5_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v5
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v5.png", plot = chrv_AT_custom1_C0_new_v5)

chrv_AT_custom1_C0_new_v6 = ggplot(data = chrv_A_T_C0_new_v6_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v6
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v6.png", plot = chrv_AT_custom1_C0_new_v6)

chrv_AT_custom1_C0_new_v7 = ggplot(data = chrv_A_T_C0_new_v7_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v7
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v7.png", plot = chrv_AT_custom1_C0_new_v7)

chrv_AT_custom1_C0_new_v8 = ggplot(data = chrv_A_T_C0_new_v8_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v8
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v8.png", plot = chrv_AT_custom1_C0_new_v8)

chrv_AT_custom1_C0_new_v9 = ggplot(data = chrv_A_T_C0_new_v9_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v9
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v9.png", plot = chrv_AT_custom1_C0_new_v9)

chrv_AT_custom1_C0_new_v10 = ggplot(data = chrv_A_T_C0_new_v10_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v10
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v10.png", plot = chrv_AT_custom1_C0_new_v10)

chrv_AT_custom1_C0_new_v11 = ggplot(data = chrv_A_T_C0_new_v11_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v11
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v11.png", plot = chrv_AT_custom1_C0_new_v11)

chrv_AT_custom1_C0_new_v12 = ggplot(data = chrv_A_T_C0_new_v12_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v12
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v12.png", plot = chrv_AT_custom1_C0_new_v12)

chrv_AT_custom1_C0_new_v13 = ggplot(data = chrv_A_T_C0_new_v13_q, aes(x = 1:50)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AT_all_C_min, chrv_AT_all_C_max)
chrv_AT_custom1_C0_new_v13
ggsave("figures/periodicity/chrv_AT_custom1_C0_new_v13.png", plot = chrv_AT_custom1_C0_new_v13)









# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile+custom/library for C26
nuc_AA_TT_AT_TA_C26_q1 = apply(nuc_C26_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C26_q1))
})
nuc_AA_TT_AT_TA_C26_q2 = apply(nuc_C26_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C26_q2))
})
nuc_AA_TT_AT_TA_C26_q3 = apply(nuc_C26_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C26_q3))
})
nuc_AA_TT_AT_TA_C26_q4 = apply(nuc_C26_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C26_q4))
})
nuc_AA_TT_AT_TA_C26_custom1 = apply(nuc_C26_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C26_custom1))
})

random_AA_TT_AT_TA_C26_q1 = apply(random_C26_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C26_q1))
})
random_AA_TT_AT_TA_C26_q2 = apply(random_C26_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C26_q2))
})
random_AA_TT_AT_TA_C26_q3 = apply(random_C26_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C26_q3))
})
random_AA_TT_AT_TA_C26_q4 = apply(random_C26_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C26_q4))
})
random_AA_TT_AT_TA_C26_custom1 = apply(random_C26_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C26_custom1))
})

tiling_AA_TT_AT_TA_C26_q1 = apply(tiling_C26_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C26_q1))
})
tiling_AA_TT_AT_TA_C26_q2 = apply(tiling_C26_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C26_q2))
})
tiling_AA_TT_AT_TA_C26_q3 = apply(tiling_C26_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C26_q3))
})
tiling_AA_TT_AT_TA_C26_q4 = apply(tiling_C26_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C26_q4))
})
tiling_AA_TT_AT_TA_C26_custom1 = apply(tiling_C26_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C26_custom1))
})

chrv_AA_TT_AT_TA_C26_q1 = apply(chrv_C26_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C26_q1))
})
chrv_AA_TT_AT_TA_C26_q2 = apply(chrv_C26_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C26_q2))
})
chrv_AA_TT_AT_TA_C26_q3 = apply(chrv_C26_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C26_q3))
})
chrv_AA_TT_AT_TA_C26_q4 = apply(chrv_C26_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C26_q4))
})
chrv_AA_TT_AT_TA_C26_custom1 = apply(chrv_C26_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C26_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C26_q = data.frame(q1 = nuc_AA_TT_AT_TA_C26_q1,
                                   q2 = nuc_AA_TT_AT_TA_C26_q2,
                                   q3 = nuc_AA_TT_AT_TA_C26_q3,
                                   q4 = nuc_AA_TT_AT_TA_C26_q4,
                                   custom1 = nuc_AA_TT_AT_TA_C26_custom1)

random_AA_TT_AT_TA_C26_q = data.frame(q1 = random_AA_TT_AT_TA_C26_q1,
                                      q2 = random_AA_TT_AT_TA_C26_q2,
                                      q3 = random_AA_TT_AT_TA_C26_q3,
                                      q4 = random_AA_TT_AT_TA_C26_q4,
                                      custom1 = random_AA_TT_AT_TA_C26_custom1)

tiling_AA_TT_AT_TA_C26_q = data.frame(q1 = tiling_AA_TT_AT_TA_C26_q1,
                                      q2 = tiling_AA_TT_AT_TA_C26_q2,
                                      q3 = tiling_AA_TT_AT_TA_C26_q3,
                                      q4 = tiling_AA_TT_AT_TA_C26_q4,
                                      custom1 = tiling_AA_TT_AT_TA_C26_custom1)

chrv_AA_TT_AT_TA_C26_q = data.frame(q1 = chrv_AA_TT_AT_TA_C26_q1,
                                    q2 = chrv_AA_TT_AT_TA_C26_q2,
                                    q3 = chrv_AA_TT_AT_TA_C26_q3,
                                    q4 = chrv_AA_TT_AT_TA_C26_q4,
                                    custom1 = chrv_AA_TT_AT_TA_C26_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C29
nuc_AA_TT_AT_TA_C29_q1 = apply(nuc_C29_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C29_q1))
})
nuc_AA_TT_AT_TA_C29_q2 = apply(nuc_C29_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C29_q2))
})
nuc_AA_TT_AT_TA_C29_q3 = apply(nuc_C29_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C29_q3))
})
nuc_AA_TT_AT_TA_C29_q4 = apply(nuc_C29_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C29_q4))
})
nuc_AA_TT_AT_TA_C29_custom1 = apply(nuc_C29_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C29_custom1))
})

random_AA_TT_AT_TA_C29_q1 = apply(random_C29_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C29_q1))
})
random_AA_TT_AT_TA_C29_q2 = apply(random_C29_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C29_q2))
})
random_AA_TT_AT_TA_C29_q3 = apply(random_C29_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C29_q3))
})
random_AA_TT_AT_TA_C29_q4 = apply(random_C29_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C29_q4))
})
random_AA_TT_AT_TA_C29_custom1 = apply(random_C29_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C29_custom1))
})

tiling_AA_TT_AT_TA_C29_q1 = apply(tiling_C29_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C29_q1))
})
tiling_AA_TT_AT_TA_C29_q2 = apply(tiling_C29_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C29_q2))
})
tiling_AA_TT_AT_TA_C29_q3 = apply(tiling_C29_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C29_q3))
})
tiling_AA_TT_AT_TA_C29_q4 = apply(tiling_C29_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C29_q4))
})
tiling_AA_TT_AT_TA_C29_custom1 = apply(tiling_C29_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C29_custom1))
})

chrv_AA_TT_AT_TA_C29_q1 = apply(chrv_C29_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C29_q1))
})
chrv_AA_TT_AT_TA_C29_q2 = apply(chrv_C29_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C29_q2))
})
chrv_AA_TT_AT_TA_C29_q3 = apply(chrv_C29_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C29_q3))
})
chrv_AA_TT_AT_TA_C29_q4 = apply(chrv_C29_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C29_q4))
})
chrv_AA_TT_AT_TA_C29_custom1 = apply(chrv_C29_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C29_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C29_q = data.frame(q1 = nuc_AA_TT_AT_TA_C29_q1,
                                   q2 = nuc_AA_TT_AT_TA_C29_q2,
                                   q3 = nuc_AA_TT_AT_TA_C29_q3,
                                   q4 = nuc_AA_TT_AT_TA_C29_q4,
                                   custom1 = nuc_AA_TT_AT_TA_C29_custom1)

random_AA_TT_AT_TA_C29_q = data.frame(q1 = random_AA_TT_AT_TA_C29_q1,
                                      q2 = random_AA_TT_AT_TA_C29_q2,
                                      q3 = random_AA_TT_AT_TA_C29_q3,
                                      q4 = random_AA_TT_AT_TA_C29_q4,
                                      custom1 = random_AA_TT_AT_TA_C29_custom1)

tiling_AA_TT_AT_TA_C29_q = data.frame(q1 = tiling_AA_TT_AT_TA_C29_q1,
                                      q2 = tiling_AA_TT_AT_TA_C29_q2,
                                      q3 = tiling_AA_TT_AT_TA_C29_q3,
                                      q4 = tiling_AA_TT_AT_TA_C29_q4,
                                      custom1 = tiling_AA_TT_AT_TA_C29_custom1)

chrv_AA_TT_AT_TA_C29_q = data.frame(q1 = chrv_AA_TT_AT_TA_C29_q1,
                                    q2 = chrv_AA_TT_AT_TA_C29_q2,
                                    q3 = chrv_AA_TT_AT_TA_C29_q3,
                                    q4 = chrv_AA_TT_AT_TA_C29_q4,
                                    custom1 = chrv_AA_TT_AT_TA_C29_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C31
nuc_AA_TT_AT_TA_C31_q1 = apply(nuc_C31_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C31_q1))
})
nuc_AA_TT_AT_TA_C31_q2 = apply(nuc_C31_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C31_q2))
})
nuc_AA_TT_AT_TA_C31_q3 = apply(nuc_C31_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C31_q3))
})
nuc_AA_TT_AT_TA_C31_q4 = apply(nuc_C31_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C31_q4))
})
nuc_AA_TT_AT_TA_C31_custom1 = apply(nuc_C31_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C31_custom1))
})

random_AA_TT_AT_TA_C31_q1 = apply(random_C31_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C31_q1))
})
random_AA_TT_AT_TA_C31_q2 = apply(random_C31_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C31_q2))
})
random_AA_TT_AT_TA_C31_q3 = apply(random_C31_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C31_q3))
})
random_AA_TT_AT_TA_C31_q4 = apply(random_C31_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C31_q4))
})
random_AA_TT_AT_TA_C31_custom1 = apply(random_C31_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C31_custom1))
})

tiling_AA_TT_AT_TA_C31_q1 = apply(tiling_C31_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C31_q1))
})
tiling_AA_TT_AT_TA_C31_q2 = apply(tiling_C31_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C31_q2))
})
tiling_AA_TT_AT_TA_C31_q3 = apply(tiling_C31_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C31_q3))
})
tiling_AA_TT_AT_TA_C31_q4 = apply(tiling_C31_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C31_q4))
})
tiling_AA_TT_AT_TA_C31_custom1 = apply(tiling_C31_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C31_custom1))
})

chrv_AA_TT_AT_TA_C31_q1 = apply(chrv_C31_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C31_q1))
})
chrv_AA_TT_AT_TA_C31_q2 = apply(chrv_C31_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C31_q2))
})
chrv_AA_TT_AT_TA_C31_q3 = apply(chrv_C31_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C31_q3))
})
chrv_AA_TT_AT_TA_C31_q4 = apply(chrv_C31_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C31_q4))
})
chrv_AA_TT_AT_TA_C31_custom1 = apply(chrv_C31_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C31_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C31_q = data.frame(q1 = nuc_AA_TT_AT_TA_C31_q1,
                                   q2 = nuc_AA_TT_AT_TA_C31_q2,
                                   q3 = nuc_AA_TT_AT_TA_C31_q3,
                                   q4 = nuc_AA_TT_AT_TA_C31_q4,
                                   custom1 = nuc_AA_TT_AT_TA_C31_custom1)

random_AA_TT_AT_TA_C31_q = data.frame(q1 = random_AA_TT_AT_TA_C31_q1,
                                      q2 = random_AA_TT_AT_TA_C31_q2,
                                      q3 = random_AA_TT_AT_TA_C31_q3,
                                      q4 = random_AA_TT_AT_TA_C31_q4,
                                      custom1 = random_AA_TT_AT_TA_C31_custom1)

tiling_AA_TT_AT_TA_C31_q = data.frame(q1 = tiling_AA_TT_AT_TA_C31_q1,
                                      q2 = tiling_AA_TT_AT_TA_C31_q2,
                                      q3 = tiling_AA_TT_AT_TA_C31_q3,
                                      q4 = tiling_AA_TT_AT_TA_C31_q4,
                                      custom1 = tiling_AA_TT_AT_TA_C31_custom1)

chrv_AA_TT_AT_TA_C31_q = data.frame(q1 = chrv_AA_TT_AT_TA_C31_q1,
                                    q2 = chrv_AA_TT_AT_TA_C31_q2,
                                    q3 = chrv_AA_TT_AT_TA_C31_q3,
                                    q4 = chrv_AA_TT_AT_TA_C31_q4,
                                    custom1 = chrv_AA_TT_AT_TA_C31_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0
nuc_AA_TT_AT_TA_C0_q1 = apply(nuc_C0_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_q1))
})
nuc_AA_TT_AT_TA_C0_q2 = apply(nuc_C0_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_q2))
})
nuc_AA_TT_AT_TA_C0_q3 = apply(nuc_C0_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_q3))
})
nuc_AA_TT_AT_TA_C0_q4 = apply(nuc_C0_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_q4))
})
nuc_AA_TT_AT_TA_C0_custom1 = apply(nuc_C0_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_custom1))
})

random_AA_TT_AT_TA_C0_q1 = apply(random_C0_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_q1))
})
random_AA_TT_AT_TA_C0_q2 = apply(random_C0_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_q2))
})
random_AA_TT_AT_TA_C0_q3 = apply(random_C0_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_q3))
})
random_AA_TT_AT_TA_C0_q4 = apply(random_C0_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_q4))
})
random_AA_TT_AT_TA_C0_custom1 = apply(random_C0_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_custom1))
})

tiling_AA_TT_AT_TA_C0_q1 = apply(tiling_C0_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_q1))
})
tiling_AA_TT_AT_TA_C0_q2 = apply(tiling_C0_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_q2))
})
tiling_AA_TT_AT_TA_C0_q3 = apply(tiling_C0_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_q3))
})
tiling_AA_TT_AT_TA_C0_q4 = apply(tiling_C0_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_q4))
})
tiling_AA_TT_AT_TA_C0_custom1 = apply(tiling_C0_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_custom1))
})

chrv_AA_TT_AT_TA_C0_q1 = apply(chrv_C0_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_q1))
})
chrv_AA_TT_AT_TA_C0_q2 = apply(chrv_C0_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_q2))
})
chrv_AA_TT_AT_TA_C0_q3 = apply(chrv_C0_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_q3))
})
chrv_AA_TT_AT_TA_C0_q4 = apply(chrv_C0_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_q4))
})
chrv_AA_TT_AT_TA_C0_custom1 = apply(chrv_C0_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_q1,
                                  q2 = nuc_AA_TT_AT_TA_C0_q2,
                                  q3 = nuc_AA_TT_AT_TA_C0_q3,
                                  q4 = nuc_AA_TT_AT_TA_C0_q4,
                                  custom1 = nuc_AA_TT_AT_TA_C0_custom1)

random_AA_TT_AT_TA_C0_q = data.frame(q1 = random_AA_TT_AT_TA_C0_q1,
                                     q2 = random_AA_TT_AT_TA_C0_q2,
                                     q3 = random_AA_TT_AT_TA_C0_q3,
                                     q4 = random_AA_TT_AT_TA_C0_q4,
                                     custom1 = random_AA_TT_AT_TA_C0_custom1)

tiling_AA_TT_AT_TA_C0_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_q1,
                                     q2 = tiling_AA_TT_AT_TA_C0_q2,
                                     q3 = tiling_AA_TT_AT_TA_C0_q3,
                                     q4 = tiling_AA_TT_AT_TA_C0_q4,
                                     custom1 = tiling_AA_TT_AT_TA_C0_custom1)

chrv_AA_TT_AT_TA_C0_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_q1,
                                   q2 = chrv_AA_TT_AT_TA_C0_q2,
                                   q3 = chrv_AA_TT_AT_TA_C0_q3,
                                   q4 = chrv_AA_TT_AT_TA_C0_q4,
                                   custom1 = chrv_AA_TT_AT_TA_C0_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v1
nuc_AA_TT_AT_TA_C0_new_v1_q1 = apply(nuc_C0_new_v1_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v1_q1))
})
nuc_AA_TT_AT_TA_C0_new_v1_q2 = apply(nuc_C0_new_v1_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v1_q2))
})
nuc_AA_TT_AT_TA_C0_new_v1_q3 = apply(nuc_C0_new_v1_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v1_q3))
})
nuc_AA_TT_AT_TA_C0_new_v1_q4 = apply(nuc_C0_new_v1_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v1_q4))
})
nuc_AA_TT_AT_TA_C0_new_v1_custom1 = apply(nuc_C0_new_v1_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v1_custom1))
})

random_AA_TT_AT_TA_C0_new_v1_q1 = apply(random_C0_new_v1_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v1_q1))
})
random_AA_TT_AT_TA_C0_new_v1_q2 = apply(random_C0_new_v1_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v1_q2))
})
random_AA_TT_AT_TA_C0_new_v1_q3 = apply(random_C0_new_v1_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v1_q3))
})
random_AA_TT_AT_TA_C0_new_v1_q4 = apply(random_C0_new_v1_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v1_q4))
})
random_AA_TT_AT_TA_C0_new_v1_custom1 = apply(random_C0_new_v1_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v1_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v1_q1 = apply(tiling_C0_new_v1_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v1_q1))
})
tiling_AA_TT_AT_TA_C0_new_v1_q2 = apply(tiling_C0_new_v1_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v1_q2))
})
tiling_AA_TT_AT_TA_C0_new_v1_q3 = apply(tiling_C0_new_v1_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v1_q3))
})
tiling_AA_TT_AT_TA_C0_new_v1_q4 = apply(tiling_C0_new_v1_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v1_q4))
})
tiling_AA_TT_AT_TA_C0_new_v1_custom1 = apply(tiling_C0_new_v1_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v1_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v1_q1 = apply(chrv_C0_new_v1_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v1_q1))
})
chrv_AA_TT_AT_TA_C0_new_v1_q2 = apply(chrv_C0_new_v1_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v1_q2))
})
chrv_AA_TT_AT_TA_C0_new_v1_q3 = apply(chrv_C0_new_v1_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v1_q3))
})
chrv_AA_TT_AT_TA_C0_new_v1_q4 = apply(chrv_C0_new_v1_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v1_q4))
})
chrv_AA_TT_AT_TA_C0_new_v1_custom1 = apply(chrv_C0_new_v1_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v1_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v1_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v1_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v1_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v1_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v1_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v1_custom1)

random_AA_TT_AT_TA_C0_new_v1_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v1_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v1_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v1_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v1_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v1_custom1)

tiling_AA_TT_AT_TA_C0_new_v1_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v1_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v1_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v1_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v1_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v1_custom1)

chrv_AA_TT_AT_TA_C0_new_v1_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v1_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v1_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v1_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v1_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v1_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v2
nuc_AA_TT_AT_TA_C0_new_v2_q1 = apply(nuc_C0_new_v2_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v2_q1))
})
nuc_AA_TT_AT_TA_C0_new_v2_q2 = apply(nuc_C0_new_v2_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v2_q2))
})
nuc_AA_TT_AT_TA_C0_new_v2_q3 = apply(nuc_C0_new_v2_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v2_q3))
})
nuc_AA_TT_AT_TA_C0_new_v2_q4 = apply(nuc_C0_new_v2_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v2_q4))
})
nuc_AA_TT_AT_TA_C0_new_v2_custom1 = apply(nuc_C0_new_v2_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v2_custom1))
})

random_AA_TT_AT_TA_C0_new_v2_q1 = apply(random_C0_new_v2_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v2_q1))
})
random_AA_TT_AT_TA_C0_new_v2_q2 = apply(random_C0_new_v2_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v2_q2))
})
random_AA_TT_AT_TA_C0_new_v2_q3 = apply(random_C0_new_v2_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v2_q3))
})
random_AA_TT_AT_TA_C0_new_v2_q4 = apply(random_C0_new_v2_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v2_q4))
})
random_AA_TT_AT_TA_C0_new_v2_custom1 = apply(random_C0_new_v2_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v2_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v2_q1 = apply(tiling_C0_new_v2_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v2_q1))
})
tiling_AA_TT_AT_TA_C0_new_v2_q2 = apply(tiling_C0_new_v2_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v2_q2))
})
tiling_AA_TT_AT_TA_C0_new_v2_q3 = apply(tiling_C0_new_v2_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v2_q3))
})
tiling_AA_TT_AT_TA_C0_new_v2_q4 = apply(tiling_C0_new_v2_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v2_q4))
})
tiling_AA_TT_AT_TA_C0_new_v2_custom1 = apply(tiling_C0_new_v2_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v2_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v2_q1 = apply(chrv_C0_new_v2_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v2_q1))
})
chrv_AA_TT_AT_TA_C0_new_v2_q2 = apply(chrv_C0_new_v2_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v2_q2))
})
chrv_AA_TT_AT_TA_C0_new_v2_q3 = apply(chrv_C0_new_v2_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v2_q3))
})
chrv_AA_TT_AT_TA_C0_new_v2_q4 = apply(chrv_C0_new_v2_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v2_q4))
})
chrv_AA_TT_AT_TA_C0_new_v2_custom1 = apply(chrv_C0_new_v2_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v2_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v2_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v2_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v2_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v2_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v2_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v2_custom1)

random_AA_TT_AT_TA_C0_new_v2_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v2_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v2_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v2_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v2_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v2_custom1)

tiling_AA_TT_AT_TA_C0_new_v2_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v2_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v2_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v2_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v2_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v2_custom1)

chrv_AA_TT_AT_TA_C0_new_v2_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v2_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v2_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v2_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v2_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v2_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v3
nuc_AA_TT_AT_TA_C0_new_v3_q1 = apply(nuc_C0_new_v3_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v3_q1))
})
nuc_AA_TT_AT_TA_C0_new_v3_q2 = apply(nuc_C0_new_v3_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v3_q2))
})
nuc_AA_TT_AT_TA_C0_new_v3_q3 = apply(nuc_C0_new_v3_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v3_q3))
})
nuc_AA_TT_AT_TA_C0_new_v3_q4 = apply(nuc_C0_new_v3_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v3_q4))
})
nuc_AA_TT_AT_TA_C0_new_v3_custom1 = apply(nuc_C0_new_v3_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v3_custom1))
})

random_AA_TT_AT_TA_C0_new_v3_q1 = apply(random_C0_new_v3_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v3_q1))
})
random_AA_TT_AT_TA_C0_new_v3_q2 = apply(random_C0_new_v3_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v3_q2))
})
random_AA_TT_AT_TA_C0_new_v3_q3 = apply(random_C0_new_v3_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v3_q3))
})
random_AA_TT_AT_TA_C0_new_v3_q4 = apply(random_C0_new_v3_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v3_q4))
})
random_AA_TT_AT_TA_C0_new_v3_custom1 = apply(random_C0_new_v3_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v3_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v3_q1 = apply(tiling_C0_new_v3_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v3_q1))
})
tiling_AA_TT_AT_TA_C0_new_v3_q2 = apply(tiling_C0_new_v3_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v3_q2))
})
tiling_AA_TT_AT_TA_C0_new_v3_q3 = apply(tiling_C0_new_v3_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v3_q3))
})
tiling_AA_TT_AT_TA_C0_new_v3_q4 = apply(tiling_C0_new_v3_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v3_q4))
})
tiling_AA_TT_AT_TA_C0_new_v3_custom1 = apply(tiling_C0_new_v3_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v3_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v3_q1 = apply(chrv_C0_new_v3_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v3_q1))
})
chrv_AA_TT_AT_TA_C0_new_v3_q2 = apply(chrv_C0_new_v3_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v3_q2))
})
chrv_AA_TT_AT_TA_C0_new_v3_q3 = apply(chrv_C0_new_v3_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v3_q3))
})
chrv_AA_TT_AT_TA_C0_new_v3_q4 = apply(chrv_C0_new_v3_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v3_q4))
})
chrv_AA_TT_AT_TA_C0_new_v3_custom1 = apply(chrv_C0_new_v3_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v3_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v3_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v3_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v3_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v3_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v3_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v3_custom1)

random_AA_TT_AT_TA_C0_new_v3_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v3_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v3_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v3_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v3_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v3_custom1)

tiling_AA_TT_AT_TA_C0_new_v3_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v3_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v3_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v3_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v3_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v3_custom1)

chrv_AA_TT_AT_TA_C0_new_v3_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v3_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v3_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v3_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v3_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v3_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v4
nuc_AA_TT_AT_TA_C0_new_v4_q1 = apply(nuc_C0_new_v4_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v4_q1))
})
nuc_AA_TT_AT_TA_C0_new_v4_q2 = apply(nuc_C0_new_v4_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v4_q2))
})
nuc_AA_TT_AT_TA_C0_new_v4_q3 = apply(nuc_C0_new_v4_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v4_q3))
})
nuc_AA_TT_AT_TA_C0_new_v4_q4 = apply(nuc_C0_new_v4_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v4_q4))
})
nuc_AA_TT_AT_TA_C0_new_v4_custom1 = apply(nuc_C0_new_v4_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v4_custom1))
})

random_AA_TT_AT_TA_C0_new_v4_q1 = apply(random_C0_new_v4_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v4_q1))
})
random_AA_TT_AT_TA_C0_new_v4_q2 = apply(random_C0_new_v4_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v4_q2))
})
random_AA_TT_AT_TA_C0_new_v4_q3 = apply(random_C0_new_v4_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v4_q3))
})
random_AA_TT_AT_TA_C0_new_v4_q4 = apply(random_C0_new_v4_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v4_q4))
})
random_AA_TT_AT_TA_C0_new_v4_custom1 = apply(random_C0_new_v4_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v4_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v4_q1 = apply(tiling_C0_new_v4_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v4_q1))
})
tiling_AA_TT_AT_TA_C0_new_v4_q2 = apply(tiling_C0_new_v4_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v4_q2))
})
tiling_AA_TT_AT_TA_C0_new_v4_q3 = apply(tiling_C0_new_v4_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v4_q3))
})
tiling_AA_TT_AT_TA_C0_new_v4_q4 = apply(tiling_C0_new_v4_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v4_q4))
})
tiling_AA_TT_AT_TA_C0_new_v4_custom1 = apply(tiling_C0_new_v4_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v4_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v4_q1 = apply(chrv_C0_new_v4_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v4_q1))
})
chrv_AA_TT_AT_TA_C0_new_v4_q2 = apply(chrv_C0_new_v4_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v4_q2))
})
chrv_AA_TT_AT_TA_C0_new_v4_q3 = apply(chrv_C0_new_v4_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v4_q3))
})
chrv_AA_TT_AT_TA_C0_new_v4_q4 = apply(chrv_C0_new_v4_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v4_q4))
})
chrv_AA_TT_AT_TA_C0_new_v4_custom1 = apply(chrv_C0_new_v4_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v4_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v4_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v4_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v4_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v4_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v4_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v4_custom1)

random_AA_TT_AT_TA_C0_new_v4_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v4_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v4_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v4_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v4_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v4_custom1)

tiling_AA_TT_AT_TA_C0_new_v4_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v4_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v4_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v4_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v4_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v4_custom1)

chrv_AA_TT_AT_TA_C0_new_v4_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v4_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v4_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v4_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v4_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v4_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v5
nuc_AA_TT_AT_TA_C0_new_v5_q1 = apply(nuc_C0_new_v5_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v5_q1))
})
nuc_AA_TT_AT_TA_C0_new_v5_q2 = apply(nuc_C0_new_v5_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v5_q2))
})
nuc_AA_TT_AT_TA_C0_new_v5_q3 = apply(nuc_C0_new_v5_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v5_q3))
})
nuc_AA_TT_AT_TA_C0_new_v5_q4 = apply(nuc_C0_new_v5_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v5_q4))
})
nuc_AA_TT_AT_TA_C0_new_v5_custom1 = apply(nuc_C0_new_v5_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v5_custom1))
})

random_AA_TT_AT_TA_C0_new_v5_q1 = apply(random_C0_new_v5_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v5_q1))
})
random_AA_TT_AT_TA_C0_new_v5_q2 = apply(random_C0_new_v5_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v5_q2))
})
random_AA_TT_AT_TA_C0_new_v5_q3 = apply(random_C0_new_v5_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v5_q3))
})
random_AA_TT_AT_TA_C0_new_v5_q4 = apply(random_C0_new_v5_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v5_q4))
})
random_AA_TT_AT_TA_C0_new_v5_custom1 = apply(random_C0_new_v5_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v5_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v5_q1 = apply(tiling_C0_new_v5_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v5_q1))
})
tiling_AA_TT_AT_TA_C0_new_v5_q2 = apply(tiling_C0_new_v5_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v5_q2))
})
tiling_AA_TT_AT_TA_C0_new_v5_q3 = apply(tiling_C0_new_v5_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v5_q3))
})
tiling_AA_TT_AT_TA_C0_new_v5_q4 = apply(tiling_C0_new_v5_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v5_q4))
})
tiling_AA_TT_AT_TA_C0_new_v5_custom1 = apply(tiling_C0_new_v5_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v5_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v5_q1 = apply(chrv_C0_new_v5_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v5_q1))
})
chrv_AA_TT_AT_TA_C0_new_v5_q2 = apply(chrv_C0_new_v5_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v5_q2))
})
chrv_AA_TT_AT_TA_C0_new_v5_q3 = apply(chrv_C0_new_v5_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v5_q3))
})
chrv_AA_TT_AT_TA_C0_new_v5_q4 = apply(chrv_C0_new_v5_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v5_q4))
})
chrv_AA_TT_AT_TA_C0_new_v5_custom1 = apply(chrv_C0_new_v5_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v5_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v5_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v5_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v5_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v5_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v5_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v5_custom1)

random_AA_TT_AT_TA_C0_new_v5_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v5_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v5_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v5_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v5_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v5_custom1)

tiling_AA_TT_AT_TA_C0_new_v5_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v5_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v5_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v5_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v5_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v5_custom1)

chrv_AA_TT_AT_TA_C0_new_v5_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v5_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v5_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v5_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v5_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v5_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v6
nuc_AA_TT_AT_TA_C0_new_v6_q1 = apply(nuc_C0_new_v6_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v6_q1))
})
nuc_AA_TT_AT_TA_C0_new_v6_q2 = apply(nuc_C0_new_v6_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v6_q2))
})
nuc_AA_TT_AT_TA_C0_new_v6_q3 = apply(nuc_C0_new_v6_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v6_q3))
})
nuc_AA_TT_AT_TA_C0_new_v6_q4 = apply(nuc_C0_new_v6_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v6_q4))
})
nuc_AA_TT_AT_TA_C0_new_v6_custom1 = apply(nuc_C0_new_v6_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v6_custom1))
})

random_AA_TT_AT_TA_C0_new_v6_q1 = apply(random_C0_new_v6_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v6_q1))
})
random_AA_TT_AT_TA_C0_new_v6_q2 = apply(random_C0_new_v6_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v6_q2))
})
random_AA_TT_AT_TA_C0_new_v6_q3 = apply(random_C0_new_v6_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v6_q3))
})
random_AA_TT_AT_TA_C0_new_v6_q4 = apply(random_C0_new_v6_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v6_q4))
})
random_AA_TT_AT_TA_C0_new_v6_custom1 = apply(random_C0_new_v6_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v6_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v6_q1 = apply(tiling_C0_new_v6_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v6_q1))
})
tiling_AA_TT_AT_TA_C0_new_v6_q2 = apply(tiling_C0_new_v6_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v6_q2))
})
tiling_AA_TT_AT_TA_C0_new_v6_q3 = apply(tiling_C0_new_v6_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v6_q3))
})
tiling_AA_TT_AT_TA_C0_new_v6_q4 = apply(tiling_C0_new_v6_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v6_q4))
})
tiling_AA_TT_AT_TA_C0_new_v6_custom1 = apply(tiling_C0_new_v6_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v6_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v6_q1 = apply(chrv_C0_new_v6_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v6_q1))
})
chrv_AA_TT_AT_TA_C0_new_v6_q2 = apply(chrv_C0_new_v6_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v6_q2))
})
chrv_AA_TT_AT_TA_C0_new_v6_q3 = apply(chrv_C0_new_v6_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v6_q3))
})
chrv_AA_TT_AT_TA_C0_new_v6_q4 = apply(chrv_C0_new_v6_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v6_q4))
})
chrv_AA_TT_AT_TA_C0_new_v6_custom1 = apply(chrv_C0_new_v6_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v6_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v6_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v6_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v6_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v6_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v6_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v6_custom1)

random_AA_TT_AT_TA_C0_new_v6_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v6_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v6_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v6_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v6_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v6_custom1)

tiling_AA_TT_AT_TA_C0_new_v6_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v6_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v6_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v6_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v6_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v6_custom1)

chrv_AA_TT_AT_TA_C0_new_v6_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v6_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v6_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v6_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v6_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v6_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v7
nuc_AA_TT_AT_TA_C0_new_v7_q1 = apply(nuc_C0_new_v7_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v7_q1))
})
nuc_AA_TT_AT_TA_C0_new_v7_q2 = apply(nuc_C0_new_v7_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v7_q2))
})
nuc_AA_TT_AT_TA_C0_new_v7_q3 = apply(nuc_C0_new_v7_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v7_q3))
})
nuc_AA_TT_AT_TA_C0_new_v7_q4 = apply(nuc_C0_new_v7_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v7_q4))
})
nuc_AA_TT_AT_TA_C0_new_v7_custom1 = apply(nuc_C0_new_v7_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v7_custom1))
})

random_AA_TT_AT_TA_C0_new_v7_q1 = apply(random_C0_new_v7_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v7_q1))
})
random_AA_TT_AT_TA_C0_new_v7_q2 = apply(random_C0_new_v7_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v7_q2))
})
random_AA_TT_AT_TA_C0_new_v7_q3 = apply(random_C0_new_v7_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v7_q3))
})
random_AA_TT_AT_TA_C0_new_v7_q4 = apply(random_C0_new_v7_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v7_q4))
})
random_AA_TT_AT_TA_C0_new_v7_custom1 = apply(random_C0_new_v7_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v7_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v7_q1 = apply(tiling_C0_new_v7_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v7_q1))
})
tiling_AA_TT_AT_TA_C0_new_v7_q2 = apply(tiling_C0_new_v7_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v7_q2))
})
tiling_AA_TT_AT_TA_C0_new_v7_q3 = apply(tiling_C0_new_v7_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v7_q3))
})
tiling_AA_TT_AT_TA_C0_new_v7_q4 = apply(tiling_C0_new_v7_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v7_q4))
})
tiling_AA_TT_AT_TA_C0_new_v7_custom1 = apply(tiling_C0_new_v7_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v7_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v7_q1 = apply(chrv_C0_new_v7_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v7_q1))
})
chrv_AA_TT_AT_TA_C0_new_v7_q2 = apply(chrv_C0_new_v7_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v7_q2))
})
chrv_AA_TT_AT_TA_C0_new_v7_q3 = apply(chrv_C0_new_v7_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v7_q3))
})
chrv_AA_TT_AT_TA_C0_new_v7_q4 = apply(chrv_C0_new_v7_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v7_q4))
})
chrv_AA_TT_AT_TA_C0_new_v7_custom1 = apply(chrv_C0_new_v7_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v7_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v7_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v7_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v7_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v7_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v7_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v7_custom1)

random_AA_TT_AT_TA_C0_new_v7_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v7_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v7_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v7_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v7_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v7_custom1)

tiling_AA_TT_AT_TA_C0_new_v7_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v7_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v7_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v7_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v7_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v7_custom1)

chrv_AA_TT_AT_TA_C0_new_v7_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v7_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v7_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v7_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v7_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v7_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v8
nuc_AA_TT_AT_TA_C0_new_v8_q1 = apply(nuc_C0_new_v8_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v8_q1))
})
nuc_AA_TT_AT_TA_C0_new_v8_q2 = apply(nuc_C0_new_v8_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v8_q2))
})
nuc_AA_TT_AT_TA_C0_new_v8_q3 = apply(nuc_C0_new_v8_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v8_q3))
})
nuc_AA_TT_AT_TA_C0_new_v8_q4 = apply(nuc_C0_new_v8_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v8_q4))
})
nuc_AA_TT_AT_TA_C0_new_v8_custom1 = apply(nuc_C0_new_v8_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v8_custom1))
})

random_AA_TT_AT_TA_C0_new_v8_q1 = apply(random_C0_new_v8_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v8_q1))
})
random_AA_TT_AT_TA_C0_new_v8_q2 = apply(random_C0_new_v8_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v8_q2))
})
random_AA_TT_AT_TA_C0_new_v8_q3 = apply(random_C0_new_v8_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v8_q3))
})
random_AA_TT_AT_TA_C0_new_v8_q4 = apply(random_C0_new_v8_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v8_q4))
})
random_AA_TT_AT_TA_C0_new_v8_custom1 = apply(random_C0_new_v8_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v8_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v8_q1 = apply(tiling_C0_new_v8_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v8_q1))
})
tiling_AA_TT_AT_TA_C0_new_v8_q2 = apply(tiling_C0_new_v8_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v8_q2))
})
tiling_AA_TT_AT_TA_C0_new_v8_q3 = apply(tiling_C0_new_v8_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v8_q3))
})
tiling_AA_TT_AT_TA_C0_new_v8_q4 = apply(tiling_C0_new_v8_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v8_q4))
})
tiling_AA_TT_AT_TA_C0_new_v8_custom1 = apply(tiling_C0_new_v8_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v8_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v8_q1 = apply(chrv_C0_new_v8_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v8_q1))
})
chrv_AA_TT_AT_TA_C0_new_v8_q2 = apply(chrv_C0_new_v8_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v8_q2))
})
chrv_AA_TT_AT_TA_C0_new_v8_q3 = apply(chrv_C0_new_v8_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v8_q3))
})
chrv_AA_TT_AT_TA_C0_new_v8_q4 = apply(chrv_C0_new_v8_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v8_q4))
})
chrv_AA_TT_AT_TA_C0_new_v8_custom1 = apply(chrv_C0_new_v8_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v8_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v8_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v8_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v8_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v8_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v8_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v8_custom1)

random_AA_TT_AT_TA_C0_new_v8_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v8_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v8_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v8_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v8_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v8_custom1)

tiling_AA_TT_AT_TA_C0_new_v8_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v8_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v8_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v8_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v8_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v8_custom1)

chrv_AA_TT_AT_TA_C0_new_v8_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v8_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v8_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v8_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v8_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v8_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v9
nuc_AA_TT_AT_TA_C0_new_v9_q1 = apply(nuc_C0_new_v9_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v9_q1))
})
nuc_AA_TT_AT_TA_C0_new_v9_q2 = apply(nuc_C0_new_v9_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v9_q2))
})
nuc_AA_TT_AT_TA_C0_new_v9_q3 = apply(nuc_C0_new_v9_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v9_q3))
})
nuc_AA_TT_AT_TA_C0_new_v9_q4 = apply(nuc_C0_new_v9_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v9_q4))
})
nuc_AA_TT_AT_TA_C0_new_v9_custom1 = apply(nuc_C0_new_v9_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v9_custom1))
})

random_AA_TT_AT_TA_C0_new_v9_q1 = apply(random_C0_new_v9_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v9_q1))
})
random_AA_TT_AT_TA_C0_new_v9_q2 = apply(random_C0_new_v9_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v9_q2))
})
random_AA_TT_AT_TA_C0_new_v9_q3 = apply(random_C0_new_v9_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v9_q3))
})
random_AA_TT_AT_TA_C0_new_v9_q4 = apply(random_C0_new_v9_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v9_q4))
})
random_AA_TT_AT_TA_C0_new_v9_custom1 = apply(random_C0_new_v9_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v9_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v9_q1 = apply(tiling_C0_new_v9_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v9_q1))
})
tiling_AA_TT_AT_TA_C0_new_v9_q2 = apply(tiling_C0_new_v9_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v9_q2))
})
tiling_AA_TT_AT_TA_C0_new_v9_q3 = apply(tiling_C0_new_v9_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v9_q3))
})
tiling_AA_TT_AT_TA_C0_new_v9_q4 = apply(tiling_C0_new_v9_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v9_q4))
})
tiling_AA_TT_AT_TA_C0_new_v9_custom1 = apply(tiling_C0_new_v9_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v9_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v9_q1 = apply(chrv_C0_new_v9_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v9_q1))
})
chrv_AA_TT_AT_TA_C0_new_v9_q2 = apply(chrv_C0_new_v9_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v9_q2))
})
chrv_AA_TT_AT_TA_C0_new_v9_q3 = apply(chrv_C0_new_v9_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v9_q3))
})
chrv_AA_TT_AT_TA_C0_new_v9_q4 = apply(chrv_C0_new_v9_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v9_q4))
})
chrv_AA_TT_AT_TA_C0_new_v9_custom1 = apply(chrv_C0_new_v9_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v9_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v9_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v9_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v9_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v9_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v9_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v9_custom1)

random_AA_TT_AT_TA_C0_new_v9_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v9_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v9_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v9_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v9_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v9_custom1)

tiling_AA_TT_AT_TA_C0_new_v9_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v9_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v9_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v9_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v9_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v9_custom1)

chrv_AA_TT_AT_TA_C0_new_v9_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v9_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v9_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v9_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v9_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v9_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v10
nuc_AA_TT_AT_TA_C0_new_v10_q1 = apply(nuc_C0_new_v10_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v10_q1))
})
nuc_AA_TT_AT_TA_C0_new_v10_q2 = apply(nuc_C0_new_v10_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v10_q2))
})
nuc_AA_TT_AT_TA_C0_new_v10_q3 = apply(nuc_C0_new_v10_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v10_q3))
})
nuc_AA_TT_AT_TA_C0_new_v10_q4 = apply(nuc_C0_new_v10_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v10_q4))
})
nuc_AA_TT_AT_TA_C0_new_v10_custom1 = apply(nuc_C0_new_v10_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v10_custom1))
})

random_AA_TT_AT_TA_C0_new_v10_q1 = apply(random_C0_new_v10_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v10_q1))
})
random_AA_TT_AT_TA_C0_new_v10_q2 = apply(random_C0_new_v10_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v10_q2))
})
random_AA_TT_AT_TA_C0_new_v10_q3 = apply(random_C0_new_v10_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v10_q3))
})
random_AA_TT_AT_TA_C0_new_v10_q4 = apply(random_C0_new_v10_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v10_q4))
})
random_AA_TT_AT_TA_C0_new_v10_custom1 = apply(random_C0_new_v10_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v10_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v10_q1 = apply(tiling_C0_new_v10_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v10_q1))
})
tiling_AA_TT_AT_TA_C0_new_v10_q2 = apply(tiling_C0_new_v10_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v10_q2))
})
tiling_AA_TT_AT_TA_C0_new_v10_q3 = apply(tiling_C0_new_v10_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v10_q3))
})
tiling_AA_TT_AT_TA_C0_new_v10_q4 = apply(tiling_C0_new_v10_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v10_q4))
})
tiling_AA_TT_AT_TA_C0_new_v10_custom1 = apply(tiling_C0_new_v10_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v10_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v10_q1 = apply(chrv_C0_new_v10_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v10_q1))
})
chrv_AA_TT_AT_TA_C0_new_v10_q2 = apply(chrv_C0_new_v10_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v10_q2))
})
chrv_AA_TT_AT_TA_C0_new_v10_q3 = apply(chrv_C0_new_v10_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v10_q3))
})
chrv_AA_TT_AT_TA_C0_new_v10_q4 = apply(chrv_C0_new_v10_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v10_q4))
})
chrv_AA_TT_AT_TA_C0_new_v10_custom1 = apply(chrv_C0_new_v10_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v10_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v10_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v10_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v10_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v10_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v10_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v10_custom1)

random_AA_TT_AT_TA_C0_new_v10_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v10_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v10_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v10_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v10_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v10_custom1)

tiling_AA_TT_AT_TA_C0_new_v10_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v10_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v10_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v10_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v10_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v10_custom1)

chrv_AA_TT_AT_TA_C0_new_v10_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v10_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v10_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v10_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v10_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v10_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v11
nuc_AA_TT_AT_TA_C0_new_v11_q1 = apply(nuc_C0_new_v11_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v11_q1))
})
nuc_AA_TT_AT_TA_C0_new_v11_q2 = apply(nuc_C0_new_v11_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v11_q2))
})
nuc_AA_TT_AT_TA_C0_new_v11_q3 = apply(nuc_C0_new_v11_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v11_q3))
})
nuc_AA_TT_AT_TA_C0_new_v11_q4 = apply(nuc_C0_new_v11_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v11_q4))
})
nuc_AA_TT_AT_TA_C0_new_v11_custom1 = apply(nuc_C0_new_v11_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v11_custom1))
})

random_AA_TT_AT_TA_C0_new_v11_q1 = apply(random_C0_new_v11_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v11_q1))
})
random_AA_TT_AT_TA_C0_new_v11_q2 = apply(random_C0_new_v11_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v11_q2))
})
random_AA_TT_AT_TA_C0_new_v11_q3 = apply(random_C0_new_v11_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v11_q3))
})
random_AA_TT_AT_TA_C0_new_v11_q4 = apply(random_C0_new_v11_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v11_q4))
})
random_AA_TT_AT_TA_C0_new_v11_custom1 = apply(random_C0_new_v11_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v11_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v11_q1 = apply(tiling_C0_new_v11_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v11_q1))
})
tiling_AA_TT_AT_TA_C0_new_v11_q2 = apply(tiling_C0_new_v11_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v11_q2))
})
tiling_AA_TT_AT_TA_C0_new_v11_q3 = apply(tiling_C0_new_v11_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v11_q3))
})
tiling_AA_TT_AT_TA_C0_new_v11_q4 = apply(tiling_C0_new_v11_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v11_q4))
})
tiling_AA_TT_AT_TA_C0_new_v11_custom1 = apply(tiling_C0_new_v11_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v11_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v11_q1 = apply(chrv_C0_new_v11_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v11_q1))
})
chrv_AA_TT_AT_TA_C0_new_v11_q2 = apply(chrv_C0_new_v11_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v11_q2))
})
chrv_AA_TT_AT_TA_C0_new_v11_q3 = apply(chrv_C0_new_v11_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v11_q3))
})
chrv_AA_TT_AT_TA_C0_new_v11_q4 = apply(chrv_C0_new_v11_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v11_q4))
})
chrv_AA_TT_AT_TA_C0_new_v11_custom1 = apply(chrv_C0_new_v11_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v11_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v11_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v11_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v11_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v11_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v11_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v11_custom1)

random_AA_TT_AT_TA_C0_new_v11_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v11_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v11_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v11_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v11_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v11_custom1)

tiling_AA_TT_AT_TA_C0_new_v11_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v11_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v11_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v11_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v11_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v11_custom1)

chrv_AA_TT_AT_TA_C0_new_v11_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v11_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v11_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v11_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v11_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v11_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v12
nuc_AA_TT_AT_TA_C0_new_v12_q1 = apply(nuc_C0_new_v12_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v12_q1))
})
nuc_AA_TT_AT_TA_C0_new_v12_q2 = apply(nuc_C0_new_v12_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v12_q2))
})
nuc_AA_TT_AT_TA_C0_new_v12_q3 = apply(nuc_C0_new_v12_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v12_q3))
})
nuc_AA_TT_AT_TA_C0_new_v12_q4 = apply(nuc_C0_new_v12_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v12_q4))
})
nuc_AA_TT_AT_TA_C0_new_v12_custom1 = apply(nuc_C0_new_v12_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v12_custom1))
})

random_AA_TT_AT_TA_C0_new_v12_q1 = apply(random_C0_new_v12_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v12_q1))
})
random_AA_TT_AT_TA_C0_new_v12_q2 = apply(random_C0_new_v12_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v12_q2))
})
random_AA_TT_AT_TA_C0_new_v12_q3 = apply(random_C0_new_v12_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v12_q3))
})
random_AA_TT_AT_TA_C0_new_v12_q4 = apply(random_C0_new_v12_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v12_q4))
})
random_AA_TT_AT_TA_C0_new_v12_custom1 = apply(random_C0_new_v12_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v12_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v12_q1 = apply(tiling_C0_new_v12_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v12_q1))
})
tiling_AA_TT_AT_TA_C0_new_v12_q2 = apply(tiling_C0_new_v12_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v12_q2))
})
tiling_AA_TT_AT_TA_C0_new_v12_q3 = apply(tiling_C0_new_v12_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v12_q3))
})
tiling_AA_TT_AT_TA_C0_new_v12_q4 = apply(tiling_C0_new_v12_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v12_q4))
})
tiling_AA_TT_AT_TA_C0_new_v12_custom1 = apply(tiling_C0_new_v12_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v12_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v12_q1 = apply(chrv_C0_new_v12_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v12_q1))
})
chrv_AA_TT_AT_TA_C0_new_v12_q2 = apply(chrv_C0_new_v12_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v12_q2))
})
chrv_AA_TT_AT_TA_C0_new_v12_q3 = apply(chrv_C0_new_v12_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v12_q3))
})
chrv_AA_TT_AT_TA_C0_new_v12_q4 = apply(chrv_C0_new_v12_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v12_q4))
})
chrv_AA_TT_AT_TA_C0_new_v12_custom1 = apply(chrv_C0_new_v12_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v12_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v12_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v12_q1,
                                         q2 = nuc_AA_TT_AT_TA_C0_new_v12_q2,
                                         q3 = nuc_AA_TT_AT_TA_C0_new_v12_q3,
                                         q4 = nuc_AA_TT_AT_TA_C0_new_v12_q4,
                                         custom1 = nuc_AA_TT_AT_TA_C0_new_v12_custom1)

random_AA_TT_AT_TA_C0_new_v12_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v12_q1,
                                            q2 = random_AA_TT_AT_TA_C0_new_v12_q2,
                                            q3 = random_AA_TT_AT_TA_C0_new_v12_q3,
                                            q4 = random_AA_TT_AT_TA_C0_new_v12_q4,
                                            custom1 = random_AA_TT_AT_TA_C0_new_v12_custom1)

tiling_AA_TT_AT_TA_C0_new_v12_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v12_q1,
                                            q2 = tiling_AA_TT_AT_TA_C0_new_v12_q2,
                                            q3 = tiling_AA_TT_AT_TA_C0_new_v12_q3,
                                            q4 = tiling_AA_TT_AT_TA_C0_new_v12_q4,
                                            custom1 = tiling_AA_TT_AT_TA_C0_new_v12_custom1)

chrv_AA_TT_AT_TA_C0_new_v12_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v12_q1,
                                          q2 = chrv_AA_TT_AT_TA_C0_new_v12_q2,
                                          q3 = chrv_AA_TT_AT_TA_C0_new_v12_q3,
                                          q4 = chrv_AA_TT_AT_TA_C0_new_v12_q4,
                                          custom1 = chrv_AA_TT_AT_TA_C0_new_v12_custom1)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile + custom/library for C0_new_v13
nuc_AA_TT_AT_TA_C0_new_v13_q1 = apply(nuc_C0_new_v13_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v13_q1))
})
nuc_AA_TT_AT_TA_C0_new_v13_q2 = apply(nuc_C0_new_v13_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v13_q2))
})
nuc_AA_TT_AT_TA_C0_new_v13_q3 = apply(nuc_C0_new_v13_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v13_q3))
})
nuc_AA_TT_AT_TA_C0_new_v13_q4 = apply(nuc_C0_new_v13_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v13_q4))
})
nuc_AA_TT_AT_TA_C0_new_v13_custom1 = apply(nuc_C0_new_v13_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(nuc_C0_new_v13_custom1))
})

random_AA_TT_AT_TA_C0_new_v13_q1 = apply(random_C0_new_v13_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v13_q1))
})
random_AA_TT_AT_TA_C0_new_v13_q2 = apply(random_C0_new_v13_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v13_q2))
})
random_AA_TT_AT_TA_C0_new_v13_q3 = apply(random_C0_new_v13_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v13_q3))
})
random_AA_TT_AT_TA_C0_new_v13_q4 = apply(random_C0_new_v13_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v13_q4))
})
random_AA_TT_AT_TA_C0_new_v13_custom1 = apply(random_C0_new_v13_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(random_C0_new_v13_custom1))
})

tiling_AA_TT_AT_TA_C0_new_v13_q1 = apply(tiling_C0_new_v13_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v13_q1))
})
tiling_AA_TT_AT_TA_C0_new_v13_q2 = apply(tiling_C0_new_v13_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v13_q2))
})
tiling_AA_TT_AT_TA_C0_new_v13_q3 = apply(tiling_C0_new_v13_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v13_q3))
})
tiling_AA_TT_AT_TA_C0_new_v13_q4 = apply(tiling_C0_new_v13_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v13_q4))
})
tiling_AA_TT_AT_TA_C0_new_v13_custom1 = apply(tiling_C0_new_v13_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(tiling_C0_new_v13_custom1))
})

chrv_AA_TT_AT_TA_C0_new_v13_q1 = apply(chrv_C0_new_v13_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v13_q1))
})
chrv_AA_TT_AT_TA_C0_new_v13_q2 = apply(chrv_C0_new_v13_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v13_q2))
})
chrv_AA_TT_AT_TA_C0_new_v13_q3 = apply(chrv_C0_new_v13_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v13_q3))
})
chrv_AA_TT_AT_TA_C0_new_v13_q4 = apply(chrv_C0_new_v13_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v13_q4))
})
chrv_AA_TT_AT_TA_C0_new_v13_custom1 = apply(chrv_C0_new_v13_custom1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(chrv_C0_new_v13_custom1))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
nuc_AA_TT_AT_TA_C0_new_v13_q = data.frame(q1 = nuc_AA_TT_AT_TA_C0_new_v13_q1,
                                          q2 = nuc_AA_TT_AT_TA_C0_new_v13_q2,
                                          q3 = nuc_AA_TT_AT_TA_C0_new_v13_q3,
                                          q4 = nuc_AA_TT_AT_TA_C0_new_v13_q4,
                                          custom1 = nuc_AA_TT_AT_TA_C0_new_v13_custom1)

random_AA_TT_AT_TA_C0_new_v13_q = data.frame(q1 = random_AA_TT_AT_TA_C0_new_v13_q1,
                                             q2 = random_AA_TT_AT_TA_C0_new_v13_q2,
                                             q3 = random_AA_TT_AT_TA_C0_new_v13_q3,
                                             q4 = random_AA_TT_AT_TA_C0_new_v13_q4,
                                             custom1 = random_AA_TT_AT_TA_C0_new_v13_custom1)

tiling_AA_TT_AT_TA_C0_new_v13_q = data.frame(q1 = tiling_AA_TT_AT_TA_C0_new_v13_q1,
                                             q2 = tiling_AA_TT_AT_TA_C0_new_v13_q2,
                                             q3 = tiling_AA_TT_AT_TA_C0_new_v13_q3,
                                             q4 = tiling_AA_TT_AT_TA_C0_new_v13_q4,
                                             custom1 = tiling_AA_TT_AT_TA_C0_new_v13_custom1)

chrv_AA_TT_AT_TA_C0_new_v13_q = data.frame(q1 = chrv_AA_TT_AT_TA_C0_new_v13_q1,
                                           q2 = chrv_AA_TT_AT_TA_C0_new_v13_q2,
                                           q3 = chrv_AA_TT_AT_TA_C0_new_v13_q3,
                                           q4 = chrv_AA_TT_AT_TA_C0_new_v13_q4,
                                           custom1 = chrv_AA_TT_AT_TA_C0_new_v13_custom1)




# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Nucleosome Library

nuc_AATT_all_C_max = max(nuc_AA_TT_AT_TA_C26_q, nuc_AA_TT_AT_TA_C29_q, nuc_AA_TT_AT_TA_C31_q, nuc_AA_TT_AT_TA_C0_q,
                         nuc_AA_TT_AT_TA_C0_new_v1_q, nuc_AA_TT_AT_TA_C0_new_v2_q, nuc_AA_TT_AT_TA_C0_new_v3_q, nuc_AA_TT_AT_TA_C0_new_v4_q,
                         nuc_AA_TT_AT_TA_C0_new_v5_q, nuc_AA_TT_AT_TA_C0_new_v6_q, nuc_AA_TT_AT_TA_C0_new_v7_q, nuc_AA_TT_AT_TA_C0_new_v8_q,
                         nuc_AA_TT_AT_TA_C0_new_v9_q, nuc_AA_TT_AT_TA_C0_new_v10_q, nuc_AA_TT_AT_TA_C0_new_v11_q, nuc_AA_TT_AT_TA_C0_new_v12_q,
                         nuc_AA_TT_AT_TA_C0_new_v13_q)
nuc_AATT_all_C_min = min(nuc_AA_TT_AT_TA_C26_q, nuc_AA_TT_AT_TA_C29_q, nuc_AA_TT_AT_TA_C31_q, nuc_AA_TT_AT_TA_C0_q,
                         nuc_AA_TT_AT_TA_C0_new_v1_q, nuc_AA_TT_AT_TA_C0_new_v2_q, nuc_AA_TT_AT_TA_C0_new_v3_q, nuc_AA_TT_AT_TA_C0_new_v4_q,
                         nuc_AA_TT_AT_TA_C0_new_v5_q, nuc_AA_TT_AT_TA_C0_new_v6_q, nuc_AA_TT_AT_TA_C0_new_v7_q, nuc_AA_TT_AT_TA_C0_new_v8_q,
                         nuc_AA_TT_AT_TA_C0_new_v9_q, nuc_AA_TT_AT_TA_C0_new_v10_q, nuc_AA_TT_AT_TA_C0_new_v11_q, nuc_AA_TT_AT_TA_C0_new_v12_q,
                         nuc_AA_TT_AT_TA_C0_new_v13_q)

nuc_AATT_q1vq4_C26 = ggplot(data = nuc_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C26
ggsave("figures/periodicity/nuc_AATT_q1vq4_C26.png", plot = nuc_AATT_q1vq4_C26)

nuc_AATT_q1vq4_C29 = ggplot(data = nuc_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C29
ggsave("figures/periodicity/nuc_AATT_q1vq4_C29.png", plot = nuc_AATT_q1vq4_C29)

nuc_AATT_q1vq4_C31 = ggplot(data = nuc_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C31
ggsave("figures/periodicity/nuc_AATT_q1vq4_C31.png", plot = nuc_AATT_q1vq4_C31)

nuc_AATT_q1vq4_C0 = ggplot(data = nuc_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0.png", plot = nuc_AATT_q1vq4_C0)

nuc_AATT_q1vq4_C0_new_v1 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v1
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v1.png", plot = nuc_AATT_q1vq4_C0_new_v1)

nuc_AATT_q1vq4_C0_new_v2 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v2
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v2.png", plot = nuc_AATT_q1vq4_C0_new_v2)

nuc_AATT_q1vq4_C0_new_v3 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v3
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v3.png", plot = nuc_AATT_q1vq4_C0_new_v3)

nuc_AATT_q1vq4_C0_new_v4 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v4
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v4.png", plot = nuc_AATT_q1vq4_C0_new_v4)

nuc_AATT_q1vq4_C0_new_v5 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v5
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v5.png", plot = nuc_AATT_q1vq4_C0_new_v5)

nuc_AATT_q1vq4_C0_new_v6 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v6
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v6.png", plot = nuc_AATT_q1vq4_C0_new_v6)

nuc_AATT_q1vq4_C0_new_v7 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v7
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v7.png", plot = nuc_AATT_q1vq4_C0_new_v7)

nuc_AATT_q1vq4_C0_new_v8 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v8
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v8.png", plot = nuc_AATT_q1vq4_C0_new_v8)

nuc_AATT_q1vq4_C0_new_v9 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v9
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v9.png", plot = nuc_AATT_q1vq4_C0_new_v9)

nuc_AATT_q1vq4_C0_new_v10 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v10
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v10.png", plot = nuc_AATT_q1vq4_C0_new_v10)

nuc_AATT_q1vq4_C0_new_v11 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v11
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v11.png", plot = nuc_AATT_q1vq4_C0_new_v11)

nuc_AATT_q1vq4_C0_new_v12 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v12
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v12.png", plot = nuc_AATT_q1vq4_C0_new_v12)

nuc_AATT_q1vq4_C0_new_v13 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q1vq4_C0_new_v13
ggsave("figures/periodicity/nuc_AATT_q1vq4_C0_new_v13.png", plot = nuc_AATT_q1vq4_C0_new_v13)


nuc_AATT_q2vq3_C26 = ggplot(data = nuc_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C26
ggsave("figures/periodicity/nuc_AATT_q2vq3_C26.png", plot = nuc_AATT_q2vq3_C26)

nuc_AATT_q2vq3_C29 = ggplot(data = nuc_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C29
ggsave("figures/periodicity/nuc_AATT_q2vq3_C29.png", plot = nuc_AATT_q2vq3_C29)

nuc_AATT_q2vq3_C31 = ggplot(data = nuc_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C31
ggsave("figures/periodicity/nuc_AATT_q2vq3_C31.png", plot = nuc_AATT_q2vq3_C31)

nuc_AATT_q2vq3_C0 = ggplot(data = nuc_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0.png", plot = nuc_AATT_q2vq3_C0)

nuc_AATT_q2vq3_C0_new_v1 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v1
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v1.png", plot = nuc_AATT_q2vq3_C0_new_v1)

nuc_AATT_q2vq3_C0_new_v2 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v2
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v2.png", plot = nuc_AATT_q2vq3_C0_new_v2)

nuc_AATT_q2vq3_C0_new_v3 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v3
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v3.png", plot = nuc_AATT_q2vq3_C0_new_v3)

nuc_AATT_q2vq3_C0_new_v4 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v4
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v4.png", plot = nuc_AATT_q2vq3_C0_new_v4)

nuc_AATT_q2vq3_C0_new_v5 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v5
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v5.png", plot = nuc_AATT_q2vq3_C0_new_v5)

nuc_AATT_q2vq3_C0_new_v6 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v6
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v6.png", plot = nuc_AATT_q2vq3_C0_new_v6)

nuc_AATT_q2vq3_C0_new_v7 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v7
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v7.png", plot = nuc_AATT_q2vq3_C0_new_v7)

nuc_AATT_q2vq3_C0_new_v8 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v8
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v8.png", plot = nuc_AATT_q2vq3_C0_new_v8)

nuc_AATT_q2vq3_C0_new_v9 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v9
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v9.png", plot = nuc_AATT_q2vq3_C0_new_v9)

nuc_AATT_q2vq3_C0_new_v10 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v10
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v10.png", plot = nuc_AATT_q2vq3_C0_new_v10)

nuc_AATT_q2vq3_C0_new_v11 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v11
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v11.png", plot = nuc_AATT_q2vq3_C0_new_v11)

nuc_AATT_q2vq3_C0_new_v12 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v12
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v12.png", plot = nuc_AATT_q2vq3_C0_new_v12)

nuc_AATT_q2vq3_C0_new_v13 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_q2vq3_C0_new_v13
ggsave("figures/periodicity/nuc_AATT_q2vq3_C0_new_v13.png", plot = nuc_AATT_q2vq3_C0_new_v13)


nuc_AATT_custom1_C26 = ggplot(data = nuc_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C26") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C26
ggsave("figures/periodicity/nuc_AATT_custom1_C26.png", plot = nuc_AATT_custom1_C26)

nuc_AATT_custom1_C29 = ggplot(data = nuc_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C29") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C29
ggsave("figures/periodicity/nuc_AATT_custom1_C29.png", plot = nuc_AATT_custom1_C29)

nuc_AATT_custom1_C31 = ggplot(data = nuc_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C31") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C31
ggsave("figures/periodicity/nuc_AATT_custom1_C31.png", plot = nuc_AATT_custom1_C31)

nuc_AATT_custom1_C0 = ggplot(data = nuc_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0
ggsave("figures/periodicity/nuc_AATT_custom1_C0.png", plot = nuc_AATT_custom1_C0)

nuc_AATT_custom1_C0_new_v1 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v1") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v1
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v1.png", plot = nuc_AATT_custom1_C0_new_v1)

nuc_AATT_custom1_C0_new_v2 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v2") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v2
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v2.png", plot = nuc_AATT_custom1_C0_new_v2)

nuc_AATT_custom1_C0_new_v3 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v3") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v3
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v3.png", plot = nuc_AATT_custom1_C0_new_v3)

nuc_AATT_custom1_C0_new_v4 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v4") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v4
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v4.png", plot = nuc_AATT_custom1_C0_new_v4)

nuc_AATT_custom1_C0_new_v5 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v5") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v5
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v5.png", plot = nuc_AATT_custom1_C0_new_v5)

nuc_AATT_custom1_C0_new_v6 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v6") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v6
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v6.png", plot = nuc_AATT_custom1_C0_new_v6)

nuc_AATT_custom1_C0_new_v7 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v7") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v7
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v7.png", plot = nuc_AATT_custom1_C0_new_v7)

nuc_AATT_custom1_C0_new_v8 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v8") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v8
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v8.png", plot = nuc_AATT_custom1_C0_new_v8)

nuc_AATT_custom1_C0_new_v9 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v9") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v9
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v9.png", plot = nuc_AATT_custom1_C0_new_v9)

nuc_AATT_custom1_C0_new_v10 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v10") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v10
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v10.png", plot = nuc_AATT_custom1_C0_new_v10)

nuc_AATT_custom1_C0_new_v11 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v11") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v11
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v11.png", plot = nuc_AATT_custom1_C0_new_v11)

nuc_AATT_custom1_C0_new_v12 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v12") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v12
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v12.png", plot = nuc_AATT_custom1_C0_new_v12)

nuc_AATT_custom1_C0_new_v13 = ggplot(data = nuc_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Nucleosome Library, C0_new_v13") +
  ylim(nuc_AATT_all_C_min, nuc_AATT_all_C_max)
nuc_AATT_custom1_C0_new_v13
ggsave("figures/periodicity/nuc_AATT_custom1_C0_new_v13.png", plot = nuc_AATT_custom1_C0_new_v13)



# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Random Library

random_AATT_all_C_max = max(random_AA_TT_AT_TA_C26_q, random_AA_TT_AT_TA_C29_q, random_AA_TT_AT_TA_C31_q, random_AA_TT_AT_TA_C0_q,
                            random_AA_TT_AT_TA_C0_new_v1_q, random_AA_TT_AT_TA_C0_new_v2_q, random_AA_TT_AT_TA_C0_new_v3_q, random_AA_TT_AT_TA_C0_new_v4_q,
                            random_AA_TT_AT_TA_C0_new_v5_q, random_AA_TT_AT_TA_C0_new_v6_q, random_AA_TT_AT_TA_C0_new_v7_q, random_AA_TT_AT_TA_C0_new_v8_q,
                            random_AA_TT_AT_TA_C0_new_v9_q, random_AA_TT_AT_TA_C0_new_v10_q, random_AA_TT_AT_TA_C0_new_v11_q, random_AA_TT_AT_TA_C0_new_v12_q,
                            random_AA_TT_AT_TA_C0_new_v13_q)
random_AATT_all_C_min = min(random_AA_TT_AT_TA_C26_q, random_AA_TT_AT_TA_C29_q, random_AA_TT_AT_TA_C31_q, random_AA_TT_AT_TA_C0_q,
                            random_AA_TT_AT_TA_C0_new_v1_q, random_AA_TT_AT_TA_C0_new_v2_q, random_AA_TT_AT_TA_C0_new_v3_q, random_AA_TT_AT_TA_C0_new_v4_q,
                            random_AA_TT_AT_TA_C0_new_v5_q, random_AA_TT_AT_TA_C0_new_v6_q, random_AA_TT_AT_TA_C0_new_v7_q, random_AA_TT_AT_TA_C0_new_v8_q,
                            random_AA_TT_AT_TA_C0_new_v9_q, random_AA_TT_AT_TA_C0_new_v10_q, random_AA_TT_AT_TA_C0_new_v11_q, random_AA_TT_AT_TA_C0_new_v12_q,
                            random_AA_TT_AT_TA_C0_new_v13_q)

random_AATT_q1vq4_C26 = ggplot(data = random_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C26") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C26
ggsave("figures/periodicity/random_AATT_q1vq4_C26.png", plot = random_AATT_q1vq4_C26)

random_AATT_q1vq4_C29 = ggplot(data = random_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C29") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C29
ggsave("figures/periodicity/random_AATT_q1vq4_C29.png", plot = random_AATT_q1vq4_C29)

random_AATT_q1vq4_C31 = ggplot(data = random_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C31") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C31
ggsave("figures/periodicity/random_AATT_q1vq4_C31.png", plot = random_AATT_q1vq4_C31)

random_AATT_q1vq4_C0 = ggplot(data = random_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0
ggsave("figures/periodicity/random_AATT_q1vq4_C0.png", plot = random_AATT_q1vq4_C0)

random_AATT_q1vq4_C0_new_v1 = ggplot(data = random_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v1") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v1
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v1.png", plot = random_AATT_q1vq4_C0_new_v1)

random_AATT_q1vq4_C0_new_v2 = ggplot(data = random_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v2") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v2
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v2.png", plot = random_AATT_q1vq4_C0_new_v2)

random_AATT_q1vq4_C0_new_v3 = ggplot(data = random_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v3") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v3
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v3.png", plot = random_AATT_q1vq4_C0_new_v3)

random_AATT_q1vq4_C0_new_v4 = ggplot(data = random_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v4") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v4
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v4.png", plot = random_AATT_q1vq4_C0_new_v4)

random_AATT_q1vq4_C0_new_v5 = ggplot(data = random_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v5") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v5
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v5.png", plot = random_AATT_q1vq4_C0_new_v5)

random_AATT_q1vq4_C0_new_v6 = ggplot(data = random_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v6") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v6
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v6.png", plot = random_AATT_q1vq4_C0_new_v6)

random_AATT_q1vq4_C0_new_v7 = ggplot(data = random_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v7") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v7
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v7.png", plot = random_AATT_q1vq4_C0_new_v7)

random_AATT_q1vq4_C0_new_v8 = ggplot(data = random_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v8") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v8
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v8.png", plot = random_AATT_q1vq4_C0_new_v8)

random_AATT_q1vq4_C0_new_v9 = ggplot(data = random_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v9") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v9
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v9.png", plot = random_AATT_q1vq4_C0_new_v9)

random_AATT_q1vq4_C0_new_v10 = ggplot(data = random_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v10") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v10
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v10.png", plot = random_AATT_q1vq4_C0_new_v10)

random_AATT_q1vq4_C0_new_v11 = ggplot(data = random_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v11") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v11
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v11.png", plot = random_AATT_q1vq4_C0_new_v11)

random_AATT_q1vq4_C0_new_v12 = ggplot(data = random_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v12") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v12
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v12.png", plot = random_AATT_q1vq4_C0_new_v12)

random_AATT_q1vq4_C0_new_v13 = ggplot(data = random_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v13") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q1vq4_C0_new_v13
ggsave("figures/periodicity/random_AATT_q1vq4_C0_new_v13.png", plot = random_AATT_q1vq4_C0_new_v13)


random_AATT_q2vq3_C26 = ggplot(data = random_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C26") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C26
ggsave("figures/periodicity/random_AATT_q2vq3_C26.png", plot = random_AATT_q2vq3_C26)

random_AATT_q2vq3_C29 = ggplot(data = random_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C29") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C29
ggsave("figures/periodicity/random_AATT_q2vq3_C29.png", plot = random_AATT_q2vq3_C29)

random_AATT_q2vq3_C31 = ggplot(data = random_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C31") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C31
ggsave("figures/periodicity/random_AATT_q2vq3_C31.png", plot = random_AATT_q2vq3_C31)

random_AATT_q2vq3_C0 = ggplot(data = random_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0
ggsave("figures/periodicity/random_AATT_q2vq3_C0.png", plot = random_AATT_q2vq3_C0)

random_AATT_q2vq3_C0_new_v1 = ggplot(data = random_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v1") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v1
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v1.png", plot = random_AATT_q2vq3_C0_new_v1)

random_AATT_q2vq3_C0_new_v2 = ggplot(data = random_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v2") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v2
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v2.png", plot = random_AATT_q2vq3_C0_new_v2)

random_AATT_q2vq3_C0_new_v3 = ggplot(data = random_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v3") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v3
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v3.png", plot = random_AATT_q2vq3_C0_new_v3)

random_AATT_q2vq3_C0_new_v4 = ggplot(data = random_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v4") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v4
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v4.png", plot = random_AATT_q2vq3_C0_new_v4)

random_AATT_q2vq3_C0_new_v5 = ggplot(data = random_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v5") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v5
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v5.png", plot = random_AATT_q2vq3_C0_new_v5)

random_AATT_q2vq3_C0_new_v6 = ggplot(data = random_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v6") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v6
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v6.png", plot = random_AATT_q2vq3_C0_new_v6)

random_AATT_q2vq3_C0_new_v7 = ggplot(data = random_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v7") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v7
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v7.png", plot = random_AATT_q2vq3_C0_new_v7)

random_AATT_q2vq3_C0_new_v8 = ggplot(data = random_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v8") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v8
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v8.png", plot = random_AATT_q2vq3_C0_new_v8)

random_AATT_q2vq3_C0_new_v9 = ggplot(data = random_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v9") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v9
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v9.png", plot = random_AATT_q2vq3_C0_new_v9)

random_AATT_q2vq3_C0_new_v10 = ggplot(data = random_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v10") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v10
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v10.png", plot = random_AATT_q2vq3_C0_new_v10)

random_AATT_q2vq3_C0_new_v11 = ggplot(data = random_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v11") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v11
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v11.png", plot = random_AATT_q2vq3_C0_new_v11)

random_AATT_q2vq3_C0_new_v12 = ggplot(data = random_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v12") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v12
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v12.png", plot = random_AATT_q2vq3_C0_new_v12)

random_AATT_q2vq3_C0_new_v13 = ggplot(data = random_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v13") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_q2vq3_C0_new_v13
ggsave("figures/periodicity/random_AATT_q2vq3_C0_new_v13.png", plot = random_AATT_q2vq3_C0_new_v13)


random_AATT_custom1_C26 = ggplot(data = random_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C26") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C26
ggsave("figures/periodicity/random_AATT_custom1_C26.png", plot = random_AATT_custom1_C26)

random_AATT_custom1_C29 = ggplot(data = random_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C29") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C29
ggsave("figures/periodicity/random_AATT_custom1_C29.png", plot = random_AATT_custom1_C29)

random_AATT_custom1_C31 = ggplot(data = random_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C31") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C31
ggsave("figures/periodicity/random_AATT_custom1_C31.png", plot = random_AATT_custom1_C31)

random_AATT_custom1_C0 = ggplot(data = random_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0
ggsave("figures/periodicity/random_AATT_custom1_C0.png", plot = random_AATT_custom1_C0)

random_AATT_custom1_C0_new_v1 = ggplot(data = random_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v1") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v1
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v1.png", plot = random_AATT_custom1_C0_new_v1)

random_AATT_custom1_C0_new_v2 = ggplot(data = random_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v2") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v2
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v2.png", plot = random_AATT_custom1_C0_new_v2)

random_AATT_custom1_C0_new_v3 = ggplot(data = random_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v3") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v3
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v3.png", plot = random_AATT_custom1_C0_new_v3)

random_AATT_custom1_C0_new_v4 = ggplot(data = random_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v4") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v4
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v4.png", plot = random_AATT_custom1_C0_new_v4)

random_AATT_custom1_C0_new_v5 = ggplot(data = random_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v5") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v5
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v5.png", plot = random_AATT_custom1_C0_new_v5)

random_AATT_custom1_C0_new_v6 = ggplot(data = random_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v6") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v6
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v6.png", plot = random_AATT_custom1_C0_new_v6)

random_AATT_custom1_C0_new_v7 = ggplot(data = random_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v7") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v7
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v7.png", plot = random_AATT_custom1_C0_new_v7)

random_AATT_custom1_C0_new_v8 = ggplot(data = random_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v8") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v8
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v8.png", plot = random_AATT_custom1_C0_new_v8)

random_AATT_custom1_C0_new_v9 = ggplot(data = random_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v9") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v9
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v9.png", plot = random_AATT_custom1_C0_new_v9)

random_AATT_custom1_C0_new_v10 = ggplot(data = random_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v10") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v10
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v10.png", plot = random_AATT_custom1_C0_new_v10)

random_AATT_custom1_C0_new_v11 = ggplot(data = random_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v11") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v11
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v11.png", plot = random_AATT_custom1_C0_new_v11)

random_AATT_custom1_C0_new_v12 = ggplot(data = random_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v12") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v12
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v12.png", plot = random_AATT_custom1_C0_new_v12)

random_AATT_custom1_C0_new_v13 = ggplot(data = random_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, C0_new_v13") +
  ylim(random_AATT_all_C_min, random_AATT_all_C_max)
random_AATT_custom1_C0_new_v13
ggsave("figures/periodicity/random_AATT_custom1_C0_new_v13.png", plot = random_AATT_custom1_C0_new_v13)




# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for Tiling Library

tiling_AATT_all_C_max = max(tiling_AA_TT_AT_TA_C26_q, tiling_AA_TT_AT_TA_C29_q, tiling_AA_TT_AT_TA_C31_q, tiling_AA_TT_AT_TA_C0_q,
                            tiling_AA_TT_AT_TA_C0_new_v1_q, tiling_AA_TT_AT_TA_C0_new_v2_q, tiling_AA_TT_AT_TA_C0_new_v3_q, tiling_AA_TT_AT_TA_C0_new_v4_q,
                            tiling_AA_TT_AT_TA_C0_new_v5_q, tiling_AA_TT_AT_TA_C0_new_v6_q, tiling_AA_TT_AT_TA_C0_new_v7_q, tiling_AA_TT_AT_TA_C0_new_v8_q,
                            tiling_AA_TT_AT_TA_C0_new_v9_q, tiling_AA_TT_AT_TA_C0_new_v10_q, tiling_AA_TT_AT_TA_C0_new_v11_q, tiling_AA_TT_AT_TA_C0_new_v12_q,
                            tiling_AA_TT_AT_TA_C0_new_v13_q)
tiling_AATT_all_C_min = min(tiling_AA_TT_AT_TA_C26_q, tiling_AA_TT_AT_TA_C29_q, tiling_AA_TT_AT_TA_C31_q, tiling_AA_TT_AT_TA_C0_q,
                            tiling_AA_TT_AT_TA_C0_new_v1_q, tiling_AA_TT_AT_TA_C0_new_v2_q, tiling_AA_TT_AT_TA_C0_new_v3_q, tiling_AA_TT_AT_TA_C0_new_v4_q,
                            tiling_AA_TT_AT_TA_C0_new_v5_q, tiling_AA_TT_AT_TA_C0_new_v6_q, tiling_AA_TT_AT_TA_C0_new_v7_q, tiling_AA_TT_AT_TA_C0_new_v8_q,
                            tiling_AA_TT_AT_TA_C0_new_v9_q, tiling_AA_TT_AT_TA_C0_new_v10_q, tiling_AA_TT_AT_TA_C0_new_v11_q, tiling_AA_TT_AT_TA_C0_new_v12_q,
                            tiling_AA_TT_AT_TA_C0_new_v13_q)

tiling_AATT_q1vq4_C26 = ggplot(data = tiling_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C26") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C26
ggsave("figures/periodicity/tiling_AATT_q1vq4_C26.png", plot = tiling_AATT_q1vq4_C26)

tiling_AATT_q1vq4_C29 = ggplot(data = tiling_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C29") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C29
ggsave("figures/periodicity/tiling_AATT_q1vq4_C29.png", plot = tiling_AATT_q1vq4_C29)

tiling_AATT_q1vq4_C31 = ggplot(data = tiling_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C31") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C31
ggsave("figures/periodicity/tiling_AATT_q1vq4_C31.png", plot = tiling_AATT_q1vq4_C31)

tiling_AATT_q1vq4_C0 = ggplot(data = tiling_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0.png", plot = tiling_AATT_q1vq4_C0)

tiling_AATT_q1vq4_C0_new_v1 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v1
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v1.png", plot = tiling_AATT_q1vq4_C0_new_v1)

tiling_AATT_q1vq4_C0_new_v2 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v2
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v2.png", plot = tiling_AATT_q1vq4_C0_new_v2)

tiling_AATT_q1vq4_C0_new_v3 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v3
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v3.png", plot = tiling_AATT_q1vq4_C0_new_v3)

tiling_AATT_q1vq4_C0_new_v4 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v4
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v4.png", plot = tiling_AATT_q1vq4_C0_new_v4)

tiling_AATT_q1vq4_C0_new_v5 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v5
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v5.png", plot = tiling_AATT_q1vq4_C0_new_v5)

tiling_AATT_q1vq4_C0_new_v6 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v6
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v6.png", plot = tiling_AATT_q1vq4_C0_new_v6)

tiling_AATT_q1vq4_C0_new_v7 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v7
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v7.png", plot = tiling_AATT_q1vq4_C0_new_v7)

tiling_AATT_q1vq4_C0_new_v8 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v8
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v8.png", plot = tiling_AATT_q1vq4_C0_new_v8)

tiling_AATT_q1vq4_C0_new_v9 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v9
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v9.png", plot = tiling_AATT_q1vq4_C0_new_v9)

tiling_AATT_q1vq4_C0_new_v10 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v10
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v10.png", plot = tiling_AATT_q1vq4_C0_new_v10)

tiling_AATT_q1vq4_C0_new_v11 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v11
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v11.png", plot = tiling_AATT_q1vq4_C0_new_v11)

tiling_AATT_q1vq4_C0_new_v12 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v12
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v12.png", plot = tiling_AATT_q1vq4_C0_new_v12)

tiling_AATT_q1vq4_C0_new_v13 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q1vq4_C0_new_v13
ggsave("figures/periodicity/tiling_AATT_q1vq4_C0_new_v13.png", plot = tiling_AATT_q1vq4_C0_new_v13)


tiling_AATT_q2vq3_C26 = ggplot(data = tiling_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C26") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C26
ggsave("figures/periodicity/tiling_AATT_q2vq3_C26.png", plot = tiling_AATT_q2vq3_C26)

tiling_AATT_q2vq3_C29 = ggplot(data = tiling_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C29") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C29
ggsave("figures/periodicity/tiling_AATT_q2vq3_C29.png", plot = tiling_AATT_q2vq3_C29)

tiling_AATT_q2vq3_C31 = ggplot(data = tiling_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C31") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C31
ggsave("figures/periodicity/tiling_AATT_q2vq3_C31.png", plot = tiling_AATT_q2vq3_C31)

tiling_AATT_q2vq3_C0 = ggplot(data = tiling_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0.png", plot = tiling_AATT_q2vq3_C0)

tiling_AATT_q2vq3_C0_new_v1 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v1
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v1.png", plot = tiling_AATT_q2vq3_C0_new_v1)

tiling_AATT_q2vq3_C0_new_v2 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v2
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v2.png", plot = tiling_AATT_q2vq3_C0_new_v2)

tiling_AATT_q2vq3_C0_new_v3 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v3
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v3.png", plot = tiling_AATT_q2vq3_C0_new_v3)

tiling_AATT_q2vq3_C0_new_v4 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v4
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v4.png", plot = tiling_AATT_q2vq3_C0_new_v4)

tiling_AATT_q2vq3_C0_new_v5 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v5
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v5.png", plot = tiling_AATT_q2vq3_C0_new_v5)

tiling_AATT_q2vq3_C0_new_v6 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v6
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v6.png", plot = tiling_AATT_q2vq3_C0_new_v6)

tiling_AATT_q2vq3_C0_new_v7 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v7
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v7.png", plot = tiling_AATT_q2vq3_C0_new_v7)

tiling_AATT_q2vq3_C0_new_v8 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v8
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v8.png", plot = tiling_AATT_q2vq3_C0_new_v8)

tiling_AATT_q2vq3_C0_new_v9 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v9
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v9.png", plot = tiling_AATT_q2vq3_C0_new_v9)

tiling_AATT_q2vq3_C0_new_v10 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v10
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v10.png", plot = tiling_AATT_q2vq3_C0_new_v10)

tiling_AATT_q2vq3_C0_new_v11 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v11
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v11.png", plot = tiling_AATT_q2vq3_C0_new_v11)

tiling_AATT_q2vq3_C0_new_v12 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v12
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v12.png", plot = tiling_AATT_q2vq3_C0_new_v12)

tiling_AATT_q2vq3_C0_new_v13 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_q2vq3_C0_new_v13
ggsave("figures/periodicity/tiling_AATT_q2vq3_C0_new_v13.png", plot = tiling_AATT_q2vq3_C0_new_v13)


tiling_AATT_custom1_C26 = ggplot(data = tiling_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C26") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C26
ggsave("figures/periodicity/tiling_AATT_custom1_C26.png", plot = tiling_AATT_custom1_C26)

tiling_AATT_custom1_C29 = ggplot(data = tiling_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C29") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C29
ggsave("figures/periodicity/tiling_AATT_custom1_C29.png", plot = tiling_AATT_custom1_C29)

tiling_AATT_custom1_C31 = ggplot(data = tiling_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C31") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C31
ggsave("figures/periodicity/tiling_AATT_custom1_C31.png", plot = tiling_AATT_custom1_C31)

tiling_AATT_custom1_C0 = ggplot(data = tiling_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0
ggsave("figures/periodicity/tiling_AATT_custom1_C0.png", plot = tiling_AATT_custom1_C0)

tiling_AATT_custom1_C0_new_v1 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v1") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v1
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v1.png", plot = tiling_AATT_custom1_C0_new_v1)

tiling_AATT_custom1_C0_new_v2 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v2") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v2
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v2.png", plot = tiling_AATT_custom1_C0_new_v2)

tiling_AATT_custom1_C0_new_v3 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v3") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v3
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v3.png", plot = tiling_AATT_custom1_C0_new_v3)

tiling_AATT_custom1_C0_new_v4 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v4") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v4
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v4.png", plot = tiling_AATT_custom1_C0_new_v4)

tiling_AATT_custom1_C0_new_v5 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v5") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v5
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v5.png", plot = tiling_AATT_custom1_C0_new_v5)

tiling_AATT_custom1_C0_new_v6 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v6") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v6
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v6.png", plot = tiling_AATT_custom1_C0_new_v6)

tiling_AATT_custom1_C0_new_v7 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v7") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v7
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v7.png", plot = tiling_AATT_custom1_C0_new_v7)

tiling_AATT_custom1_C0_new_v8 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v8") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v8
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v8.png", plot = tiling_AATT_custom1_C0_new_v8)

tiling_AATT_custom1_C0_new_v9 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v9") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v9
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v9.png", plot = tiling_AATT_custom1_C0_new_v9)

tiling_AATT_custom1_C0_new_v10 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v10") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v10
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v10.png", plot = tiling_AATT_custom1_C0_new_v10)

tiling_AATT_custom1_C0_new_v11 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v11") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v11
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v11.png", plot = tiling_AATT_custom1_C0_new_v11)

tiling_AATT_custom1_C0_new_v12 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v12") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v12
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v12.png", plot = tiling_AATT_custom1_C0_new_v12)

tiling_AATT_custom1_C0_new_v13 = ggplot(data = tiling_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, C0_new_v13") +
  ylim(tiling_AATT_all_C_min, tiling_AATT_all_C_max)
tiling_AATT_custom1_C0_new_v13
ggsave("figures/periodicity/tiling_AATT_custom1_C0_new_v13.png", plot = tiling_AATT_custom1_C0_new_v13)





# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
# q2 v q3, and custom) and for each C26/C29/C31/C0/C0_news for ChrV Library

chrv_AATT_all_C_max = max(chrv_AA_TT_AT_TA_C26_q, chrv_AA_TT_AT_TA_C29_q, chrv_AA_TT_AT_TA_C31_q, chrv_AA_TT_AT_TA_C0_q,
                          chrv_AA_TT_AT_TA_C0_new_v1_q, chrv_AA_TT_AT_TA_C0_new_v2_q, chrv_AA_TT_AT_TA_C0_new_v3_q, chrv_AA_TT_AT_TA_C0_new_v4_q,
                          chrv_AA_TT_AT_TA_C0_new_v5_q, chrv_AA_TT_AT_TA_C0_new_v6_q, chrv_AA_TT_AT_TA_C0_new_v7_q, chrv_AA_TT_AT_TA_C0_new_v8_q,
                          chrv_AA_TT_AT_TA_C0_new_v9_q, chrv_AA_TT_AT_TA_C0_new_v10_q, chrv_AA_TT_AT_TA_C0_new_v11_q, chrv_AA_TT_AT_TA_C0_new_v12_q,
                          chrv_AA_TT_AT_TA_C0_new_v13_q)
chrv_AATT_all_C_min = min(chrv_AA_TT_AT_TA_C26_q, chrv_AA_TT_AT_TA_C29_q, chrv_AA_TT_AT_TA_C31_q, chrv_AA_TT_AT_TA_C0_q,
                          chrv_AA_TT_AT_TA_C0_new_v1_q, chrv_AA_TT_AT_TA_C0_new_v2_q, chrv_AA_TT_AT_TA_C0_new_v3_q, chrv_AA_TT_AT_TA_C0_new_v4_q,
                          chrv_AA_TT_AT_TA_C0_new_v5_q, chrv_AA_TT_AT_TA_C0_new_v6_q, chrv_AA_TT_AT_TA_C0_new_v7_q, chrv_AA_TT_AT_TA_C0_new_v8_q,
                          chrv_AA_TT_AT_TA_C0_new_v9_q, chrv_AA_TT_AT_TA_C0_new_v10_q, chrv_AA_TT_AT_TA_C0_new_v11_q, chrv_AA_TT_AT_TA_C0_new_v12_q,
                          chrv_AA_TT_AT_TA_C0_new_v13_q)

chrv_AATT_q1vq4_C26 = ggplot(data = chrv_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C26") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C26
ggsave("figures/periodicity/chrv_AATT_q1vq4_C26.png", plot = chrv_AATT_q1vq4_C26)

chrv_AATT_q1vq4_C29 = ggplot(data = chrv_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C29") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C29
ggsave("figures/periodicity/chrv_AATT_q1vq4_C29.png", plot = chrv_AATT_q1vq4_C29)

chrv_AATT_q1vq4_C31 = ggplot(data = chrv_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C31") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C31
ggsave("figures/periodicity/chrv_AATT_q1vq4_C31.png", plot = chrv_AATT_q1vq4_C31)

chrv_AATT_q1vq4_C0 = ggplot(data = chrv_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0.png", plot = chrv_AATT_q1vq4_C0)

chrv_AATT_q1vq4_C0_new_v1 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v1
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v1.png", plot = chrv_AATT_q1vq4_C0_new_v1)

chrv_AATT_q1vq4_C0_new_v2 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v2
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v2.png", plot = chrv_AATT_q1vq4_C0_new_v2)

chrv_AATT_q1vq4_C0_new_v3 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v3
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v3.png", plot = chrv_AATT_q1vq4_C0_new_v3)

chrv_AATT_q1vq4_C0_new_v4 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v4
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v4.png", plot = chrv_AATT_q1vq4_C0_new_v4)

chrv_AATT_q1vq4_C0_new_v5 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v5
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v5.png", plot = chrv_AATT_q1vq4_C0_new_v5)

chrv_AATT_q1vq4_C0_new_v6 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v6
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v6.png", plot = chrv_AATT_q1vq4_C0_new_v6)

chrv_AATT_q1vq4_C0_new_v7 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v7
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v7.png", plot = chrv_AATT_q1vq4_C0_new_v7)

chrv_AATT_q1vq4_C0_new_v8 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v8
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v8.png", plot = chrv_AATT_q1vq4_C0_new_v8)

chrv_AATT_q1vq4_C0_new_v9 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v9
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v9.png", plot = chrv_AATT_q1vq4_C0_new_v9)

chrv_AATT_q1vq4_C0_new_v10 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v10
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v10.png", plot = chrv_AATT_q1vq4_C0_new_v10)

chrv_AATT_q1vq4_C0_new_v11 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v11
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v11.png", plot = chrv_AATT_q1vq4_C0_new_v11)

chrv_AATT_q1vq4_C0_new_v12 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v12
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v12.png", plot = chrv_AATT_q1vq4_C0_new_v12)

chrv_AATT_q1vq4_C0_new_v13 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q1vq4_C0_new_v13
ggsave("figures/periodicity/chrv_AATT_q1vq4_C0_new_v13.png", plot = chrv_AATT_q1vq4_C0_new_v13)


chrv_AATT_q2vq3_C26 = ggplot(data = chrv_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C26") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C26
ggsave("figures/periodicity/chrv_AATT_q2vq3_C26.png", plot = chrv_AATT_q2vq3_C26)

chrv_AATT_q2vq3_C29 = ggplot(data = chrv_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C29") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C29
ggsave("figures/periodicity/chrv_AATT_q2vq3_C29.png", plot = chrv_AATT_q2vq3_C29)

chrv_AATT_q2vq3_C31 = ggplot(data = chrv_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C31") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C31
ggsave("figures/periodicity/chrv_AATT_q2vq3_C31.png", plot = chrv_AATT_q2vq3_C31)

chrv_AATT_q2vq3_C0 = ggplot(data = chrv_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0.png", plot = chrv_AATT_q2vq3_C0)

chrv_AATT_q2vq3_C0_new_v1 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v1
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v1.png", plot = chrv_AATT_q2vq3_C0_new_v1)

chrv_AATT_q2vq3_C0_new_v2 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v2
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v2.png", plot = chrv_AATT_q2vq3_C0_new_v2)

chrv_AATT_q2vq3_C0_new_v3 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v3
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v3.png", plot = chrv_AATT_q2vq3_C0_new_v3)

chrv_AATT_q2vq3_C0_new_v4 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v4
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v4.png", plot = chrv_AATT_q2vq3_C0_new_v4)

chrv_AATT_q2vq3_C0_new_v5 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v5
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v5.png", plot = chrv_AATT_q2vq3_C0_new_v5)

chrv_AATT_q2vq3_C0_new_v6 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v6
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v6.png", plot = chrv_AATT_q2vq3_C0_new_v6)

chrv_AATT_q2vq3_C0_new_v7 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v7
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v7.png", plot = chrv_AATT_q2vq3_C0_new_v7)

chrv_AATT_q2vq3_C0_new_v8 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v8
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v8.png", plot = chrv_AATT_q2vq3_C0_new_v8)

chrv_AATT_q2vq3_C0_new_v9 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v9
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v9.png", plot = chrv_AATT_q2vq3_C0_new_v9)

chrv_AATT_q2vq3_C0_new_v10 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v10
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v10.png", plot = chrv_AATT_q2vq3_C0_new_v10)

chrv_AATT_q2vq3_C0_new_v11 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v11
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v11.png", plot = chrv_AATT_q2vq3_C0_new_v11)

chrv_AATT_q2vq3_C0_new_v12 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v12
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v12.png", plot = chrv_AATT_q2vq3_C0_new_v12)

chrv_AATT_q2vq3_C0_new_v13 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_q2vq3_C0_new_v13
ggsave("figures/periodicity/chrv_AATT_q2vq3_C0_new_v13.png", plot = chrv_AATT_q2vq3_C0_new_v13)


chrv_AATT_custom1_C26 = ggplot(data = chrv_AA_TT_AT_TA_C26_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C26 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C26") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C26
ggsave("figures/periodicity/chrv_AATT_custom1_C26.png", plot = chrv_AATT_custom1_C26)

chrv_AATT_custom1_C29 = ggplot(data = chrv_AA_TT_AT_TA_C29_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C29 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C29") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C29
ggsave("figures/periodicity/chrv_AATT_custom1_C29.png", plot = chrv_AATT_custom1_C29)

chrv_AATT_custom1_C31 = ggplot(data = chrv_AA_TT_AT_TA_C31_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C31 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C31") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C31
ggsave("figures/periodicity/chrv_AATT_custom1_C31.png", plot = chrv_AATT_custom1_C31)

chrv_AATT_custom1_C0 = ggplot(data = chrv_AA_TT_AT_TA_C0_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0
ggsave("figures/periodicity/chrv_AATT_custom1_C0.png", plot = chrv_AATT_custom1_C0)

chrv_AATT_custom1_C0_new_v1 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v1_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v1") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v1
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v1.png", plot = chrv_AATT_custom1_C0_new_v1)

chrv_AATT_custom1_C0_new_v2 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v2_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v2") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v2
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v2.png", plot = chrv_AATT_custom1_C0_new_v2)

chrv_AATT_custom1_C0_new_v3 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v3_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v3") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v3
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v3.png", plot = chrv_AATT_custom1_C0_new_v3)

chrv_AATT_custom1_C0_new_v4 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v4_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v4") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v4
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v4.png", plot = chrv_AATT_custom1_C0_new_v4)

chrv_AATT_custom1_C0_new_v5 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v5_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v5") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v5
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v5.png", plot = chrv_AATT_custom1_C0_new_v5)

chrv_AATT_custom1_C0_new_v6 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v6_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v6") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v6
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v6.png", plot = chrv_AATT_custom1_C0_new_v6)

chrv_AATT_custom1_C0_new_v7 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v7_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v7") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v7
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v7.png", plot = chrv_AATT_custom1_C0_new_v7)

chrv_AATT_custom1_C0_new_v8 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v8_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v8") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v8
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v8.png", plot = chrv_AATT_custom1_C0_new_v8)

chrv_AATT_custom1_C0_new_v9 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v9_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v9") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v9
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v9.png", plot = chrv_AATT_custom1_C0_new_v9)

chrv_AATT_custom1_C0_new_v10 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v10_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v10") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v10
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v10.png", plot = chrv_AATT_custom1_C0_new_v10)

chrv_AATT_custom1_C0_new_v11 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v11_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v11") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v11
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v11.png", plot = chrv_AATT_custom1_C0_new_v11)

chrv_AATT_custom1_C0_new_v12 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v12_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v12") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v12
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v12.png", plot = chrv_AATT_custom1_C0_new_v12)

chrv_AATT_custom1_C0_new_v13 = ggplot(data = chrv_AA_TT_AT_TA_C0_new_v13_q, aes(x = 1:49)) +
  geom_line(aes(y=custom1, color="70-100%")) +
  scale_color_manual(breaks = c("70-100%"),
                     values = c("dodgerblue"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in ChrV Library, C0_new_v13") +
  ylim(chrv_AATT_all_C_min, chrv_AATT_all_C_max)
chrv_AATT_custom1_C0_new_v13
ggsave("figures/periodicity/chrv_AATT_custom1_C0_new_v13.png", plot = chrv_AATT_custom1_C0_new_v13)




