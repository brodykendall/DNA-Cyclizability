library(tidyverse)
library(gridExtra)

source("scripts/functions/sequence-functions.R")
source("scripts/functions/plotting-functions.R")

dat_tiling = read.csv("data/Created/tiling_ir_lstm_cn_tiling_post_smoothed.csv")
dat_chrv = read.csv("data/Created/chrv_ir_lstm_cn_tiling_post_smoothed_matched.csv")

# Trained on Random:
dat_chrI = read.csv("data/Created/ir_lstm_cn_random_yeast_chrI_1bpresolution_subsequence50_smoothC0.csv")

# Trained on Tiling:
dat_chrI = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrI_1bpresolution_subsequence50_smoothC0.csv")
# dat_chrV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0.csv")
# dat_chrV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_adj.csv")
dat_chrVI = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrVI_1bpresolution_subsequence50_smoothC0.csv")
# dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_full_reconstructed_smoothC0.csv")

dat_chrV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11.csv")
# dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_full_reconstructed_smoothC0_10_11.csv")
dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_contracted_smoothC0_10_11.csv")

dat_exp_random1 = read.csv("data/Created/ir_lstm_cn_tiling_expanded_random_library_smoothC0_10_11.csv")
dat_exp_random2 = read.csv("data/Created/ir_lstm_cn_tiling_expanded_random_library2_smoothC0_10_11.csv")

# Trained on ChrV:
# dat_chrVI = read.csv("data/Created/ir_lstm_cn_chrv_yeast_chrVI_1bpresolution_subsequence50_smoothC0.csv")


# dat_tiling = dat_tiling[which(!is.na(dat_tiling$C26)),]
# dat_chrv = dat_chrv[which(!is.na(dat_chrv$C26)),]
# dat_tiling = dat_tiling[which(!is.na(dat_tiling$smooth_10.4bp_cn_mean2)),]
dat_tiling = dat_tiling[which(!is.na(dat_tiling$smoothC0)),]
dat_chrv = dat_chrv[which(!is.na(dat_chrv$smooth_10.4bp_cn_mean2)),]

dat_chrI = dat_chrI[which(!is.na(dat_chrI$smoothC0)),]
dat_chrV = dat_chrV[which(!is.na(dat_chrV$smoothC0)),]
dat_chrVI = dat_chrVI[which(!is.na(dat_chrVI$smoothC0)),]

dat_exp_random1 = dat_exp_random1[which(!is.na(dat_exp_random1$smoothC0)),]
dat_exp_random1 = dat_exp_random1[which(!is.na(dat_exp_random1$orig_idx)),]
dat_exp_random2 = dat_exp_random2[which(!is.na(dat_exp_random2$smoothC0)),]
dat_exp_random2 = dat_exp_random2[which(!is.na(dat_exp_random2$orig_idx)),]

dat_tiling = dat_tiling %>%
  rename(x50mer = sequence)
dat_chrv = dat_chrv %>%
  rename(x50mer = sequence)

dat_chrI = dat_chrI %>%
  rename(x50mer = sequence)
dat_chrV = dat_chrV %>%
  rename(x50mer = sequence)
dat_chrVI = dat_chrVI %>%
  rename(x50mer = sequence)

dat_exp_random1 = dat_exp_random1 %>%
  rename(x50mer = sequence)
dat_exp_random2 = dat_exp_random2 %>%
  rename(x50mer = sequence)

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")

sequence_1_df_tiling <- dat_tiling %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_tiling <- dat_tiling %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_tiling), byrow = TRUE)} %>%
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


sequence_1_df_chrv <- dat_chrv %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrv <- dat_chrv %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv), byrow = TRUE)} %>%
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

sequence_1_df_chrI <- dat_chrI %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrI <- dat_chrI %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrI), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_chrI <- sequence_1_df_chrI %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_chrI <- sequence_2_df_chrI %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_1_df_chrV <- dat_chrV %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrV <- dat_chrV %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrV), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_chrV <- sequence_1_df_chrV %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_chrV <- sequence_2_df_chrV %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_1_df_chrVI <- dat_chrVI %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrVI <- dat_chrVI %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrVI), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_chrVI <- sequence_1_df_chrVI %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_chrVI <- sequence_2_df_chrVI %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_1_df_exp_random1 <- dat_exp_random1 %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_exp_random1 <- dat_exp_random1 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_exp_random1), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_exp_random1 <- sequence_1_df_exp_random1 %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_exp_random1 <- sequence_2_df_exp_random1 %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

sequence_1_df_exp_random2 <- dat_exp_random2 %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_exp_random2 <- dat_exp_random2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_exp_random2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_exp_random2 <- sequence_1_df_exp_random2 %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_exp_random2 <- sequence_2_df_exp_random2 %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X50di)

dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)

dat_chrI = cbind(dat_chrI, sequence_1_factor_chrI, sequence_2_factor_chrI)
# dat_chrI = cbind(dat_chrI, sequence_1_factor_chrI)
# dat_chrV = dat_chrV[2:(nrow(dat_chrV)-1),]
dat_chrV = cbind(dat_chrV, sequence_1_factor_chrV, sequence_2_factor_chrV)
# dat_chrV = cbind(dat_chrV, sequence_1_factor_chrV)
dat_chrVI = cbind(dat_chrVI, sequence_1_factor_chrVI, sequence_2_factor_chrVI)

dat_exp_random1 = cbind(dat_exp_random1, sequence_1_factor_exp_random1, sequence_2_factor_exp_random1)
dat_exp_random2 = cbind(dat_exp_random2, sequence_1_factor_exp_random2, sequence_2_factor_exp_random2)


# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.6, 0.65)
tiling_ylims_1000s = c(0.55, 0.775)

# tiling_AT_smoothC0_plot_list = construct_periodicity_plots_AT(
#   dat_tiling, smooth_10.4bp_cn_mean2, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
#   include_title=FALSE)
tiling_AT_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_tiling, smoothC0, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_smoothC0_plot_list[[1]]
tiling_AT_smoothC0_plot_list[[2]]
tiling_AT_smoothC0_plot_list[[3]]


# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.65)
chrv_ylims_1000s = c(0.525, 0.85)

chrv_AT_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_chrv, smooth_10.4bp_cn_mean2, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_smoothC0_plot_list[[1]]
chrv_AT_smoothC0_plot_list[[2]]
chrv_AT_smoothC0_plot_list[[3]]


# ChrI (A/T):
chrI_ylims_quartiles = c(0.575, 0.65)
chrI_ylims_1000s = c(0.525, 0.85)

chrI_AT_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_chrI, smoothC0, ylims_quartiles=chrI_ylims_quartiles, ylims_1000s=chrI_ylims_1000s,
  include_title=FALSE)

chrI_AT_smoothC0_plot_list[[1]]
chrI_AT_smoothC0_plot_list[[2]]
chrI_AT_smoothC0_plot_list[[3]]


# ChrV (A/T):
chrV_ylims_quartiles = c(0.575, 0.675)
chrV_ylims_1000s = c(0.525, 0.85)

chrV_AT_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_chrV, smoothC0, ylims_quartiles=chrV_ylims_quartiles, ylims_1000s=chrV_ylims_1000s,
  include_title=FALSE)

chrV_AT_smoothC0_plot_list[[1]]
chrV_AT_smoothC0_plot_list[[2]]
chrV_AT_smoothC0_plot_list[[3]]



# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

# tiling_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
#   dat_tiling, smooth_10.4bp_cn_mean2, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
#   include_title=FALSE)
tiling_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, smoothC0, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

tiling_AATT_smoothC0_plot_list[[1]]
tiling_AATT_smoothC0_plot_list[[2]]
tiling_AATT_smoothC0_plot_list[[3]]



# ChrI (AA/AT/TA/TT):
chrI_ylims_di_quartiles = c(0.35, 0.425)
chrI_ylims_di_1000s = c(0.3, 0.6)

chrI_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrI, smoothC0, ylims_quartiles=chrI_ylims_di_quartiles, ylims_1000s=chrI_ylims_di_1000s,
  include_title=FALSE)

chrI_AATT_smoothC0_plot_list[[1]]
chrI_AATT_smoothC0_plot_list[[2]]
chrI_AATT_smoothC0_plot_list[[3]]


# ChrV (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)
# 
# chrV_ylims_di_quartiles = c(0.1, 0.15)
# chrV_ylims_di_1000s = c(0.3, 0.7)
# 
# chrV_ylims_di_quartiles = c(0.2, 0.3)
# 
# chrV_ylims_di_quartiles = c(0.075, 0.1)
# 
# chrV_ylims_di_quartiles = c(0.05, 0.1)

chrV_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, smoothC0, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_smoothC0_plot_list[[1]]
chrV_AATT_smoothC0_plot_list[[2]]
chrV_AATT_smoothC0_plot_list[[3]]


# ChrVI (AA/AT/TA/TT):
chrVI_ylims_di_quartiles = c(0.35, 0.425)
chrVI_ylims_di_1000s = c(0.3, 0.6)

chrVI_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrVI, smoothC0, ylims_quartiles=chrVI_ylims_di_quartiles, ylims_1000s=chrVI_ylims_di_1000s,
  include_title=FALSE)

chrVI_AATT_smoothC0_plot_list[[1]]
chrVI_AATT_smoothC0_plot_list[[2]]
chrVI_AATT_smoothC0_plot_list[[3]]


# Expanded Random Library 1 (AA/AT/TA/TT):
exp_random_ylims_di_quartiles = c(0.2, 0.3)
exp_random_ylims_di_1000s = c(0.325, 0.65)

exp_random1_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_exp_random1, smoothC0, ylims_quartiles=exp_random_ylims_di_quartiles, ylims_1000s=exp_random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

exp_random1_AATT_smoothC0_plot_list[[1]]
exp_random1_AATT_smoothC0_plot_list[[2]]
exp_random1_AATT_smoothC0_plot_list[[3]]


# Expanded Random Library 2 (AA/AT/TA/TT):
exp_random2_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_exp_random2, smoothC0, ylims_quartiles=exp_random_ylims_di_quartiles, ylims_1000s=exp_random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

exp_random2_AATT_smoothC0_plot_list[[1]]
exp_random2_AATT_smoothC0_plot_list[[2]]
exp_random2_AATT_smoothC0_plot_list[[3]]


# Contracted Random Library 1 (AA/AT/TA/TT):
dat_contracted_random1 = dat_exp_random1[which(!is.na(dat_exp_random1$C0)),]

random_ylims_di_quartiles = c(0.2, 0.3)

contracted_random1_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_contracted_random1, smoothC0, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=exp_random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

contracted_random1_AATT_smoothC0_plot_list[[1]]
contracted_random1_AATT_smoothC0_plot_list[[2]]
contracted_random1_AATT_smoothC0_plot_list[[3]]


# Contracted Random Library 2 (AA/AT/TA/TT):
dat_contracted_random2 = dat_exp_random2[which(!is.na(dat_exp_random2$C0)),]

contracted_random2_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_contracted_random2, smoothC0, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=exp_random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

contracted_random2_AATT_smoothC0_plot_list[[1]]
contracted_random2_AATT_smoothC0_plot_list[[2]]
contracted_random2_AATT_smoothC0_plot_list[[3]]

