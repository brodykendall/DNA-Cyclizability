library(tidyverse)
library(gridExtra)

source("scripts/functions/sequence-functions.R")
source("scripts/functions/plotting-functions.R")

dat_nuc = read.csv("data/predictions/ir_lstm_smoothC0_10_11_tiling_cycle1_best_fold_predictions.csv")
dat_random = read.csv("data/predictions/ir_lstm_smoothC0_10_11_tiling_cycle3_best_fold_predictions.csv")
dat_tiling = read.csv("data/predictions/ir_lstm_smoothC0_10_11_tiling_cycle5_best_fold_predictions.csv")
dat_chrv = read.csv("data/predictions/ir_lstm_smoothC0_10_11_tiling_cycle6_best_fold_predictions.csv")

# dat_nuc = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle1_predictions.csv")
# dat_random = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle3_predictions.csv")
# dat_tiling = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle5_predictions.csv")
# dat_chrv = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle6_predictions.csv")

# dat_chrI = read.csv("data/predictions/ir_lstm_smoothC0_tiling_ir_lstm_cn_tiling_yeast_chrI_1bpresolution_subsequence50_smoothC0_predictions.csv")
dat_chrV = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11_best_fold_predictions.csv")
# dat_nuc = read.csv("data/predictions/ir_lstm_smoothC0_tiling_cycle1_predictions.csv")
# dat_random = read.csv("data/predictions/ir_lstm_smoothC0_tiling_cycle3_predictions.csv")
# dat_tiling = read.csv("data/predictions/ir_lstm_smoothC0_tiling_cycle5_predictions.csv")
# dat_chrv = read.csv("data/predictions/ir_lstm_smoothC0_tiling_cycle6_predictions.csv")

dat_chrI = dat_chrI %>%
  rename(x50mer = sequence,
         Predicted_SmoothC0 = smoothC0_predictions)
dat_chrV = dat_chrV %>%
  rename(x50mer = sequence,
         Predicted_SmoothC0 = smoothC0_predictions)

colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")


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


sequence_1_df_nuc <- dat_nuc %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_nuc <- dat_nuc %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_nuc), byrow = TRUE)} %>%
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


sequence_1_df_random <- dat_random %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_random <- dat_random %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_random), byrow = TRUE)} %>%
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


dat_chrI = cbind(dat_chrI, sequence_1_factor_chrI, sequence_2_factor_chrI)
dat_chrV = cbind(dat_chrV, sequence_1_factor_chrV, sequence_2_factor_chrV)

dat_nuc = cbind(dat_nuc, sequence_1_factor_nuc, sequence_2_factor_nuc)
# dat_nuc = cbind(dat_nuc, sequence_1_factor_nuc)
dat_random = cbind(dat_random, sequence_1_factor_random, sequence_2_factor_random)
# dat_random = cbind(dat_random, sequence_1_factor_random)
dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
# dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling)
dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)
# dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv)


# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_predicted_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_nuc, Predicted_SmoothC0, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_predicted_smoothC0_plot_list[[1]]
nuc_AT_predicted_smoothC0_plot_list[[2]]
nuc_AT_predicted_smoothC0_plot_list[[3]]


# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_predicted_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_random, Predicted_SmoothC0, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_predicted_smoothC0_plot_list[[1]]
random_AT_predicted_smoothC0_plot_list[[2]]
random_AT_predicted_smoothC0_plot_list[[3]]

# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.6, 0.65)
tiling_ylims_1000s = c(0.55, 0.775)

tiling_AT_predicted_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_tiling, Predicted_SmoothC0, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_predicted_smoothC0_plot_list[[1]]
tiling_AT_predicted_smoothC0_plot_list[[2]]
tiling_AT_predicted_smoothC0_plot_list[[3]]


# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.65)
chrv_ylims_1000s = c(0.525, 0.675)

chrv_AT_predicted_smoothC0_plot_list = construct_periodicity_plots_AT(
  dat_chrv, Predicted_SmoothC0, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_predicted_smoothC0_plot_list[[1]]
chrv_AT_predicted_smoothC0_plot_list[[2]]
chrv_AT_predicted_smoothC0_plot_list[[3]]





# ChrI (AA/AT/TA/TT):
chrI_ylims_di_quartiles = c(0.35, 0.425)
chrI_ylims_di_1000s = c(0.3, 0.6)

chrI_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrI, Predicted_SmoothC0, ylims_quartiles=chrI_ylims_di_quartiles, ylims_1000s=chrI_ylims_di_1000s,
  include_title=FALSE)

chrI_AATT_smoothC0_plot_list[[1]]
chrI_AATT_smoothC0_plot_list[[2]]
chrI_AATT_smoothC0_plot_list[[3]]


# ChrV (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)

chrV_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, Predicted_SmoothC0, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_smoothC0_plot_list[[1]]
chrV_AATT_smoothC0_plot_list[[2]]
chrV_AATT_smoothC0_plot_list[[3]]


# Nucleosome Library (AA/AT/TA/TT):
# nuc_ylims_di_quartiles = c(0.225, 0.575)
nuc_ylims_di_quartiles = c(0.225, 0.625)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, Predicted_SmoothC0, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

nuc_AATT_smoothC0_plot_list[[1]]
nuc_AATT_smoothC0_plot_list[[2]]
nuc_AATT_smoothC0_plot_list[[3]]


# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.2, 0.3)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, Predicted_SmoothC0, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE)

random_AATT_smoothC0_plot_list[[1]]
random_AATT_smoothC0_plot_list[[2]]
random_AATT_smoothC0_plot_list[[3]]


# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

tiling_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, Predicted_SmoothC0, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE)

tiling_AATT_smoothC0_plot_list[[1]]
tiling_AATT_smoothC0_plot_list[[2]]
tiling_AATT_smoothC0_plot_list[[3]]


# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.35, 0.425)
chrv_ylims_di_1000s = c(0.3, 0.6)

chrv_AATT_smoothC0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, Predicted_SmoothC0, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE)

chrv_AATT_smoothC0_plot_list[[1]]
chrv_AATT_smoothC0_plot_list[[2]]
chrv_AATT_smoothC0_plot_list[[3]]
