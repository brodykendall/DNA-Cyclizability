library(tidyverse)
library(gridExtra)

source("scripts/functions/sequence-functions.R")
source("scripts/functions/plotting-functions.R")


dat_nuc = read.csv("cycle1.txt")
dat_random = read.csv("cycle3.txt")
dat_tiling = read.csv("cycle5.txt")
dat_chrv = read.csv("cycle6.txt")

colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")


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


dat_nuc = cbind(dat_nuc, sequence_1_factor_nuc, sequence_2_factor_nuc)
dat_random = cbind(dat_random, sequence_1_factor_random, sequence_2_factor_random)
dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)


# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_original_C26_plot_list = construct_periodicity_plots_AT(
  dat_nuc, C26, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_original_C26_plot_list[[1]]
nuc_AT_original_C26_plot_list[[3]]


# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_original_C26_plot_list = construct_periodicity_plots_AT(
  dat_random, C26, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_original_C26_plot_list[[1]]
random_AT_original_C26_plot_list[[3]]



# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.6, 0.65)
tiling_ylims_1000s = c(0.55, 0.775)

tiling_AT_original_C26_plot_list = construct_periodicity_plots_AT(
  dat_tiling, C26, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_C26_plot_list[[1]]
tiling_AT_original_C26_plot_list[[3]]


tiling_AT_original_C26_plot_list = construct_periodicity_plots_AT(
  dat_tiling, C26, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_C26_plot_list[[1]]
tiling_AT_original_C26_plot_list[[3]]




# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.65)
chrv_ylims_1000s = c(0.525, 0.675)

chrv_AT_original_C26_plot_list = construct_periodicity_plots_AT(
  dat_chrv, C26, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AT_original_C26_plot_list[[1]]
chrv_AT_original_C26_plot_list[[3]]



grid.arrange(nuc_AT_original_C26_plot_list[[1]], random_AT_original_C26_plot_list[[1]],
             tiling_AT_original_C26_plot_list[[1]], chrv_AT_original_C26_plot_list[[1]])



# Nucleosome Library (AA/AT/TA/TT):
nuc_ylims_di_quartiles = c(0.225, 0.625)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_original_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, C26, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

nuc_AATT_original_C26_plot_list[[1]]

nuc_AATT_original_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, C29, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

nuc_AATT_original_C29_plot_list[[1]]

nuc_AATT_original_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, C31, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

nuc_AATT_original_C31_plot_list[[1]]

nuc_AATT_original_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, C0, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

nuc_AATT_original_C0_plot_list[[1]]



# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.15, 0.375)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_original_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, C26, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

random_AATT_original_C26_plot_list[[1]]

random_AATT_original_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, C29, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

random_AATT_original_C29_plot_list[[1]]

random_AATT_original_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, C31, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

random_AATT_original_C31_plot_list[[1]]

random_AATT_original_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, C0, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

random_AATT_original_C0_plot_list[[1]]


# Tiling Library (AA/AT/TA/TT):

tiling_ylims_di_1000s = c(0.325, 0.6)
tiling_ylims_di_quartiles = c(0.275, 0.5)

tiling_AATT_original_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, C26, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

tiling_AATT_original_C26_plot_list[[1]]

tiling_AATT_original_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, C29, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

tiling_AATT_original_C29_plot_list[[1]]

tiling_AATT_original_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, C31, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

tiling_AATT_original_C31_plot_list[[1]]

tiling_ylims_di_quartiles = c(0.35, 0.425)

tiling_AATT_original_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, C0, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

tiling_AATT_original_C0_plot_list[[1]]





# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.275, 0.475)
chrv_ylims_di_1000s = c(0.3, 0.45)

chrv_AATT_original_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, C26, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_C26_plot_list[[1]]

chrv_AATT_original_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, C29, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_C29_plot_list[[1]]

chrv_AATT_original_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, C31, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_C31_plot_list[[1]]

chrv_ylims_di_quartiles = c(0.35, 0.41)

chrv_AATT_original_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, C0, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_C0_plot_list[[1]]







grid.arrange(nuc_AATT_original_C26_plot_list[[1]], random_AATT_original_C26_plot_list[[1]],
             tiling_AATT_original_C26_plot_list[[1]], chrv_AATT_original_C26_plot_list[[1]])

