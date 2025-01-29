library(tidyverse)
library(gridExtra)

source("scripts/functions/sequence-functions.R")
source("scripts/functions/plotting-functions.R")


dat_nuc = read.csv("cycle1.txt")
dat_random = read.csv("cycle3.txt")
dat_tiling = read.csv("cycle5.txt")
dat_chrv = read.csv("cycle6.txt")

dat_chrV = read.csv("data/Created/yeast_chrV_1bpresolution_subsequence50.csv")

colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

colnames(dat_chrV) = c("x50mer", "chrID", "position")

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


dat_nuc = cbind(dat_nuc, sequence_1_factor_nuc, sequence_2_factor_nuc)
dat_random = cbind(dat_random, sequence_1_factor_random, sequence_2_factor_random)
dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)
# dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv)

dat_chrV = cbind(dat_chrV, sequence_1_factor_chrV, sequence_2_factor_chrV)


############## DNAcycP:
DNAcycP_pred_nuc = read.csv("data/predictions/cycle1_cycle_norm_DNAcycP.txt", 
                            header=FALSE, col.names="DNAcycP_pred")
DNAcycP_pred_random = read.csv("data/predictions/cycle3_cycle_norm_DNAcycP.txt",
                               header=FALSE, col.names="DNAcycP_pred")
DNAcycP_pred_tiling = read.csv("data/predictions/cycle5_cycle_norm_DNAcycP.txt",
                               header=FALSE, col.names="DNAcycP_pred")
DNAcycP_pred_chrv = read.csv("data/predictions/cycle6_cycle_norm_DNAcycP.txt",
                             header=FALSE, col.names="DNAcycP_pred")

DNAcycP_pred_chrV = read.csv("data/predictions/chrv_1bp_cycle_norm_DNAcycP.txt",
                             header=FALSE, col.names="DNAcycP_pred")

dat_nuc$DNAcycP_pred = unlist(DNAcycP_pred_nuc)
dat_random$DNAcycP_pred = unlist(DNAcycP_pred_random)
dat_tiling$DNAcycP_pred = unlist(DNAcycP_pred_tiling)
dat_chrv$DNAcycP_pred = unlist(DNAcycP_pred_chrv)

dat_chrV$DNAcycP_pred = unlist(DNAcycP_pred_chrV)

# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AT(
  dat_nuc, DNAcycP_pred, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_original_DNAcycP_pred_plot_list[[1]]
nuc_AT_original_DNAcycP_pred_plot_list[[3]]

# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AT(
  dat_random, DNAcycP_pred, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_original_DNAcycP_pred_plot_list[[1]]
random_AT_original_DNAcycP_pred_plot_list[[3]]

# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.6, 0.65)
tiling_ylims_1000s = c(0.55, 0.8)

tiling_AT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AT(
  dat_tiling, DNAcycP_pred, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_DNAcycP_pred_plot_list[[1]]
tiling_AT_original_DNAcycP_pred_plot_list[[3]]

# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.675)
chrv_ylims_1000s = c(0.55, 0.8)

chrv_AT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AT(
  dat_chrv, DNAcycP_pred, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_original_DNAcycP_pred_plot_list[[1]]
chrv_AT_original_DNAcycP_pred_plot_list[[2]]
chrv_AT_original_DNAcycP_pred_plot_list[[3]]


grid.arrange(nuc_AT_original_DNAcycP_pred_plot_list[[1]], random_AT_original_DNAcycP_pred_plot_list[[1]],
             tiling_AT_original_DNAcycP_pred_plot_list[[1]], chrv_AT_original_DNAcycP_pred_plot_list[[1]])
grid.arrange(tiling_AT_original_DNAcycP_pred_plot_list[[1]], chrv_AT_original_DNAcycP_pred_plot_list[[1]])


# Nucleosome Library (AA/AT/TA/TT):
nuc_ylims_di_quartiles = c(0.225, 0.575)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, DNAcycP_pred, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE)

nuc_AATT_original_DNAcycP_pred_plot_list[[1]]
nuc_AATT_original_DNAcycP_pred_plot_list[[3]]

# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.2, 0.3)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, DNAcycP_pred, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE)

random_AATT_original_DNAcycP_pred_plot_list[[1]]
random_AATT_original_DNAcycP_pred_plot_list[[3]]

# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

tiling_AATT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, DNAcycP_pred, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE)

tiling_AATT_original_DNAcycP_pred_plot_list[[1]]
tiling_AATT_original_DNAcycP_pred_plot_list[[3]]

# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.35, 0.425)
chrv_ylims_di_1000s = c(0.3, 0.6)

chrv_AATT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, DNAcycP_pred, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_DNAcycP_pred_plot_list[[1]]
chrv_AATT_original_DNAcycP_pred_plot_list[[3]]

grid.arrange(nuc_AATT_original_DNAcycP_pred_plot_list[[1]], random_AATT_original_DNAcycP_pred_plot_list[[1]],
             tiling_AATT_original_DNAcycP_pred_plot_list[[1]], chrv_AATT_original_DNAcycP_pred_plot_list[[1]])

grid.arrange(tiling_AATT_original_DNAcycP_pred_plot_list[[1]], chrv_AATT_original_DNAcycP_pred_plot_list[[1]])


# ChrV 1bp (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)

chrV_AATT_original_DNAcycP_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, DNAcycP_pred, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_original_DNAcycP_pred_plot_list[[1]]
chrV_AATT_original_DNAcycP_pred_plot_list[[3]]




############## CycPred:
CycPred_pred_nuc = read.csv("data/predictions/cycpred_cycle1.txt", 
                            header=FALSE, col.names="CycPred_pred")
CycPred_pred_random = read.csv("data/predictions/cycpred_cycle3.txt",
                               header=FALSE, col.names="CycPred_pred")
CycPred_pred_tiling = read.csv("data/predictions/cycpred_cycle5.txt",
                               header=FALSE, col.names="CycPred_pred")
CycPred_pred_chrv = read.csv("data/predictions/cycpred_cycle6.txt",
                             header=FALSE, col.names="CycPred_pred")

CycPred_pred_chrV = read.csv("data/predictions/cycpred_chrv_1bp.txt",
                             header=FALSE, col.names="CycPred_pred")

dat_nuc$CycPred_pred = unlist(CycPred_pred_nuc)
dat_random$CycPred_pred = unlist(CycPred_pred_random)
dat_tiling$CycPred_pred = unlist(CycPred_pred_tiling)
dat_chrv$CycPred_pred = unlist(CycPred_pred_chrv)

dat_chrV$CycPred_pred = unlist(CycPred_pred_chrV)

# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_original_CycPred_pred_plot_list = construct_periodicity_plots_AT(
  dat_nuc, CycPred_pred, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_original_CycPred_pred_plot_list[[1]]
nuc_AT_original_CycPred_pred_plot_list[[3]]

# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_original_CycPred_pred_plot_list = construct_periodicity_plots_AT(
  dat_random, CycPred_pred, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_original_CycPred_pred_plot_list[[1]]
random_AT_original_CycPred_pred_plot_list[[3]]

# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.575, 0.65)
tiling_ylims_1000s = c(0.55, 0.8)

tiling_AT_original_CycPred_pred_plot_list = construct_periodicity_plots_AT(
  dat_tiling, CycPred_pred, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_CycPred_pred_plot_list[[1]]
tiling_AT_original_CycPred_pred_plot_list[[3]]

# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.675)
chrv_ylims_1000s = c(0.55, 0.8)

chrv_AT_original_CycPred_pred_plot_list = construct_periodicity_plots_AT(
  dat_chrv, CycPred_pred, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_original_CycPred_pred_plot_list[[1]]
chrv_AT_original_CycPred_pred_plot_list[[3]]


grid.arrange(nuc_AT_original_CycPred_pred_plot_list[[1]], random_AT_original_CycPred_pred_plot_list[[1]],
             tiling_AT_original_CycPred_pred_plot_list[[1]], chrv_AT_original_CycPred_pred_plot_list[[1]])
grid.arrange(tiling_AT_original_CycPred_pred_plot_list[[1]], chrv_AT_original_CycPred_pred_plot_list[[1]])


# Nucleosome Library (AA/AT/TA/TT):
nuc_ylims_di_quartiles = c(0.225, 0.575)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_original_CycPred_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, CycPred_pred, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE)

nuc_AATT_original_CycPred_pred_plot_list[[1]]
nuc_AATT_original_CycPred_pred_plot_list[[3]]

# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.2, 0.3)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_original_CycPred_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, CycPred_pred, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE)

random_AATT_original_CycPred_pred_plot_list[[1]]
random_AATT_original_CycPred_pred_plot_list[[3]]

# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

tiling_AATT_original_CycPred_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, CycPred_pred, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE)

tiling_AATT_original_CycPred_pred_plot_list[[1]]
tiling_AATT_original_CycPred_pred_plot_list[[3]]

# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.35, 0.425)
chrv_ylims_di_1000s = c(0.3, 0.6)

chrv_AATT_original_CycPred_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, CycPred_pred, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE)

chrv_AATT_original_CycPred_pred_plot_list[[1]]
chrv_AATT_original_CycPred_pred_plot_list[[3]]

grid.arrange(nuc_AATT_original_CycPred_pred_plot_list[[1]], random_AATT_original_CycPred_pred_plot_list[[1]],
             tiling_AATT_original_CycPred_pred_plot_list[[1]], chrv_AATT_original_CycPred_pred_plot_list[[1]])

grid.arrange(tiling_AATT_original_CycPred_pred_plot_list[[1]], chrv_AATT_original_CycPred_pred_plot_list[[1]])

# ChrV 1bp (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)

chrV_AATT_original_CycPred_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, CycPred_pred, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_original_CycPred_pred_plot_list[[1]]
chrV_AATT_original_CycPred_pred_plot_list[[3]]






############## DeepBend:
DeepBend_pred_nuc = read.csv("data/predictions/deepbend_cycle1.txt", 
                            header=FALSE, col.names="DeepBend_pred")
DeepBend_pred_random = read.csv("data/predictions/deepbend_cycle3.txt",
                               header=FALSE, col.names="DeepBend_pred")
DeepBend_pred_tiling = read.csv("data/predictions/deepbend_cycle5.txt",
                               header=FALSE, col.names="DeepBend_pred")
DeepBend_pred_chrv = read.csv("data/predictions/deepbend_cycle6.txt",
                             header=FALSE, col.names="DeepBend_pred")

DeepBend_pred_chrV = read.csv("data/predictions/deepbend_chrv_1bp.txt",
                             header=FALSE, col.names="DeepBend_pred")

dat_nuc$DeepBend_pred = unlist(DeepBend_pred_nuc)
dat_random$DeepBend_pred = unlist(DeepBend_pred_random)
dat_tiling$DeepBend_pred = unlist(DeepBend_pred_tiling)
dat_chrv$DeepBend_pred = unlist(DeepBend_pred_chrv)

dat_chrV$DeepBend_pred = unlist(DeepBend_pred_chrV)

# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_nuc, DeepBend_pred, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_original_DeepBend_pred_plot_list[[1]]
nuc_AT_original_DeepBend_pred_plot_list[[3]]

# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_random, DeepBend_pred, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_original_DeepBend_pred_plot_list[[1]]
random_AT_original_DeepBend_pred_plot_list[[3]]

# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.575, 0.65)
tiling_ylims_1000s = c(0.55, 0.8)

tiling_AT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_tiling, DeepBend_pred, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_DeepBend_pred_plot_list[[1]]
tiling_AT_original_DeepBend_pred_plot_list[[3]]

# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.675)
chrv_ylims_1000s = c(0.55, 0.8)

chrv_AT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_chrv, DeepBend_pred, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_original_DeepBend_pred_plot_list[[1]]
chrv_AT_original_DeepBend_pred_plot_list[[3]]


grid.arrange(nuc_AT_original_DeepBend_pred_plot_list[[1]], random_AT_original_DeepBend_pred_plot_list[[1]],
             tiling_AT_original_DeepBend_pred_plot_list[[1]], chrv_AT_original_DeepBend_pred_plot_list[[1]])
grid.arrange(tiling_AT_original_DeepBend_pred_plot_list[[1]], chrv_AT_original_DeepBend_pred_plot_list[[1]])


# Nucleosome Library (AA/AT/TA/TT):
nuc_ylims_di_quartiles = c(0.225, 0.575)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, DeepBend_pred, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE)

nuc_AATT_original_DeepBend_pred_plot_list[[1]]
nuc_AATT_original_DeepBend_pred_plot_list[[3]]

# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.2, 0.3)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, DeepBend_pred, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE)

random_AATT_original_DeepBend_pred_plot_list[[1]]
random_AATT_original_DeepBend_pred_plot_list[[3]]

# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

tiling_AATT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, DeepBend_pred, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE)

tiling_AATT_original_DeepBend_pred_plot_list[[1]]
tiling_AATT_original_DeepBend_pred_plot_list[[3]]

# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.35, 0.425)
chrv_ylims_di_1000s = c(0.3, 0.6)

chrv_AATT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, DeepBend_pred, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrv_AATT_original_DeepBend_pred_plot_list[[1]]
chrv_AATT_original_DeepBend_pred_plot_list[[3]]

grid.arrange(nuc_AATT_original_DeepBend_pred_plot_list[[1]], random_AATT_original_DeepBend_pred_plot_list[[1]],
             tiling_AATT_original_DeepBend_pred_plot_list[[1]], chrv_AATT_original_DeepBend_pred_plot_list[[1]])
grid.arrange(tiling_AATT_original_DeepBend_pred_plot_list[[1]], chrv_AATT_original_DeepBend_pred_plot_list[[1]])

# ChrV 1bp (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)

chrV_AATT_original_DeepBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, DeepBend_pred, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_original_DeepBend_pred_plot_list[[1]]
chrV_AATT_original_DeepBend_pred_plot_list[[3]]







############## MixBend:
MixBend_pred_nuc = read.csv("data/predictions/mixbend_cycle1.txt", sep="\t")
MixBend_pred_random = read.csv("data/predictions/mixbend_cycle3.txt", sep="\t")
MixBend_pred_tiling = read.csv("data/predictions/mixbend_cycle5.txt", sep="\t")
MixBend_pred_chrv = read.csv("data/predictions/mixbend_cycle6.txt", sep="\t")

MixBend_pred_chrV = read.csv("data/predictions/mixbend_chrv_1bp.txt",
                             header=TRUE, col.names="MixBend_pred")

dat_nuc$MixBend_pred = MixBend_pred_nuc$pred_c0
dat_random$MixBend_pred = MixBend_pred_random$pred_c0
dat_tiling$MixBend_pred = MixBend_pred_tiling$pred_c0
dat_chrv$MixBend_pred = MixBend_pred_chrv$pred_c0

dat_chrV$MixBend_pred = unlist(MixBend_pred_chrV)

# Nucleosome Library (A/T):
nuc_ylims_quartiles = c(0.5, 0.8)
nuc_ylims_1000s = c(0.45, 0.85)

nuc_AT_original_MixBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_nuc, MixBend_pred, ylims_quartiles=nuc_ylims_quartiles, ylims_1000s=nuc_ylims_1000s,
  include_title=FALSE)

nuc_AT_original_MixBend_pred_plot_list[[1]]
nuc_AT_original_MixBend_pred_plot_list[[3]]

# Random Library (A/T):
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.45, 0.575)

random_AT_original_MixBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_random, MixBend_pred, ylims_quartiles=random_ylims_quartiles, ylims_1000s=random_ylims_1000s,
  include_title=FALSE)

random_AT_original_MixBend_pred_plot_list[[1]]
random_AT_original_MixBend_pred_plot_list[[3]]

# Tiling Library (A/T):
tiling_ylims_quartiles = c(0.6, 0.65)
tiling_ylims_1000s = c(0.55, 0.8)

tiling_AT_original_MixBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_tiling, MixBend_pred, ylims_quartiles=tiling_ylims_quartiles, ylims_1000s=tiling_ylims_1000s,
  include_title=FALSE)

tiling_AT_original_MixBend_pred_plot_list[[1]]
tiling_AT_original_MixBend_pred_plot_list[[3]]

# ChrV Library (A/T):
chrv_ylims_quartiles = c(0.575, 0.675)
chrv_ylims_1000s = c(0.55, 0.8)

chrv_AT_original_MixBend_pred_plot_list = construct_periodicity_plots_AT(
  dat_chrv, MixBend_pred, ylims_quartiles=chrv_ylims_quartiles, ylims_1000s=chrv_ylims_1000s,
  include_title=FALSE)

chrv_AT_original_MixBend_pred_plot_list[[1]]
chrv_AT_original_MixBend_pred_plot_list[[3]]


grid.arrange(nuc_AT_original_MixBend_pred_plot_list[[1]], random_AT_original_MixBend_pred_plot_list[[1]],
             tiling_AT_original_MixBend_pred_plot_list[[1]], chrv_AT_original_MixBend_pred_plot_list[[1]])
grid.arrange(tiling_AT_original_MixBend_pred_plot_list[[1]], chrv_AT_original_MixBend_pred_plot_list[[1]])


# Nucleosome Library (AA/AT/TA/TT):
nuc_ylims_di_quartiles = c(0.225, 0.575)
nuc_ylims_di_1000s = c(0.2, 0.625)

nuc_AATT_original_MixBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_nuc, MixBend_pred, ylims_quartiles=nuc_ylims_di_quartiles, ylims_1000s=nuc_ylims_di_1000s,
  include_title=FALSE)

nuc_AATT_original_MixBend_pred_plot_list[[1]]
nuc_AATT_original_MixBend_pred_plot_list[[3]]

# Random Library (AA/AT/TA/TT):
random_ylims_di_quartiles = c(0.2, 0.3)
random_ylims_di_1000s = c(0.2, 0.325)

random_AATT_original_MixBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_random, MixBend_pred, ylims_quartiles=random_ylims_di_quartiles, ylims_1000s=random_ylims_di_1000s,
  include_title=FALSE)

random_AATT_original_MixBend_pred_plot_list[[1]]
random_AATT_original_MixBend_pred_plot_list[[3]]

# Tiling Library (AA/AT/TA/TT):
tiling_ylims_di_quartiles = c(0.35, 0.425)
tiling_ylims_di_1000s = c(0.325, 0.65)

tiling_AATT_original_MixBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_tiling, MixBend_pred, ylims_quartiles=tiling_ylims_di_quartiles, ylims_1000s=tiling_ylims_di_1000s,
  include_title=FALSE)

tiling_AATT_original_MixBend_pred_plot_list[[1]]
tiling_AATT_original_MixBend_pred_plot_list[[3]]

# ChrV Library (AA/AT/TA/TT):
chrv_ylims_di_quartiles = c(0.35, 0.425)
chrv_ylims_di_1000s = c(0.3, 0.6)

chrv_AATT_original_MixBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrv, MixBend_pred, ylims_quartiles=chrv_ylims_di_quartiles, ylims_1000s=chrv_ylims_di_1000s,
  include_title=FALSE)

chrv_AATT_original_MixBend_pred_plot_list[[1]]
chrv_AATT_original_MixBend_pred_plot_list[[3]]

grid.arrange(nuc_AATT_original_MixBend_pred_plot_list[[1]], random_AATT_original_MixBend_pred_plot_list[[1]],
             tiling_AATT_original_MixBend_pred_plot_list[[1]], chrv_AATT_original_MixBend_pred_plot_list[[1]])

grid.arrange(tiling_AATT_original_MixBend_pred_plot_list[[1]], chrv_AATT_original_MixBend_pred_plot_list[[1]])

# ChrV 1bp (AA/AT/TA/TT):
chrV_ylims_di_quartiles = c(0.35, 0.425)
chrV_ylims_di_1000s = c(0.3, 0.7)

chrV_AATT_original_MixBend_pred_plot_list = construct_periodicity_plots_AA_TT_AT_TA(
  dat_chrV, MixBend_pred, ylims_quartiles=chrV_ylims_di_quartiles, ylims_1000s=chrV_ylims_di_1000s,
  include_title=FALSE, include_labels=FALSE, include_legend=FALSE)

chrV_AATT_original_MixBend_pred_plot_list[[1]]
chrV_AATT_original_MixBend_pred_plot_list[[3]]


