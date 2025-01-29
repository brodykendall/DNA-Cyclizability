library(tidyverse)

# dat_tiling = read.csv("data/Created/tiling_smoothC0.csv") %>% select(-X)
dat_tiling = read.csv("cycle5.txt")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

# dat_chrv = read.csv("data/Created/chrv_smoothC0.csv") %>% select(-X)
dat_chrv = read.csv("cycle6.txt")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

# 
dat_random <- read_csv("cycle3.txt")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
# 
dat_random_temp <- read.csv("data/Created/cycle3_predicted_smooth_C0.csv")
dat_tiling_temp <- read.csv("data/Created/cycle5_predicted_smooth_C0.csv")
dat_chrv_temp <- read.csv("data/Created/cycle6_predicted_smooth_C0.csv")
colnames(dat_random_temp) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "predicted_smooth_C0")
colnames(dat_tiling_temp) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "predicted_smooth_C0")
colnames(dat_chrv_temp) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "predicted_smooth_C0")
# 
# # Remove first sequence from each gene:
# dat_tiling_temp = dat_tiling_temp[-seq(1, nrow(dat_tiling_temp), by=143),]
# # Remove final sequence from each gene:
# dat_tiling_temp = dat_tiling_temp[-seq(142, nrow(dat_tiling_temp), by=142),]
# 
# # Remove first sequence:
# dat_chrv_temp = dat_chrv_temp[-1,]
# # Remove final sequence:
# dat_chrv_temp = dat_chrv_temp[-nrow(dat_chrv_temp),]
# 
dat_random$predicted_smooth_C0 = dat_random_temp$predicted_smooth_C0
dat_tiling$predicted_smooth_C0 = dat_tiling_temp$predicted_smooth_C0
dat_chrv$predicted_smooth_C0 = dat_chrv_temp$predicted_smooth_C0

# dat_chrv_2 <- read.csv("data/Created/yeast_chrV_ir_lstm_tiling_post_smoothed_matched.csv")
# # # Remove irrelevant columns:
# # dat_chrv_2 = dat_chrv_2[,c(1:19, 25:32)]
# dat_chrv_2 = dat_chrv_2[,c(1:19)]
# # # Remove first sequence:
# dat_chrv_2 = dat_chrv_2[-1,]
# # # Remove final sequence:
# dat_chrv_2 = dat_chrv_2[-nrow(dat_chrv_2),]
# 
# dat_chrv_2 = dat_chrv_2 %>%
#   rename(x50mer = Sequence,
#          C26 = n.26,
#          C29 = n.29,
#          C31 = n.31
#   )

dat_chrv_2 = read.csv("data/Created/cycle6_ir_lstm_cn_tiling_post_smoothed_matched.csv")

dat_chrv_2 = dat_chrv_2 %>%
  rename(x50mer = sequence,
         C26 = n.26,
         C29 = n.29,
         C31 = n.31
  )

dat_chrv_3 = read.csv("data/predictions/cycle6_ir_lstm_cn_tiling_post_smoothed_matched_predictions.csv")

dat_chrv_3 = dat_chrv_3 %>%
  select(smooth_10.4bp_C02_predictions, smooth_11bp_C02_predictions)

dat_chrv_2 = cbind(dat_chrv_2, dat_chrv_3)

# smooth_C0_v2_tiling_pred_tiling = read.csv("data/Created/ir_lstm_smoothC0_v2_tiling_pred_tiling.csv", header=FALSE)
# smooth_C0_v2_chrv_pred_tiling = read.csv("data/Created/ir_lstm_smoothC0_v2_chrv_pred_tiling.csv", header=FALSE)
# smooth_C0_v2_random_pred_tiling = read.csv("data/Created/ir_lstm_smoothC0_v2_random_pred_tiling.csv", header=FALSE)

# dat_tiling$smooth_C0_v2_pred_tiling = smooth_C0_v2_tiling_pred_tiling %>% unlist()
# dat_chrv$smooth_C0_v2_pred_tiling = smooth_C0_v2_chrv_pred_tiling %>% unlist()
# dat_random$smooth_C0_v2_pred_tiling = smooth_C0_v2_random_pred_tiling %>% unlist()

nucleotides <- c("A", "C", "G", "T")
dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")

source("scripts/functions/sequence-functions.R")

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

sequence_1_df_chrv_2 <- dat_chrv_2 %>%
  pull(x50mer) %>%
  sequence_df()
sequence_2_df_chrv_2 <- dat_chrv_2 %>%
  pull(x50mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat_chrv_2), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)
sequence_1_factor_chrv_2 <- sequence_1_df_chrv_2 %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))
sequence_2_factor_chrv_2 <- sequence_2_df_chrv_2 %>%
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

dat_tiling = cbind(dat_tiling, sequence_1_factor_tiling, sequence_2_factor_tiling)
dat_chrv = cbind(dat_chrv, sequence_1_factor_chrv, sequence_2_factor_chrv)
dat_chrv_2 = cbind(dat_chrv_2, sequence_1_factor_chrv_2, sequence_2_factor_chrv_2)
dat_random = cbind(dat_random, sequence_1_factor_random, sequence_2_factor_random)

construct_periodicity_plots_AT = function(data, var, ylims_quartiles=c(0.575, 0.65), 
                                          ylims_1000s=c(0.57, 0.77)) {
  
  var_string = deparse(substitute(var))
  data_string = deparse(substitute(data))
  var_of_interest = data %>% select({{var}})
  
  # Find cutoffs for first-fourth quartiles + custom
  cutoffs_0_25 = quantile(var_of_interest[,1], c(0, 0.25))
  cutoffs_25_50 = quantile(var_of_interest[,1], c(0.25, 0.5))
  cutoffs_50_75 = quantile(var_of_interest[,1], c(0.5, 0.75))
  cutoffs_75_100 = quantile(var_of_interest[,1], c(0.75, 1))
  cutoffs_custom1 = c(sort(var_of_interest[,1], TRUE)[1000], cutoffs_75_100[2])
  cutoffs_custom2 = c(cutoffs_0_25[1], sort(var_of_interest[,1], FALSE)[1000])
  
  # Divide each library into quartiles + custom
  var_q1 = data %>%
    filter({{var}} >= cutoffs_0_25[1] & {{var}} <= cutoffs_0_25[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q2 = data %>%
    filter({{var}} >= cutoffs_25_50[1] & {{var}} <= cutoffs_25_50[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q3 = data %>%
    filter({{var}} >= cutoffs_50_75[1] & {{var}} <= cutoffs_50_75[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q4 = data %>%
    filter({{var}} >= cutoffs_75_100[1] & {{var}} <= cutoffs_75_100[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom1 = data %>%
    filter({{var}} >= cutoffs_custom1[1] & {{var}} <= cutoffs_custom1[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom2 = data %>%
    filter({{var}} >= cutoffs_custom2[1] & {{var}} <= cutoffs_custom2[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  
  # Find the relative frequencies of A and T at each position (1-50) for each quartile + custom
  A_T_var_q1 = apply(var_q1 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q1))
  })
  A_T_var_q2 = apply(var_q2 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q2))
  })
  A_T_var_q3 = apply(var_q3 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q3))
  })
  A_T_var_q4 = apply(var_q4 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_q4))
  })
  A_T_var_custom1 = apply(var_custom1 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_custom1))
  })
  A_T_var_custom2 = apply(var_custom2 %>% select(all_of(ps1)), 2, function(col) {
    A_freq = sum(col == "A")
    T_freq = sum(col == "T")
    return((A_freq + T_freq)/nrow(var_custom2))
  })
  
  # Construct dataframe containing the relative frequences of A and T by position for all four quartiles + custom
  A_T_var_q = data.frame(q1 = A_T_var_q1,
                         q2 = A_T_var_q2,
                         q3 = A_T_var_q3,
                         q4 = A_T_var_q4, 
                         custom1 = A_T_var_custom1,
                         custom2 = A_T_var_custom2)
  
  # Plot the relative frequencies of A and T by position for each quartile (q1 v q4, 
  # q2 v q3, and custom)
  AT_q1vq4_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=q1, color="0-25%")) +
    geom_line(aes(y=q4, color="75-100%")) +
    scale_color_manual(breaks = c("0-25%", "75-100%"),
                       values = c("blue", "red"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_quartiles)
  
  AT_q1vq4_var
  
  AT_q2vq3_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=q2, color="25-50%")) +
    geom_line(aes(y=q3, color="50-75%")) +
    scale_color_manual(breaks = c("25-50%", "50-75%"),
                       values = c("green", "orange"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_quartiles)
  
  AT_custom1_var = ggplot(data = A_T_var_q, aes(x = 1:50)) +
    geom_line(aes(y=custom1, color="Top 1000")) +
    geom_line(aes(y=custom2, color="Bottom 1000")) +
    scale_color_manual(breaks = c("Top 1000", "Bottom 1000"),
                       values = c("firebrick1", "dodgerblue"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average A/T Content") +
    labs(title = paste0("Average A/T Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_1000s)
  
  plot_list = list(AT_q1vq4_var, AT_q2vq3_var, AT_custom1_var)
  
  return(plot_list)
}

library(gridExtra)

dat_random$actual_C0 = find_c0Aphi(dat_random%>%select(C26,C29,C31), aa=c(1,0.82,0.7))[,1]
dat_random_actual_C0_plot_list = construct_periodicity_plots_AT(dat_random, actual_C0, ylims_quartiles=c(0.4,0.6))

dat_random_actual_C0_plot_list[[1]]
# ChrV Predicted Cn
ylims_orig = c(.4, .85)

chrv_2_AT_predicted_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C26, ylims_1000s = ylims_orig)
chrv_2_AT_predicted_C29_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C29, ylims_1000s = ylims_orig)
chrv_2_AT_predicted_C31_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C31, ylims_1000s = ylims_orig)

chrv_2_AT_predicted_C26_plot_list[[3]]

chrv_2_AT_predicted_C29_plot_list[[3]]

chrv_2_AT_predicted_C31_plot_list[[3]]


# ChrV Raw Cn
chrv_2_AT_raw_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C26_raw, ylims_1000s = ylims_orig)
chrv_2_AT_raw_C29_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C29_raw, ylims_1000s = ylims_orig)
chrv_2_AT_raw_C31_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C31_raw, ylims_1000s = ylims_orig)

chrv_2_AT_raw_C26_plot_list[[3]]

chrv_2_AT_raw_C29_plot_list[[3]]

chrv_2_AT_raw_C31_plot_list[[3]]


# ChrV k=3
chrv_2_AT_smooth_3bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_3bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_3bp_C26_plot_list[[1]]


# ChrV k=5
chrv_2_AT_smooth_5bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_5bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_5bp_C26_plot_list[[1]]


# ChrV k=7
chrv_2_AT_smooth_7bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_7bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_7bp_C26_plot_list[[1]]


# ChrV k=9
chrv_2_AT_smooth_9bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_9bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_9bp_C26_plot_list[[1]]


# ChrV k=11
chrv_2_AT_smooth_11bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_11bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_11bp_C26_plot_list[[1]]

chrv_2_AT_smooth_11bp_C26_plot_list[[3]]


# ChrV k=10.7
chrv_2_AT_smooth_10.7bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.7bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_10.7bp_C26_plot_list[[1]]

chrv_2_AT_smooth_10.7bp_C26_plot_list[[3]]


# ChrV k=10.4
chrv_2_AT_smooth_10.4bp_C26_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_C26, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)
chrv_2_AT_smooth_10.4bp_C29_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_C29, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)
chrv_2_AT_smooth_10.4bp_C31_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_C31, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_10.4bp_cn_mean_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_cn_mean, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)
chrv_2_AT_smooth_10.4bp_C02_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_C02, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)
chrv_2_AT_smooth_10.4bp_C02_predictions_plot_list = construct_periodicity_plots_AT(dat_chrv_2, smooth_10.4bp_C02_predictions, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_smooth_10.4bp_C26_plot_list[[1]]
chrv_2_AT_smooth_10.4bp_C29_plot_list[[1]]
chrv_2_AT_smooth_10.4bp_C31_plot_list[[1]]

chrv_2_AT_smooth_10.4bp_cn_mean_plot_list[[1]]
chrv_2_AT_smooth_10.4bp_C02_plot_list[[1]]
chrv_2_AT_smooth_10.4bp_C02_predictions_plot_list[[1]]

chrv_2_AT_smooth_10.4bp_C26_plot_list[[3]]
chrv_2_AT_smooth_10.4bp_C29_plot_list[[3]]
chrv_2_AT_smooth_10.4bp_C31_plot_list[[3]]

chrv_2_AT_smooth_10.4bp_C02_plot_list[[3]]


# ChrV DNACycP
chrv_2_AT_DNACycP_plot_list = construct_periodicity_plots_AT(dat_chrv_2, DNAcycP_pred_chrV, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_DNACycP_plot_list[[1]]


# ChrV CycPred
chrv_2_AT_CycPred_plot_list = construct_periodicity_plots_AT(dat_chrv_2, CycPred_pred_chrV, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_CycPred_plot_list[[1]]


# ChrV DeepBend
chrv_2_AT_DeepBend_plot_list = construct_periodicity_plots_AT(dat_chrv_2, DeepBend_pred_chrV, ylims_quartiles = c(0.51, .71), ylims_1000s = ylims_orig)

chrv_2_AT_DeepBend_plot_list[[1]]








# Tiling Original C0:
# tiling_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_tiling, C0, ylims_quartiles = c(.51, .71))
tiling_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_tiling, C0)

tiling_AT_original_C0_plot_list[[1]]

tiling_AT_original_C0_plot_list[[3]]


# Tiling Predicted (From Cn) Post-Smoothed C0 Comparison:
tiling_AT_predicted_smooth_C0_plot_list = construct_periodicity_plots_AT(dat_tiling, predicted_smooth_C0)

tiling_AT_predicted_smooth_C0_plot_list[[1]]

grid.arrange(tiling_AT_original_C0_plot_list[[1]], tiling_AT_predicted_smooth_C0_plot_list[[1]])

tiling_AT_predicted_smooth_C0_plot_list[[3]]

grid.arrange(tiling_AT_original_C0_plot_list[[3]], tiling_AT_predicted_smooth_C0_plot_list[[3]])

tiling_smooth_C0_lm = lm(predicted_smooth_C0~C26+C29+C31, data=dat_tiling)
tiling_C0_lm = lm(C0~C26+C29+C31, data=dat_tiling)


# Tiling Smooth C0 Comparison:
tiling_AT_smooth_C0_v1_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C0_v1)
tiling_AT_smooth_C0_v2_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C0_v2)
tiling_AT_smooth_C0_v3_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C0_v3)
tiling_AT_C0_lag1_plot_list = construct_periodicity_plots_AT(dat_tiling, C0_lag1)

tiling_AT_smooth_C0_v1_plot_list[[1]]
tiling_AT_smooth_C0_v2_plot_list[[1]]
tiling_AT_smooth_C0_v3_plot_list[[1]]
tiling_AT_C0_lag1_plot_list[[1]]

grid.arrange(tiling_AT_original_C0_plot_list[[1]], tiling_AT_smooth_C0_v1_plot_list[[1]])
grid.arrange(tiling_AT_original_C0_plot_list[[1]], tiling_AT_smooth_C0_v2_plot_list[[1]])
grid.arrange(tiling_AT_original_C0_plot_list[[1]], tiling_AT_smooth_C0_v3_plot_list[[1]])
grid.arrange(tiling_AT_original_C0_plot_list[[1]], tiling_AT_C0_lag1_plot_list[[1]])

tiling_AT_smooth_C0_v1_plot_list[[3]]
tiling_AT_smooth_C0_v2_plot_list[[3]]
tiling_AT_smooth_C0_v3_plot_list[[3]]

grid.arrange(tiling_AT_original_C0_plot_list[[3]], tiling_AT_smooth_C0_v1_plot_list[[3]])
grid.arrange(tiling_AT_original_C0_plot_list[[3]], tiling_AT_smooth_C0_v2_plot_list[[3]])
grid.arrange(tiling_AT_original_C0_plot_list[[3]], tiling_AT_smooth_C0_v3_plot_list[[3]])

grid.arrange(tiling_AT_smooth_C0_v2_plot_list[[3]], tiling_AT_smooth_C0_v3_plot_list[[3]])



# ChrV Original C0:
chrv_ylims_1000s = c(0.55, 0.75)

# chrv_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, C0, ylims_quartiles = c(.51, .71), ylims_1000s = chrv_ylims_1000s)
chrv_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, C0, ylims_1000s = chrv_ylims_1000s)

chrv_2_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_chrv_2, C0, ylims_1000s = chrv_ylims_1000s)

chrv_AT_original_C0_plot_list[[1]]
chrv_2_AT_original_C0_plot_list[[1]]

chrv_AT_original_C0_plot_list[[3]]
chrv_2_AT_original_C0_plot_list[[3]]


# ChrV DNAcycP pred:
chrv_2_AT_DNAcycP_pred_plot_list = construct_periodicity_plots_AT(dat_chrv_2, DNAcycP_pred_chrV, ylims_1000s = chrv_ylims_1000s)

chrv_2_AT_DNAcycP_pred_plot_list[[1]]

grid.arrange(chrv_2_AT_original_C0_plot_list[[1]], chrv_2_AT_DNAcycP_pred_plot_list[[1]])


# ChrV DNAcycP Post-Smoothed C0 Comparison:
chrv_2_AT_DNAcycP_post_smooth_5bp_plot_list = construct_periodicity_plots_AT(dat_chrv_2, DNAcycP_post_smooth_5bp, ylims_1000s = chrv_ylims_1000s)

chrv_2_AT_DNAcycP_post_smooth_5bp_plot_list[[1]]

grid.arrange(chrv_2_AT_original_C0_plot_list[[1]], chrv_2_AT_DNAcycP_post_smooth_5bp_plot_list[[1]])

chrv_2_AT_DNAcycP_post_smooth_5bp_plot_list[[3]]

grid.arrange(chrv_2_AT_original_C0_plot_list[[3]], chrv_2_AT_DNAcycP_post_smooth_5bp_plot_list[[3]])


# ChrV Predicted (From Sequence) Post-Smoothed C0 Comparison:
chrv_AT_pred_tiling_post_smooth_plot_list = construct_periodicity_plots_AT(dat_chrv_2, pred_tiling_post_smooth, ylims_1000s = chrv_ylims_1000s)

chrv_AT_pred_tiling_post_smooth_plot_list[[1]]

grid.arrange(chrv_AT_original_C0_plot_list[[1]], chrv_AT_pred_tiling_post_smooth_plot_list[[1]])

chrv_AT_pred_tiling_post_smooth_plot_list[[3]]

grid.arrange(chrv_AT_original_C0_plot_list[[3]], chrv_AT_pred_tiling_post_smooth_plot_list[[3]])




# ChrV Post-Smoothed C0 Comparison:
chrv_2_AT_post_smooth_5bp_plot_list = construct_periodicity_plots_AT(dat_chrv_2, post_smooth_5bp, ylims_1000s = chrv_ylims_1000s)

chrv_2_AT_post_smooth_5bp_plot_list[[1]]

grid.arrange(chrv_2_AT_original_C0_plot_list[[1]], chrv_2_AT_post_smooth_5bp_plot_list[[1]])
grid.arrange(chrv_2_AT_post_smooth_5bp_plot_list[[1]], chrv_AT_pred_tiling_post_smooth_plot_list[[1]])

chrv_2_AT_post_smooth_5bp_plot_list[[3]]

grid.arrange(chrv_2_AT_original_C0_plot_list[[3]], chrv_2_AT_post_smooth_5bp_plot_list[[3]])
grid.arrange(chrv_2_AT_post_smooth_5bp_plot_list[[3]], chrv_AT_pred_tiling_post_smooth_plot_list[[3]])


# ChrV Predicted (From Cn) Post-Smoothed C0 Comparison:
chrv_AT_predicted_smooth_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, predicted_smooth_C0, ylims_1000s = chrv_ylims_1000s)
# chrv_AT_predicted_smooth_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, lm_predicted_smoothC0, ylims_1000s = chrv_ylims_1000s)

chrv_AT_predicted_smooth_C0_plot_list[[1]]

grid.arrange(chrv_AT_original_C0_plot_list[[1]], chrv_AT_predicted_smooth_C0_plot_list[[1]])

chrv_AT_predicted_smooth_C0_plot_list[[3]]

grid.arrange(chrv_AT_original_C0_plot_list[[3]], chrv_AT_predicted_smooth_C0_plot_list[[3]])

chrv_smooth_C0_lm = lm(predicted_smooth_C0~C26+C29+C31, data=dat_chrv)
cor(chrv_smooth_C0_lm$fitted.values, dat_chrv$predicted_smooth_C0)
cor(chrv_smooth_C0_lm$fitted.values, dat_chrv$C0)
dat_chrv$lm_predicted_smoothC0 = chrv_smooth_C0_lm$fitted.values


# ChrV Smooth C0 Comparison:
chrv_AT_smooth_C0_v1_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v1, ylims_1000s = chrv_ylims_1000s)
chrv_AT_smooth_C0_v2_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v2, ylims_1000s = chrv_ylims_1000s)
chrv_AT_smooth_C0_v3_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v3, ylims_1000s = chrv_ylims_1000s)

chrv_AT_smooth_C0_v1_plot_list[[1]]
chrv_AT_smooth_C0_v2_plot_list[[1]]
chrv_AT_smooth_C0_v3_plot_list[[1]]

grid.arrange(chrv_AT_original_C0_plot_list[[1]], chrv_AT_smooth_C0_v1_plot_list[[1]])
grid.arrange(chrv_AT_original_C0_plot_list[[1]], chrv_AT_smooth_C0_v2_plot_list[[1]])
grid.arrange(chrv_AT_original_C0_plot_list[[1]], chrv_AT_smooth_C0_v3_plot_list[[1]])

grid.arrange(chrv_AT_smooth_C0_v1_plot_list[[1]], chrv_2_AT_post_smooth_5bp_plot_list[[1]])
grid.arrange(chrv_AT_smooth_C0_v2_plot_list[[1]], chrv_2_AT_post_smooth_5bp_plot_list[[1]])

grid.arrange(chrv_AT_smooth_C0_v2_plot_list[[1]], chrv_AT_smooth_C0_v3_plot_list[[1]])

chrv_AT_smooth_C0_v1_plot_list[[3]]
chrv_AT_smooth_C0_v2_plot_list[[3]]
chrv_AT_smooth_C0_v3_plot_list[[3]]

grid.arrange(chrv_AT_original_C0_plot_list[[3]], chrv_AT_smooth_C0_v1_plot_list[[3]])
grid.arrange(chrv_AT_original_C0_plot_list[[3]], chrv_AT_smooth_C0_v2_plot_list[[3]])
grid.arrange(chrv_AT_original_C0_plot_list[[3]], chrv_AT_smooth_C0_v3_plot_list[[3]])

grid.arrange(chrv_AT_smooth_C0_v1_plot_list[[3]], chrv_2_AT_post_smooth_5bp_plot_list[[3]])
grid.arrange(chrv_AT_smooth_C0_v2_plot_list[[3]], chrv_2_AT_post_smooth_5bp_plot_list[[3]])

grid.arrange(chrv_AT_smooth_C0_v2_plot_list[[3]], chrv_AT_smooth_C0_v3_plot_list[[3]])



# ChrV C0 Comparison:
# chrv_AT_smooth_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v2, ylims_1000s = c(0.55, 0.71))
chrv_AT_smooth_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0, ylims_1000s = )
chrv_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_chrv, C0, ylims_1000s = c(0.55, 0.71))

chrv_AT_smooth_C0_plot_list[[1]]
chrv_AT_original_C0_plot_list[[1]]

grid.arrange(chrv_AT_smooth_C0_plot_list[[1]], chrv_AT_original_C0_plot_list[[1]])

chrv_AT_smooth_C0_plot_list[[3]]
chrv_AT_original_C0_plot_list[[3]]

grid.arrange(chrv_AT_smooth_C0_plot_list[[3]], chrv_AT_original_C0_plot_list[[3]])




# Random Original C0:
random_ylims_quartiles = c(0.45, 0.55)
random_ylims_1000s = c(0.4, 0.6)

random_AT_original_C0_plot_list = construct_periodicity_plots_AT(dat_random, C0, 
                                                                 ylims_quartiles = random_ylims_quartiles,
                                                                 ylims_1000s = random_ylims_1000s)

random_AT_original_C0_plot_list[[1]]

random_AT_original_C0_plot_list[[3]]

random_AT_original_C0_plot_list[[2]]





# Tiling smooth_C0_v2 prediction (trained on tiling) Comparison:
# tiling_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C0_v2_pred_tiling_v2)
tiling_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C0_v2_pred_tiling)

tiling_AT_smooth_C0_plot_list[[1]]
tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[1]]
tiling_AT_original_C0_plot_list[[1]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[1]], tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[1]])
grid.arrange(tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[1]], tiling_AT_original_C0_plot_list[[1]])


tiling_AT_smooth_C0_plot_list[[3]]
tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[3]]
tiling_AT_original_C0_plot_list[[3]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[3]], tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[3]])
grid.arrange(tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[3]], tiling_AT_original_C0_plot_list[[3]])

tiling_AT_smooth_C0_v2_pred_tiling_plot_list[[2]]



# ChrV smooth_C0_v2 prediction (trained on tiling) Comparison:
# chrv_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v2_pred_tiling_v2)
chrv_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C0_v2_pred_tiling)

chrv_AT_smooth_C0_plot_list[[1]]
chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[1]]
chrv_AT_original_C0_plot_list[[1]]

grid.arrange(chrv_AT_smooth_C0_plot_list[[1]], chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[1]])
grid.arrange(chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[1]], chrv_AT_original_C0_plot_list[[1]])


chrv_AT_smooth_C0_plot_list[[3]]
chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[3]]
chrv_AT_original_C0_plot_list[[3]]

grid.arrange(chrv_AT_smooth_C0_plot_list[[3]], chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[3]])
grid.arrange(chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[3]], chrv_AT_original_C0_plot_list[[3]])

chrv_AT_smooth_C0_v2_pred_tiling_plot_list[[2]]





# Random smooth_C0_v2 prediction (trained on tiling) Comparison:
# random_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_random, smooth_C0_v2_pred_tiling_v2)
random_AT_smooth_C0_v2_pred_tiling_plot_list = construct_periodicity_plots_AT(dat_random, smooth_C0_v2_pred_tiling, 
                                                                              ylims_quartiles = random_ylims_quartiles,
                                                                              ylims_1000s = random_ylims_1000s)

random_AT_smooth_C0_v2_pred_tiling_plot_list[[1]]
random_AT_original_C0_plot_list[[1]]

grid.arrange(random_AT_smooth_C0_v2_pred_tiling_plot_list[[1]], random_AT_original_C0_plot_list[[1]])


random_AT_smooth_C0_v2_pred_tiling_plot_list[[3]]
random_AT_original_C0_plot_list[[3]]

grid.arrange(random_AT_smooth_C0_v2_pred_tiling_plot_list[[3]], random_AT_original_C0_plot_list[[3]])





# Tiling grouped phi C0 Comparison:
tiling_AT_C0_grouped_phi_plot_list = construct_periodicity_plots_AT(grouped_phi_df_tiling, C0_grouped_phi)

tiling_AT_smooth_C0_plot_list[[1]]
tiling_AT_C0_grouped_phi_plot_list[[1]]
tiling_AT_original_C0_plot_list[[1]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[1]], tiling_AT_C0_grouped_phi_plot_list[[1]])
grid.arrange(tiling_AT_C0_grouped_phi_plot_list[[1]], tiling_AT_original_C0_plot_list[[1]])


tiling_AT_smooth_C0_plot_list[[3]]
tiling_AT_C0_grouped_phi_plot_list[[3]]
tiling_AT_original_C0_plot_list[[3]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[3]], tiling_AT_C0_grouped_phi_plot_list[[3]])
grid.arrange(tiling_AT_C0_grouped_phi_plot_list[[3]], tiling_AT_original_C0_plot_list[[3]])

tiling_AT_C0_grouped_phi_plot_list[[2]]



# Tiling grouped a C0 Comparison:
tiling_AT_C0_grouped_a_plot_list = construct_periodicity_plots_AT(grouped_a_df_tiling, C0_grouped_a)

tiling_AT_smooth_C0_plot_list[[1]]
tiling_AT_C0_grouped_a_plot_list[[1]]
tiling_AT_original_C0_plot_list[[1]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[1]], tiling_AT_C0_grouped_a_plot_list[[1]])
grid.arrange(tiling_AT_C0_grouped_a_plot_list[[1]], tiling_AT_original_C0_plot_list[[1]])
grid.arrange(tiling_AT_C0_grouped_a_plot_list[[1]], tiling_AT_C0_grouped_phi_plot_list[[1]])


tiling_AT_smooth_C0_plot_list[[3]]
tiling_AT_C0_grouped_a_plot_list[[3]]
tiling_AT_original_C0_plot_list[[3]]

grid.arrange(tiling_AT_smooth_C0_plot_list[[3]], tiling_AT_C0_grouped_a_plot_list[[3]])
grid.arrange(tiling_AT_C0_grouped_a_plot_list[[3]], tiling_AT_original_C0_plot_list[[3]])

tiling_AT_C0_grouped_a_plot_list[[2]]










Cn_ylims_quartiles = c(0.525, 0.7)
Cn_ylims_1000s = c(0.45, 0.825)

# Tiling C26 Comparison:
tiling_AT_smooth_C26_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C26, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)
tiling_AT_original_C26_plot_list = construct_periodicity_plots_AT(dat_tiling, C26, 
                                                                  ylims_quartiles=Cn_ylims_quartiles,
                                                                  ylims_1000s=Cn_ylims_1000s)

tiling_AT_smooth_C26_plot_list[[1]]
tiling_AT_original_C26_plot_list[[1]]

tiling_AT_smooth_C26_plot_list[[3]]
tiling_AT_original_C26_plot_list[[3]]


# Tiling C29 Comparison:
tiling_AT_smooth_C29_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C29, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)
tiling_AT_original_C29_plot_list = construct_periodicity_plots_AT(dat_tiling, C29, 
                                                                  ylims_quartiles=Cn_ylims_quartiles,
                                                                  ylims_1000s=Cn_ylims_1000s)

tiling_AT_smooth_C29_plot_list[[1]]
tiling_AT_original_C29_plot_list[[1]]

tiling_AT_smooth_C29_plot_list[[3]]
tiling_AT_original_C29_plot_list[[3]]


# Tiling C31 Comparison:
tiling_AT_smooth_C31_plot_list = construct_periodicity_plots_AT(dat_tiling, smooth_C31, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)
tiling_AT_original_C31_plot_list = construct_periodicity_plots_AT(dat_tiling, C31, 
                                                                  ylims_quartiles=Cn_ylims_quartiles,
                                                                  ylims_1000s=Cn_ylims_1000s)

tiling_AT_smooth_C31_plot_list[[1]]
tiling_AT_original_C31_plot_list[[1]]

tiling_AT_smooth_C31_plot_list[[3]]
tiling_AT_original_C31_plot_list[[3]]


# ChrV C26 Comparison:
chrv_AT_smooth_C26_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C26, 
                                                              ylims_quartiles=Cn_ylims_quartiles,
                                                              ylims_1000s=Cn_ylims_1000s)
chrv_AT_original_C26_plot_list = construct_periodicity_plots_AT(dat_chrv, C26, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)

chrv_AT_smooth_C26_plot_list[[1]]
chrv_AT_original_C26_plot_list[[1]]

chrv_AT_smooth_C26_plot_list[[3]]
chrv_AT_original_C26_plot_list[[3]]


# ChrV C29 Comparison:
chrv_AT_smooth_C29_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C29, 
                                                              ylims_quartiles=Cn_ylims_quartiles,
                                                              ylims_1000s=Cn_ylims_1000s)
chrv_AT_original_C29_plot_list = construct_periodicity_plots_AT(dat_chrv, C29, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)

chrv_AT_smooth_C29_plot_list[[1]]
chrv_AT_original_C29_plot_list[[1]]

chrv_AT_smooth_C29_plot_list[[3]]
chrv_AT_original_C29_plot_list[[3]]


# ChrV C31 Comparison:
chrv_AT_smooth_C31_plot_list = construct_periodicity_plots_AT(dat_chrv, smooth_C31, 
                                                              ylims_quartiles=Cn_ylims_quartiles,
                                                              ylims_1000s=Cn_ylims_1000s)
chrv_AT_original_C31_plot_list = construct_periodicity_plots_AT(dat_chrv, C31, 
                                                                ylims_quartiles=Cn_ylims_quartiles,
                                                                ylims_1000s=Cn_ylims_1000s)

chrv_AT_smooth_C31_plot_list[[1]]
chrv_AT_original_C31_plot_list[[1]]

chrv_AT_smooth_C31_plot_list[[3]]
chrv_AT_original_C31_plot_list[[3]]

# ggsave("figures/periodicity/nuc_AT_q1vq4_C26.png", plot = nuc_AT_q1vq4_C26)
# 
# ggsave("figures/periodicity/nuc_AT_q2vq3_C26.png", plot = nuc_AT_q2vq3_C26)
# 
# ggsave("figures/periodicity/nuc_AT_custom1_C26.png", plot = nuc_AT_custom1_C26)






## AA, TT, AT, TA:

construct_periodicity_plots_AA_TT_AT_TA = function(data, var, ylims_quartiles=c(0.35, 0.425), 
                                                   ylims_1000s=c(0.325, 0.65)) {
  
  var_string = deparse(substitute(var))
  data_string = deparse(substitute(data))
  var_of_interest = data %>% select({{var}})
  
  # Find cutoffs for first-fourth quartiles + custom
  cutoffs_0_25 = quantile(var_of_interest[,1], c(0, 0.25))
  cutoffs_25_50 = quantile(var_of_interest[,1], c(0.25, 0.5))
  cutoffs_50_75 = quantile(var_of_interest[,1], c(0.5, 0.75))
  cutoffs_75_100 = quantile(var_of_interest[,1], c(0.75, 1))
  cutoffs_custom1 = c(sort(var_of_interest[,1], TRUE)[1000], cutoffs_75_100[2])
  cutoffs_custom2 = c(cutoffs_0_25[1], sort(var_of_interest[,1], FALSE)[1000])
  
  # Divide each library into quartiles + custom
  var_q1 = data %>%
    filter({{var}} >= cutoffs_0_25[1] & {{var}} <= cutoffs_0_25[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q2 = data %>%
    filter({{var}} >= cutoffs_25_50[1] & {{var}} <= cutoffs_25_50[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q3 = data %>%
    filter({{var}} >= cutoffs_50_75[1] & {{var}} <= cutoffs_50_75[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_q4 = data %>%
    filter({{var}} >= cutoffs_75_100[1] & {{var}} <= cutoffs_75_100[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom1 = data %>%
    filter({{var}} >= cutoffs_custom1[1] & {{var}} <= cutoffs_custom1[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  var_custom2 = data %>%
    filter({{var}} >= cutoffs_custom2[1] & {{var}} <= cutoffs_custom2[2]) %>%
    select(all_of(ps1), all_of(ps2), {{var}})
  
  # Find the relative frequencies of AA, TT, AT and TA at each position (1-49) for each quartile + custom
  AA_TT_AT_TA_var_q1 = apply(var_q1 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q1))
  })
  AA_TT_AT_TA_var_q2 = apply(var_q2 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q2))
  })
  AA_TT_AT_TA_var_q3 = apply(var_q3 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q3))
  })
  AA_TT_AT_TA_var_q4 = apply(var_q4 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_q4))
  })
  AA_TT_AT_TA_var_custom1 = apply(var_custom1 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_custom1))
  })
  AA_TT_AT_TA_var_custom2 = apply(var_custom2 %>% select(all_of(ps2)), 2, function(col) {
    AA_freq = sum(col == "AA")
    TT_freq = sum(col == "TT")
    AT_freq = sum(col == "AT")
    TA_freq = sum(col == "TA")
    return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(var_custom2))
  })
  
  # Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles + custom
  AA_TT_AT_TA_var_q = data.frame(q1 = AA_TT_AT_TA_var_q1,
                                 q2 = AA_TT_AT_TA_var_q2,
                                 q3 = AA_TT_AT_TA_var_q3,
                                 q4 = AA_TT_AT_TA_var_q4, 
                                 custom1 = AA_TT_AT_TA_var_custom1,
                                 custom2 = AA_TT_AT_TA_var_custom2)
  
  # Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4, 
  # q2 v q3, and custom)
  AATT_q1vq4_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=q1, color="0-25%")) +
    geom_line(aes(y=q4, color="75-100%")) +
    scale_color_manual(breaks = c("0-25%", "75-100%"),
                       values = c("blue", "red"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") +
    labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_quartiles)
  
  # AATT_q1vq4_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
  #   geom_line(aes(y=q1, color="0-25%")) + 
  #   geom_line(aes(y=q4, color="75-100%"))
  
  AATT_q2vq3_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=q2, color="25-50%")) +
    geom_line(aes(y=q3, color="50-75%")) +
    scale_color_manual(breaks = c("25-50%", "50-75%"),
                       values = c("green", "orange"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") +
    labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_quartiles)
  
  AATT_custom1_var = ggplot(data = AA_TT_AT_TA_var_q, aes(x = 1:49)) +
    geom_line(aes(y=custom1, color="Top 1000")) +
    geom_line(aes(y=custom2, color="Bottom 1000")) +
    scale_color_manual(breaks = c("Top 1000", "Bottom 1000"),
                       values = c("firebrick1", "dodgerblue"),
                       name="Cutoffs") +
    xlab("Nucleotide Position") +
    ylab("Average AA/TT/AT/TA Content") +
    labs(title = paste0("Average AA/TT/AT/TA Content by Position in ", data_string, ", ", var_string)) +
    ylim(ylims_1000s)
  
  plot_list = list(AATT_q1vq4_var, AATT_q2vq3_var, AATT_custom1_var)
  
  return(plot_list)
}


# ChrV Predicted Cn
ylims_orig = c(.4, .85)

chrv_2_AATT_predicted_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C26, ylims_1000s = c(0.2, 0.7))
chrv_2_AATT_predicted_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C29, ylims_1000s = c(0.2, 0.7))
chrv_2_AATT_predicted_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C31, ylims_1000s = c(0.2, 0.7))

chrv_2_AATT_predicted_C26_plot_list[[3]]

chrv_2_AATT_predicted_C29_plot_list[[3]]

chrv_2_AATT_predicted_C31_plot_list[[3]]


# ChrV Raw Cn
chrv_2_AATT_raw_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C26_raw, ylims_1000s = c(0.2, 0.6))
chrv_2_AATT_raw_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C29_raw, ylims_1000s = c(0.2, 0.6))
chrv_2_AATT_raw_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, C31_raw, ylims_1000s = c(0.2, 0.6))

chrv_2_AATT_raw_C26_plot_list[[3]]

chrv_2_AATT_raw_C29_plot_list[[3]]

chrv_2_AATT_raw_C31_plot_list[[3]]


# ChrV k=3
chrv_2_AATT_smooth_3bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_3bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = ylims_orig)

chrv_2_AATT_smooth_3bp_C26_plot_list[[1]]


# ChrV k=5
chrv_2_AATT_smooth_5bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_5bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = ylims_orig)

chrv_2_AATT_smooth_5bp_C26_plot_list[[1]]


# ChrV k=7
chrv_2_AATT_smooth_7bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_7bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = ylims_orig)

chrv_2_AATT_smooth_7bp_C26_plot_list[[1]]


# ChrV k=9
chrv_2_AATT_smooth_9bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_9bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = ylims_orig)

chrv_2_AATT_smooth_9bp_C26_plot_list[[1]]


# ChrV k=11
chrv_2_AATT_smooth_11bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_11bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))

chrv_2_AATT_smooth_11bp_C26_plot_list[[1]]

chrv_2_AATT_smooth_11bp_C26_plot_list[[3]]


# ChrV k=10.7
chrv_2_AATT_smooth_10.7bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_10.7bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))

chrv_2_AATT_smooth_10.7bp_C26_plot_list[[1]]

chrv_2_AATT_smooth_10.7bp_C26_plot_list[[3]]


# ChrV k=10.4
chrv_2_AATT_smooth_10.4bp_C26_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_10.4bp_C26, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))
chrv_2_AATT_smooth_10.4bp_C29_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_10.4bp_C29, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))
chrv_2_AATT_smooth_10.4bp_C31_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_10.4bp_C31, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))

chrv_2_AATT_smooth_10.4bp_C02_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_chrv_2, smooth_10.4bp_C02, ylims_quartiles = c(0.3, 0.5), ylims_1000s = c(0.2, 0.6))

chrv_2_AATT_smooth_10.4bp_C26_plot_list[[1]]
chrv_2_AATT_smooth_10.4bp_C29_plot_list[[1]]
chrv_2_AATT_smooth_10.4bp_C31_plot_list[[1]]

chrv_2_AATT_smooth_10.4bp_C02_plot_list[[1]]

chrv_2_AATT_smooth_10.4bp_C26_plot_list[[3]]
chrv_2_AATT_smooth_10.4bp_C29_plot_list[[3]]
chrv_2_AATT_smooth_10.4bp_C31_plot_list[[3]]

chrv_2_AATT_smooth_10.4bp_C02_plot_list[[3]]


# Tiling C0 Comparison:
# tiling_AATT_smooth_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_tiling, smooth_C0_v2)
tiling_AATT_smooth_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_tiling, smooth_C0)
tiling_AATT_original_C0_plot_list = construct_periodicity_plots_AA_TT_AT_TA(dat_tiling, C0)

tiling_AATT_smooth_C0_plot_list[[1]]
tiling_AATT_original_C0_plot_list[[1]]

grid.arrange(tiling_AATT_smooth_C0_plot_list[[1]], tiling_AATT_original_C0_plot_list[[1]])

tiling_AATT_smooth_C0_plot_list[[3]]
tiling_AATT_original_C0_plot_list[[3]]

grid.arrange(tiling_AATT_smooth_C0_plot_list[[3]], tiling_AATT_original_C0_plot_list[[3]])









dat_tiling %>%
  arrange(smooth_C0) %>%
  mutate(smooth_vs_orig_difference = smooth_C0 - C0) %>%
  select(x50mer, smooth_vs_orig_difference)

