library(tidyverse)
library(gridExtra)

source("scripts/functions/sequence-functions.R")
source("scripts/functions/plotting-functions.R")

dat_nuc = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle1_predictions.csv")
dat_random = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle3_predictions.csv")
dat_tiling = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle5_predictions.csv")
dat_chrv = read.csv("data/predictions/ir_lstm_smooth_104bp_cn_mean2_tiling_cycle6_predictions.csv")

colnames(dat_nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")
colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase", "Predicted_SmoothC0")


nuc_cutoffs_0_25 = quantile(dat_nuc$Predicted_SmoothC0, c(0, 0.25))
nuc_cutoffs_75_100 = quantile(dat_nuc$Predicted_SmoothC0, c(0.75, 1))

random_cutoffs_0_25 = quantile(dat_random$Predicted_SmoothC0, c(0, 0.25))
random_cutoffs_75_100 = quantile(dat_random$Predicted_SmoothC0, c(0.75, 1))

tiling_cutoffs_0_25 = quantile(dat_tiling$Predicted_SmoothC0, c(0, 0.25))
tiling_cutoffs_75_100 = quantile(dat_tiling$Predicted_SmoothC0, c(0.75, 1))

chrv_cutoffs_0_25 = quantile(dat_chrv$Predicted_SmoothC0, c(0, 0.25))
chrv_cutoffs_75_100 = quantile(dat_chrv$Predicted_SmoothC0, c(0.75, 1))


dat_nuc_q1 = dat_nuc %>%
  filter(Predicted_SmoothC0 >= nuc_cutoffs_0_25[1] & Predicted_SmoothC0 <= nuc_cutoffs_0_25[2])
dat_nuc_q4 = dat_nuc %>%
  filter(Predicted_SmoothC0 >= nuc_cutoffs_75_100[1] & Predicted_SmoothC0 <= nuc_cutoffs_75_100[2])

dat_random_q1 = dat_random %>%
  filter(Predicted_SmoothC0 >= random_cutoffs_0_25[1] & Predicted_SmoothC0 <= random_cutoffs_0_25[2])
dat_random_q4 = dat_random %>%
  filter(Predicted_SmoothC0 >= random_cutoffs_75_100[1] & Predicted_SmoothC0 <= random_cutoffs_75_100[2])

dat_tiling_q1 = dat_tiling %>%
  filter(Predicted_SmoothC0 >= tiling_cutoffs_0_25[1] & Predicted_SmoothC0 <= tiling_cutoffs_0_25[2])
dat_tiling_q4 = dat_tiling %>%
  filter(Predicted_SmoothC0 >= tiling_cutoffs_75_100[1] & Predicted_SmoothC0 <= tiling_cutoffs_75_100[2])

dat_chrv_q1 = dat_chrv %>%
  filter(Predicted_SmoothC0 >= chrv_cutoffs_0_25[1] & Predicted_SmoothC0 <= chrv_cutoffs_0_25[2])
dat_chrv_q4 = dat_chrv %>%
  filter(Predicted_SmoothC0 >= chrv_cutoffs_75_100[1] & Predicted_SmoothC0 <= chrv_cutoffs_75_100[2])


# (1): # Proportion of sequences containing at least one polyA/polyT/polyAandT:
prop_seqs_with_polyAT <- function(df, column_name, min_length) {
  pattern_AorT <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern_AandT <- paste0("(A|T){", min_length, ",}")
  
  has_polyAorT <- grepl(pattern_AorT, df[[column_name]])
  has_polyAandT <- grepl(pattern_AandT, df[[column_name]])
  
  count_polyAorT <- sum(has_polyAorT)
  count_polyAandT <- sum(has_polyAandT)
  
  return(c(count_polyAorT/dim(df)[1], count_polyAandT/dim(df)[1]))
}


prop_seqs_with_polyAT(dat_nuc_q1, "x50mer", 4)
# [1] 0.6580269 0.9600161
prop_seqs_with_polyAT(dat_nuc_q4, "x50mer", 4)
# [1] 0.6311031 0.9863372
prop_seqs_with_polyAT(dat_nuc, "x50mer", 4)
# [1] 0.5611092 0.9650877

prop_seqs_with_polyAT(dat_random_q1, "x50mer", 4)
# [1] 0.2905709 0.7950609
prop_seqs_with_polyAT(dat_random_q4, "x50mer", 4)
# [1] 0.3159076 0.8948044
prop_seqs_with_polyAT(dat_random, "x50mer", 4)
# [1] 0.2463117 0.8260103

prop_seqs_with_polyAT(dat_tiling_q1, "x50mer", 4)
# [1] 0.7149378 0.9355089
prop_seqs_with_polyAT(dat_tiling_q4, "x50mer", 4)
# [1] 0.6372863 0.9698427
prop_seqs_with_polyAT(dat_tiling, "x50mer", 4)
# [1] 0.5809659 0.9388233

prop_seqs_with_polyAT(dat_chrv_q1, "x50mer", 4)
# [1] 0.7114218 0.9582059
prop_seqs_with_polyAT(dat_chrv_q4, "x50mer", 4)
# [1] 0.6184651 0.9776224
prop_seqs_with_polyAT(dat_chrv, "x50mer", 4)
# [1] 0.5695840 0.9554755

# Q1 Higher proportion of AorT (except random)
# Q4 Higher proportion of AandT across the board (interesting result)




# (2): # Average length of any given polyA/polyT/polyAandT:
avg_polyAT_length <- function(df, column_name, min_length) {
  pattern_AorT <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern_AandT <- paste0("(A|T){", min_length, ",}")
  
  matches_AorT <- gregexpr(pattern_AorT, df[[column_name]], perl = TRUE)
  matches_AandT <- gregexpr(pattern_AandT, df[[column_name]], perl = TRUE)
  
  AorT_lengths = lapply(matches_AorT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(match_lengths)
    # for(len in attr(match_obj, "match.length")) {
    #   if(len != -1) {return(len)}
    # }
  }) %>% unlist()
  AandT_lengths = lapply(matches_AandT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(match_lengths)
    # for(len in attr(match_obj, "match.length")) {
    #   if(len != -1) {return(len)}
    # }
  }) %>% unlist()
  
  return(c(mean(AorT_lengths), mean(AandT_lengths)))
}

avg_polyAT_length(dat_nuc_q1, "x50mer", 4)
# [1] 4.619028 5.528519
avg_polyAT_length(dat_nuc_q4, "x50mer", 4)
# [1] 4.483248 5.519916
avg_polyAT_length(dat_nuc, "x50mer", 4)
# [1] 4.501303 5.430173

avg_polyAT_length(dat_random_q1, "x50mer", 4)
# [1] 4.351834 4.989982
avg_polyAT_length(dat_random_q4, "x50mer", 4)
# [1] 4.306104 5.000544
avg_polyAT_length(dat_random, "x50mer", 4)
# [1] 4.306069 4.961757

avg_polyAT_length(dat_tiling_q1, "x50mer", 4)
# [1] 5.085949 5.966561
avg_polyAT_length(dat_tiling_q4, "x50mer", 4)
# [1] 4.717143 5.673527
avg_polyAT_length(dat_tiling, "x50mer", 4)
# [1] 4.808349 5.641722

avg_polyAT_length(dat_chrv_q1, "x50mer", 4)
# [1] 4.840084 5.771236
avg_polyAT_length(dat_chrv_q4, "x50mer", 4)
# [1] 4.566924 5.593728
avg_polyAT_length(dat_chrv, "x50mer", 4)
# [1] 4.640883 5.547685

# Q1 longer on average across the board except random AandT






# (3): # Average position of polyA/polyT/polyAandT:
avg_polyAT_position <- function(df, column_name, min_length) {
  pattern_AorT <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern_AandT <- paste0("(A|T){", min_length, ",}")
  
  matches_AorT <- gregexpr(pattern_AorT, df[[column_name]], perl = TRUE)
  matches_AandT <- gregexpr(pattern_AandT, df[[column_name]], perl = TRUE)
  
  AorT_lengths = lapply(matches_AorT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(match_lengths)
  }) %>% unlist()
  AandT_lengths = lapply(matches_AandT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(match_lengths)
  }) %>% unlist()
  
  AorT_starts = lapply(matches_AorT, function(match_obj) {
    if(match_obj[1] == -1) {return()}
    return(match_obj)
  }) %>% unlist()
  AandT_starts = lapply(matches_AandT, function(match_obj) {
    if(match_obj[1] == -1) {return()}
    return(match_obj)
  }) %>% unlist()
  
  AorT_spanned = numeric(50)
  for(i in 1:length(AorT_starts)) {
    for(j in AorT_starts[i]:(AorT_starts[i] + AorT_lengths[i]-1)) {
      AorT_spanned[j] = AorT_spanned[j] + 1
    }
  }
  
  AandT_spanned = numeric(50)
  for(i in 1:length(AandT_starts)) {
    for(j in AandT_starts[i]:(AandT_starts[i] + AandT_lengths[i]-1)) {
      AandT_spanned[j] = AandT_spanned[j] + 1
      if(j>50){
        print(c(AandT_starts[i], AandT_lengths[i]))
      }
    }
  }
  
  ret = list()
  ret[[1]] = AorT_spanned/dim(df)[1]
  ret[[2]] = AandT_spanned/dim(df)[1]
  return(ret)
}

nuc_list1_ylims = c(0, 0.2)
dat_nuc_q1_avg_polyAT_position = avg_polyAT_position(dat_nuc_q1, "x50mer", 4)
# plot(dat_nuc_q1_avg_polyAT_position[[1]], ylim=nuc_list1_ylims)
# plot(dat_nuc_q1_avg_polyAT_position[[2]])
dat_nuc_q4_avg_polyAT_position = avg_polyAT_position(dat_nuc_q4, "x50mer", 4)
# plot(dat_nuc_q4_avg_polyAT_position[[1]], ylim=nuc_list1_ylims)
# plot(dat_nuc_q4_avg_polyAT_position[[2]])
dat_nuc_avg_polyAT_position = avg_polyAT_position(dat_nuc, "x50mer", 4)
# plot(dat_nuc_avg_polyAT_position[[1]], ylim=nuc_list1_ylims)
# plot(dat_nuc_avg_polyAT_position[[2]])

dat_nuc_avg_polyAT_position_df = data.frame("q1" = dat_nuc_q1_avg_polyAT_position[[1]],
                                            "q4" = dat_nuc_q4_avg_polyAT_position[[1]],
                                            position=1:50)
ggplot(data=dat_nuc_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red")

# plot(dat_nuc_q1_avg_polyAT_position[[1]] - dat_nuc_avg_polyAT_position[[1]])
# plot(dat_nuc_q4_avg_polyAT_position[[1]] - dat_nuc_avg_polyAT_position[[1]])

random_list1_ylims = c(0, 0.05)
dat_random_q1_avg_polyAT_position = avg_polyAT_position(dat_random_q1, "x50mer", 4)
# plot(dat_random_q1_avg_polyAT_position[[1]], ylim=random_list1_ylims)
# plot(dat_random_q1_avg_polyAT_position[[2]])
dat_random_q4_avg_polyAT_position = avg_polyAT_position(dat_random_q4, "x50mer", 4)
# plot(dat_random_q4_avg_polyAT_position[[1]], ylim=random_list1_ylims)
# plot(dat_random_q4_avg_polyAT_position[[2]])
dat_random_avg_polyAT_position = avg_polyAT_position(dat_random, "x50mer", 4)
# plot(dat_random_avg_polyAT_position[[1]], ylim=random_list1_ylims)
# plot(dat_random_avg_polyAT_position[[2]])

dat_random_avg_polyAT_position_df = data.frame("q1" = dat_random_q1_avg_polyAT_position[[1]],
                                               "q4" = dat_random_q4_avg_polyAT_position[[1]],
                                               position=1:50)
ggplot(data=dat_random_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red")

# plot(dat_random_q1_avg_polyAT_position[[1]] - dat_random_avg_polyAT_position[[1]])
# plot(dat_random_q4_avg_polyAT_position[[1]] - dat_random_avg_polyAT_position[[1]])

tiling_list1_ylims = c(0, 0.2)
# tiling_list2_ylims = c(0, 0.4)
dat_tiling_q1_avg_polyAT_position = avg_polyAT_position(dat_tiling_q1, "x50mer", 4)
# plot(dat_tiling_q1_avg_polyAT_position[[1]], ylim=tiling_list1_ylims)
# plot(dat_tiling_q1_avg_polyAT_position[[2]], ylim=tiling_list2_ylims)
dat_tiling_q4_avg_polyAT_position = avg_polyAT_position(dat_tiling_q4, "x50mer", 4)
# plot(dat_tiling_q4_avg_polyAT_position[[1]], ylim=tiling_list1_ylims)
# plot(dat_tiling_q4_avg_polyAT_position[[2]], ylim=tiling_list2_ylims)
dat_tiling_avg_polyAT_position = avg_polyAT_position(dat_tiling, "x50mer", 4)
# plot(dat_tiling_avg_polyAT_position[[1]], ylim=tiling_list1_ylims)
# plot(dat_tiling_avg_polyAT_position[[2]], ylim=tiling_list2_ylims)

dat_tiling_avg_polyAT_position_df = data.frame("q1" = dat_tiling_q1_avg_polyAT_position[[1]],
                                               "q4" = dat_tiling_q4_avg_polyAT_position[[1]],
                                               position=1:50)
ggplot(data=dat_tiling_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red")

# plot(dat_tiling_q1_avg_polyAT_position[[1]] - dat_tiling_avg_polyAT_position[[1]])
# plot(dat_tiling_q4_avg_polyAT_position[[1]] - dat_tiling_avg_polyAT_position[[1]])

chrv_list1_ylims = c(0, 0.2)
# chrv_list2_ylims = c(0, 0.4)
dat_chrv_q1_avg_polyAT_position = avg_polyAT_position(dat_chrv_q1, "x50mer", 4)
# plot(dat_chrv_q1_avg_polyAT_position[[1]], ylim=chrv_list1_ylims)
# plot(dat_chrv_q1_avg_polyAT_position[[2]], ylim=chrv_list2_ylims)
dat_chrv_q4_avg_polyAT_position = avg_polyAT_position(dat_chrv_q4, "x50mer", 4)
# plot(dat_chrv_q4_avg_polyAT_position[[1]], ylim=chrv_list1_ylims)
# plot(dat_chrv_q4_avg_polyAT_position[[2]], ylim=chrv_list2_ylims)
dat_chrv_avg_polyAT_position = avg_polyAT_position(dat_chrv, "x50mer", 4)
# plot(dat_chrv_avg_polyAT_position[[1]], ylim=chrv_list1_ylims)
# plot(dat_chrv_avg_polyAT_position[[2]], ylim=chrv_list2_ylims)

dat_chrv_avg_polyAT_position_df = data.frame("q1" = dat_chrv_q1_avg_polyAT_position[[1]],
                                             "q4" = dat_chrv_q4_avg_polyAT_position[[1]],
                                             position=1:50)
ggplot(data=dat_chrv_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red")

# plot(dat_chrv_q1_avg_polyAT_position[[1]] - dat_chrv_avg_polyAT_position[[1]])
# plot(dat_chrv_q4_avg_polyAT_position[[1]] - dat_chrv_avg_polyAT_position[[1]])







# (4): # Average length of longest polyA/polyT/polyAandT in a sequence 
# (given that it has at least one polyA/polyT/polyAandT:
avg_polyAT_position <- function(df, column_name, min_length) {
  pattern_AorT <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern_AandT <- paste0("(A|T){", min_length, ",}")
  
  matches_AorT <- gregexpr(pattern_AorT, df[[column_name]], perl = TRUE)
  matches_AandT <- gregexpr(pattern_AandT, df[[column_name]], perl = TRUE)
  
  AorT_lengths = lapply(matches_AorT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(max(match_lengths))
  }) %>% unlist()
  AandT_lengths = lapply(matches_AandT, function(match_obj) {
    match_lengths = attr(match_obj, "match.length")
    if(match_lengths[1] == -1) {return()}
    return(max(match_lengths))
  }) %>% unlist()
  
  return(c(mean(AorT_lengths), mean(AandT_lengths)))
}

avg_polyAT_position(dat_nuc_q1, "x50mer", 4)
# [1] 4.863511 6.906865
avg_polyAT_position(dat_nuc_q4, "x50mer", 4)
# [1] 4.636422 6.926665
avg_polyAT_position(dat_nuc, "x50mer", 4)
# [1] 4.655237 6.698001

avg_polyAT_position(dat_random_q1, "x50mer", 4)
# [1] 4.397351 5.479629
avg_polyAT_position(dat_random_q4, "x50mer", 4)
# [1] 4.341117 5.568817
avg_polyAT_position(dat_random, "x50mer", 4)
# [1] 4.336589 5.454766

avg_polyAT_position(dat_tiling_q1, "x50mer", 4)
# [1] 5.627836 7.720619
avg_polyAT_position(dat_tiling_q4, "x50mer", 4)
# [1] 4.947040 7.183065
avg_polyAT_position(dat_tiling, "x50mer", 4)
# [1] 5.080413 7.030881

avg_polyAT_position(dat_chrv_q1, "x50mer", 4)
# [1] 5.202852 7.424823
avg_polyAT_position(dat_chrv_q4, "x50mer", 4)
# [1] 4.731418 7.023337
avg_polyAT_position(dat_chrv, "x50mer", 4)
# [1] 4.835734 6.891382





