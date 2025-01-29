library(tidyverse)

# dat_tiling = read.csv("data/Created/tiling_ir_lstm_cn_tiling_post_smoothed.csv")
# dat_chrv = read.csv("data/Created/chrv_ir_lstm_cn_tiling_post_smoothed_matched.csv")

dat_chrI = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrI_1bpresolution_subsequence50_smoothC0_10_11.csv")
dat_chrV = read.csv("data/Created/ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11.csv")

dat_chrV = read.csv("data/predictions/ir_lstm_smoothC0_10_11_contracted_tiling_ir_lstm_cn_tiling_yeast_chrV_1bpresolution_subsequence50_smoothC0_10_11_best_fold_predictions.csv")


# dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_full_reconstructed_smoothC0.csv")
dat_tiling = read.csv("data/Created/ir_lstm_cn_tiling_tiling_contracted_smoothC0_10_11.csv")

# dat_tiling = dat_tiling[which(!is.na(dat_tiling$smooth_10.4bp_cn_mean2)),]
# dat_chrv = dat_chrv[which(!is.na(dat_chrv$smooth_10.4bp_cn_mean2)),]

dat_chrI = dat_chrI[which(!is.na(dat_chrI$smoothC0)),]
dat_chrV = dat_chrV[which(!is.na(dat_chrV$smoothC0)),]

dat_tiling = dat_tiling[which(!is.na(dat_tiling$smoothC0)),]


# dat_tiling = dat_tiling %>%
#   rename(x50mer = sequence)
# dat_chrv = dat_chrv %>%
#   rename(x50mer = sequence)

dat_chrI = dat_chrI %>%
  rename(x50mer = sequence)
dat_chrV = dat_chrV %>%
  rename(x50mer = sequence)

dat_tiling = dat_tiling %>%
  rename(x50mer = sequence)


chrI_cutoffs_0_25 = quantile(dat_chrI$smoothC0, c(0, 0.25))
chrI_cutoffs_75_100 = quantile(dat_chrI$smoothC0, c(0.75, 1))

chrV_cutoffs_0_25 = quantile(dat_chrV$smoothC0, c(0, 0.25))
chrV_cutoffs_75_100 = quantile(dat_chrV$smoothC0, c(0.75, 1))

chrV_cutoffs_0_25 = quantile(dat_chrV$smoothC0_predictions, c(0, 0.25))
chrV_cutoffs_75_100 = quantile(dat_chrV$smoothC0_predictions, c(0.75, 1))

tiling_cutoffs_0_25 = quantile(dat_tiling$smoothC0, c(0, 0.25))
tiling_cutoffs_75_100 = quantile(dat_tiling$smoothC0, c(0.75, 1))


dat_chrI_q1 = dat_chrI %>%
  filter(smoothC0 >= chrI_cutoffs_0_25[1] & smoothC0 <= chrI_cutoffs_0_25[2])
dat_chrI_q4 = dat_chrI %>%
  filter(smoothC0 >= chrI_cutoffs_75_100[1] & smoothC0 <= chrI_cutoffs_75_100[2])

# Original:
dat_chrV_q1 = dat_chrV %>%
  filter(smoothC0 >= chrV_cutoffs_0_25[1] & smoothC0 <= chrV_cutoffs_0_25[2])
dat_chrV_q4 = dat_chrV %>%
  filter(smoothC0 >= chrV_cutoffs_75_100[1] & smoothC0 <= chrV_cutoffs_75_100[2])

# Predictions:
dat_chrV_q1 = dat_chrV %>%
  filter(smoothC0_predictions >= chrV_cutoffs_0_25[1] & smoothC0_predictions <= chrV_cutoffs_0_25[2])
dat_chrV_q4 = dat_chrV %>%
  filter(smoothC0_predictions >= chrV_cutoffs_75_100[1] & smoothC0_predictions <= chrV_cutoffs_75_100[2])

dat_tiling_q1 = dat_tiling %>%
  filter(smoothC0 >= tiling_cutoffs_0_25[1] & smoothC0 <= tiling_cutoffs_0_25[2])
dat_tiling_q4 = dat_tiling %>%
  filter(smoothC0 >= tiling_cutoffs_75_100[1] & smoothC0 <= tiling_cutoffs_75_100[2])


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


chrV_list1_ylims = c(0, 0.2)
dat_chrV_q1_avg_polyAT_position = avg_polyAT_position(dat_chrV_q1, "x50mer", 4)
dat_chrV_q4_avg_polyAT_position = avg_polyAT_position(dat_chrV_q4, "x50mer", 4)
dat_chrV_avg_polyAT_position = avg_polyAT_position(dat_chrV, "x50mer", 4)

dat_chrV_avg_polyAT_position_df = data.frame("q1" = dat_chrV_q1_avg_polyAT_position[[1]][4:47],
                                               "q4" = dat_chrV_q4_avg_polyAT_position[[1]][4:47],
                                               position=4:47)

ggplot(data=dat_chrV_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


tiling_list1_ylims = c(0, 0.2)
dat_tiling_q1_avg_polyAT_position = avg_polyAT_position(dat_tiling_q1, "x50mer", 4)
dat_tiling_q4_avg_polyAT_position = avg_polyAT_position(dat_tiling_q4, "x50mer", 4)
dat_tiling_avg_polyAT_position = avg_polyAT_position(dat_tiling, "x50mer", 4)

dat_tiling_avg_polyAT_position_df = data.frame("q1" = dat_tiling_q1_avg_polyAT_position[[1]][4:47],
                                               "q4" = dat_tiling_q4_avg_polyAT_position[[1]][4:47],
                                               position=4:47)
ggplot(data=dat_tiling_avg_polyAT_position_df, aes(x=position)) +
  geom_line(aes(y=q1), color="blue") + 
  geom_line(aes(y=q4), color="red") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


count_polyAT_seqs <- function(df, column_name, min_length) {
  # Create regex pattern for polyA/T of specified length
  # pattern <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern <- paste0("(A|T){", min_length, ",}")
  
  # Check for the pattern in each sequence
  has_polyA_T <- grepl(pattern, df[[column_name]])
  
  # Count the number of sequences containing the pattern
  count <- sum(has_polyA_T)
  
  return(count/dim(dat_tiling_q1)[1])
}
count_polyAT_subseqs <- function(sequence, min_length) {
  # Create regex pattern for polyA/T of specified length
  # pattern <- paste0("A{", min_length, ",}|T{", min_length, ",}")
  pattern <- paste0("(A|T){", min_length, ",}")
  
  # Find all matches of the pattern
  matches <- gregexpr(pattern, sequence, perl = TRUE)
  
  # Count the number of matches
  count <- sum(sapply(matches, function(x) ifelse(x[1] == -1, 0, length(x))))
  
  return(count)
}
find_ATcontent <- function(sequence) {
  count_A = sum(strsplit(sequence, "")[[1]] == "A")
  count_T = sum(strsplit(sequence, "")[[1]] == "T")
  
  return(count_A + count_T)
}

count_polyAT_seqs(dat_tiling_q1, "x50mer", 4)
# 0.9325098
count_polyAT_seqs(dat_tiling_q4, "x50mer", 4)
# 0.9693049

count_polyAT_seqs(dat_tiling_q1, "x50mer", 5)
# 0.7896427
count_polyAT_seqs(dat_tiling_q4, "x50mer", 5)
# 0.8664015

count_polyAT_seqs(dat_tiling_q1, "x50mer", 6)
# 0.631763
count_polyAT_seqs(dat_tiling_q4, "x50mer", 6)
# 0.6986949

count_polyAT_seqs(dat_tiling_q1, "x50mer", 7)
# 0.5019333
count_polyAT_seqs(dat_tiling_q4, "x50mer", 7)
# 0.5071329

count_polyAT_seqs(dat_tiling_q1, "x50mer", 8)
# 0.3867253
count_polyAT_seqs(dat_tiling_q4, "x50mer", 8)
# 0.3362158



count_polyAT_seqs(dat_chrv_q1, "x50mer", 4)
# 0.961488
count_polyAT_seqs(dat_chrv_q4, "x50mer", 4)
# 0.9840941

count_polyAT_seqs(dat_chrv_q1, "x50mer", 5)
# 0.8380025
count_polyAT_seqs(dat_chrv_q4, "x50mer", 5)
# 0.8933626


dat_tiling_q1$polyAT3_count = sapply(dat_tiling_q1$x50mer, count_polyAT_subseqs, min_length=3)
dat_tiling_q4$polyAT3_count = sapply(dat_tiling_q4$x50mer, count_polyAT_subseqs, min_length=3)

mean(dat_tiling_q1$polyAT3_count)
# 4.189084
mean(dat_tiling_q4$polyAT3_count)
# 4.427436



dat_chrv_q1$polyAT3_count = sapply(dat_chrv_q1$x50mer, count_polyAT_subseqs, min_length=3)
dat_chrv_q4$polyAT3_count = sapply(dat_chrv_q4$x50mer, count_polyAT_subseqs, min_length=3)

mean(dat_chrv_q1$polyAT3_count)
# 
mean(dat_chrv_q4$polyAT3_count)
# 


dat_tiling_q1$polyAT4_count = sapply(dat_tiling_q1$x50mer, count_polyAT_subseqs, min_length=4)
dat_tiling_q4$polyAT4_count = sapply(dat_tiling_q4$x50mer, count_polyAT_subseqs, min_length=4)

mean(dat_tiling_q1$polyAT4_count)
# 
mean(dat_tiling_q4$polyAT4_count)
# 


dat_chrv_q1$polyAT4_count = sapply(dat_chrv_q1$x50mer, count_polyAT_subseqs, min_length=4)
dat_chrv_q4$polyAT4_count = sapply(dat_chrv_q4$x50mer, count_polyAT_subseqs, min_length=4)

mean(dat_chrv_q1$polyAT4_count)
# 
mean(dat_chrv_q4$polyAT4_count)
# 



count_exact_polyAT_starts <- function(sequence, length) {
  # Initialize a vector to hold counts for each position
  starts <- numeric(50 - length + 1)
  
  # Create regex pattern for polyA/T of exact specified length
  # pattern <- paste0("A{", length, "}|T{", length, "}")
  # pattern = paste0("(?<!A)A{", length, "}(?!A)|(?<!T)T{", length, "}(?!T)")
  pattern = paste0("(?<!(A|T))(A|T){",length,"}(?!(A|T))")
  # Find all matches of the pattern
  matches <- gregexpr(pattern, sequence, perl = TRUE)
  
  # Process the matches to count starts at each position
  for (match_start in matches[[1]]) {
    if (match_start != -1 && match_start <= (50 - length + 1)) {
      starts[match_start] <- starts[match_start] + 1
    }
  }
  
  return(starts)
}
calculate_average_polyAT_starts <- function(df, column_name, length) {
  # Apply the counting function to each sequence
  all_starts <- t(sapply(df[[column_name]], count_exact_polyAT_starts, length = length))
  
  # Sum the starts at each position across all sequences
  total_starts <- colSums(all_starts)
  
  # Calculate the average number of starts at each position
  average_starts <- total_starts / nrow(df)
  
  return(average_starts)
}






count_exact_polyAT_starts <- function(sequence, length) {
  # Initialize a vector to hold counts for each position
  starts <- numeric(50 - length + 1)
  
  # Create regex pattern for polyA/T of exact specified length
  pattern <- paste0("(?<!A)A{", length, "}(?!A)|(?<!T)T{", length, "}(?!T)")
  # pattern <- paste0("(?<!(A|T))(A|T){", length, "}(?!(A|T))")
  
  # Find all matches of the pattern
  matches <- gregexpr(pattern, sequence, perl = TRUE)
  
  # Process the matches to count starts at each position
  for (match_start in matches[[1]]) {
    if (match_start != -1) {
      # Adjust the start position to ensure it is within valid range
      valid_start <- ifelse(match_start == 1, 1, match_start)
      if (valid_start <= (50 - length + 1)) {
        starts[valid_start] <- starts[valid_start] + 1
      }
    }
  }
  
  return(starts)
}

# Function to calculate average starts across all sequences
calculate_average_polyAT_starts <- function(df, column_name, length) {
  # Apply the counting function to each sequence
  all_starts <- t(sapply(df[[column_name]], count_exact_polyAT_starts, length = length))
  
  # Sum the starts at each position across all sequences
  total_starts <- colSums(all_starts)
  
  # Calculate the average number of starts at each position
  average_starts <- total_starts / nrow(df)
  
  return(average_starts)
}


# Example usage
average_starts <- calculate_average_polyA_T_starts(dat_tiling_q1, "x50mer", 5)
print(average_starts)
plot(average_starts)




calculate_average_polyAT_starts <- function(df, column_name, length) {
  # Apply the counting function to each sequence
  # pattern = paste0("(?=((A|T){",length,"}))")
  # pattern = paste0("(^|[^AT])[AT]{",length,"}($|[^AT])")
  pattern = paste0("(?<![AT])[AT]{",length,"}(?![AT])")
  start_indices = unlist(gregexpr(pattern, df[[column_name]], perl=TRUE))
  
  total_starts = numeric(50-length-1)
  for(idx in start_indices) {
    if(idx == -1){
      next
    }
    total_starts[idx] = total_starts[idx]+1
  }
  
  # Calculate the average number of starts at each position
  average_starts <- total_starts / nrow(df)
  
  return(average_starts)
}


dat_tiling_q1_polyAT2_starts = calculate_average_polyAT_starts(dat_tiling_q1, "x50mer", 2)
dat_tiling_q4_polyAT2_starts = calculate_average_polyAT_starts(dat_tiling_q4, "x50mer", 2)

plot(dat_tiling_q1_polyAT2_starts)
plot(dat_tiling_q4_polyAT2_starts)



dat_tiling_q1_polyAT3_starts = calculate_average_polyAT_starts(dat_tiling_q1, "x50mer", 3)
dat_tiling_q4_polyAT3_starts = calculate_average_polyAT_starts(dat_tiling_q4, "x50mer", 3)

plot(dat_tiling_q1_polyAT3_starts)
plot(dat_tiling_q4_polyAT3_starts)



dat_tiling_q1_polyAT4_starts = calculate_average_polyAT_starts(dat_tiling_q1, "x50mer", 4)
dat_tiling_q4_polyAT4_starts = calculate_average_polyAT_starts(dat_tiling_q4, "x50mer", 4)
dat_tiling_polyAT4_starts = calculate_average_polyAT_starts(dat_tiling, "x50mer", 4)


plot(dat_tiling_q1_polyAT4_starts)
plot(dat_tiling_q4_polyAT4_starts)
plot(dat_tiling_polyAT4_starts)



dat_tiling_q1_polyAT5_starts = calculate_average_polyAT_starts(dat_tiling_q1, "x50mer", 5)
dat_tiling_q4_polyAT5_starts = calculate_average_polyAT_starts(dat_tiling_q4, "x50mer", 5)

plot(dat_tiling_q1_polyAT5_starts)
plot(dat_tiling_q4_polyAT5_starts)





dat_tiling_q1$AT_count = sapply(dat_tiling_q1$x50mer, find_ATcontent)
dat_tiling_q4$AT_count = sapply(dat_tiling_q4$x50mer, find_ATcontent)

