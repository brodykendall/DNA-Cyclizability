library(tidyverse)

# Position specific variable names 
# (e.g. X1di represents the dinucleotide in the 1st position of the sequence):
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")

################################################################################
# Tiling Library Quartiles
################################################################################

# Load tiling library training dataset
dat = readRDS("data/Created/processed_tiling_newC0_v2.rds")

# Find cutoffs for first-fourth quartiles of intrinsic cyclizability (C0)
cutoffs_0_25 = quantile(dat$C0_new, c(0, 0.25)) 
cutoffs_25_50 = quantile(dat$C0_new, c(0.25, 0.5)) 
cutoffs_50_75 = quantile(dat$C0_new, c(0.5, 0.75)) 
cutoffs_75_100 = quantile(dat$C0_new, c(0.75, 1)) 

# Divide tiling library dataset into quartiles
dat_q1 = dat %>%
  filter(C0_new >= cutoffs_0_25[1] & C0_new <= cutoffs_0_25[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q2 = dat %>%
  filter(C0_new >= cutoffs_25_50[1] & C0_new <= cutoffs_25_50[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q3 = dat %>%
  filter(C0_new >= cutoffs_50_75[1] & C0_new <= cutoffs_50_75[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q4 = dat %>%
  filter(C0_new >= cutoffs_75_100[1] & C0_new <= cutoffs_75_100[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each quartile
A_T_q1 = apply(dat_q1 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q1))
})
A_T_q2 = apply(dat_q2 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q2))
})
A_T_q3 = apply(dat_q3 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q3))
})
A_T_q4 = apply(dat_q4 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q4))
})

# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
A_T_q = data.frame(q1 = A_T_q1,
                   q2 = A_T_q2,
                   q3 = A_T_q3,
                   q4 = A_T_q4)

# Plot the relative frequencies of A and T by position for each quartile (q1 v q4 and q2 v q3)
tiling_AT_q1vq4_newC0 = ggplot(data = A_T_q, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, New C0")
tiling_AT_q1vq4_newC0
ggsave("figures/periodicity/tiling_AT_q1vq4_newC0_v2.png", plot = tiling_AT_q1vq4_newC0)
tiling_AT_q2vq3_newC0 = ggplot(data = A_T_q, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, New C0")
tiling_AT_q2vq3_newC0
ggsave("figures/periodicity/tiling_AT_q2vq3_newC0_v2.png", plot = tiling_AT_q2vq3_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile
AA_TT_AT_TA_q1 = apply(dat_q1 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q1))
})
AA_TT_AT_TA_q2 = apply(dat_q2 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q2))
})
AA_TT_AT_TA_q3 = apply(dat_q3 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q3))
})
AA_TT_AT_TA_q4 = apply(dat_q4 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q4))
})

# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles
AA_TT_AT_TA_q = data.frame(q1 = AA_TT_AT_TA_q1,
                           q2 = AA_TT_AT_TA_q2,
                           q3 = AA_TT_AT_TA_q3,
                           q4 = AA_TT_AT_TA_q4)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
tiling_AATT_q1vq4_newC0 = ggplot(data = AA_TT_AT_TA_q, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, New C0")
tiling_AATT_q1vq4_newC0
ggsave("figures/periodicity/tiling_AATT_q1vq4_newC0_v2.png", plot = tiling_AATT_q1vq4_newC0)
tiling_AATT_q2vq3_newC0 = ggplot(data = AA_TT_AT_TA_q, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, New C0")
tiling_AATT_q2vq3_newC0
ggsave("figures/periodicity/tiling_AATT_q2vq3_newC0_v2.png", plot = tiling_AATT_q2vq3_newC0)










################################################################################
# Tiling Library Top v Bottom 1000
################################################################################

# Divide tiling library dataset into top and bottom 1000
dat_top1000 = dat %>%
  slice_max(C0_new, n=1000)
dat_bot1000 = dat %>%
  slice_min(C0_new, n=1000)

# Find the relative frequencies of A and T at each position (1-50) for top and bottom 1000
A_T_top1000 = apply(dat_top1000 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
A_T_bot1000 = apply(dat_bot1000 %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
# Construct dataframe containing the relative frequences of A and T by position for top and bottom 1000
A_T_1000 = data.frame(top = A_T_top1000,
                      bot = A_T_bot1000)
# Plot the relative frequencies of A and T by position
tiling_AT_topvbot1000_newC0 = ggplot(data = A_T_1000, aes(x = 1:50)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Tiling Library, New C0")
tiling_AT_topvbot1000_newC0
ggsave("figures/periodicity/tiling_AT_topvbot1000_newC0_v2.png", plot = tiling_AT_topvbot1000_newC0)

# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for top and bottom 1000
AA_TT_AT_TA_top1000 = apply(dat_top1000 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
AA_TT_AT_TA_bot1000 = apply(dat_bot1000 %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position
# for top and bottom 1000
AA_TT_AT_TA_1000 = data.frame(top = AA_TT_AT_TA_top1000,
                              bot = AA_TT_AT_TA_bot1000)
# Plot the relative frequencies of AA, TT, AT, and TA by position
tiling_AATT_topvbot1000_newC0 = ggplot(data = AA_TT_AT_TA_1000, aes(x = 1:49)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Tiling Library, New C0")
tiling_AATT_topvbot1000_newC0
ggsave("figures/periodicity/tiling_AATT_topvbot1000_newC0_v2.png", plot = tiling_AATT_topvbot1000_newC0)









################################################################################
# Random Library Quartiles
################################################################################

# Load random library training dataset
dat_random_library = readRDS("data/Created/processed_random_newC0_v2.rds")

# Find cutoffs for first-fourth quartiles of intrinsic cyclizability (C0)
cutoffs_0_25_random_library = quantile(dat_random_library$C0_new, c(0, 0.25)) 
cutoffs_25_50_random_library = quantile(dat_random_library$C0_new, c(0.25, 0.5)) 
cutoffs_50_75_random_library = quantile(dat_random_library$C0_new, c(0.5, 0.75)) 
cutoffs_75_100_random_library = quantile(dat_random_library$C0_new, c(0.75, 1)) 

# Divide random library dataset into quartiles
dat_q1_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_0_25_random_library[1] & C0_new <= cutoffs_0_25_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q2_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_25_50_random_library[1] & C0_new <= cutoffs_25_50_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q3_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_50_75_random_library[1] & C0_new <= cutoffs_50_75_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q4_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_75_100_random_library[1] & C0_new <= cutoffs_75_100_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each quartile
A_T_q1_random_library = apply(dat_q1_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q1_random_library))
})
A_T_q2_random_library = apply(dat_q2_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q2_random_library))
})
A_T_q3_random_library = apply(dat_q3_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q3_random_library))
})
A_T_q4_random_library = apply(dat_q4_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q4_random_library))
})
# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
A_T_q_random_library = data.frame(q1 = A_T_q1_random_library,
                                  q2 = A_T_q2_random_library,
                                  q3 = A_T_q3_random_library,
                                  q4 = A_T_q4_random_library)
# Plot the relative frequencies of A and T by position for each quartile (q1 v q4 and q2 v q3)
random_AT_q1vq4_newC0 = ggplot(data = A_T_q_random_library, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, New C0")
random_AT_q1vq4_newC0
ggsave("figures/periodicity/random_AT_q1vq4_newC0_v2.png", plot = random_AT_q1vq4_newC0)
random_AT_q2vq3_newC0 = ggplot(data = A_T_q_random_library, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, New C0")
random_AT_q2vq3_newC0
ggsave("figures/periodicity/random_AT_q2vq3_newC0_v2.png", plot = random_AT_q2vq3_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile
AA_TT_AT_TA_q1_random_library = apply(dat_q1_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q1_random_library))
})
AA_TT_AT_TA_q2_random_library = apply(dat_q2_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q2_random_library))
})
AA_TT_AT_TA_q3_random_library = apply(dat_q3_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q3_random_library))
})
AA_TT_AT_TA_q4_random_library = apply(dat_q4_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q4_random_library))
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles
AA_TT_AT_TA_q_random_library = data.frame(q1 = AA_TT_AT_TA_q1_random_library,
                                          q2 = AA_TT_AT_TA_q2_random_library,
                                          q3 = AA_TT_AT_TA_q3_random_library,
                                          q4 = AA_TT_AT_TA_q4_random_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
random_AATT_q1vq4_newC0 = ggplot(data = AA_TT_AT_TA_q_random_library, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, New C0")
random_AATT_q1vq4_newC0
ggsave("figures/periodicity/random_AATT_q1vq4_newC0_v2.png", plot = random_AATT_q1vq4_newC0)
random_AATT_q2vq3_newC0 = ggplot(data = AA_TT_AT_TA_q_random_library, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, New C0")
random_AATT_q2vq3_newC0
ggsave("figures/periodicity/random_AATT_q2vq3_newC0_v2.png", plot = random_AATT_q2vq3_newC0)
















################################################################################
# Random Library Top 1000
################################################################################

# Divide random library dataset into top and bottom 1000
dat_top1000_random_library = dat_random_library %>%
  slice_max(C0_new, n=1000)
dat_bot1000_random_library = dat_random_library %>%
  slice_min(C0_new, n=1000)

# Find the relative frequencies of A and T at each position (1-50) for top and bottom 1000
A_T_top1000_random_library = apply(dat_top1000_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
A_T_bot1000_random_library = apply(dat_bot1000_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
# Construct dataframe containing the relative frequences of A and T by position for top and bottom 1000
A_T_1000_random_library = data.frame(top = A_T_top1000_random_library,
                                     bot = A_T_bot1000_random_library)
# Plot the relative frequencies of A and T by position
random_AT_topvbot1000_newC0 = ggplot(data = A_T_1000_random_library, aes(x = 1:50)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, New C0")
random_AT_topvbot1000_newC0
ggsave("figures/periodicity/random_AT_topvbot1000_newC0_v2.png", plot = random_AT_topvbot1000_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for top and bottom 1000
AA_TT_AT_TA_top1000_random_library = apply(dat_top1000_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
AA_TT_AT_TA_bot1000_random_library = apply(dat_bot1000_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position 
# for top and bottom 1000
AA_TT_AT_TA_1000_random_library = data.frame(top = AA_TT_AT_TA_top1000_random_library,
                                             bot = AA_TT_AT_TA_bot1000_random_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position
random_AATT_topvbot1000_newC0 = ggplot(data = AA_TT_AT_TA_1000_random_library, aes(x = 1:49)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, New C0")
random_AATT_topvbot1000_newC0
ggsave("figures/periodicity/random_AATT_topvbot1000_newC0_v2.png", plot = random_AATT_topvbot1000_newC0)








################################################################################
# Random Library Halves
################################################################################

# Find cutoffs for first and second halves of intrinsic cyclizability (C0)
cutoffs_0_50_random_library = quantile(dat_random_library$C0_new, c(0, 0.5)) 
cutoffs_50_100_random_library = quantile(dat_random_library$C0_new, c(0.5, 1)) 

# Divide random library dataset into halves
dat_h1_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_0_50_random_library[1] & C0_new <= cutoffs_0_50_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_h2_random_library = dat_random_library %>%
  filter(C0_new >= cutoffs_50_100_random_library[1] & C0_new <= cutoffs_50_100_random_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each half
A_T_h1_random_library = apply(dat_h1_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_h1_random_library))
})
A_T_h2_random_library = apply(dat_h2_random_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_h2_random_library))
})
# Construct dataframe containing the relative frequences of A and T by position for both halves
A_T_h_random_library = data.frame(h1 = A_T_h1_random_library,
                                  h2 = A_T_h2_random_library)
# Plot the relative frequencies of A and T by position for each half
random_AT_topvbothalf_newC0 = ggplot(data = A_T_h_random_library, aes(x = 1:50)) +
  geom_line(aes(y=h1, color="0-50%")) +
  geom_line(aes(y=h2, color="50-100%")) +
  scale_color_manual(breaks = c("0-50%", "50-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in Random Library, New C0")
random_AT_topvbothalf_newC0
ggsave("figures/periodicity/random_AT_topvbothalf_newC0_v2.png", plot = random_AT_topvbothalf_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each half
AA_TT_AT_TA_h1_random_library = apply(dat_h1_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_h1_random_library))
})
AA_TT_AT_TA_h2_random_library = apply(dat_h2_random_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_h2_random_library))
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for both halves
AA_TT_AT_TA_h_random_library = data.frame(h1 = AA_TT_AT_TA_h1_random_library,
                                          h2 = AA_TT_AT_TA_h2_random_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
random_AATT_topvbothalf_newC0 = ggplot(data = AA_TT_AT_TA_h_random_library, aes(x = 1:49)) +
  geom_line(aes(y=h1, color="0-50%")) +
  geom_line(aes(y=h2, color="50-100%")) +
  scale_color_manual(breaks = c("0-50%", "50-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in Random Library, New C0")
random_AATT_topvbothalf_newC0
ggsave("figures/periodicity/random_AATT_topvbothalf_newC0_v2.png", plot = random_AATT_topvbothalf_newC0)












################################################################################
# ChrV Library Quartiles
################################################################################

# Load chromosome V library training dataset
dat_chrV_library = readRDS("data/Created/processed_chrV_newC0_v2.rds")

# Find cutoffs for first-fourth quartiles of intrinsic cyclizability (C0)
cutoffs_0_25_chrV_library = quantile(dat_chrV_library$C0_new, c(0, 0.25)) 
cutoffs_25_50_chrV_library = quantile(dat_chrV_library$C0_new, c(0.25, 0.5)) 
cutoffs_50_75_chrV_library = quantile(dat_chrV_library$C0_new, c(0.5, 0.75)) 
cutoffs_75_100_chrV_library = quantile(dat_chrV_library$C0_new, c(0.75, 1)) 

# Divide chromosome V library dataset into quartiles
dat_q1_chrV_library = dat_chrV_library %>%
  filter(C0_new >= cutoffs_0_25_chrV_library[1] & C0_new <= cutoffs_0_25_chrV_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q2_chrV_library = dat_chrV_library %>%
  filter(C0_new >= cutoffs_25_50_chrV_library[1] & C0_new <= cutoffs_25_50_chrV_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q3_chrV_library = dat_chrV_library %>%
  filter(C0_new >= cutoffs_50_75_chrV_library[1] & C0_new <= cutoffs_50_75_chrV_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q4_chrV_library = dat_chrV_library %>%
  filter(C0_new >= cutoffs_75_100_chrV_library[1] & C0_new <= cutoffs_75_100_chrV_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each quartile
A_T_q1_chrV_library = apply(dat_q1_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q1_chrV_library))
})
A_T_q2_chrV_library = apply(dat_q2_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q2_chrV_library))
})
A_T_q3_chrV_library = apply(dat_q3_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q3_chrV_library))
})
A_T_q4_chrV_library = apply(dat_q4_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q4_chrV_library))
})
# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
A_T_q_chrV_library = data.frame(q1 = A_T_q1_chrV_library,
                                q2 = A_T_q2_chrV_library,
                                q3 = A_T_q3_chrV_library,
                                q4 = A_T_q4_chrV_library)
# Plot the relative frequencies of A and T by position for each quartile (q1 v q4 and q2 v q3)
chrV_AT_q1vq4_newC0 = ggplot(data = A_T_q_chrV_library, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in chrV Library, New C0")
chrV_AT_q1vq4_newC0
ggsave("figures/periodicity/chrV_AT_q1vq4_newC0_v2.png", plot = chrV_AT_q1vq4_newC0)
chrV_AT_q2vq3_newC0 = ggplot(data = A_T_q_chrV_library, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in chrV Library, New C0")
chrV_AT_q2vq3_newC0
ggsave("figures/periodicity/chrV_AT_q2vq3_newC0_v2.png", plot = chrV_AT_q2vq3_newC0)



# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile
AA_TT_AT_TA_q1_chrV_library = apply(dat_q1_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q1_chrV_library))
})
AA_TT_AT_TA_q2_chrV_library = apply(dat_q2_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q2_chrV_library))
})
AA_TT_AT_TA_q3_chrV_library = apply(dat_q3_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q3_chrV_library))
})
AA_TT_AT_TA_q4_chrV_library = apply(dat_q4_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q4_chrV_library))
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles
AA_TT_AT_TA_q_chrV_library = data.frame(q1 = AA_TT_AT_TA_q1_chrV_library,
                                        q2 = AA_TT_AT_TA_q2_chrV_library,
                                        q3 = AA_TT_AT_TA_q3_chrV_library,
                                        q4 = AA_TT_AT_TA_q4_chrV_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
chrV_AATT_q1vq4_newC0 = ggplot(data = AA_TT_AT_TA_q_chrV_library, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in chrV Library, New C0")
chrV_AATT_q1vq4_newC0
ggsave("figures/periodicity/chrV_AATT_q1vq4_newC0_v2.png", plot = chrV_AATT_q1vq4_newC0)
chrV_AATT_q2vq3_newC0 = ggplot(data = AA_TT_AT_TA_q_chrV_library, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in chrV Library, New C0")
chrV_AATT_q2vq3_newC0
ggsave("figures/periodicity/chrV_AATT_q2vq3_newC0_v2.png", plot = chrV_AATT_q2vq3_newC0)
















################################################################################
# chrV Library Top 1000
################################################################################

# Divide chromosome V library dataset into top and bottom 1000
dat_top1000_chrV_library = dat_chrV_library %>%
  slice_max(C0_new, n=1000)
dat_bot1000_chrV_library = dat_chrV_library %>%
  slice_min(C0_new, n=1000)

# Find the relative frequencies of A and T at each position (1-50) for top and bottom 1000
A_T_top1000_chrV_library = apply(dat_top1000_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
A_T_bot1000_chrV_library = apply(dat_bot1000_chrV_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
# Construct dataframe containing the relative frequences of A and T by position for top and bottom 1000
A_T_1000_chrV_library = data.frame(top = A_T_top1000_chrV_library,
                                   bot = A_T_bot1000_chrV_library)
# Plot the relative frequencies of A and T by position
chrV_AT_topvbot1000_newC0 = ggplot(data = A_T_1000_chrV_library, aes(x = 1:50)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in chrV Library, New C0")
chrV_AT_topvbot1000_newC0
ggsave("figures/periodicity/chrV_AT_topvbot1000_newC0_v2.png", plot = chrV_AT_topvbot1000_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for top and bottom 1000
AA_TT_AT_TA_top1000_chrV_library = apply(dat_top1000_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
AA_TT_AT_TA_bot1000_chrV_library = apply(dat_bot1000_chrV_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position 
# for top and bottom 1000
AA_TT_AT_TA_1000_chrV_library = data.frame(top = AA_TT_AT_TA_top1000_chrV_library,
                                           bot = AA_TT_AT_TA_bot1000_chrV_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position
chrV_AATT_topvbot1000_newC0 = ggplot(data = AA_TT_AT_TA_1000_chrV_library, aes(x = 1:49)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in chrV Library, New C0")
chrV_AATT_topvbot1000_newC0
ggsave("figures/periodicity/chrV_AATT_topvbot1000_newC0_v2.png", plot = chrV_AATT_topvbot1000_newC0)


















################################################################################
# yeast Library Quartiles
################################################################################

# Load yeast library training dataset
dat_yeast_library = readRDS("data/Created/processed_yeast_newC0_v2.rds")

# Find cutoffs for first-fourth quartiles of intrinsic cyclizability (C0)
cutoffs_0_25_yeast_library = quantile(dat_yeast_library$C0_new, c(0, 0.25)) 
cutoffs_25_50_yeast_library = quantile(dat_yeast_library$C0_new, c(0.25, 0.5)) 
cutoffs_50_75_yeast_library = quantile(dat_yeast_library$C0_new, c(0.5, 0.75)) 
cutoffs_75_100_yeast_library = quantile(dat_yeast_library$C0_new, c(0.75, 1)) 

# Divide yeast library dataset into quartiles
dat_q1_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_0_25_yeast_library[1] & C0_new <= cutoffs_0_25_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q2_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_25_50_yeast_library[1] & C0_new <= cutoffs_25_50_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q3_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_50_75_yeast_library[1] & C0_new <= cutoffs_50_75_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_q4_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_75_100_yeast_library[1] & C0_new <= cutoffs_75_100_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each quartile
A_T_q1_yeast_library = apply(dat_q1_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q1_yeast_library))
})
A_T_q2_yeast_library = apply(dat_q2_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q2_yeast_library))
})
A_T_q3_yeast_library = apply(dat_q3_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q3_yeast_library))
})
A_T_q4_yeast_library = apply(dat_q4_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_q4_yeast_library))
})
# Construct dataframe containing the relative frequences of A and T by position for all four quartiles
A_T_q_yeast_library = data.frame(q1 = A_T_q1_yeast_library,
                                 q2 = A_T_q2_yeast_library,
                                 q3 = A_T_q3_yeast_library,
                                 q4 = A_T_q4_yeast_library)
# Plot the relative frequencies of A and T by position for each quartile (q1 v q4 and q2 v q3)
yeast_AT_q1vq4_newC0 = ggplot(data = A_T_q_yeast_library, aes(x = 1:50)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in yeast Library, New C0")
yeast_AT_q1vq4_newC0
ggsave("figures/periodicity/yeast_AT_q1vq4_newC0_v2.png", plot = yeast_AT_q1vq4_newC0)
yeast_AT_q2vq3_newC0 = ggplot(data = A_T_q_yeast_library, aes(x = 1:50)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in yeast Library, New C0")
yeast_AT_q2vq3_newC0
ggsave("figures/periodicity/yeast_AT_q2vq3_newC0_v2.png", plot = yeast_AT_q2vq3_newC0)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each quartile
AA_TT_AT_TA_q1_yeast_library = apply(dat_q1_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q1_yeast_library))
})
AA_TT_AT_TA_q2_yeast_library = apply(dat_q2_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q2_yeast_library))
})
AA_TT_AT_TA_q3_yeast_library = apply(dat_q3_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q3_yeast_library))
})
AA_TT_AT_TA_q4_yeast_library = apply(dat_q4_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_q4_yeast_library))
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for all four quartiles
AA_TT_AT_TA_q_yeast_library = data.frame(q1 = AA_TT_AT_TA_q1_yeast_library,
                                         q2 = AA_TT_AT_TA_q2_yeast_library,
                                         q3 = AA_TT_AT_TA_q3_yeast_library,
                                         q4 = AA_TT_AT_TA_q4_yeast_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
yeast_AATT_q1vq4_newC0 = ggplot(data = AA_TT_AT_TA_q_yeast_library, aes(x = 1:49)) +
  geom_line(aes(y=q1, color="0-25%")) +
  geom_line(aes(y=q4, color="75-100%")) +
  scale_color_manual(breaks = c("0-25%", "75-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in yeast Library, New C0")
yeast_AATT_q1vq4_newC0
ggsave("figures/periodicity/yeast_AATT_q1vq4_newC0_v2.png", plot = yeast_AATT_q1vq4_newC0)
yeast_AATT_q2vq3_newC0 = ggplot(data = AA_TT_AT_TA_q_yeast_library, aes(x = 1:49)) +
  geom_line(aes(y=q2, color="25-50%")) +
  geom_line(aes(y=q3, color="50-75%")) +
  scale_color_manual(breaks = c("25-50%", "50-75%"),
                     values = c("green", "orange"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in yeast Library, New C0")
yeast_AATT_q2vq3_newC0
ggsave("figures/periodicity/yeast_AATT_q2vq3_newC0_v2.png", plot = yeast_AATT_q2vq3_newC0)















################################################################################
# yeast Library Top 1000
################################################################################

# Divide yeast library dataset into top and bottom 1000
dat_top1000_yeast_library = dat_yeast_library %>%
  slice_max(C0_new, n=1000)
dat_bot1000_yeast_library = dat_yeast_library %>%
  slice_min(C0_new, n=1000)

# Find the relative frequencies of A and T at each position (1-50) for top and bottom 1000
A_T_top1000_yeast_library = apply(dat_top1000_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
A_T_bot1000_yeast_library = apply(dat_bot1000_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/1000)
})
# Construct dataframe containing the relative frequences of A and T by position for top and bottom 1000
A_T_1000_yeast_library = data.frame(top = A_T_top1000_yeast_library,
                                    bot = A_T_bot1000_yeast_library)
# Plot the relative frequencies of A and T by position
yeast_AT_topvbot1000_newC0 = ggplot(data = A_T_1000_yeast_library, aes(x = 1:50)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in yeast Library, New C0")
yeast_AT_topvbot1000_newC0
ggsave("figures/periodicity/yeast_AT_topvbot1000_newC0_v2.png", plot = yeast_AT_topvbot1000_newC0)



# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for top and bottom 1000
AA_TT_AT_TA_top1000_yeast_library = apply(dat_top1000_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
AA_TT_AT_TA_bot1000_yeast_library = apply(dat_bot1000_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/1000)
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position
# for top and bottom 1000
AA_TT_AT_TA_1000_yeast_library = data.frame(top = AA_TT_AT_TA_top1000_yeast_library,
                                            bot = AA_TT_AT_TA_bot1000_yeast_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position
yeast_AATT_topvbot1000_newC0 = ggplot(data = AA_TT_AT_TA_1000_yeast_library, aes(x = 1:49)) +
  geom_line(aes(y=top, color="Top 1000")) +
  geom_line(aes(y=bot, color="Bottom 1000")) +
  scale_color_manual(breaks = c("Bottom 1000", "Top 1000"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in yeast Library, New C0")
yeast_AATT_topvbot1000_newC0
ggsave("figures/periodicity/yeast_AATT_topvbot1000_newC0_v2.png", plot = yeast_AATT_topvbot1000_newC0)



















################################################################################
# yeast Library Halves
################################################################################

# Find cutoffs for first and second halves of intrinsic cyclizability (C0)
cutoffs_0_50_yeast_library = quantile(dat_yeast_library$C0_new, c(0, 0.5)) 
cutoffs_50_100_yeast_library = quantile(dat_yeast_library$C0_new, c(0.5, 1)) 

# Divide yeast library dataset into halves
dat_h1_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_0_50_yeast_library[1] & C0_new <= cutoffs_0_50_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)
dat_h2_yeast_library = dat_yeast_library %>%
  filter(C0_new >= cutoffs_50_100_yeast_library[1] & C0_new <= cutoffs_50_100_yeast_library[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

# Find the relative frequencies of A and T at each position (1-50) for each half
A_T_h1_yeast_library = apply(dat_h1_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_h1_yeast_library))
})
A_T_h2_yeast_library = apply(dat_h2_yeast_library %>% select(all_of(ps1)), 2, function(col) {
  A_freq = sum(col == "A")
  T_freq = sum(col == "T")
  return((A_freq + T_freq)/nrow(dat_h2_yeast_library))
})
# Construct dataframe containing the relative frequences of A and T by position for both halves
A_T_h_yeast_library = data.frame(h1 = A_T_h1_yeast_library,
                                 h2 = A_T_h2_yeast_library)
# Plot the relative frequencies of A and T by position for each half
yeast_AT_topvbothalf = ggplot(data = A_T_h_yeast_library, aes(x = 1:50)) +
  geom_line(aes(y=h1, color="0-50%")) +
  geom_line(aes(y=h2, color="50-100%")) +
  scale_color_manual(breaks = c("0-50%", "50-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Nucleotide Position") +
  ylab("Average A/T Content") +
  labs(title = "Average A/T Content by Position in yeast Library")
yeast_AT_topvbothalf
ggsave("figures/periodicity/yeast_AT_topvbothalf_newC0_v2.png", plot = yeast_AT_topvbothalf)


# Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each half
AA_TT_AT_TA_h1_yeast_library = apply(dat_h1_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_h1_yeast_library))
})
AA_TT_AT_TA_h2_yeast_library = apply(dat_h2_yeast_library %>% select(all_of(ps2)), 2, function(col) {
  AA_freq = sum(col == "AA")
  TT_freq = sum(col == "TT")
  AT_freq = sum(col == "AT")
  TA_freq = sum(col == "TA")
  return((AA_freq + TT_freq + AT_freq + TA_freq)/nrow(dat_h2_yeast_library))
})
# Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position for both halves
AA_TT_AT_TA_h_yeast_library = data.frame(h1 = AA_TT_AT_TA_h1_yeast_library,
                                         h2 = AA_TT_AT_TA_h2_yeast_library)
# Plot the relative frequencies of AA, TT, AT, and TA by position for each quartile (q1 v q4 and q2 v q3)
yeast_AATT_topvbothalf = ggplot(data = AA_TT_AT_TA_h_yeast_library, aes(x = 1:49)) +
  geom_line(aes(y=h1, color="0-50%")) +
  geom_line(aes(y=h2, color="50-100%")) +
  scale_color_manual(breaks = c("0-50%", "50-100%"),
                     values = c("blue", "red"),
                     name="C0 Cutoffs") +
  xlab("Dinucleotide Position") +
  ylab("Average AA/TT/AT/TA Content") +
  labs(title = "Average AA/TT/AT/TA Content by Position in yeast Library")
yeast_AATT_topvbothalf
ggsave("figures/periodicity/yeast_AATT_topvbothalf_newC0_v2.png", plot = yeast_AATT_topvbothalf)

