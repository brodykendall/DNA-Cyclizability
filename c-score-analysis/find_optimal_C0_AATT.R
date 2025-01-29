library(tidyverse)
# Default Tether positions:
# n <- c(26, 29, 31)

# Potential new phase-shift positions (specifically correcting for n=29):
# n <- c(26, 33, 31)
# or
# n <- c(26, 24, 31)
# or 
# n <- c(26, 29+(10.4/2), 31)
# or
# n <- c(26, 29-(10.4/2), 31)

# Ratios provided by Basu:
# aa <- c(1, 0.82, 0.7)

# Ratios actually used by Basu:
# aa <- c(1, 1/0.82, 1/0.7)

# Period
# k <- 2*pi/10.4

# dat = readRDS("data/Created/processed_tiling_newC0.rds")
# dat_test = readRDS("data/Created/processed_tiling_test_newC0.rds")
# random = readRDS("data/Created/processed_random_newC0.rds")
# random_test = readRDS("data/Created/processed_random_test_newC0.rds")

##### Based on AA/TT/AT/TA:

ps2 <- paste0("X", 1:49, "di")

find_topm_fourier_power_AATT = function(aa, m, n, data) {
  temp_dat_train = data %>%
    select("x50mer","C26","C29","C31","C0","Amplitude","Phase", all_of(ps2))
  
  temp_dat_train[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(temp_dat_train[,c("C26","C29","C31")], n=n, aa=aa)
  
  ### Want to find Fourier power for AA/TT/AT/TA 9.8bp periodicity for a grid of values of f
  # for the m sequences with the highest C0_new values
  temp_dat_train_topm = temp_dat_train %>%
    top_n(m, C0_new)
  
  Xtwo_topm = temp_dat_train_topm %>% select(all_of(ps2))
  Xtwo_AorT_topm = matrix(nrow=nrow(Xtwo_topm), ncol=ncol(Xtwo_topm))
  colnames(Xtwo_AorT_topm) = colnames(Xtwo_topm)
  Xtwo_AorT_topm[] = ((Xtwo_topm == "AA") | (Xtwo_topm == "TT") | (Xtwo_topm == "AT") | (Xtwo_topm == "TA")) %>% 
    as.matrix() %>% as.numeric()
  
  # Collapse all m sequences (This is what Basu does on Page 12 of Supplementary 
  # Information - "Refinements of the model"):
  Xtwo_AorT_topm = apply(Xtwo_AorT_topm, 2, mean)
  
  temp_fourier = fft(Xtwo_AorT_topm)[6]
  return(Re(temp_fourier)^2 + Im(temp_fourier)^2)
}

# Example using values from the paper, 
# (with the exception of fourier power based on AA/TT/AT/TA instead of A/T)
find_topm_fourier_power_AATT(aa=c(1, 1/0.82, 1/0.7), m=1000, n=c(26,29,31), data=random)



find_optimal_aa_AATT = function(data, library_name, m=1000, n=c(26,29,31),
                           a2_values=seq(from=0.1,to=3,by=0.1), 
                           a3_values=seq(from=0.1,to=2.5,by=0.1)) {
  
  aa_grid = map2(rep(a2_values, length(a3_values)), rep(a3_values, each=length(a2_values)),~c(1, .x, .y))
  
  power_vector = aa_grid %>%
    map(~find_topm_fourier_power_AATT(aa=.x, m=m, n=n, data=data)) %>%
    unlist() 
  
  power_grid = power_vector %>%
    matrix(nrow=length(a2_values), ncol=length(a3_values))
  
  # Rows (a2) correspond to (n=29)/(n=26)
  # Columns (a3) correspond to (n=31)/(n=26)
  rownames(power_grid) = a2_values
  colnames(power_grid) = a3_values
  
  # Heatmap:
  gg_format = expand.grid(a2=a2_values, a3=a3_values)
  gg_format$weight = power_vector
  heat = ggplot(gg_format, aes(a2, a3)) + 
    geom_tile(aes(fill=weight)) +
    scale_fill_gradient2("Fourier Power", low="blue", high="red", midpoint=mean(power_vector)) +
    ggtitle(paste0("Optimized Using top ", m, " sequences from ", library_name, " Library"))
  
  min_row = gg_format[min(which(gg_format$weight == min(gg_format$weight))),]
  
  aa_optimal = c(1, min_row$a2, min_row$a3)
  
  to_return = list()
  to_return$power_grid = power_grid
  to_return$heatmap = heat
  to_return$aa_optimal = aa_optimal
  
  return(to_return)
}

# Optimize Based on Top Half of Random Library:
optimal_random_tophalf_prelim_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2)
optimal_random_tophalf_prelim_AATT$aa_optimal
# 1.0 2.3 1.5
optimal_random_tophalf_prelim_AATT$heatmap
optimal_random_tophalf_refined_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2,
                                                 a2_values=seq(from=2.1,to=2.5,by=0.01),
                                                 a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_random_tophalf_refined_AATT$aa_optimal
# 1.00 2.36 1.47
optimal_random_tophalf_refined_AATT$heatmap



# Optimize Based on Top Quartile of Random Library:
optimal_random_topquartile_prelim_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4)
optimal_random_topquartile_prelim_AATT$aa_optimal
# 1.0 1.8 1.3
optimal_random_topquartile_prelim_AATT$heatmap
optimal_random_topquartile_refined_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4,
                                                     a2_values=seq(from=1.6,to=2.0,by=0.01),
                                                     a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_topquartile_refined_AATT$aa_optimal
# 1.00 1.94 1.34
optimal_random_topquartile_refined_AATT$heatmap



# Optimize Based on Top 1000 of Random Library:
optimal_random_top1000_prelim_AATT = find_optimal_aa_AATT(random, "Random")
optimal_random_top1000_prelim_AATT$aa_optimal
# 1.0 1.1 1.3
optimal_random_top1000_prelim_AATT$heatmap
optimal_random_top1000_refined_AATT = find_optimal_aa_AATT(random, "Random",
                                                 a2_values=seq(from=0.9,to=1.3,by=0.01),
                                                 a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_top1000_refined_AATT$aa_optimal
# 1.00 1.04 1.26
optimal_random_top1000_refined_AATT$heatmap



# Optimize Based on Top Half of Random Library, new n:
# n=c(26,33,31)
n=c(26, 29-(10.4/2), 31)
optimal_random_tophalf_prelim_new_n_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2, n=n,
                                                      a2_values=seq(from=5,to=10,by=0.1))
optimal_random_tophalf_prelim_new_n_AATT$aa_optimal
# 1.0 9.1 1.5
optimal_random_tophalf_prelim_new_n_AATT$heatmap
optimal_random_tophalf_refined_new_n_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2, n=n,
                                                       a2_values=seq(from=8.8,to=9.2,by=0.01),
                                                       a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_random_tophalf_refined_new_n_AATT$aa_optimal
# 1.00 8.88 1.53
optimal_random_tophalf_refined_new_n_AATT$heatmap



# Optimize Based on Top Quartile of Random Library, new n:
optimal_random_topquartile_prelim_new_n_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4, n=n,
                                                          a2_values=seq(from=5,to=10,by=0.1))
optimal_random_topquartile_prelim_new_n_AATT$aa_optimal
# 1.0 9.4 1.4
optimal_random_topquartile_prelim_new_n_AATT$heatmap
optimal_random_topquartile_refined_new_n_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4, n=n,
                                                           a2_values=seq(from=9.2,to=9.6,by=0.01),
                                                           a3_values=seq(from=1.2,to=1.6,by=0.01))
optimal_random_topquartile_refined_new_n_AATT$aa_optimal
# 1.0 9.4 1.4
optimal_random_topquartile_refined_new_n_AATT$heatmap



# Optimize Based on Top 1000 of Random Library, new n:
optimal_random_top1000_prelim_new_n_AATT = find_optimal_aa_AATT(random, "Random", n=n,
                                                      a2_values=seq(from=5,to=10,by=0.1))
optimal_random_top1000_prelim_new_n_AATT$aa_optimal
# 1.0 7.8 1.5
optimal_random_top1000_prelim_new_n_AATT$heatmap
optimal_random_top1000_refined_new_n_AATT = find_optimal_aa_AATT(random, "Random", n=n,
                                                       a2_values=seq(from=7.6,to=8.0,by=0.01),
                                                       a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_random_top1000_refined_new_n_AATT$aa_optimal
# 1.00 7.76 1.50
optimal_random_top1000_refined_new_n_AATT$heatmap




# Optimize Based on Top Half of Random Library, original n, negative A26:
optimal_random_tophalf_prelim_neg_a2_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2,
                                                       a2_values=seq(from=-10,to=-5,by=0.1))
optimal_random_tophalf_prelim_neg_a2_AATT$aa_optimal
# 1.0 -9.3  1.5
optimal_random_tophalf_prelim_neg_a2_AATT$heatmap
optimal_random_tophalf_refined_neg_a2_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/2,
                                                        a2_values=seq(from=2.1,to=2.5,by=0.01),
                                                        a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_random_tophalf_refined_neg_a2_AATT$aa_optimal
# 
optimal_random_tophalf_refined_neg_a2_AATT$heatmap



# Optimize Based on Top Quartile of Random Library, original n, negative A26:
optimal_random_topquartile_prelim_neg_a2_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4,
                                                           a2_values=seq(from=-2.5,to=-0.1,by=0.1))
optimal_random_topquartile_prelim_neg_a2_AATT$aa_optimal
# 
optimal_random_topquartile_prelim_neg_a2_AATT$heatmap
optimal_random_topquartile_refined_neg_a2_AATT = find_optimal_aa_AATT(random, "Random", nrow(random)/4,
                                                            a2_values=seq(from=1.6,to=2.0,by=0.01),
                                                            a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_topquartile_refined_neg_a2_AATT$aa_optimal
# 
optimal_random_topquartile_refined_neg_a2_AATT$heatmap



# Optimize Based on Top 1000 of Random Library, original n, negative A26:
optimal_random_top1000_prelim_neg_a2_AATT = find_optimal_aa_AATT(random, "Random",
                                                       a2_values=seq(from=-2.5,to=-0.1,by=0.1))
optimal_random_top1000_prelim_neg_a2_AATT$aa_optimal
# 
optimal_random_top1000_prelim_neg_a2_AATT$heatmap
optimal_random_top1000_refined_neg_a2_AATT = find_optimal_aa_AATT(random, "Random",
                                                        a2_values=seq(from=0.9,to=1.3,by=0.01),
                                                        a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_top1000_refined_neg_a2_AATT$aa_optimal
# 
optimal_random_top1000_refined_neg_a2_AATT$heatmap


# TILING:

# Optimize Based on Top Half of Tiling Library:
optimal_tiling_tophalf_prelim_AATT = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/2))
optimal_tiling_tophalf_prelim_AATT$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_tophalf_prelim_AATT$heatmap
optimal_tiling_tophalf_refined_AATT = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/2),
                                                 a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                 a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_tophalf_refined_AATT$aa_optimal
# 1.00 0.36 1.19
optimal_tiling_tophalf_refined_AATT$heatmap


# Optimize Based on Top Quartile of Tiling Library:
optimal_tiling_topquartile_prelim_AATT = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/4))
optimal_tiling_topquartile_prelim_AATT$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_topquartile_prelim_AATT$heatmap
optimal_tiling_topquartile_refined_AATT = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/4),
                                                     a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                     a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_topquartile_refined_AATT$aa_optimal
# 1.00 0.37 1.15
optimal_tiling_topquartile_refined_AATT$heatmap


# Optimize Based on Top 1000 of Tiling Library:
optimal_tiling_top1000_prelim = find_optimal_aa_AATT(tiling, "Tiling")
optimal_tiling_top1000_prelim$aa_optimal
# 1.0 0.6 1.2
optimal_tiling_top1000_prelim$heatmap
optimal_tiling_top1000_refined = find_optimal_aa_AATT(tiling, "Tiling", 
                                                 a2_values=seq(from=0.4,to=0.8,by=0.01),
                                                 a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_top1000_refined$aa_optimal
# 1.00 0.57 1.14
optimal_tiling_top1000_refined$heatmap


# Optimize Based on Top Half of Tiling Library, new n:
optimal_tiling_tophalf_prelim_new_n = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/2), n=n)
optimal_tiling_tophalf_prelim_new_n$aa_optimal
# 1.0 2.5 1.5
optimal_tiling_tophalf_prelim_new_n$heatmap
optimal_tiling_tophalf_refined_new_n = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/2), n=n,
                                                       a2_values=seq(from=2.3,to=2.7,by=0.01),
                                                       a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_tiling_tophalf_refined_new_n$aa_optimal
# 1.00 2.69 1.50
optimal_tiling_tophalf_refined_new_n$heatmap


# Optimize Based on Top Quartile of Tiling Library, new n:
optimal_tiling_topquartile_prelim_new_n = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/4), n=n,
                                                          a2_values=seq(from=5,to=10,by=0.1))
optimal_tiling_topquartile_prelim_new_n$aa_optimal
# 1.0 9.7 1.4
optimal_tiling_topquartile_prelim_new_n$heatmap
optimal_tiling_topquartile_refined_new_n = find_optimal_aa_AATT(tiling, "Tiling", floor(nrow(tiling)/4), n=n,
                                                           a2_values=seq(from=9.5,to=9.9,by=0.01),
                                                           a3_values=seq(from=1.2,to=1.6,by=0.01))
optimal_tiling_topquartile_refined_new_n$aa_optimal
# 1.00 9.78 1.39
optimal_tiling_topquartile_refined_new_n$heatmap


# Optimize Based on Top 1000 of Tiling Library, new n:
optimal_tiling_top1000_prelim_new_n = find_optimal_aa_AATT(tiling, "Tiling", n=n)
optimal_tiling_top1000_prelim_new_n$aa_optimal
# 1.0 2.5 1.4
optimal_tiling_top1000_prelim_new_n$heatmap
optimal_tiling_top1000_refined_new_n = find_optimal_aa_AATT(tiling, "Tiling", n=n,
                                                       a2_values=seq(from=2.3,to=2.7,by=0.01),
                                                       a3_values=seq(from=1.2,to=1.6,by=0.01))
optimal_tiling_top1000_refined_new_n$aa_optimal
# 1.00 2.70 1.39
optimal_tiling_top1000_refined_new_n$heatmap




# a3 (The optimized ratio of A(n=31) to A(n=26)) is generally quite constant 
# across all the different sets we optimize over, usually around 1.2-1.3
# However, a2 (the optimized ratio of A(n=29) to A(n=26)) is much more variable,
# but most noticeably between libraries. The value ranges from 1.04-2.36 for the Random Library,
# and from 0.36-0.57 for the Tiling Library. Interestingly, with larger sample sizes
# for the Random Library we see increasing a2 values (>1), but with larger sample sizes for the
# Tiling Library we see decreasing a2 values (<1).

# What if we base analysis/optimization on a cutoff value for C0 instead? I think there is
# a higher proportion and magnitude of highly cyclizable sequences in the Tiling Library
# compared to the Random Library

# Could also potentially try keeping the original values of n, but allowing for a negative A29.

# Let's look into the Fourier power at various phases for n=29


# Combined Tiling and Random Libraries:
tiling_and_random = rbind(tiling %>%
                                select(all_of(intersect(colnames(tiling), colnames(random)))), 
                              random %>%
                                select(all_of(intersect(colnames(tiling), colnames(random)))))

optimal_tiling_and_random_tophalf_prelim = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", floor(nrow(tiling_and_random)/2))
optimal_tiling_and_random_tophalf_prelim$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_and_random_tophalf_prelim$heatmap
optimal_tiling_and_random_tophalf_refined = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", floor(nrow(tiling_and_random)/2),
                                                            a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                            a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_and_random_tophalf_refined$aa_optimal
# 1.00 0.39 1.19
optimal_tiling_and_random_tophalf_refined$heatmap




optimal_tiling_and_random_topquartile_prelim = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", floor(nrow(tiling_and_random)/4))
optimal_tiling_and_random_topquartile_prelim$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_and_random_topquartile_prelim$heatmap
optimal_tiling_and_random_topquartile_refined = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", floor(nrow(tiling_and_random)/4),
                                                                a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                                a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_and_random_topquartile_refined$aa_optimal
# 1.00 0.41 1.16
optimal_tiling_and_random_topquartile_refined$heatmap




optimal_tiling_and_random_top2000_prelim = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", 2000)
optimal_tiling_and_random_top2000_prelim$aa_optimal
# 1.0 0.5 1.1
optimal_tiling_and_random_top2000_prelim$heatmap
optimal_tiling_and_random_top2000_refined = find_optimal_aa_AATT(tiling_and_random, "Tiling and Random", 2000,
                                                            a2_values=seq(from=0.3,to=0.7,by=0.01),
                                                            a3_values=seq(from=0.9,to=1.3,by=0.01))
optimal_tiling_and_random_top2000_refined$aa_optimal
# 1.00 0.51 1.09
optimal_tiling_and_random_top2000_refined$heatmap




# TODO: What if we minimize the sum of the fourier powers of each library?




# plot_periodicity = function(top_v_bottom, data, library_name, aa_optimal) {
#   
#   temp_dat = data %>%
#     select("x50mer","C26","C29","C31","C0","Amplitude","Phase", all_of(ps2))
#   
#   temp_dat[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(temp_dat, aa_optimal)
#   
#   temp_dat_top = temp_dat %>%
#     top_n(top_v_bottom, C0_new)
#   temp_dat_bot = temp_dat %>%
#     top_n(-top_v_bottom, C0_new)
#   
#   # Find the relative frequencies of AA, TT, AT, and TA at each position (1-49) for each subset
#   AA_TT_AT_TA_top = apply(temp_dat_top %>% select(all_of(ps2)), 2, function(col) {
#     AA_freq = sum(col == "AA")
#     TT_freq = sum(col == "TT")
#     AT_freq = sum(col == "AT")
#     TA_freq = sum(col == "TA")
#     return((AA_freq + TT_freq + AT_freq + TA_freq)/top_v_bottom)
#   })
#   AA_TT_AT_TA_bot = apply(temp_dat_bot %>% select(all_of(ps2)), 2, function(col) {
#     AA_freq = sum(col == "AA")
#     TT_freq = sum(col == "TT")
#     AT_freq = sum(col == "AT")
#     TA_freq = sum(col == "TA")
#     return((AA_freq + TT_freq + AT_freq + TA_freq)/top_v_bottom)
#   })
#   
#   # Construct dataframe containing the relative frequences of AA, TT, AT, and TA by position 
#   AA_TT_AT_TA = data.frame(top = AA_TT_AT_TA_top,
#                            bot = AA_TT_AT_TA_bot)
#   # Plot the relative frequencies of AA, TT, AT, and TA by position
#   return(ggplot(data = AA_TT_AT_TA, aes(x = 1:49)) +
#            geom_line(aes(y=top, color=paste0("Top ", top_v_bottom))) +
#            geom_line(aes(y=bot, color=paste0("Bottom ", top_v_bottom))) +
#            scale_color_manual(breaks = c(paste0("Bottom ", top_v_bottom), paste0("Top ", top_v_bottom)),
#                               values = c("blue", "red"),
#                               name="C0 Cutoffs") +
#            xlab("Dinucleotide Position") +
#            ylab("Average AA/TT/AT/TA Content") +
#            labs(title = paste0("AA/TT/AT/TA Content in ", library_name, " Library, Optimized C0"))
#   )
# }
# 
# 
# # C0 optimized using the top half of the Random Library:
# plot_periodicity(1000, random, "Random", optimal_random_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_random_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_random_tophalf_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_random_tophalf_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_random_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_random_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_random_tophalf_refined$aa_optimal)
# 
# # plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_random_tophalf_refined$aa_optimal)
# 
# 
# 
# 
# 
# # C0 optimized using the top quartile of the Random Library:
# plot_periodicity(1000, random, "Random", optimal_random_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_random_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_random_topquartile_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_random_topquartile_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_random_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_random_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_random_topquartile_refined$aa_optimal)
# 
# # plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_random_topquartile_refined$aa_optimal)
# 
# 
# 
# 
# 
# # C0 optimized using the top 1000 of the Random Library:
# plot_periodicity(1000, random, "Random", optimal_random_top1000_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_random_top1000_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_random_top1000_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_random_top1000_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_random_top1000_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_random_top1000_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_random_top1000_refined$aa_optimal)
# 
# # plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_random_top1000_refined$aa_optimal)
# 
# 
# 
# 
# 
# # C0 optimized using the top half of the Tiling Library:
# plot_periodicity(1000, random, "Random", optimal_tiling_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_tiling_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_tiling_tophalf_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_tiling_tophalf_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_tiling_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_tiling_tophalf_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_tiling_tophalf_refined$aa_optimal)
# 
# # plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_tiling_tophalf_refined$aa_optimal)
# 
# 
# 
# 
# 
# # C0 optimized using the top quartile of the Tiling Library:
# plot_periodicity(1000, random, "Random", optimal_tiling_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_tiling_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_tiling_topquartile_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_tiling_topquartile_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_tiling_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_tiling_topquartile_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_tiling_topquartile_refined$aa_optimal)
# 
# # plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_tiling_topquartile_refined$aa_optimal)
# 
# 
# 
# 
# 
# # C0 optimized using the top 1000 of the Tiling Library:
# plot_periodicity(1000, random, "Random", optimal_tiling_top1000_refined$aa_optimal)
# plot_periodicity(nrow(random)/4, random, "Random", optimal_tiling_top1000_refined$aa_optimal)
# plot_periodicity(nrow(random)/2, random, "Random", optimal_tiling_top1000_refined$aa_optimal)
# 
# plot_periodicity(nrow(random_test)/4, random_test, "Random Test", optimal_tiling_top1000_refined$aa_optimal)
# 
# plot_periodicity(1000, tiling, "Tiling", optimal_tiling_top1000_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/4, tiling, "Tiling", optimal_tiling_top1000_refined$aa_optimal)
# plot_periodicity(nrow(tiling)/2, tiling, "Tiling", optimal_tiling_top1000_refined$aa_optimal)

# plot_periodicity(nrow(dat_test)/4, dat_test, "Tiling Test", optimal_tiling_top1000_refined$aa_optimal)




# General qualitative conclusions:
# It seems like optimizing on either library is successfully at collapsing the phase on itself,
# however, it does not successfully do so on the other library.
# I would argue that this effect is more noticeable when optimizing on the Random Library and
# applying the it to the Tiling Library - leading me to conclude that optimizing on 
# the Tiling Library is more effective.
# Of the optimizations we tried on the Tiling Library, the one based on the top quartile
# of sequences seems to be the best.






