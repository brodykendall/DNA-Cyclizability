ps1 <- paste0("X", 1:50, "mono")

find_ps1_fourier = function(var, data) {
  data = data %>% select(all_of(ps1))
  Xone_AorT_data = matrix(nrow=nrow(data), ncol=ncol(data))
  colnames(Xone_AorT_data) = colnames(data)
  Xone_AorT_data[] = ((data == "A") | (data == "T")) %>% 
    as.matrix() %>% as.numeric()
  
  # Collapse all sequences (This is what Basu does on Page 12 of Supplementary 
  # Information - "Refinements of the model"):
  Xone_AorT_data = apply(Xone_AorT_data, 2, mean)
  
  temp_fourier = fft(Xone_AorT_data)[6]
  return(temp_fourier)
}

find_iterative_fourier_ps1 = function(var, group_size=100, data) {
  temp_dat_train = data %>%
    select("x50mer", {{var}}, all_of(ps1))
  
  temp_dat_train = temp_dat_train %>%
    arrange(desc({{var}}))
  
  group = rep(1:(nrow(temp_dat_train)/group_size), each=group_size)
  group = c(group, rep(max(group)+1, times=nrow(temp_dat_train)%%group_size))
  temp_dat_train$group = group
  
  fourier_list = temp_dat_train %>%
    group_by(group) %>%
    group_map(~find_ps1_fourier(var, .x)) %>%
    unlist()
  
  mean_var_df = temp_dat_train %>%
    group_by(group) %>%
    summarize(mean_var = mean({{var}}))
  
  out = data.frame(fourier = fourier_list, mean_var = mean_var_df$mean_var)
  
  return(out)
}


group_size = 100
# group_size = 1

## Nucleosome:
nuc_C26_iterative_fourier_AT = find_iterative_fourier_ps1(C26, group_size, nuc)
nuc_C29_iterative_fourier_AT = find_iterative_fourier_ps1(C29, group_size, nuc)
nuc_C31_iterative_fourier_AT = find_iterative_fourier_ps1(C31, group_size, nuc)

# plot(Mod(nuc_C26_iterative_fourier_AT$fourier))
# plot(Mod(nuc_C29_iterative_fourier_AT$fourier))
# plot(Mod(nuc_C31_iterative_fourier_AT$fourier))


## Random:
random_C26_iterative_fourier_AT = find_iterative_fourier_ps1(C26, group_size, random)
random_C29_iterative_fourier_AT = find_iterative_fourier_ps1(C29, group_size, random)
random_C31_iterative_fourier_AT = find_iterative_fourier_ps1(C31, group_size, random)

png(filename="figures/c-score-analysis/random_C26_fourier_amplitude_group_size_100_AT.png")
plot(Mod(random_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Amplitude of AT")
title("Random Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C29_fourier_amplitude_group_size_100_AT.png")
plot(Mod(random_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Amplitude of AT")
title("Random Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C31_fourier_amplitude_group_size_100_AT.png")
plot(Mod(random_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Amplitude of AT")
title("Random Library Grouped by Increasing C31, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/random_C26_fourier_phase_group_size_100_AT.png")
plot(Arg(random_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C29_fourier_phase_group_size_100_AT.png")
plot(Arg(random_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C31_fourier_phase_group_size_100_AT.png")
plot(Arg(random_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C31, Group Size 100")
dev.off()

## Tiling:
tiling_C26_iterative_fourier_AT = find_iterative_fourier_ps1(C26, group_size, tiling)
tiling_C29_iterative_fourier_AT = find_iterative_fourier_ps1(C29, group_size, tiling)
tiling_C31_iterative_fourier_AT = find_iterative_fourier_ps1(C31, group_size, tiling)

png(filename="figures/c-score-analysis/tiling_C26_fourier_amplitude_group_size_100_AT.png")
plot(Mod(tiling_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Amplitude of AT")
title("Tiling Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C29_fourier_amplitude_group_size_100_AT.png")
plot(Mod(tiling_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Amplitude of AT")
title("Tiling Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C31_fourier_amplitude_group_size_100_AT.png")
plot(Mod(tiling_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Amplitude of AT")
title("Tiling Library Grouped by Increasing C31, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/tiling_C26_fourier_phase_group_size_100_AT.png")
plot(Arg(tiling_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Phase of AT")
title("Tiling Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C31_fourier_phase_group_size_100_AT.png")
plot(Arg(tiling_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Phase of AT")
title("Tiling Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C31_fourier_phase_group_size_100_AT.png")
plot(Arg(tiling_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Phase of AT")
title("Tiling Library Grouped by Increasing C31, Group Size 100")
dev.off()

## ChrV:
chrv_C26_iterative_fourier_AT = find_iterative_fourier_ps1(C26, group_size, chrv)
chrv_C29_iterative_fourier_AT = find_iterative_fourier_ps1(C29, group_size, chrv)
chrv_C31_iterative_fourier_AT = find_iterative_fourier_ps1(C31, group_size, chrv)

png(filename="figures/c-score-analysis/chrv_C26_fourier_amplitude_group_size_100_AT.png")
plot(Mod(chrv_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Amplitude of AT")
title("ChrV Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C29_fourier_amplitude_group_size_100_AT.png")
plot(Mod(chrv_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Amplitude of AT")
title("ChrV Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C31_fourier_amplitude_group_size_100_AT.png")
plot(Mod(chrv_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Amplitude of AT")
title("ChrV Library Grouped by Increasing C31, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/chrv_C26_fourier_phase_group_size_100_AT.png")
plot(Arg(chrv_C26_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C26", 
     ylab="Phase of AT")
title("ChrV Library Grouped by Increasing C26, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C29_fourier_phase_group_size_100_AT.png")
plot(Arg(chrv_C29_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C29", 
     ylab="Phase of AT")
title("ChrV Library Grouped by Increasing C29, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C31_fourier_phase_group_size_100_AT.png")
plot(Arg(chrv_C31_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C31", 
     ylab="Phase of AT")
title("ChrV Library Grouped by Increasing C31, Group Size 100")
dev.off()









##### Sort by Amplitude rather than by C26/C29/C31:

find_iterative_fourier_ps1_sorted_amp = function(var, data, group_size=1) {
  temp_dat_train = data %>%
    select("x50mer", {{var}}, all_of(ps1))
  
  group = rep(1:(nrow(temp_dat_train)/group_size), each=group_size)
  group = c(group, rep(max(group)+1, times=nrow(temp_dat_train)%%group_size))
  temp_dat_train$group = group
  
  fourier_list = temp_dat_train %>%
    group_by(group) %>%
    group_map(~find_ps1_fourier(var, .x)) %>%
    unlist()
  
  temp_dat_train = temp_dat_train %>%
    mutate(fourier_amp = sapply(fourier_list, Mod),
           fourier_phase = sapply(fourier_list, Arg)) %>%
    arrange(desc(fourier_amp)) %>%
    select(-"x50mer", -all_of(ps1), -group)
  
  # temp_dat_train = temp_dat_train %>%
  #   arrange(desc({{var}}))
  
  # group = rep(1:(nrow(temp_dat_train)/group_size), each=group_size)
  # group = c(group, rep(max(group)+1, times=nrow(temp_dat_train)%%group_size))
  # temp_dat_train$group = group
  
  # fourier_list = temp_dat_train %>%
  #   group_by(group) %>%
  #   group_map(~find_ps1_fourier(var, .x)) %>%
  #   unlist()
  # 
  # mean_var_df = temp_dat_train %>%
  #   group_by(group) %>%
  #   summarize(mean_var = mean({{var}}))
  
  # out = data.frame(fourier = fourier_list, mean_var = mean_var_df$mean_var)
  
  return(temp_dat_train)
}

# ## Nucleosome:
# nuc_C26_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C26, nuc)
# nuc_C29_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C29, nuc)
# nuc_C31_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C31, nuc)
# 
# plot(nuc_C26_iterative_fourier_sorted_amp_AT$fourier_amp, 
#      xlab="Index, Sorted by Increasing Fourier Amplitude", 
#      ylab="Amplitude of AT")
# title("Nucleosome Library Indexed by Increasing Fourier Amplitude")
# 
# 
# plot(x=nuc_C26_iterative_fourier_sorted_amp_AT$fourier_amp, 
#      y=nuc_C26_iterative_fourier_sorted_amp_AT$C26,
#      xlab="Index, Sorted by Increasing Fourier Amplitude", 
#      ylab="Amplitude of AT")
# title("Nucleosome Library Indexed by Increasing Fourier Amplitude")
# 
# plot(nuc_C26_iterative_fourier_sorted_amp_AT$C26)
# plot(nuc_C29_iterative_fourier_sorted_amp_AT$C29)
# plot(nuc_C31_iterative_fourier_sorted_amp_AT$C31)
# 
# 
# 
# plot(x=nuc_C26_iterative_fourier_sorted_amp_AT$fourier_phase, 
#      y=nuc_C26_iterative_fourier_sorted_amp_AT$C26,
#      xlab="Index, Sorted by Increasing Fourier Amplitude", 
#      ylab="Phase of AT")
# title("Random Library Grouped by Increasing Fourier Amplitude, Group Size 100")
# 
# plot(nuc_C29_iterative_fourier_sorted_amp_AT$fourier_phase, 
#      nuc_C29_iterative_fourier_sorted_amp_AT$C29,
#      xlab="Index, Sorted by Increasing Fourier Amplitude", 
#      ylab="Phase of AT")
# title("Random Library Grouped by Increasing C29, Group Size 100")
# 
# plot(nuc_C31_iterative_fourier_sorted_amp_AT$fourier_phase, 
#      nuc_C31_iterative_fourier_sorted_amp_AT$C31,
#      xlab="Index, Sorted by Increasing Fourier Amplitude", 
#      ylab="Phase of AT")
# title("Random Library Grouped by Increasing C31, Group Size 100")


## Random:
random_C26_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C26, random)
random_C29_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C29, random)
random_C31_iterative_fourier_sorted_amp_AT = find_iterative_fourier_ps1_sorted_amp(C31, random)

plot(random_C26_iterative_fourier_sorted_amp_AT$fourier_amp, 
     xlab="Index, Sorted by Increasing Fourier Amplitude", 
     ylab="Amplitude of AT")
title("Random Library Indexed by Increasing Fourier Amplitude")


plot(x=random_C26_iterative_fourier_sorted_amp_AT$fourier_amp, 
     y=random_C26_iterative_fourier_sorted_amp_AT$C26,
     xlab="Index, Sorted by Increasing Fourier Amplitude", 
     ylab="Amplitude of AT")
title("Random Library Indexed by Increasing Fourier Amplitude")

plot(random_C26_iterative_fourier_sorted_amp_AT$C26)
plot(random_C29_iterative_fourier_sorted_amp_AT$C29)
plot(random_C31_iterative_fourier_sorted_amp_AT$C31)



plot(x=random_C26_iterative_fourier_sorted_amp_AT$fourier_phase, 
     y=random_C26_iterative_fourier_sorted_amp_AT$C26,
     xlab="Index, Sorted by Increasing Fourier Amplitude", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing Fourier Amplitude, Group Size 100")

plot(random_C29_iterative_fourier_sorted_amp_AT$fourier_phase, 
     random_C29_iterative_fourier_sorted_amp_AT$C29,
     xlab="Index, Sorted by Increasing Fourier Amplitude", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C29, Group Size 100")

plot(random_C31_iterative_fourier_sorted_amp_AT$fourier_phase, 
     random_C31_iterative_fourier_sorted_amp_AT$C31,
     xlab="Index, Sorted by Increasing Fourier Amplitude", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C31, Group Size 100")
