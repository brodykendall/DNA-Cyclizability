ps1 <- paste0("X", 1:50, "mono")

find_ps1_fourier = function(data) {
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
    group_map(~find_ps1_fourier(.x)) %>%
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
nuc_C0_iterative_fourier_AT = find_iterative_fourier_ps1(C0, group_size, nuc)
nuc_C0_new_v12_iterative_fourier_AT = find_iterative_fourier_ps1(C0_new_v12, group_size, nuc)

# # plot(Mod(nuc_C0_iterative_fourier_AT$fourier))
# plot(x=nuc_C0_new_v12_iterative_fourier_AT$mean_var,
#      y=Mod(nuc_C0_new_v12_iterative_fourier_AT$fourier))
# 
# # plot(x=nuc_C0_iterative_fourier_AT$mean_var,
#      # y=Arg(nuc_C0_iterative_fourier_AT$fourier))
# plot(x=nuc_C0_new_v12_iterative_fourier_AT$mean_var,
#      y=Arg(nuc_C0_new_v12_iterative_fourier_AT$fourier))
# 
# plot(Mod(nuc_C0_new_v12_iterative_fourier_AT$fourier))
# 
# plot(Arg(nuc_C0_new_v12_iterative_fourier_AT$fourier))




## Random:
random_C0_iterative_fourier_AT = find_iterative_fourier_ps1(C0, group_size, random)
random_C0_new_v12_iterative_fourier_AT = find_iterative_fourier_ps1(C0_new_v12, group_size, random)

png(filename="figures/c-score-analysis/random_C0_fourier_amplitude_group_size_100_AT.png")
plot(Mod(random_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Amplitude of AT")
title("Random Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C0_new_v12_fourier_amplitude_group_size_100_AT.png")
plot(Mod(random_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Amplitude of AT")
title("Random Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/random_C0_fourier_phase_group_size_100_AT.png")
plot(Arg(random_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/random_C0_new_v12_fourier_phase_group_size_100_AT.png")
plot(Arg(random_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Phase of AT")
title("Random Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()



# plot(Mod(random_C0_iterative_fourier_AT$fourier))
plot(x=random_C0_new_v12_iterative_fourier_AT$mean_var,
     y=Mod(random_C0_new_v12_iterative_fourier_AT$fourier))

# plot(Arg(random_C0_iterative_fourier_AT$fourier))
plot(x=random_C0_new_v12_iterative_fourier_AT$mean_var,
     y=Arg(random_C0_new_v12_iterative_fourier_AT$fourier))

plot(Mod(random_C0_new_v12_iterative_fourier_AT$fourier))



random_C0_iterative_fourier_AT = find_iterative_fourier_ps1(C0, group_size=1, random)

plot(Arg(random_C0_iterative_fourier_AT$fourier))

plot(Mod(random_C0_iterative_fourier_AT$fourier))



## Tiling:
tiling_C0_iterative_fourier_AT = find_iterative_fourier_ps1(C0, group_size, tiling)
tiling_C0_new_v12_iterative_fourier_AT = find_iterative_fourier_ps1(C0_new_v12, group_size, tiling)


png(filename="figures/c-score-analysis/tiling_C0_fourier_amplitude_group_size_100_AT.png")
plot(Mod(tiling_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Amplitude of AT")
title("Tiling Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C0_new_v12_fourier_amplitude_group_size_100_AT.png")
plot(Mod(tiling_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Amplitude of AT")
title("Tiling Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/tiling_C0_fourier_phase_group_size_100_AT.png")
plot(Arg(tiling_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Phase of AT")
title("Tiling Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/tiling_C0_new_v12_fourier_phase_group_size_100_AT.png")
plot(Arg(tiling_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Phase of AT")
title("Tiling Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()


# plot(Mod(tiling_C0_iterative_fourier_AT$fourier))
plot(x=tiling_C0_new_v12_iterative_fourier_AT$mean_var,
     y=Mod(tiling_C0_new_v12_iterative_fourier_AT$fourier))

# plot(Arg(tiling_C0_iterative_fourier_AT$fourier))
plot(x=tiling_C0_new_v12_iterative_fourier_AT$mean_var,
     y=Arg(tiling_C0_new_v12_iterative_fourier_AT$fourier))


plot(Mod(tiling_C0_new_v12_iterative_fourier_AT$fourier))

plot(Arg(tiling_C0_new_v12_iterative_fourier_AT$fourier))




## ChrV:
chrv_C0_iterative_fourier_AT = find_iterative_fourier_ps1(C0, group_size, chrv)
chrv_C0_new_v12_iterative_fourier_AT = find_iterative_fourier_ps1(C0_new_v12, group_size, chrv)

png(filename="figures/c-score-analysis/chrv_C0_fourier_amplitude_group_size_100_AT.png")
plot(Mod(chrv_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Amplitude of AT")
title("ChrV Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C0_new_v12_fourier_amplitude_group_size_100_AT.png")
plot(Mod(chrv_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Amplitude of AT")
title("ChrV Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()

png(filename="figures/c-score-analysis/chrv_C0_fourier_phase_group_size_100_AT.png")
plot(Arg(chrv_C0_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0", 
     ylab="Phase of AT")
title("ChrV Library Grouped by Increasing C0, Group Size 100")
dev.off()
png(filename="figures/c-score-analysis/chrv_C0_new_v12_fourier_phase_group_size_100_AT.png")
plot(Arg(chrv_C0_new_v12_iterative_fourier_AT$fourier), 
     xlab="Group Index, Sorted by Increasing C0_new_v12", 
     ylab="Phase of AT")
title("ChrV Library Grouped by Increasing C0_new_v12, Group Size 100")
dev.off()
