ps2 <- paste0("X", 1:49, "di")

find_ps2_fourier = function(var, data) {
  data = data %>% select(all_of(ps2))
  Xtwo_AorT_data = matrix(nrow=nrow(data), ncol=ncol(data))
  colnames(Xtwo_AorT_data) = colnames(data)
  Xtwo_AorT_data[] = ((data == "AA") | (data == "TT") | (data == "AT") | (data == "TA")) %>% 
    as.matrix() %>% as.numeric()
  
  # Collapse all sequences (This is what Basu does on Page 12 of Supplementary 
  # Information - "Refinements of the model"):
  Xtwo_AorT_data = apply(Xtwo_AorT_data, 2, mean)
  
  temp_fourier = fft(Xtwo_AorT_data)[6]
  return(temp_fourier)
}

find_iterative_fourier_ps2 = function(var, group_size=100, data) {
  temp_dat_train = data %>%
    select("x50mer", {{var}}, all_of(ps2))
  
  temp_dat_train = temp_dat_train %>%
    arrange(desc({{var}}))
  
  group = rep(1:(nrow(temp_dat_train)/group_size), each=group_size)
  group = c(group, rep(max(group)+1, times=nrow(temp_dat_train)%%group_size))
  temp_dat_train$group = group
  
  fourier_list = temp_dat_train %>%
    group_by(group) %>%
    group_map(~find_ps2_fourier(var, .x)) %>%
    unlist()
  
  mean_var_df = temp_dat_train %>%
    group_by(group) %>%
    summarize(mean_var = mean({{var}}))
  
  out = data.frame(fourier = fourier_list, mean_var = mean_var_df$mean_var)
  
  return(out)
}

group_size = 100
# group_size = 1

### First, Look at Raw Data:

## Nucleosome:
nuc_C26_iterative_fourier_AATT = find_iterative_fourier_ps2(C26, group_size, nuc)
nuc_C29_iterative_fourier_AATT = find_iterative_fourier_ps2(C29, group_size, nuc)
nuc_C31_iterative_fourier_AATT = find_iterative_fourier_ps2(C31, group_size, nuc)

plot(Mod(nuc_C26_iterative_fourier_AATT$fourier))
plot(Mod(nuc_C29_iterative_fourier_AATT$fourier))
plot(Mod(nuc_C31_iterative_fourier_AATT$fourier))

png(filename="figures/c-score-analysis/nuc_C26_fourier_amplitude_group_size_100_AATT.png")
plot(x=nuc_C26_iterative_fourier_AATT$mean_var, 
     y=Mod(nuc_C26_iterative_fourier_AATT$fourier))
title("nuc C26 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/nuc_C29_fourier_amplitude_group_size_100_AATT.png")
plot(x=nuc_C29_iterative_fourier_AATT$mean_var, 
     y=Mod(nuc_C29_iterative_fourier_AATT$fourier))
title("nuc C29 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/nuc_C31_fourier_amplitude_group_size_100_AATT.png")
plot(x=nuc_C31_iterative_fourier_AATT$mean_var, 
     y=Mod(nuc_C31_iterative_fourier_AATT$fourier))
title("nuc C31 fourier amplitude")
dev.off()

png(filename="figures/c-score-analysis/nuc_C26_fourier_phase_group_size_100_AATT.png")
plot(x=nuc_C26_iterative_fourier_AATT$mean_var, 
     y=Arg(nuc_C26_iterative_fourier_AATT$fourier)^2)
title("nuc C26 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/nuc_C29_fourier_phase_group_size_100_AATT.png")
plot(x=nuc_C29_iterative_fourier_AATT$mean_var, 
     y=Arg(nuc_C29_iterative_fourier_AATT$fourier)^2)
title("nuc C29 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/nuc_C31_fourier_phase_group_size_100_AATT.png")
plot(x=nuc_C31_iterative_fourier_AATT$mean_var, 
     y=Arg(nuc_C31_iterative_fourier_AATT$fourier)^2)
title("nuc C31 fourier phase")
dev.off()

## Random:
random_C26_iterative_fourier_AATT = find_iterative_fourier_ps2(C26, group_size, random)
random_C29_iterative_fourier_AATT = find_iterative_fourier_ps2(C29, group_size, random)
random_C31_iterative_fourier_AATT = find_iterative_fourier_ps2(C31, group_size, random)

plot(Mod(random_C26_iterative_fourier_AATT$fourier))
plot(Mod(random_C29_iterative_fourier_AATT$fourier))
plot(Mod(random_C31_iterative_fourier_AATT$fourier))

png(filename="figures/c-score-analysis/random_C26_fourier_amplitude_group_size_100_AATT.png")
plot(x=random_C26_iterative_fourier_AATT$mean_var, 
     y=Mod(random_C26_iterative_fourier_AATT$fourier))
title("random C26 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/random_C29_fourier_amplitude_group_size_100_AATT.png")
plot(x=random_C29_iterative_fourier_AATT$mean_var, 
     y=Mod(random_C29_iterative_fourier_AATT$fourier))
title("random C29 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/random_C31_fourier_amplitude_group_size_100_AATT.png")
plot(x=random_C31_iterative_fourier_AATT$mean_var, 
     y=Mod(random_C31_iterative_fourie_AATTr$fourier))
title("random C31 fourier amplitude")
dev.off()

png(filename="figures/c-score-analysis/random_C26_fourier_phase_group_size_100_AATT.png")
plot(x=random_C26_iterative_fourier_AATT$mean_var, 
     y=Arg(random_C26_iterative_fourier_AATT$fourier))
title("random C26 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/random_C29_fourier_phase_group_size_100_AATT.png")
plot(x=random_C29_iterative_fourier_AATT$mean_var, 
     y=Arg(random_C29_iterative_fourier_AATT$fourier))
title("random C29 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/random_C31_fourier_phase_group_size_100_AATT.png")
plot(x=random_C31_iterative_fourier_AATT$mean_var, 
     y=Arg(random_C31_iterative_fourier_AATT$fourier))
title("random C31 fourier phase")
dev.off()


## Tiling:
tiling_C26_iterative_fourier_AATT = find_iterative_fourier_ps2(C26, group_size, tiling)
tiling_C29_iterative_fourier_AATT = find_iterative_fourier_ps2(C29, group_size, tiling)
tiling_C31_iterative_fourier_AATT = find_iterative_fourier_ps2(C31, group_size, tiling)

plot(Mod(tiling_C26_iterative_fourier_AATT$fourier))
plot(Mod(tiling_C29_iterative_fourier_AATT$fourier))
plot(Mod(tiling_C31_iterative_fourier_AATT$fourier))

png(filename="figures/c-score-analysis/tiling_C26_fourier_amplitude_group_size_100_AATT.png")
plot(x=tiling_C26_iterative_fourier_AATT$mean_var, 
     y=Mod(tiling_C26_iterative_fourier_AATT$fourier))
title("tiling C26 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/tiling_C29_fourier_amplitude_group_size_100.png")
plot(x=tiling_C29_iterative_fourier$mean_var, 
     y=Mod(tiling_C29_iterative_fourier$fourier))
title("tiling C29 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/tiling_C31_fourier_amplitude_group_size_100.png")
plot(x=tiling_C31_iterative_fourier$mean_var, 
     y=Mod(tiling_C31_iterative_fourier$fourier))
title("tiling C31 fourier amplitude")
dev.off()

png(filename="figures/c-score-analysis/tiling_C26_fourier_phase_group_size_100.png")
plot(x=tiling_C26_iterative_fourier$mean_var, 
     y=Arg(tiling_C26_iterative_fourier$fourier)^2)
title("tiling C26 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/tiling_C29_fourier_phase_group_size_100.png")
plot(x=tiling_C29_iterative_fourier$mean_var, 
     y=Arg(tiling_C29_iterative_fourier$fourier)^2)
title("tiling C29 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/tiling_C31_fourier_phase_group_size_100.png")
plot(x=tiling_C31_iterative_fourier$mean_var, 
     y=Arg(tiling_C31_iterative_fourier$fourier)^2)
title("tiling C31 fourier phase")
dev.off()


## ChrV:
chrv_C26_iterative_fourier = find_iterative_fourier_ps2(C26, group_size, chrv)
chrv_C29_iterative_fourier = find_iterative_fourier_ps2(C29, group_size, chrv)
chrv_C31_iterative_fourier = find_iterative_fourier_ps2(C31, group_size, chrv)

plot(Mod(chrv_C26_iterative_fourier$fourier))
plot(Mod(chrv_C29_iterative_fourier$fourier))
plot(Mod(chrv_C31_iterative_fourier$fourier))

png(filename="figures/c-score-analysis/chrv_C26_fourier_amplitude_group_size_100.png")
plot(x=chrv_C26_iterative_fourier$mean_var, 
     y=Mod(chrv_C26_iterative_fourier$fourier))
title("chrv C26 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/chrv_C29_fourier_amplitude_group_size_100.png")
plot(x=chrv_C29_iterative_fourier$mean_var, 
     y=Mod(chrv_C29_iterative_fourier$fourier))
title("chrv C29 fourier amplitude")
dev.off()
png(filename="figures/c-score-analysis/chrv_C31_fourier_amplitude_group_size_100.png")
plot(x=chrv_C31_iterative_fourier$mean_var, 
     y=Mod(chrv_C31_iterative_fourier$fourier))
title("chrv C31 fourier amplitude")
dev.off()

png(filename="figures/c-score-analysis/chrv_C26_fourier_phase_group_size_100.png")
plot(x=chrv_C26_iterative_fourier$mean_var, 
     y=Arg(chrv_C26_iterative_fourier$fourier)^2)
title("chrv C26 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/chrv_C29_fourier_phase_group_size_100.png")
plot(x=chrv_C29_iterative_fourier$mean_var, 
     y=Arg(chrv_C29_iterative_fourier$fourier)^2)
title("chrv C29 fourier phase")
dev.off()
png(filename="figures/c-score-analysis/chrv_C31_fourier_phase_group_size_100.png")
plot(x=chrv_C31_iterative_fourier$mean_var, 
     y=Arg(chrv_C31_iterative_fourier$fourier)^2)
title("chrv C31 fourier phase")
dev.off()





### Now, Look at Their C0:

new_group_size = 100
# new_group_size = 1

## Nucleosome:
nuc_C0_iterative_fourier = find_iterative_fourier_ps2(C0, new_group_size, nuc)

png(filename="figures/c-score-analysis/nuc_C0_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_iterative_fourier$fourier))
title("nuc C0 fourier amplitude")
dev.off()
cor(nuc_C0_iterative_fourier$mean_var,
    Mod(nuc_C0_iterative_fourier$fourier))
# group_size=100: 0.8691267
# group_size=1: 0.3485084

png(filename="figures/c-score-analysis/nuc_C0_fourier_phase_group_size_100.png")
plot(x=nuc_C0_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_iterative_fourier$fourier))
title("nuc C0 fourier phase")
dev.off()
cor(nuc_C0_iterative_fourier$mean_var,
    Arg(nuc_C0_iterative_fourier$fourier)^2)
# group_size=100: -0.08659521
# group_size=1: 0.1265932

## Random:
random_C0_iterative_fourier = find_iterative_fourier_ps2(C0, new_group_size, random)

png(filename="figures/c-score-analysis/random_C0_fourier_amplitude_group_size_100.png")
plot(x=random_C0_iterative_fourier$mean_var, 
     y=Mod(random_C0_iterative_fourier$fourier))
title("random C0 fourier amplitude")
dev.off()
cor(random_C0_iterative_fourier$mean_var,
    Mod(random_C0_iterative_fourier$fourier))
# group_size=100: 0.3804952
# group_size=1: 0.2826443

png(filename="figures/c-score-analysis/random_C0_fourier_phase_group_size_100.png")
plot(x=random_C0_iterative_fourier$mean_var, 
     y=Arg(random_C0_iterative_fourier$fourier))
title("random C0 fourier phase")
dev.off()
cor(random_C0_iterative_fourier$mean_var,
    Arg(random_C0_iterative_fourier$fourier)^2)
# group_size=100: 0.1793607
# group_size=1: -0.006438992

## Tiling:
tiling_C0_iterative_fourier = find_iterative_fourier_ps2(C0, new_group_size, tiling)

png(filename="figures/c-score-analysis/tiling_C0_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_iterative_fourier$fourier))
title("tiling C0 fourier amplitude")
dev.off()
cor(tiling_C0_iterative_fourier$mean_var,
    Mod(tiling_C0_iterative_fourier$fourier))
# group_size=100: 0.3970173
# group_size=1: 0.3405582

png(filename="figures/c-score-analysis/tiling_C0_fourier_phase_group_size_100.png")
plot(x=tiling_C0_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_iterative_fourier$fourier))
title("tiling C0 fourier phase")
dev.off()
cor(tiling_C0_iterative_fourier$mean_var,
    Arg(tiling_C0_iterative_fourier$fourier)^2)
# group_size=100: 0.2333639
# group_size=1: 0.007517642

## ChrV:
chrv_C0_iterative_fourier = find_iterative_fourier_ps2(C0, new_group_size, chrv)

png(filename="figures/c-score-analysis/chrv_C0_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_iterative_fourier$fourier))
title("chrv C0 fourier amplitude")
dev.off()
cor(chrv_C0_iterative_fourier$mean_var,
    Mod(chrv_C0_iterative_fourier$fourier))
# group_size=100: 0.301645
# group_size=1: 0.2814013

png(filename="figures/c-score-analysis/chrv_C0_fourier_phase_group_size_100.png")
plot(x=chrv_C0_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_iterative_fourier$fourier)^2)
title("chrv C0 fourier phase")
dev.off()
cor(chrv_C0_iterative_fourier$mean_var,
    Arg(chrv_C0_iterative_fourier$fourier)^2)
# group_size=100: 0.19392
# group_size=1: 0.01781341



### Now, Look at New C0s:

## Nucleosome:
nuc_C0_new_v1_iterative_fourier = find_iterative_fourier_ps2(C0_new_v1, new_group_size, nuc)
nuc_C0_new_v2_iterative_fourier = find_iterative_fourier_ps2(C0_new_v2, new_group_size, nuc)
nuc_C0_new_v3_iterative_fourier = find_iterative_fourier_ps2(C0_new_v3, new_group_size, nuc)
nuc_C0_new_v4_iterative_fourier = find_iterative_fourier_ps2(C0_new_v4, new_group_size, nuc)
nuc_C0_new_v5_iterative_fourier = find_iterative_fourier_ps2(C0_new_v5, new_group_size, nuc)
nuc_C0_new_v6_iterative_fourier = find_iterative_fourier_ps2(C0_new_v6, new_group_size, nuc)
nuc_C0_new_v7_iterative_fourier = find_iterative_fourier_ps2(C0_new_v7, new_group_size, nuc)
nuc_C0_new_v8_iterative_fourier = find_iterative_fourier_ps2(C0_new_v8, new_group_size, nuc)

png(filename="figures/c-score-analysis/nuc_C0_new_v1_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v1_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v1_iterative_fourier$fourier))
title("nuc C0_new_v1 fourier amplitude")
dev.off()
cor(nuc_C0_new_v1_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v1_iterative_fourier$fourier))
# group_size=100: 0.7866793 ()
# group_size=1: 0.3272889

png(filename="figures/c-score-analysis/nuc_C0_new_v1_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v1_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v1_iterative_fourier$fourier)^2)
title("nuc C0_new_v1 fourier phase")
dev.off()
cor(nuc_C0_new_v1_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v1_iterative_fourier$fourier)^2)
# group_size=100: -0.4400623 ()
# group_size=1: 0.0568588


png(filename="figures/c-score-analysis/nuc_C0_new_v2_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v2_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v2_iterative_fourier$fourier))
title("nuc C0_new_v2 fourier amplitude")
dev.off()
cor(nuc_C0_new_v2_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v2_iterative_fourier$fourier))
# group_size=100: 0.664045 ()
# group_size=1: 0.2986372

png(filename="figures/c-score-analysis/nuc_C0_new_v2_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v2_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v2_iterative_fourier$fourier)^2)
title("nuc C0_new_v2 fourier phase")
dev.off()
cor(nuc_C0_new_v2_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v2_iterative_fourier$fourier)^2)
# group_size=100: -0.329724 ()
# group_size=1: -0.0148448


png(filename="figures/c-score-analysis/nuc_C0_new_v3_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v3_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v3_iterative_fourier$fourier))
title("nuc C0_new_v3 fourier amplitude")
dev.off()
cor(nuc_C0_new_v3_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v3_iterative_fourier$fourier))
# group_size=100: 0.8201253 ()
# group_size=1: 0.3393222

png(filename="figures/c-score-analysis/nuc_C0_new_v3_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v3_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v3_iterative_fourier$fourier)^2)
title("nuc C0_new_v3 fourier phase")
dev.off()
cor(nuc_C0_new_v3_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v3_iterative_fourier$fourier)^2)
# group_size=100: -0.4812385 ()
# group_size=1: 0.07295388


png(filename="figures/c-score-analysis/nuc_C0_new_v4_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v4_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v4_iterative_fourier$fourier))
title("nuc C0_new_v4 fourier amplitude")
dev.off()
cor(nuc_C0_new_v4_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v4_iterative_fourier$fourier))
# group_size=100: 0.8952093 ()
# group_size=1: 0.3567956

png(filename="figures/c-score-analysis/nuc_C0_new_v4_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v4_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v4_iterative_fourier$fourier)^2)
title("nuc C0_new_v4 fourier phase")
dev.off()
cor(nuc_C0_new_v4_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v4_iterative_fourier$fourier)^2)
# group_size=100: -0.2974169 ()
# group_size=1: 0.1339683


png(filename="figures/c-score-analysis/nuc_C0_new_v5_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v5_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v5_iterative_fourier$fourier))
title("nuc C0_new_v5 fourier amplitude")
dev.off()
cor(nuc_C0_new_v5_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v5_iterative_fourier$fourier))
# group_size=100: 0.8797394 ()
# group_size=1: 0.3539867

png(filename="figures/c-score-analysis/nuc_C0_new_v5_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v5_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v5_iterative_fourier$fourier)^2)
title("nuc C0_new_v5 fourier phase")
dev.off()
cor(nuc_C0_new_v5_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v5_iterative_fourier$fourier)^2)
# group_size=100: -0.190022 ()
# group_size=1: 0.1306574


png(filename="figures/c-score-analysis/nuc_C0_new_v6_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v6_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v6_iterative_fourier$fourier))
title("nuc C0_new_v6 fourier amplitude")
dev.off()
cor(nuc_C0_new_v6_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v6_iterative_fourier$fourier))
# group_size=100: 0.7875781 ()
# group_size=1: 0.3229214

png(filename="figures/c-score-analysis/nuc_C0_new_v6_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v6_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v6_iterative_fourier$fourier)^2)
title("nuc C0_new_v6 fourier phase")
dev.off()
cor(nuc_C0_new_v6_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v6_iterative_fourier$fourier)^2)
# group_size=100: 0.09506803 ()
# group_size=1: 0.08721717

png(filename="figures/c-score-analysis/nuc_C0_new_v7_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v7_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v7_iterative_fourier$fourier))
title("nuc C0_new_v7 fourier amplitude")
dev.off()
cor(nuc_C0_new_v7_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v7_iterative_fourier$fourier))
# group_size=100: 0.8691585 ()
# group_size=1: 0.350529

png(filename="figures/c-score-analysis/nuc_C0_new_v7_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v7_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v7_iterative_fourier$fourier)^2)
title("nuc C0_new_v7 fourier phase")
dev.off()
cor(nuc_C0_new_v7_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v7_iterative_fourier$fourier)^2)
# group_size=100: -0.1770919 ()
# group_size=1: 0.1175745


png(filename="figures/c-score-analysis/nuc_C0_new_v8_fourier_amplitude_group_size_100.png")
plot(x=nuc_C0_new_v8_iterative_fourier$mean_var, 
     y=Mod(nuc_C0_new_v8_iterative_fourier$fourier))
title("nuc C0_new_v8 fourier amplitude")
dev.off()
cor(nuc_C0_new_v8_iterative_fourier$mean_var,
    Mod(nuc_C0_new_v8_iterative_fourier$fourier))
# group_size=100: 0.7898344 ()
# group_size=1: 0.3162723

png(filename="figures/c-score-analysis/nuc_C0_new_v8_fourier_phase_group_size_100.png")
plot(x=nuc_C0_new_v8_iterative_fourier$mean_var, 
     y=Arg(nuc_C0_new_v8_iterative_fourier$fourier)^2)
title("nuc C0_new_v8 fourier phase")
dev.off()
cor(nuc_C0_new_v8_iterative_fourier$mean_var,
    Arg(nuc_C0_new_v8_iterative_fourier$fourier)^2)
# group_size=100: 0.1149842 ()
# group_size=1: 0.07933138

## Random:
random_C0_new_v1_iterative_fourier = find_iterative_fourier_ps2(C0_new_v1, new_group_size, random)
random_C0_new_v2_iterative_fourier = find_iterative_fourier_ps2(C0_new_v2, new_group_size, random)
random_C0_new_v3_iterative_fourier = find_iterative_fourier_ps2(C0_new_v3, new_group_size, random)
random_C0_new_v4_iterative_fourier = find_iterative_fourier_ps2(C0_new_v4, new_group_size, random)
random_C0_new_v5_iterative_fourier = find_iterative_fourier_ps2(C0_new_v5, new_group_size, random)
random_C0_new_v6_iterative_fourier = find_iterative_fourier_ps2(C0_new_v6, new_group_size, random)
random_C0_new_v7_iterative_fourier = find_iterative_fourier_ps2(C0_new_v7, new_group_size, random)
random_C0_new_v8_iterative_fourier = find_iterative_fourier_ps2(C0_new_v8, new_group_size, random)

png(filename="figures/c-score-analysis/random_C0_new_v1_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v1_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v1_iterative_fourier$fourier))
title("random C0_new_v1 fourier amplitude")
dev.off()
cor(random_C0_new_v1_iterative_fourier$mean_var,
    Mod(random_C0_new_v1_iterative_fourier$fourier))
# group_size=100: 0.2545164 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v1_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v1_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v1_iterative_fourier$fourier))
title("random C0_new_v1 fourier phase")
dev.off()
cor(random_C0_new_v1_iterative_fourier$mean_var,
    Arg(random_C0_new_v1_iterative_fourier$fourier)^2)
# group_size=100: -0.5670716 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v2_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v2_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v2_iterative_fourier$fourier))
title("random C0_new_v2 fourier amplitude")
dev.off()
cor(random_C0_new_v2_iterative_fourier$mean_var,
    Mod(random_C0_new_v2_iterative_fourier$fourier))
# group_size=100: 0.3714517 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v2_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v2_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v2_iterative_fourier$fourier)^2)
title("random C0_new_v2 fourier phase")
dev.off()
cor(random_C0_new_v2_iterative_fourier$mean_var,
    Arg(random_C0_new_v2_iterative_fourier$fourier)^2)
# group_size=100: -0.7197404 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v3_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v3_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v3_iterative_fourier$fourier))
title("random C0_new_v3 fourier amplitude")
dev.off()
cor(random_C0_new_v3_iterative_fourier$mean_var,
    Mod(random_C0_new_v3_iterative_fourier$fourier))
# group_size=100: 0.3997128 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v3_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v3_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v3_iterative_fourier$fourier)^2)
title("random C0_new_v3 fourier phase")
dev.off()
cor(random_C0_new_v3_iterative_fourier$mean_var,
    Arg(random_C0_new_v3_iterative_fourier$fourier)^2)
# group_size=100: -0.4482491 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v4_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v4_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v4_iterative_fourier$fourier))
title("random C0_new_v4 fourier amplitude")
dev.off()
cor(random_C0_new_v4_iterative_fourier$mean_var,
    Mod(random_C0_new_v4_iterative_fourier$fourier))
# group_size=100: 0.3972147 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v4_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v4_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v4_iterative_fourier$fourier)^2)
title("random C0_new_v4 fourier phase")
dev.off()
cor(random_C0_new_v4_iterative_fourier$mean_var,
    Arg(random_C0_new_v4_iterative_fourier$fourier)^2)
# group_size=100: 0.01510235 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v5_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v5_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v5_iterative_fourier$fourier))
title("random C0_new_v5 fourier amplitude")
dev.off()
cor(random_C0_new_v5_iterative_fourier$mean_var,
    Mod(random_C0_new_v5_iterative_fourier$fourier))
# group_size=100: 0.311049 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v5_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v5_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v5_iterative_fourier$fourier)^2)
title("random C0_new_v5 fourier phase")
dev.off()
cor(random_C0_new_v5_iterative_fourier$mean_var,
    Arg(random_C0_new_v5_iterative_fourier$fourier)^2)
# group_size=100: 0.0679625 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v6_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v6_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v6_iterative_fourier$fourier))
title("random C0_new_v6 fourier amplitude")
dev.off()
cor(random_C0_new_v6_iterative_fourier$mean_var,
    Mod(random_C0_new_v6_iterative_fourier$fourier))
# group_size=100: 0.4575471 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v6_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v6_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v6_iterative_fourier$fourier)^2)
title("random C0_new_v6 fourier phase")
dev.off()
cor(random_C0_new_v6_iterative_fourier$mean_var,
    Arg(random_C0_new_v6_iterative_fourier$fourier)^2)
# group_size=100: -0.1713985 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v7_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v7_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v7_iterative_fourier$fourier))
title("random C0_new_v7 fourier amplitude")
dev.off()
cor(random_C0_new_v7_iterative_fourier$mean_var,
    Mod(random_C0_new_v7_iterative_fourier$fourier))
# group_size=100: 0.3574101 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v7_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v7_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v7_iterative_fourier$fourier)^2)
title("random C0_new_v7 fourier phase")
dev.off()
cor(random_C0_new_v7_iterative_fourier$mean_var,
    Arg(random_C0_new_v7_iterative_fourier$fourier)^2)
# group_size=100: -0.01056585 ()
# group_size=1: 


png(filename="figures/c-score-analysis/random_C0_new_v8_fourier_amplitude_group_size_100.png")
plot(x=random_C0_new_v8_iterative_fourier$mean_var, 
     y=Mod(random_C0_new_v8_iterative_fourier$fourier))
title("random C0_new_v8 fourier amplitude")
dev.off()
cor(random_C0_new_v8_iterative_fourier$mean_var,
    Mod(random_C0_new_v8_iterative_fourier$fourier))
# group_size=100: 0.5247695 ()
# group_size=1: 

png(filename="figures/c-score-analysis/random_C0_new_v8_fourier_phase_group_size_100.png")
plot(x=random_C0_new_v8_iterative_fourier$mean_var, 
     y=Arg(random_C0_new_v8_iterative_fourier$fourier)^2)
title("random C0_new_v8 fourier phase")
dev.off()
cor(random_C0_new_v8_iterative_fourier$mean_var,
    Arg(random_C0_new_v8_iterative_fourier$fourier)^2)
# group_size=100: -0.2674374 ()
# group_size=1: 



## Tiling:
tiling_C0_new_v1_iterative_fourier = find_iterative_fourier_ps2(C0_new_v1, new_group_size, tiling)
tiling_C0_new_v2_iterative_fourier = find_iterative_fourier_ps2(C0_new_v2, new_group_size, tiling)
tiling_C0_new_v3_iterative_fourier = find_iterative_fourier_ps2(C0_new_v3, new_group_size, tiling)
tiling_C0_new_v4_iterative_fourier = find_iterative_fourier_ps2(C0_new_v4, new_group_size, tiling)
tiling_C0_new_v5_iterative_fourier = find_iterative_fourier_ps2(C0_new_v5, new_group_size, tiling)
tiling_C0_new_v6_iterative_fourier = find_iterative_fourier_ps2(C0_new_v6, new_group_size, tiling)
tiling_C0_new_v7_iterative_fourier = find_iterative_fourier_ps2(C0_new_v7, new_group_size, tiling)
tiling_C0_new_v8_iterative_fourier = find_iterative_fourier_ps2(C0_new_v8, new_group_size, tiling)

png(filename="figures/c-score-analysis/tiling_C0_new_v1_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v1_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v1_iterative_fourier$fourier))
title("tiling C0_new_v1 fourier amplitude")
dev.off()
cor(tiling_C0_new_v1_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v1_iterative_fourier$fourier))
# group_size=100: 0.3047383 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v1_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v1_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v1_iterative_fourier$fourier)^2)
title("tiling C0_new_v1 fourier phase")
dev.off()
cor(tiling_C0_new_v1_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v1_iterative_fourier$fourier)^2)
# group_size=100: -0.4098702 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v2_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v2_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v2_iterative_fourier$fourier))
title("tiling C0_new_v2 fourier amplitude")
dev.off()
cor(tiling_C0_new_v2_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v2_iterative_fourier$fourier))
# group_size=100: 0.3832398 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v2_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v2_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v2_iterative_fourier$fourier)^2)
title("tiling C0_new_v2 fourier phase")
dev.off()
cor(tiling_C0_new_v2_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v2_iterative_fourier$fourier)^2)
# group_size=100: -0.6566997 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v3_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v3_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v3_iterative_fourier$fourier))
title("tiling C0_new_v3 fourier amplitude")
dev.off()
cor(tiling_C0_new_v3_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v3_iterative_fourier$fourier))
# group_size=100: 0.4546317 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v3_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v3_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v3_iterative_fourier$fourier)^2)
title("tiling C0_new_v3 fourier phase")
dev.off()
cor(tiling_C0_new_v3_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v3_iterative_fourier$fourier)^2)
# group_size=100: -0.301346 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v4_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v4_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v4_iterative_fourier$fourier))
title("tiling C0_new_v4 fourier amplitude")
dev.off()
cor(tiling_C0_new_v4_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v4_iterative_fourier$fourier))
# group_size=100: 0.4527697 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v4_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v4_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v4_iterative_fourier$fourier)^2)
title("tiling C0_new_v4 fourier phase")
dev.off()
cor(tiling_C0_new_v4_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v4_iterative_fourier$fourier)^2)
# group_size=100: 0.1984834 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v5_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v5_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v5_iterative_fourier$fourier))
title("tiling C0_new_v5 fourier amplitude")
dev.off()
cor(tiling_C0_new_v5_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v5_iterative_fourier$fourier))
# group_size=100: 0.42934 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v5_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v5_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v5_iterative_fourier$fourier)^2)
title("tiling C0_new_v5 fourier phase")
dev.off()
cor(tiling_C0_new_v5_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v5_iterative_fourier$fourier)^2)
# group_size=100: 0.1489829 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v6_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v6_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v6_iterative_fourier$fourier))
title("tiling C0_new_v6 fourier amplitude")
dev.off()
cor(tiling_C0_new_v6_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v6_iterative_fourier$fourier))
# group_size=100: 0.3072414 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v6_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v6_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v6_iterative_fourier$fourier)^2)
title("tiling C0_new_v6 fourier phase")
dev.off()
cor(tiling_C0_new_v6_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v6_iterative_fourier$fourier)^2)
# group_size=100: -0.08049581 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v7_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v7_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v7_iterative_fourier$fourier))
title("tiling C0_new_v7 fourier amplitude")
dev.off()
cor(tiling_C0_new_v7_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v7_iterative_fourier$fourier))
# group_size=100: 0.4162789 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v7_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v7_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v7_iterative_fourier$fourier)^2)
title("tiling C0_new_v7 fourier phase")
dev.off()
cor(tiling_C0_new_v7_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v7_iterative_fourier$fourier)^2)
# group_size=100: 0.07647496 ()
# group_size=1: 


png(filename="figures/c-score-analysis/tiling_C0_new_v8_fourier_amplitude_group_size_100.png")
plot(x=tiling_C0_new_v8_iterative_fourier$mean_var, 
     y=Mod(tiling_C0_new_v8_iterative_fourier$fourier))
title("tiling C0_new_v8 fourier amplitude")
dev.off()
cor(tiling_C0_new_v8_iterative_fourier$mean_var,
    Mod(tiling_C0_new_v8_iterative_fourier$fourier))
# group_size=100: 0.3513792 ()
# group_size=1: 

png(filename="figures/c-score-analysis/tiling_C0_new_v8_fourier_phase_group_size_100.png")
plot(x=tiling_C0_new_v8_iterative_fourier$mean_var, 
     y=Arg(tiling_C0_new_v8_iterative_fourier$fourier)^2)
title("tiling C0_new_v8 fourier phase")
dev.off()
cor(tiling_C0_new_v8_iterative_fourier$mean_var,
    Arg(tiling_C0_new_v8_iterative_fourier$fourier)^2)
# group_size=100: -0.1505237 ()
# group_size=1: 



## ChrV:
chrv_C0_new_v1_iterative_fourier = find_iterative_fourier_ps2(C0_new_v1, new_group_size, chrv)
chrv_C0_new_v2_iterative_fourier = find_iterative_fourier_ps2(C0_new_v2, new_group_size, chrv)
chrv_C0_new_v3_iterative_fourier = find_iterative_fourier_ps2(C0_new_v3, new_group_size, chrv)
chrv_C0_new_v4_iterative_fourier = find_iterative_fourier_ps2(C0_new_v4, new_group_size, chrv)
chrv_C0_new_v5_iterative_fourier = find_iterative_fourier_ps2(C0_new_v5, new_group_size, chrv)
chrv_C0_new_v6_iterative_fourier = find_iterative_fourier_ps2(C0_new_v6, new_group_size, chrv)
chrv_C0_new_v7_iterative_fourier = find_iterative_fourier_ps2(C0_new_v7, new_group_size, chrv)
chrv_C0_new_v8_iterative_fourier = find_iterative_fourier_ps2(C0_new_v8, new_group_size, chrv)

png(filename="figures/c-score-analysis/chrv_C0_new_v1_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v1_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v1_iterative_fourier$fourier))
title("chrv C0_new_v1 fourier amplitude")
dev.off()
cor(chrv_C0_new_v1_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v1_iterative_fourier$fourier))
# group_size=100: 0.2841569 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v1_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v1_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v1_iterative_fourier$fourier)^2)
title("chrv C0_new_v1 fourier phase")
dev.off()
cor(chrv_C0_new_v1_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v1_iterative_fourier$fourier)^2)
# group_size=100: -0.3080718 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v2_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v2_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v2_iterative_fourier$fourier))
title("chrv C0_new_v2 fourier amplitude")
dev.off()
cor(chrv_C0_new_v2_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v2_iterative_fourier$fourier))
# group_size=100: 0.3168109 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v2_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v2_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v2_iterative_fourier$fourier)^2)
title("chrv C0_new_v2 fourier phase")
dev.off()
cor(chrv_C0_new_v2_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v2_iterative_fourier$fourier)^2)
# group_size=100: -0.6004006 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v3_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v3_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v3_iterative_fourier$fourier))
title("chrv C0_new_v3 fourier amplitude")
dev.off()
cor(chrv_C0_new_v3_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v3_iterative_fourier$fourier))
# group_size=100: 0.4284334 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v3_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v3_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v3_iterative_fourier$fourier)^2)
title("chrv C0_new_v3 fourier phase")
dev.off()
cor(chrv_C0_new_v3_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v3_iterative_fourier$fourier)^2)
# group_size=100: -0.2437397 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v4_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v4_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v4_iterative_fourier$fourier))
title("chrv C0_new_v4 fourier amplitude")
dev.off()
cor(chrv_C0_new_v4_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v4_iterative_fourier$fourier))
# group_size=100: 0.4341921 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v4_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v4_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v4_iterative_fourier$fourier)^2)
title("chrv C0_new_v4 fourier phase")
dev.off()
cor(chrv_C0_new_v4_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v4_iterative_fourier$fourier)^2)
# group_size=100: 0.1846693 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v5_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v5_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v5_iterative_fourier$fourier))
title("chrv C0_new_v5 fourier amplitude")
dev.off()
cor(chrv_C0_new_v5_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v5_iterative_fourier$fourier))
# group_size=100: 0.3401086 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v5_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v5_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v5_iterative_fourier$fourier)^2)
title("chrv C0_new_v5 fourier phase")
dev.off()
cor(chrv_C0_new_v5_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v5_iterative_fourier$fourier)^2)
# group_size=100: 0.1944237 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v6_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v6_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v6_iterative_fourier$fourier))
title("chrv C0_new_v6 fourier amplitude")
dev.off()
cor(chrv_C0_new_v6_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v6_iterative_fourier$fourier))
# group_size=100: 0.1824184 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v6_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v6_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v6_iterative_fourier$fourier)^2)
title("chrv C0_new_v6 fourier phase")
dev.off()
cor(chrv_C0_new_v6_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v6_iterative_fourier$fourier)^2)
# group_size=100: 0.02251557 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v7_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v7_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v7_iterative_fourier$fourier))
title("chrv C0_new_v7 fourier amplitude")
dev.off()
cor(chrv_C0_new_v7_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v7_iterative_fourier$fourier))
# group_size=100: 0.3488099 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v7_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v7_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v7_iterative_fourier$fourier)^2)
title("chrv C0_new_v7 fourier phase")
dev.off()
cor(chrv_C0_new_v7_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v7_iterative_fourier$fourier)^2)
# group_size=100: 0.119158 ()
# group_size=1: 


png(filename="figures/c-score-analysis/chrv_C0_new_v8_fourier_amplitude_group_size_100.png")
plot(x=chrv_C0_new_v8_iterative_fourier$mean_var, 
     y=Mod(chrv_C0_new_v8_iterative_fourier$fourier))
title("chrv C0_new_v8 fourier amplitude")
dev.off()
cor(chrv_C0_new_v8_iterative_fourier$mean_var,
    Mod(chrv_C0_new_v8_iterative_fourier$fourier))
# group_size=100: 0.2526444 ()
# group_size=1: 

png(filename="figures/c-score-analysis/chrv_C0_new_v8_fourier_phase_group_size_100.png")
plot(x=chrv_C0_new_v8_iterative_fourier$mean_var, 
     y=Arg(chrv_C0_new_v8_iterative_fourier$fourier)^2)
title("chrv C0_new_v8 fourier phase")
dev.off()
cor(chrv_C0_new_v8_iterative_fourier$mean_var,
    Arg(chrv_C0_new_v8_iterative_fourier$fourier)^2)
# group_size=100: -0.01780233 ()
# group_size=1: 

