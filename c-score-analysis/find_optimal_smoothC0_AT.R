library(tidyverse)

##### Based on A/T:

ps1 <- paste0("X", 1:50, "mono")

find_topm_fourier_power_AT_original = function(f, m, n, data) {
  # temp_dat_train = data %>%
  #   select("x50mer","C26","C29","C31","C0","Amplitude","Phase", all_of(ps1))
  temp_dat_train = data %>%
    select("x50mer","C26","C29","C31","C0", all_of(ps1))
  
  temp_dat_train = temp_dat_train %>%
    mutate(Q_f = C26 + f*C31)
  # temp_dat_train[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(temp_dat_train[,c("C26","C29","C31")], n=n, aa=aa)
  
  ### Want to find Fourier power for A/T 10bp periodicity for a sequence of values of f
  # for the m sequences with the highest Q(f) values
  temp_dat_train_topm = temp_dat_train %>%
    top_n(m, Q_f)
  
  Xone_topm = temp_dat_train_topm %>% select(all_of(ps1))
  Xone_AorT_topm = matrix(nrow=nrow(Xone_topm), ncol=ncol(Xone_topm))
  colnames(Xone_AorT_topm) = colnames(Xone_topm)
  Xone_AorT_topm[] = ((Xone_topm == "A") | (Xone_topm == "T")) %>% 
    as.matrix() %>% as.numeric()
  
  # Collapse all m sequences (This is what Basu does on Page 12 of Supplementary 
  # Information - "Refinements of the model"):
  Xone_AorT_topm = apply(Xone_AorT_topm, 2, mean)
  
  temp_fourier = fft(Xone_AorT_topm)[6]
  return(Re(temp_fourier)^2 + Im(temp_fourier)^2)
}

find_topm_fourier_power_AT_original(f=0.7, m=1000, n=c(26,29,31), data=dat_chrv_2)
find_topm_fourier_power_AT_original(f=0.7, m=1000, n=c(26,29,31), data=dat_random)


find_optimal_f_AT = function(data, library_name, m=1000, n=c(26,29,31),
                             f_values=seq(from=0.1,to=2,by=0.1)) {
  
  # aa_grid = map2(rep(a2_values, length(a3_values)), rep(a3_values, each=length(a2_values)),~c(1, .x, .y))
  power_vector = f_values %>%
    map(~find_topm_fourier_power_AT_original(f=.x, m=m, n=n, data=data)) %>%
    unlist()
  # power_vector = aa_grid %>%
  #   map(~find_topm_fourier_power_AT(aa=.x, m=m, n=n, data=data)) %>%
  #   unlist() 
  
  # power_grid = power_vector %>%
  #   matrix(nrow=length(a2_values), ncol=length(a3_values))
  
  # Rows (a2) correspond to (n=29)/(n=26)
  # Columns (a3) correspond to (n=31)/(n=26)
  
  names(power_vector) = f_values
  
  # rownames(power_grid) = a2_values
  # colnames(power_grid) = a3_values
  
  f_plot = ggplot() +
    geom_line(aes(x=f_values, y=power_vector))
  
  # Heatmap:
  # gg_format = expand.grid(a2=a2_values, a3=a3_values)
  # gg_format$weight = power_vector
  # heat = ggplot(gg_format, aes(a2, a3)) + 
  #   geom_tile(aes(fill=weight)) +
  #   scale_fill_gradient2("Fourier Power", low="blue", high="red", midpoint=mean(power_vector)) +
  #   ggtitle(paste0("Optimized Using top ", m, " sequences from ", library_name, " Library"))
  
  # min_row = gg_format[min(which(gg_format$weight == min(gg_format$weight))),]
  
  f_optimal = as.numeric(names(which(power_vector == min(power_vector))))
  
  # aa_optimal = c(1, min_row$a2, min_row$a3)
  
  to_return = list()
  # to_return$power_grid = power_grid
  to_return$power_vector = power_vector
  to_return$f_plot = f_plot
  to_return$f_optimal = f_optimal
  
  return(to_return)
}

optimal_f_chrv_2_tophalf_prelim_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/2),
                                                       f_values=seq(from=0.1,to=1,by=0.01))
optimal_f_chrv_2_tophalf_prelim_AT$f_optimal
# 0.8
optimal_f_chrv_2_tophalf_prelim_AT$f_plot

optimal_f_chrv_2_tophalf_refined_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/2),
                                                        f_values=seq(from=0.6,to=1,by=0.01))
optimal_f_chrv_2_tophalf_refined_AT$f_optimal
# 0.78
optimal_f_chrv_2_tophalf_refined_AT$f_plot



optimal_f_chrv_2_topquartile_prelim_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/4),
                                                           f_values=seq(from=0.1,to=2,by=0.1))
optimal_f_chrv_2_topquartile_prelim_AT$f_optimal
# 0.9
optimal_f_chrv_2_topquartile_prelim_AT$f_plot

optimal_f_chrv_2_topquartile_refined_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/4),
                                                            f_values=seq(from=0.6,to=1,by=0.01))
optimal_f_chrv_2_topquartile_refined_AT$f_optimal
# 0.85
optimal_f_chrv_2_topquartile_refined_AT$f_plot



optimal_f_chrv_2_top1000_prelim_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", 1000,
                                                       f_values=seq(from=0.1,to=2,by=0.1))
optimal_f_chrv_2_top1000_prelim_AT$f_optimal
# 0.8
optimal_f_chrv_2_top1000_prelim_AT$f_plot

optimal_f_chrv_2_top1000_refined_AT = find_optimal_f_AT(dat_chrv_2, "ChrV", 1000,
                                                        f_values=seq(from=0.6,to=1,by=0.01))
optimal_f_chrv_2_top1000_refined_AT$f_optimal
# 0.81
optimal_f_chrv_2_top1000_refined_AT$f_plot



optimal_f_dat_random_top1000_prelim_AT = find_optimal_f_AT(dat_random, "Random", 1000,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_random_top1000_prelim_AT$f_optimal
# 0.69
optimal_f_dat_random_top1000_prelim_AT$f_plot


optimal_f_dat_random_top2000_prelim_AT = find_optimal_f_AT(dat_random, "Random", 2000,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_random_top2000_prelim_AT$f_optimal
# 0.71
optimal_f_dat_random_top2000_prelim_AT$f_plot


optimal_f_dat_random_topquartile_prelim_AT = find_optimal_f_AT(dat_random, "Random", nrow(dat_random)/4,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_random_topquartile_prelim_AT$f_optimal
# 0.73
optimal_f_dat_random_topquartile_prelim_AT$f_plot


optimal_f_dat_tiling_top500_prelim_AT = find_optimal_f_AT(dat_tiling, "Tiling", 500,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_tiling_top500_prelim_AT$f_optimal
# 0.70, 0.71
optimal_f_dat_tiling_top500_prelim_AT$f_plot



optimal_f_dat_tiling_top1000_prelim_AT = find_optimal_f_AT(dat_tiling, "Tiling", 1000,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_tiling_top1000_prelim_AT$f_optimal
# 0.73
optimal_f_dat_tiling_top1000_prelim_AT$f_plot


optimal_f_dat_tiling_top2000_prelim_AT = find_optimal_f_AT(dat_tiling, "Tiling", 2000,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_tiling_top2000_prelim_AT$f_optimal
# 0.75
optimal_f_dat_tiling_top2000_prelim_AT$f_plot


optimal_f_dat_tiling_topquartile_prelim_AT = find_optimal_f_AT(dat_tiling, "Tiling", nrow(dat_tiling)/4,
                                                           f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_tiling_topquartile_prelim_AT$f_optimal
# 0.77
optimal_f_dat_tiling_topquartile_prelim_AT$f_plot


optimal_f_dat_tiling_tophalf_prelim_AT = find_optimal_f_AT(dat_tiling, "Tiling", nrow(dat_tiling)/2,
                                                               f_values=seq(from=0.3,to=1,by=0.01))
optimal_f_dat_tiling_tophalf_prelim_AT$f_optimal
# 0.71
optimal_f_dat_tiling_tophalf_prelim_AT$f_plot


find_c0Aphi <- function(dat, n=c(26,29,31), aa=c(1,1,1)){
  # Find C0, Amplitude, and Phase for each sequence
  # dat: Each row corresponds to a sequence.
  #   First three columns need to be C26, C29, and C31, in that order
  # n: Phase shift (bp) in C26, C29, and C31, respectively, by default c(26,29,31)
  # aa: Ratio A26/A26, A29/A26, A31/A26, respectively, by default c(1,1,1) 
  #   Note: In Basu, they claim to use c(1,0.82,0.7), but they actually use c(1, 1/0.82, 1/0.7)
  # 
  # Returns: columns C0, A26, Phi
  mat <- matrix(0, nrow = 3, ncol = 3)
  k <- 2*pi/10.4
  mat[1:3,1] <- 1
  mat[1:3,2] <- sin(n*k)
  mat[1:3,3] <- cos(n*k)
  mat[1,2:3] <- mat[1,2:3]*aa[1]
  mat[2,2:3] <- mat[2,2:3]*aa[2]
  mat[3,2:3] <- mat[3,2:3]*aa[3]
  inv_mat <- solve(mat)
  c0A1A2 <- as.matrix(dat[,1:3]) %*% t(inv_mat)
  c0Aphi <- c0A1A2
  c0Aphi[,1] <- c0A1A2[,1]
  c0Aphi[,2] <- sqrt(c0A1A2[,2]^2 + c0A1A2[,3]^2)
  c0Aphi[,3] <- sign(c0A1A2[,3]) * acos(c0A1A2[,2]/c0Aphi[,2])
  return(c0Aphi)
}

find_topm_fourier_power_AT = function(aa, m, n, data) {
  # temp_dat_train = data %>%
  #   select("x50mer","C26","C29","C31","C0","Amplitude","Phase", all_of(ps1))
  temp_dat_train = data %>%
    select("x50mer","C26","C29","C31","C0", all_of(ps1))
  
  temp_dat_train[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(temp_dat_train[,c("C26","C29","C31")], n=n, aa=aa)
  
  ### Want to find Fourier power for A/T 10bp periodicity for a grid of values of f
  # for the m sequences with the highest C0_new values
  temp_dat_train_topm = temp_dat_train %>%
    top_n(m, C0_new)
  
  Xone_topm = temp_dat_train_topm %>% select(all_of(ps1))
  Xone_AorT_topm = matrix(nrow=nrow(Xone_topm), ncol=ncol(Xone_topm))
  colnames(Xone_AorT_topm) = colnames(Xone_topm)
  Xone_AorT_topm[] = ((Xone_topm == "A") | (Xone_topm == "T")) %>% 
    as.matrix() %>% as.numeric()
  
  # Collapse all m sequences (This is what Basu does on Page 12 of Supplementary 
  # Information - "Refinements of the model"):
  Xone_AorT_topm = apply(Xone_AorT_topm, 2, mean)
  
  temp_fourier = fft(Xone_AorT_topm)[6]
  return(Re(temp_fourier)^2 + Im(temp_fourier)^2)
}

find_topm_fourier_power_AT(aa=c(1, 1/0.82, 1/0.7), m=1000, n=c(26,29,31), data=dat_chrv_2)


find_optimal_aa_AT = function(data, library_name, m=1000, n=c(26,29,31),
                              a2_values=seq(from=0.1,to=3,by=0.1), 
                              a3_values=seq(from=0.1,to=2.5,by=0.1)) {
  
  aa_grid = map2(rep(a2_values, length(a3_values)), rep(a3_values, each=length(a2_values)),~c(1, .x, .y))
  
  power_vector = aa_grid %>%
    map(~find_topm_fourier_power_AT(aa=.x, m=m, n=n, data=data)) %>%
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



# Optimize Based on Top Quartile of ChrV Library:
optimal_chrv_2_topquartile_prelim_AT = find_optimal_aa_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/4),
                                                          a2_values=seq(from=0.1,to=2,by=0.1),
                                                          a3_values=seq(from=0.5,to=2,by=0.1))
optimal_chrv_2_topquartile_prelim_AT$aa_optimal
# 1.0 0.4 1.1
optimal_chrv_2_topquartile_prelim_AT$heatmap
optimal_chrv_2_topquartile_refined_AT = find_optimal_aa_AT(dat_chrv_2, "ChrV", floor(nrow(dat_chrv_2)/4),
                                                           a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                           a3_values=seq(from=0.9,to=1.3,by=0.01))
optimal_chrv_2_topquartile_refined_AT$aa_optimal
# 1.00 0.37 1.05
optimal_chrv_2_topquartile_refined_AT$heatmap


# Optimize Based on Top 1000 of ChrV Library:
optimal_chrv_2_top1000_prelim_AT = find_optimal_aa_AT(dat_chrv_2, "ChrV", 1000,
                                                      a2_values=seq(from=0.1,to=2,by=0.1),
                                                      a3_values=seq(from=0.5,to=2,by=0.1))
optimal_chrv_2_top1000_prelim_AT$aa_optimal
# 1.0 0.8 1.1
optimal_chrv_2_top1000_prelim_AT$heatmap
optimal_chrv_2_top1000_refined_AT = find_optimal_aa_AT(dat_chrv_2, "ChrV", 1000,
                                                       a2_values=seq(from=0.6,to=1.0,by=0.01),
                                                       a3_values=seq(from=0.9,to=1.3,by=0.01))
optimal_chrv_2_top1000_refined_AT$aa_optimal
# 1.00 0.72 1.06
optimal_chrv_2_top1000_refined_AT$heatmap
