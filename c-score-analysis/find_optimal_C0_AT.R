library(tidyverse)

##### Based on A/T:

ps1 <- paste0("X", 1:50, "mono")

find_topm_fourier_power_AT = function(aa, m, n, data) {
  temp_dat_train = data %>%
    select("x50mer","C26","C29","C31","C0","Amplitude","Phase", all_of(ps1))
  
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


# Example using values from the paper, 
find_topm_fourier_power_AT(aa=c(1, 1/0.82, 1/0.7), m=1000, n=c(26,29,31), data=random)


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


# Optimize Based on Top Half of Random Library:
optimal_random_tophalf_prelim_AT = find_optimal_aa_AT(random, "Random", nrow(random)/2)
optimal_random_tophalf_prelim_AT$aa_optimal
# 
optimal_random_tophalf_prelim_AT$heatmap
optimal_random_tophalf_refined_AT = find_optimal_aa_AT(random, "Random", nrow(random)/2,
                                                       a2_values=seq(from=2.6,to=3,by=0.01),
                                                       a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_random_tophalf_refined_AT$aa_optimal
# 
optimal_random_tophalf_refined_AT$heatmap


# Optimize Based on Top Quartile of Random Library:
optimal_random_topquartile_prelim_AT = find_optimal_aa_AT(random, "Random", nrow(random)/4)
optimal_random_topquartile_prelim_AT$aa_optimal
# 1.0 2.2 1.4
optimal_random_topquartile_prelim_AT$heatmap
optimal_random_topquartile_refined_AT = find_optimal_aa_AT(random, "Random", nrow(random)/4,
                                                           a2_values=seq(from=2.0,to=2.4,by=0.01),
                                                           a3_values=seq(from=1.2,to=1.6,by=0.01))
optimal_random_topquartile_refined_AT$aa_optimal
# 1.00 2.02 1.38
optimal_random_topquartile_refined_AT$heatmap


# Optimize Based on Top 1000 of Random Library:
optimal_random_top1000_prelim_AT = find_optimal_aa_AT(random, "Random")
optimal_random_top1000_prelim_AT$aa_optimal
# 1.0 1.1 1.3
optimal_random_top1000_prelim_AT$heatmap
min(optimal_random_top1000_prelim_AT$power_grid)
# 0.00287718
optimal_random_top1000_refined_AT = find_optimal_aa_AT(random, "Random",
                                                       a2_values=seq(from=0.9,to=1.3,by=0.01),
                                                       a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_top1000_refined_AT$aa_optimal
# 1.00 1.27 1.35
optimal_random_top1000_refined_AT$heatmap
min(optimal_random_top1000_refined_AT$power_grid)
# 0.0002080229



# Optimize Based on Top 4000 of Random Library:
optimal_random_top4000_prelim_AT = find_optimal_aa_AT(random, "Random", 4000)
optimal_random_top4000_prelim_AT$aa_optimal
# 1.0 1.7 1.4
optimal_random_top4000_prelim_AT$heatmap
min(optimal_random_top4000_prelim_AT$power_grid)
# 6.495431e-06
optimal_random_top4000_refined_AT = find_optimal_aa_AT(random, "Random", 4000,
                                                       a2_values=seq(from=1.5,to=1.9,by=0.01),
                                                       a3_values=seq(from=1.2,to=1.6,by=0.01))
optimal_random_top4000_refined_AT$aa_optimal
# 1.00 1.66 1.40
optimal_random_top4000_refined_AT$heatmap
min(optimal_random_top4000_refined_AT$power_grid)
# 3.43112e-06



# Optimize Based on Top 2000 of Random Library:
optimal_random_top2000_prelim_AT = find_optimal_aa_AT(random, "Random", 2000)
optimal_random_top2000_prelim_AT$aa_optimal
# 1.0 1.2 1.3
optimal_random_top2000_prelim_AT$heatmap
min(optimal_random_top2000_prelim_AT$power_grid)
# 0.0004660397
optimal_random_top2000_refined_AT = find_optimal_aa_AT(random, "Random", 2000,
                                                       a2_values=seq(from=1,to=1.4,by=0.01),
                                                       a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_top2000_refined_AT$aa_optimal
# 1.00 1.36 1.34
optimal_random_top2000_refined_AT$heatmap
min(optimal_random_top2000_refined_AT$power_grid)
# 4.539719e-07




# TILING:

# Optimize Based on Top Half of Tiling Library:
optimal_tiling_tophalf_prelim_AT = find_optimal_aa_AT(tiling, "Tiling", nrow(tiling)/2)
optimal_tiling_tophalf_prelim_AT$aa_optimal
# 
optimal_tiling_tophalf_prelim_AT$heatmap
optimal_tiling_tophalf_refined_AT = find_optimal_aa_AT(tiling, "Tiling", nrow(tiling)/2,
                                                       a2_values=seq(from=2.6,to=3,by=0.01),
                                                       a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_tiling_tophalf_refined_AT$aa_optimal
# 
optimal_tiling_tophalf_refined_AT$heatmap


# Optimize Based on Top Quartile of Tiling Library:
optimal_tiling_topquartile_prelim_AT = find_optimal_aa_AT(tiling, "Tiling", nrow(tiling)/4)
optimal_tiling_topquartile_prelim_AT$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_topquartile_prelim_AT$heatmap
optimal_tiling_topquartile_refined_AT = find_optimal_aa_AT(tiling, "Tiling", nrow(tiling)/4,
                                                           a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                           a3_values=seq(from=1.0,to=1.4,by=0.01))
optimal_tiling_topquartile_refined_AT$aa_optimal
# 1.00 0.39 1.16
optimal_tiling_topquartile_refined_AT$heatmap


# Optimize Based on Top 1000 of Tiling Library:
optimal_tiling_top1000_prelim_AT = find_optimal_aa_AT(tiling, "Tiling")
optimal_tiling_top1000_prelim_AT$aa_optimal
# 
optimal_tiling_top1000_prelim_AT$heatmap
min(optimal_tiling_top1000_prelim_AT$power_grid)
# 
optimal_tiling_top1000_refined_AT = find_optimal_aa_AT(tiling, "Tiling",
                                                       a2_values=seq(from=1.1,to=1.5,by=0.01),
                                                       a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_tiling_top1000_refined_AT$aa_optimal
# 
optimal_tiling_top1000_refined_AT$heatmap
min(optimal_tiling_top1000_refined_AT$power_grid)
# 



# Optimize Based on Top 25000 of Tiling Library:
optimal_tiling_top25000_prelim_AT = find_optimal_aa_AT(tiling, "Tiling", 25000)
optimal_tiling_top25000_prelim_AT$aa_optimal
# 1.0 0.4 1.2
optimal_tiling_top25000_prelim_AT$heatmap
min(optimal_tiling_top25000_prelim_AT$power_grid)
# 0.0002396757
optimal_tiling_top25000_refined_AT = find_optimal_aa_AT(tiling, "Tiling", 25000,
                                                        a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                        a3_values=seq(from=1,to=1.4,by=0.01))
optimal_tiling_top25000_refined_AT$aa_optimal
# 1.00 0.39 1.18
optimal_tiling_top25000_refined_AT$heatmap
min(optimal_tiling_top25000_refined_AT$power_grid)
# 1.217288e-05


# Optimize Based on Top 15000 of Tiling Library:
optimal_tiling_top15000_prelim_AT = find_optimal_aa_AT(tiling, "Tiling", 15000)
optimal_tiling_top15000_prelim_AT$aa_optimal
# 1.0 0.4 1.1
optimal_tiling_top15000_prelim_AT$heatmap
min(optimal_tiling_top15000_prelim_AT$power_grid)
# 
optimal_tiling_top15000_refined_AT = find_optimal_aa_AT(tiling, "Tiling", 15000,
                                                        a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                        a3_values=seq(from=1,to=1.4,by=0.01))
optimal_tiling_top15000_refined_AT$aa_optimal
# 
optimal_tiling_top15000_refined_AT$heatmap
min(optimal_tiling_top15000_refined_AT$power_grid)
# 




# CHR V:

# Optimize Based on Top Half of ChrV Library:
optimal_chrv_tophalf_prelim_AT = find_optimal_aa_AT(chrv, "ChrV", nrow(chrv)/2)
optimal_chrv_tophalf_prelim_AT$aa_optimal
# 
optimal_chrv_tophalf_prelim_AT$heatmap
optimal_chrv_tophalf_refined_AT = find_optimal_aa_AT(chrv, "ChrV", nrow(chrv)/2,
                                                     a2_values=seq(from=2.6,to=3,by=0.01),
                                                     a3_values=seq(from=1.3,to=1.7,by=0.01))
optimal_chrv_tophalf_refined_AT$aa_optimal
# 
optimal_chrv_tophalf_refined_AT$heatmap


# Optimize Based on Top Quartile of ChrV Library:
optimal_chrv_topquartile_prelim_AT = find_optimal_aa_AT(chrv, "ChrV", nrow(chrv)/4,
                                                        a2_values=seq(from=2,to=4,by=0.1))
optimal_chrv_topquartile_prelim_AT$aa_optimal
# 
optimal_chrv_topquartile_prelim_AT$heatmap
optimal_chrv_topquartile_refined_AT = find_optimal_aa_AT(chrv, "ChrV", nrow(chrv)/4,
                                                         a2_values=seq(from=1.6,to=2.0,by=0.01),
                                                         a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_chrv_topquartile_refined_AT$aa_optimal
# 
optimal_chrv_topquartile_refined_AT$heatmap


# Optimize Based on Top 1000 of ChrV Library:
optimal_chrv_top1000_prelim_AT = find_optimal_aa_AT(chrv, "ChrV")
optimal_chrv_top1000_prelim_AT$aa_optimal
# 
optimal_chrv_top1000_prelim_AT$heatmap
min(optimal_chrv_top1000_prelim_AT$power_grid)
# 
optimal_chrv_top1000_refined_AT = find_optimal_aa_AT(chrv, "ChrV",
                                                     a2_values=seq(from=1.1,to=1.5,by=0.01),
                                                     a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_chrv_top1000_refined_AT$aa_optimal
# 
optimal_chrv_top1000_refined_AT$heatmap
min(optimal_chrv_top1000_refined_AT$power_grid)
# 



# Optimize Based on Top 25000 of ChrV Library:
optimal_chrv_top25000_prelim_AT = find_optimal_aa_AT(chrv, "ChrV", 25000)
optimal_chrv_top25000_prelim_AT$aa_optimal
# 1.0 0.4 1.1
optimal_chrv_top25000_prelim_AT$heatmap
min(optimal_chrv_top25000_prelim_AT$power_grid)
# 0.002613723
optimal_chrv_top25000_refined_AT = find_optimal_aa_AT(chrv, "ChrV", 25000,
                                                      a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                      a3_values=seq(from=0.9,to=1.3,by=0.01))
optimal_chrv_top25000_refined_AT$aa_optimal
# 1.00 0.35 1.11
optimal_chrv_top25000_refined_AT$heatmap
min(optimal_chrv_top25000_refined_AT$power_grid)
# 1.808002e-06



# Optimize Based on Top 15000 of ChrV Library:
optimal_chrv_top15000_prelim_AT = find_optimal_aa_AT(chrv, "ChrV", 15000)
optimal_chrv_top15000_prelim_AT$aa_optimal
# 1.0 0.4 1.1
optimal_chrv_top15000_prelim_AT$heatmap
min(optimal_chrv_top15000_prelim_AT$power_grid)
# 0.005825466
optimal_chrv_top15000_refined_AT = find_optimal_aa_AT(chrv, "ChrV", 15000,
                                                      a2_values=seq(from=0.2,to=0.6,by=0.01),
                                                      a3_values=seq(from=0.9,to=1.3,by=0.01))
optimal_chrv_top15000_refined_AT$aa_optimal
# 1.00 0.35 1.11
optimal_chrv_top15000_refined_AT$heatmap
min(optimal_chrv_top15000_refined_AT$power_grid)
# 3.521615e-05





# COMBINE RESULTS FROM RANDOM AND TILING:

# Combine Random top 4000 and Tiling top 25000:
tiling_top25000_random_top4000_min_indicies = which(optimal_tiling_top25000_prelim_AT$power_grid + 
                                                      optimal_random_top4000_prelim_AT$power_grid == 
                                                      min(optimal_tiling_top25000_prelim_AT$power_grid + 
                                                            optimal_random_top4000_prelim_AT$power_grid), 
                                                    arr.ind=TRUE)
c(rownames(optimal_tiling_top25000_prelim_AT$power_grid)[tiling_top25000_random_top4000_min_indicies[1,1]], 
  colnames(optimal_tiling_top25000_prelim_AT$power_grid)[tiling_top25000_random_top4000_min_indicies[1,2]])
# "0.7" "1.3"

# Refine Random:
optimal_random_top4000_refined_combined_AT = find_optimal_aa_AT(random, "Random", 4000,
                                                                a2_values=seq(from=0.5,to=0.9,by=0.01),
                                                                a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_random_top4000_refined_combined_AT$aa_optimal
# 1.00 0.90 1.32
optimal_random_top4000_refined_combined_AT$heatmap
min(optimal_random_top4000_refined_combined_AT$power_grid)
# 0.00326381

# Refine Tiling:
optimal_tiling_top25000_refined_combined_AT = find_optimal_aa_AT(tiling, "Tiling", 25000,
                                                                 a2_values=seq(from=0.5,to=0.9,by=0.01),
                                                                 a3_values=seq(from=1.1,to=1.5,by=0.01))
optimal_tiling_top25000_refined_combined_AT$aa_optimal
# 1.00 0.50 1.23
optimal_tiling_top25000_refined_combined_AT$heatmap
min(optimal_tiling_top25000_refined_combined_AT$power_grid)
# 0.004161046

tiling_top25000_random_top4000_min_indicies_refined = which(optimal_tiling_top25000_refined_combined_AT$power_grid + 
                                                              optimal_random_top4000_refined_combined_AT$power_grid == 
                                                              min(optimal_tiling_top25000_refined_combined_AT$power_grid + 
                                                                    optimal_random_top4000_refined_combined_AT$power_grid), 
                                                            arr.ind=TRUE)
c(rownames(optimal_tiling_top25000_refined_combined_AT$power_grid)[tiling_top25000_random_top4000_min_indicies_refined[1,1]], 
  colnames(optimal_tiling_top25000_refined_combined_AT$power_grid)[tiling_top25000_random_top4000_min_indicies_refined[1,2]])
# "0.64" "1.29"




# Combine Random top 2000 and Tiling top 15000:
tiling_top15000_random_top2000_min_indicies = which(optimal_tiling_top15000_prelim_AT$power_grid + 
                                                      optimal_random_top2000_prelim_AT$power_grid == 
                                                      min(optimal_tiling_top15000_prelim_AT$power_grid + 
                                                            optimal_random_top2000_prelim_AT$power_grid), 
                                                    arr.ind=TRUE)
c(rownames(optimal_tiling_top15000_prelim_AT$power_grid)[tiling_top15000_random_top2000_min_indicies[1,1]], 
  colnames(optimal_tiling_top15000_prelim_AT$power_grid)[tiling_top15000_random_top2000_min_indicies[1,2]])
# "0.7" "1.2"

# Refine Random:
optimal_random_top2000_refined_combined_AT = find_optimal_aa_AT(random, "Random", 2000,
                                                                a2_values=seq(from=0.5,to=0.9,by=0.01),
                                                                a3_values=seq(from=1,to=1.4,by=0.01))
optimal_random_top2000_refined_combined_AT$aa_optimal
# 1.00 0.89 1.28
optimal_random_top2000_refined_combined_AT$heatmap
min(optimal_random_top2000_refined_combined_AT$power_grid)
# 0.005717874

# Refine Tiling:
optimal_tiling_top15000_refined_combined_AT = find_optimal_aa_AT(tiling, "Tiling", 15000,
                                                                 a2_values=seq(from=0.5,to=0.9,by=0.01),
                                                                 a3_values=seq(from=1,to=1.4,by=0.01))
optimal_tiling_top15000_refined_combined_AT$aa_optimal
# 1.00 0.50 1.19
optimal_tiling_top15000_refined_combined_AT$heatmap
min(optimal_tiling_top15000_refined_combined_AT$power_grid)
# 0.007916504

tiling_top15000_random_top2000_min_indicies_refined = which(optimal_tiling_top15000_refined_combined_AT$power_grid + 
                                                              optimal_random_top2000_refined_combined_AT$power_grid == 
                                                              min(optimal_tiling_top15000_refined_combined_AT$power_grid + 
                                                                    optimal_random_top2000_refined_combined_AT$power_grid), 
                                                            arr.ind=TRUE)
c(rownames(optimal_tiling_top15000_refined_combined_AT$power_grid)[tiling_top15000_random_top2000_min_indicies_refined[1,1]], 
  colnames(optimal_tiling_top15000_refined_combined_AT$power_grid)[tiling_top15000_random_top2000_min_indicies_refined[1,2]])
# "0.68" "1.24"





# COMBINE RESULTS FROM TILING, RANDOM, AND CHR V:

# Combine Random top 4000, Tiling top 25000, and ChrV top 25000:
tiling_top25000_chrv_top25000_random_top4000_min_indicies = which(optimal_tiling_top25000_prelim_AT$power_grid +
                                                                    optimal_chrv_top25000_prelim_AT$power_grid + 
                                                                    optimal_random_top4000_prelim_AT$power_grid == 
                                                                    min(optimal_tiling_top25000_prelim_AT$power_grid +
                                                                          optimal_chrv_top25000_prelim_AT$power_grid + 
                                                                          optimal_random_top4000_prelim_AT$power_grid), 
                                                                  arr.ind=TRUE)
c(rownames(optimal_chrv_top25000_prelim_AT$power_grid)[tiling_top25000_chrv_top25000_random_top4000_min_indicies[1,1]], 
  colnames(optimal_chrv_top25000_prelim_AT$power_grid)[tiling_top25000_chrv_top25000_random_top4000_min_indicies[1,2]])
# "0.5" "1.2"

# Refine Random:
optimal_random_top4000_refined_combined2_AT = find_optimal_aa_AT(random, "Random", 4000,
                                                                 a2_values=seq(from=0.3,to=0.7,by=0.01),
                                                                 a3_values=seq(from=1,to=1.4,by=0.01))
optimal_random_top4000_refined_combined2_AT$aa_optimal
# 1.00 0.70 1.31
optimal_random_top4000_refined_combined2_AT$heatmap
min(optimal_random_top4000_refined_combined2_AT$power_grid)
# 0.009954119

# Refine Tiling:
optimal_tiling_top25000_refined_combined2_AT = find_optimal_aa_AT(tiling, "Tiling", 25000,
                                                                  a2_values=seq(from=0.3,to=0.7,by=0.01),
                                                                  a3_values=seq(from=1,to=1.4,by=0.01))
optimal_tiling_top25000_refined_combined2_AT$aa_optimal
# 1.00 0.39 1.18
optimal_tiling_top25000_refined_combined2_AT$heatmap
min(optimal_tiling_top25000_refined_combined2_AT$power_grid)
# 1.217288e-05

# Refine ChrV:
optimal_chrv_top25000_refined_combined2_AT = find_optimal_aa_AT(chrv, "ChrV", 25000,
                                                                a2_values=seq(from=0.3,to=0.7,by=0.01),
                                                                a3_values=seq(from=1,to=1.4,by=0.01))
optimal_chrv_top25000_refined_combined2_AT$aa_optimal
# 1.00 0.35 1.11
optimal_chrv_top25000_refined_combined2_AT$heatmap
min(optimal_chrv_top25000_refined_combined2_AT$power_grid)
# 1.808002e-06

tiling_top25000_chrv_top25000_random_top4000_min_indicies_refined = 
  which(optimal_tiling_top25000_refined_combined2_AT$power_grid +
          optimal_chrv_top25000_refined_combined2_AT$power_grid + 
          optimal_random_top4000_refined_combined2_AT$power_grid == 
          min(optimal_chrv_top25000_refined_combined2_AT$power_grid + 
                optimal_random_top4000_refined_combined2_AT$power_grid +
                optimal_tiling_top25000_refined_combined2_AT$power_grid), 
        arr.ind=TRUE)
c(rownames(optimal_chrv_top25000_refined_combined2_AT$power_grid)[tiling_top25000_chrv_top25000_random_top4000_min_indicies_refined[1,1]], 
  colnames(optimal_chrv_top25000_refined_combined2_AT$power_grid)[tiling_top25000_chrv_top25000_random_top4000_min_indicies_refined[1,2]])
# "0.55" "1.22"

tiling_top25000_chrv_top25000_random_top4000_gg_format = expand.grid(
  a2=as.numeric(rownames(optimal_tiling_top25000_refined_combined2_AT$power_grid)),
  a3=as.numeric(colnames(optimal_tiling_top25000_refined_combined2_AT$power_grid))
)
tiling_top25000_chrv_top25000_random_top4000_gg_format$weight = 
  c(optimal_tiling_top25000_refined_combined2_AT$power_grid +
      optimal_chrv_top25000_refined_combined2_AT$power_grid + 
      optimal_random_top4000_refined_combined2_AT$power_grid)
ggplot(tiling_top25000_chrv_top25000_random_top4000_gg_format, aes(a2, a3)) + 
  geom_tile(aes(fill=weight)) +
  scale_fill_gradient2("Sum Fourier Power", low="blue", high="red", 
                       midpoint=mean(tiling_top25000_chrv_top25000_random_top4000_gg_format$weight)) +
  ggtitle(paste0("Top 25000 ChrV + Top 25000 Tiling + Top 4000 Random"))


# Combine Random top 2000, Tiling top 15000, and ChrV top 15000:
tiling_top15000_chrv_top15000_random_top2000_min_indicies = which(optimal_tiling_top15000_prelim_AT$power_grid +
                                                                    optimal_chrv_top15000_prelim_AT$power_grid + 
                                                                    optimal_random_top2000_prelim_AT$power_grid == 
                                                                    min(optimal_tiling_top15000_prelim_AT$power_grid +
                                                                          optimal_chrv_top15000_prelim_AT$power_grid + 
                                                                          optimal_random_top2000_prelim_AT$power_grid), 
                                                                  arr.ind=TRUE)
c(rownames(optimal_chrv_top15000_prelim_AT$power_grid)[tiling_top15000_chrv_top15000_random_top2000_min_indicies[1,1]], 
  colnames(optimal_chrv_top15000_prelim_AT$power_grid)[tiling_top15000_chrv_top15000_random_top2000_min_indicies[1,2]])
# "0.6" "1.2"

# Refine Random:
optimal_random_top2000_refined_combined2_AT = find_optimal_aa_AT(random, "Random", 2000,
                                                                 a2_values=seq(from=0.4,to=0.8,by=0.01),
                                                                 a3_values=seq(from=1,to=1.4,by=0.01))
optimal_random_top2000_refined_combined2_AT$aa_optimal
# 1.00 0.80 1.25
optimal_random_top2000_refined_combined2_AT$heatmap
min(optimal_random_top2000_refined_combined2_AT$power_grid)
# 0.01082304

# Refine Tiling:
optimal_tiling_top15000_refined_combined2_AT = find_optimal_aa_AT(tiling, "Tiling", 15000,
                                                                  a2_values=seq(from=0.4,to=0.8,by=0.01),
                                                                  a3_values=seq(from=1,to=1.4,by=0.01))
optimal_tiling_top15000_refined_combined2_AT$aa_optimal
# 1.00 0.40 1.14
optimal_tiling_top15000_refined_combined2_AT$heatmap
min(optimal_tiling_top15000_refined_combined2_AT$power_grid)
# 2.192235e-05

# Refine ChrV:
optimal_chrv_top15000_refined_combined2_AT = find_optimal_aa_AT(chrv, "ChrV", 15000,
                                                                a2_values=seq(from=0.4,to=0.8,by=0.01),
                                                                a3_values=seq(from=1,to=1.4,by=0.01))
optimal_chrv_top15000_refined_combined2_AT$aa_optimal
# 1.00 0.40 1.15
optimal_chrv_top15000_refined_combined2_AT$heatmap
min(optimal_chrv_top15000_refined_combined2_AT$power_grid)
# 0.003047803

tiling_top15000_chrv_top15000_random_top2000_min_indicies_refined = 
  which(optimal_tiling_top15000_refined_combined2_AT$power_grid +
          optimal_chrv_top15000_refined_combined2_AT$power_grid + 
          optimal_random_top2000_refined_combined2_AT$power_grid == 
          min(optimal_tiling_top15000_refined_combined2_AT$power_grid + 
                optimal_chrv_top15000_refined_combined2_AT$power_grid +
                optimal_random_top2000_refined_combined2_AT$power_grid), 
        arr.ind=TRUE)
c(rownames(optimal_chrv_top15000_refined_combined2_AT$power_grid)[tiling_top15000_chrv_top15000_random_top2000_min_indicies_refined[1,1]], 
  colnames(optimal_chrv_top15000_refined_combined2_AT$power_grid)[tiling_top15000_chrv_top15000_random_top2000_min_indicies_refined[1,2]])
# "0.59" "1.2" 

tiling_top15000_chrv_top15000_random_top2000_gg_format = expand.grid(
  a2=as.numeric(rownames(optimal_tiling_top15000_refined_combined2_AT$power_grid)),
  a3=as.numeric(colnames(optimal_tiling_top15000_refined_combined2_AT$power_grid))
)
tiling_top15000_chrv_top15000_random_top2000_gg_format$weight = 
  c(optimal_tiling_top15000_refined_combined2_AT$power_grid +
      optimal_chrv_top15000_refined_combined2_AT$power_grid + 
      optimal_random_top2000_refined_combined2_AT$power_grid)
ggplot(tiling_top15000_chrv_top15000_random_top2000_gg_format, aes(a2, a3)) + 
  geom_tile(aes(fill=weight)) +
  scale_fill_gradient2("Sum Fourier Power", low="blue", high="red", 
                       midpoint=mean(tiling_top15000_chrv_top15000_random_top2000_gg_format$weight)) +
  ggtitle(paste0("Top 15000 ChrV + Top 15000 Tiling + Top 2000 Random"))




# Combine Random top 2000 and ChrV top 15000:
chrv_top15000_random_top2000_min_indicies = which(optimal_chrv_top15000_prelim_AT$power_grid + 
                                                    optimal_random_top2000_prelim_AT$power_grid == 
                                                    min(optimal_chrv_top15000_prelim_AT$power_grid + 
                                                          optimal_random_top2000_prelim_AT$power_grid), 
                                                  arr.ind=TRUE)
c(rownames(optimal_chrv_top15000_prelim_AT$power_grid)[chrv_top15000_random_top2000_min_indicies[1,1]], 
  colnames(optimal_chrv_top15000_prelim_AT$power_grid)[chrv_top15000_random_top2000_min_indicies[1,2]])
# 

# Refine Random:
# optimal_random_top2000_refined_combined_AT = find_optimal_aa_AT(random, "Random", 2000,
#                                                                 a2_values=seq(from=0.5,to=0.9,by=0.01),
#                                                                 a3_values=seq(from=1,to=1.4,by=0.01))
# optimal_random_top2000_refined_combined_AT$aa_optimal
# # 1.00 0.89 1.28
# optimal_random_top2000_refined_combined_AT$heatmap
# min(optimal_random_top2000_refined_combined_AT$power_grid)
# # 0.005717874

# Refine ChrV:
optimal_chrv_top15000_refined_combined_AT = find_optimal_aa_AT(chrv, "ChrV", 15000,
                                                               a2_values=seq(from=0.5,to=0.9,by=0.01),
                                                               a3_values=seq(from=1,to=1.4,by=0.01))
optimal_chrv_top15000_refined_combined_AT$aa_optimal
# 
optimal_chrv_top15000_refined_combined_AT$heatmap
min(optimal_chrv_top15000_refined_combined_AT$power_grid)
# 

chrv_top15000_random_top2000_min_indicies_refined = which(optimal_chrv_top15000_refined_combined_AT$power_grid + 
                                                            optimal_random_top2000_refined_combined_AT$power_grid == 
                                                            min(optimal_chrv_top15000_refined_combined_AT$power_grid + 
                                                                  optimal_random_top2000_refined_combined_AT$power_grid), 
                                                          arr.ind=TRUE)
c(rownames(optimal_chrv_top15000_refined_combined_AT$power_grid)[chrv_top15000_random_top2000_min_indicies_refined[1,1]], 
  colnames(optimal_chrv_top15000_refined_combined_AT$power_grid)[chrv_top15000_random_top2000_min_indicies_refined[1,2]])
# 



