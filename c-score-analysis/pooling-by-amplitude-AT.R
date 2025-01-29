ps1 <- paste0("X", 1:50, "mono")

# Fourier for each sequence individually:
find_all_ps1_fourier = function(data) {
  data = data %>% select(all_of(ps1))
  Xone_AorT_data = matrix(nrow=nrow(data), ncol=ncol(data))
  colnames(Xone_AorT_data) = colnames(data)
  Xone_AorT_data[] = ((data == "A") | (data == "T")) %>% 
    as.matrix() %>% as.numeric()
  
  temp_fourier = t(apply(Xone_AorT_data, 1, fft))[,6]
  return(temp_fourier)
}

# Fourier for sequences collapsed:
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


all_ps1_fourier_random = find_all_ps1_fourier(random)
all_ps1_fourier_tiling = find_all_ps1_fourier(tiling)
all_ps1_fourier_chrv = find_all_ps1_fourier(chrv)

num_groups = 20
# amp_lower_bounds = seq(from=0, to=8, length.out=num_groups+1)[-(num_groups+1)]
# amp_upper_bounds = seq(from=0, to=8, length.out=num_groups+1)[-1]
# library(emdbook)
# amp_lower_bounds = lseq(from=1e-2, to=5, length.out=num_groups+1)[-(num_groups+1)]
# amp_upper_bounds = lseq(from=1e-2, to=5, length.out=num_groups+1)[-1]
amp_lower_bounds = map((0:num_groups)/20, ~quantile(Mod(all_ps1_fourier_tiling), .x))[-(num_groups+1)] %>% unlist()
amp_upper_bounds = map((0:num_groups)/20, ~quantile(Mod(all_ps1_fourier_tiling), .x))[-1] %>% unlist()

amp_centers = map2(amp_lower_bounds, amp_upper_bounds, ~mean(c(.x, .y))) %>%
  unlist()

amp_lower_bounds[1] = -Inf
amp_upper_bounds[num_groups] = Inf


group_assignment_list_random = map2(amp_lower_bounds, amp_upper_bounds, ~(Mod(all_ps1_fourier_random) >= .x) & 
                                      (Mod(all_ps1_fourier_random) < .y))
group_assignment_list_tiling = map2(amp_lower_bounds, amp_upper_bounds, ~(Mod(all_ps1_fourier_tiling) >= .x) & 
                                      (Mod(all_ps1_fourier_tiling) < .y))
group_assignment_list_chrv = map2(amp_lower_bounds, amp_upper_bounds, ~(Mod(all_ps1_fourier_chrv) >= .x) & 
                                    (Mod(all_ps1_fourier_chrv) < .y))

# Size of each group:
plot(sapply(group_assignment_list_random, sum))
plot(sapply(group_assignment_list_tiling, sum))
plot(sapply(group_assignment_list_chrv, sum))

grouped_df_list_random = map(group_assignment_list_random, ~random[.x,])
grouped_df_list_tiling = map(group_assignment_list_tiling, ~tiling[.x,])
grouped_df_list_chrv = map(group_assignment_list_chrv, ~chrv[.x,])

grouped_all_ps1_fourier_random = sapply(grouped_df_list_random, find_all_ps1_fourier)
grouped_all_ps1_fourier_tiling = sapply(grouped_df_list_tiling, find_all_ps1_fourier)
grouped_all_ps1_fourier_chrv = sapply(grouped_df_list_chrv, find_all_ps1_fourier)

# Plot mean phase of each sequence in group:
plot(map(grouped_all_ps1_fourier_random, ~mean(Arg(.x)))%>%unlist())
plot(map(grouped_all_ps1_fourier_tiling, ~mean(Arg(.x)))%>%unlist())
plot(map(grouped_all_ps1_fourier_chrv, ~mean(Arg(.x)))%>%unlist())

# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
#      map(grouped_all_ps1_fourier_random, ~mean(Mod(.x)))%>%unlist(),
#      xlab = "Mean C26 Value by group",
#      ylab = "Mean Amplitude")
# title("random C26")
# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
#      map(grouped_all_ps1_fourier_random, ~mean(Mod(.x)))%>%unlist(),
#      xlab = "Mean C29 Value by group",
#      ylab = "Mean Amplitude")
# title("random C29")
# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
#      map(grouped_all_ps1_fourier_random, ~mean(Mod(.x)))%>%unlist(),
#      xlab = "Mean C31 Value by group",
#      ylab = "Mean Amplitude")
# title("random C31")

plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("random C26")
plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("random C29")
plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("random C31")
plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("random C0")

plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C26")
plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C29")
plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C31")
plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C0")

plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C26")
plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C29")
plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C31")
plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C0")

grouped_ps1_fourier_random = sapply(grouped_df_list_random, find_ps1_fourier)
grouped_ps1_fourier_tiling = sapply(grouped_df_list_tiling, find_ps1_fourier)
grouped_ps1_fourier_chrv = sapply(grouped_df_list_chrv, find_ps1_fourier)

# Plot phase of collapsed sequences by group:
plot(Arg(grouped_ps1_fourier_random))
plot(Arg(grouped_ps1_fourier_tiling))
plot(Arg(grouped_ps1_fourier_chrv))






# Given an A, find C0 and phi for each sequence:
find_c0phi <- function(dat, A, B=c(1,1,1)){
  # Find C0 and Amplitude for each sequence
  # dat: Each row corresponds to a sequence (only for one group).
  #   First three columns need to be C26, C29, and C31, in that order
  # 
  # Returns: columns C0, A
  
  # mat <- matrix(0, nrow = 3, ncol = 2)
  k <- 2*pi/10.4
  n <- c(26, 29, 31)
  
  # LEAST-SQUARES SOLUTION:
  mat <- matrix(0, nrow=3, ncol=3)
  mat[1:3,1] <- 1
  mat[1:3,2] <- B*A*sin(n*k)
  mat[1:3,3] <- B*A*cos(n*k)

  
  mat2 <- solve(t(mat)%*%mat)%*%t(mat)
  
  c0cosphisinphi <-as.matrix(dat[,1:3]) %*% t(mat2)
  
  c0phi <- matrix(0, nrow=nrow(c0cosphisinphi), ncol=2)
  c0phi[,1] <- c0cosphisinphi[,1]
  c0phi[,2] <- atan2(c0cosphisinphi[,3], c0cosphisinphi[,2])
  # c0phiphi <- c0cosphisinphi
  # c0phiphi[,2] <- acos(c0cosphisinphi[,2])
  # c0phiphi[,3] <- asin(c0cosphisinphi[,3])
  # print(c0phiphi)
  # print(c0phi)
  
  # Check:
  # if (c0phiphi[,2] != c0phiphi[,3]) {
  #   print(c0cosphisinphi)
  #   print(c0phiphi)
  # }
  
  # c0phi = c0phiphi[,1:2]
  
  # CLOSED-FORM dS0/dC0=0 and dS0/dA=0 SOLUTION (equivalent)
  # cnbar = apply(dat[,1:3], 1, mean)
  # sinsum = sum(sin(n*k + phi))
  # sin2sum = sum(sin(n*k + phi)^2)
  # sin_dev = t(apply(dat[,1:3]-cnbar, 1, function(row) {row*sin(n*k+phi)}))
  # sin_dev_sum = apply(sin_dev, 1, sum)
  # c0 = cnbar - (sinsum*sin_dev_sum)/(3*sin2sum - sinsum^2)
  # A = 3*sin_dev_sum/(3*sin2sum - sinsum^2)
  # 
  # # c0A = matrix(nrow=length(c0), ncol=2)
  # # c0A[,1] = c0
  # # c0A[,2] = A
  
  # If A < 0, then A <- 0 and C0 <- mean(Cn)
  # negA_rowidx = c0A[,2] < 0
  # c0A[negA_rowidx,1] <- apply(dat[negA_rowidx, 1:3], 1, mean)
  # c0A[negA_rowidx,2] <- 0
  return(c0phi)
}



# Try brute-forcing our way to optimal a_k (for a given group):
find_optimal_a_AT = function(data, a_values=seq(from=0.01,to=2,by=0.01)) {
  
  c0phi_list = map(a_values, ~find_c0phi(data %>% select(C26, C29, C31), .x))
  
  k = 2*pi/10.4
  loss_by_component_list = map2(c0phi_list, a_values, ~c(sum((data %>% select(C26) - (.x[,1] + .y*sin(k*26+.x[,2])))^2),
                                                         sum((data %>% select(C29) - (.x[,1] + .y*sin(k*29+.x[,2])))^2),
                                                         sum((data %>% select(C31) - (.x[,1] + .y*sin(k*31+.x[,2])))^2)))
  
  total_loss_list = map(loss_by_component_list, ~mean(.x)/dim(data)[1]) %>% unlist()
  
  to_return = list()
  
  optimal_index = which(total_loss_list == min(total_loss_list))
  
  if(length(optimal_index) > 1) {
    optimal_index=optimal_index[1]
  }
  
  to_return$optimal_a = a_values[optimal_index]
  to_return$c0phi = c0phi_list[[optimal_index]]
  to_return$a_values = a_values
  to_return$total_loss_list = total_loss_list
  to_return$loss_by_componenet_list = loss_by_component_list
  
  return(to_return)
}

optimize_a_list_random = lapply(grouped_df_list_random, find_optimal_a_AT)
optimize_a_list_tiling = lapply(grouped_df_list_tiling, find_optimal_a_AT)
optimize_a_list_chrv = lapply(grouped_df_list_chrv, find_optimal_a_AT)


# Optimal a for each group:
plot(map(optimize_a_list_random, ~.x$optimal_a) %>% unlist(),
     ylab="optimal a by group, random",
     ylim=c(-pi,pi))
plot(map(optimize_a_list_tiling, ~.x$optimal_a) %>% unlist(),
     ylab="optimal a by group, tiling",
     ylim=c(-pi,pi))
plot(map(optimize_a_list_chrv, ~.x$optimal_a) %>% unlist(),
     ylab="optimal a by group, chrv",
     ylim=c(-pi,pi))


ylims=c(0,0.15)
# ylims=c(0,0.5)
# Plot loss vs a for each group
map(1:length(optimize_a_list_random), ~plot(optimize_a_list_random[[.x]]$a_values, 
                                              optimize_a_list_random[[.x]]$total_loss_list,
                                              ylab=paste(c("Loss for group ", .x)),
                                              xlab="A Values, random",
                                              ylim=ylims))
map(1:length(optimize_a_list_tiling), ~plot(optimize_a_list_tiling[[.x]]$a_values, 
                                              optimize_a_list_tiling[[.x]]$total_loss_list,
                                              ylab=paste(c("Loss for group ", .x)),
                                              xlab="A Values, tiling",
                                              ylim=ylims))
map(1:length(optimize_a_list_chrv), ~plot(optimize_a_list_chrv[[.x]]$a_values, 
                                            optimize_a_list_chrv[[.x]]$total_loss_list,
                                            ylab=paste(c("Loss for group ", .x)),
                                            xlab="A Values, chrv",
                                            ylim=ylims))


# Plot loss by component vs a for each group:
plot_loss_by_component = function(optimize_a_obj, group_number, library_name, ylims) {
  plot(optimize_a_obj$a_values,
       matrix(optimize_a_obj$loss_by_componenet_list %>% unlist(), nrow=3)[1,]/length(optimize_a_obj$a_values),
       col="blue",
       ylab=paste(c("Loss for group ", group_number)),
       xlab=paste(c("A Values"), library_name),
       ylim=ylims)
  points(optimize_a_obj$a_values,
         matrix(optimize_a_obj$loss_by_componenet_list %>% unlist(), nrow=3)[2,]/length(optimize_a_obj$a_values),
         col="red")
  points(optimize_a_obj$a_values,
         matrix(optimize_a_obj$loss_by_componenet_list %>% unlist(), nrow=3)[3,]/length(optimize_a_obj$a_values),
         col="green")
  legend("topright", c("26", "29", "31"), fill=c("blue", "red", "green"))
}
map(1:length(optimize_a_list_random), ~plot_loss_by_component(optimize_a_list_random[[.x]], 
                                                                .x, "random", c(0,3)))
map(1:length(optimize_a_list_tiling), ~plot_loss_by_component(optimize_a_list_tiling[[.x]], 
                                                                .x, "tiling", c(0,30)))
map(1:length(optimize_a_list_chrv), ~plot_loss_by_component(optimize_a_list_chrv[[.x]], 
                                                              .x, "chrv", c(0, 40)))



plot(map(optimize_a_list_random, ~mean(.x$c0phi[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("random new C0")

plot(map(optimize_a_list_tiling, ~mean(.x$c0phi[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling new C0")

plot(map(optimize_a_list_chrv, ~mean(.x$c0phi[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv new C0")



for(i in 1:length(grouped_df_list_random)) {
  grouped_df_list_random[[i]]$C0_grouped_a = optimize_a_list_random[[i]]$c0phi[,1]
  grouped_df_list_random[[i]]$phi_grouped_a = optimize_a_list_random[[i]]$c0phi[,2]
}
grouped_a_df_random = bind_rows(grouped_df_list_random)

for(i in 1:length(grouped_df_list_tiling)) {
  grouped_df_list_tiling[[i]]$C0_grouped_a = optimize_a_list_tiling[[i]]$c0phi[,1]
  grouped_df_list_tiling[[i]]$phi_grouped_a = optimize_a_list_tiling[[i]]$c0phi[,2]
}
grouped_a_df_tiling = bind_rows(grouped_df_list_tiling)

for(i in 1:length(grouped_df_list_chrv)) {
  grouped_df_list_chrv[[i]]$C0_grouped_a = optimize_a_list_chrv[[i]]$c0phi[,1]
  grouped_df_list_chrv[[i]]$phi_grouped_a = optimize_a_list_chrv[[i]]$c0phi[,2]
}
grouped_a_df_chrv = bind_rows(grouped_df_list_chrv)






# ALLOW FOR DEPENDENCE OF A ON n:

# Given a multiplier B for C31 amplitude, find the correlation between 
# approximate C0 and Fourier amplitude
corr_C0_fourier_amp = function(dat, B, fourier_amp) {
  approxC0 = (dat$C31/B + dat$C26)/(1+1/B)
  return(cor(approxC0, fourier_amp))
}

library(emdbook)
# B_values = seq(.25, 4, by=0.1)
# B_values = lseq(.1, 10, length.out=49)
B_values = lseq(.1, 10, length.out=100)
corr_C0_fourier_amp_all_random = map(B_values, 
                                     ~corr_C0_fourier_amp(random, .x, 
                                                          Mod(all_ps1_fourier_random))) %>%
  unlist()
corr_C0_fourier_amp_all_tiling = map(B_values, 
                                     ~corr_C0_fourier_amp(tiling, .x, 
                                                          Mod(all_ps1_fourier_tiling))) %>%
  unlist()
corr_C0_fourier_amp_all_chrv = map(B_values, 
                                   ~corr_C0_fourier_amp(chrv, .x, 
                                                        Mod(all_ps1_fourier_chrv))) %>%
  unlist()



plot(B_values, corr_C0_fourier_amp_all_random,
     ylab = "Correlation", log='x')
abline(v=B_values[which.max(corr_C0_fourier_amp_all_random)])
title("random")
plot(B_values, corr_C0_fourier_amp_all_tiling,
     ylab = "Correlation", log='x')
abline(v=B_values[which.max(corr_C0_fourier_amp_all_tiling)])
title("tiling")
plot(B_values, corr_C0_fourier_amp_all_chrv,
     ylab = "Correlation", log='x')
abline(v=B_values[which.max(corr_C0_fourier_amp_all_chrv)])
title("chrv")

optimal_B_val_random = B_values[which.max(corr_C0_fourier_amp_all_random)]
optimal_B_val_tiling = B_values[which.max(corr_C0_fourier_amp_all_tiling)]
optimal_B_val_chrv = B_values[which.max(corr_C0_fourier_amp_all_chrv)]

plot()

# At first glance, it seems like the optimal B is very library-dependent
# What if the optimal B is phase-dependent, though?
# Let's try finding the optimal B for each phase-based subgroup

B_values = lseq(.1, 10, length.out=100)
# B_values = lseq(.01, 100, length.out=100)
grouped_corr_C0_fourier_amp_list_random = map2(grouped_df_list_random, 
                                               grouped_all_ps1_fourier_random, 
                                               function(x, y) map(
                                                 B_values, ~corr_C0_fourier_amp(
                                                   x, .x, Mod(y))) %>% unlist())
grouped_corr_C0_fourier_amp_list_tiling = map2(grouped_df_list_tiling, 
                                               grouped_all_ps1_fourier_tiling, 
                                               function(x, y) map(
                                                 B_values, ~corr_C0_fourier_amp(
                                                   x, .x, Mod(y))) %>% unlist())
grouped_corr_C0_fourier_amp_list_chrv = map2(grouped_df_list_chrv, 
                                             grouped_all_ps1_fourier_chrv, 
                                             function(x, y) map(
                                               B_values, ~corr_C0_fourier_amp(
                                                 x, .x, Mod(y))) %>% unlist())


optimal_B_values_random = map(grouped_corr_C0_fourier_amp_list_random, ~B_values[which.max(.x)]) %>% 
  unlist()
optimal_B_values_tiling = map(grouped_corr_C0_fourier_amp_list_tiling, ~B_values[which.max(.x)]) %>% 
  unlist()
optimal_B_values_chrv = map(grouped_corr_C0_fourier_amp_list_chrv, ~B_values[which.max(.x)]) %>% 
  unlist()

plot(optimal_B_values_random, ylab="Optimal B Value", xlab="Group Index",
     log='y')
title("random")
plot(optimal_B_values_tiling, ylab="Optimal B Value", xlab="Group Index",
     log='y')
title("tiling")
plot(optimal_B_values_chrv, ylab="Optimal B Value", xlab="Group Index",
     log='y')
title("chrv")

optimal_B_values_random
optimal_B_values_tiling
optimal_B_values_chrv

# Plot correlation vs B for each group (by library first)
map(1:length(grouped_corr_C0_fourier_amp_list_random), ~plot(
  B_values, grouped_corr_C0_fourier_amp_list_random[[.x]],
  ylab=paste(c("Correlation for group ", .x)),
  xlab="B Values, random", log='x'))
map(1:length(grouped_corr_C0_fourier_amp_list_tiling), ~plot(
  B_values, grouped_corr_C0_fourier_amp_list_tiling[[.x]],
  ylab=paste(c("Correlation for group ", .x)),
  xlab="B Values, tiling", log='x'))
map(1:length(grouped_corr_C0_fourier_amp_list_chrv), ~plot(
  B_values, grouped_corr_C0_fourier_amp_list_chrv[[.x]],
  ylab=paste(c("Correlation for group ", .x)),
  xlab="B Values, chrv", log='x'))

# Plot correlation vs B for each group (by group first)
map(1:length(grouped_corr_C0_fourier_amp_list_random), function(x) {
  par(mfrow=c(1,3))
  plot(
    B_values, grouped_corr_C0_fourier_amp_list_random[[x]],
    ylab=paste("Correlation for group ", x),
    xlab="B Values, random", log='x')
  plot(
    B_values, grouped_corr_C0_fourier_amp_list_tiling[[x]],
    ylab=paste("Correlation for group ", x),
    xlab="B Values, tiling", log='x')
  plot(
    B_values, grouped_corr_C0_fourier_amp_list_chrv[[x]],
    ylab=paste("Correlation for group ", x),
    xlab="B Values, chrv", log='x')
})

group_cur=7
B_cur=B_values[which.max(grouped_corr_C0_fourier_amp_list_random[[group_cur]])]
plot(Mod(grouped_all_ps1_fourier_random[[group_cur]]),
     (grouped_df_list_random[[group_cur]]$C31/B_cur + grouped_df_list_random[[group_cur]]$C26)/
       (1+1/B_cur),
     xlab="Fourier Amplitude",
     ylab="Approx C0")
title(paste0("random, group ", group_cur, " with B=", round(B_cur, 2)))


group_cur=15
ylims=c(-3, 3)
map(round(lseq(1, length(B_values), length.out=10)), 
    ~plot(Mod(grouped_all_ps1_fourier_random[[group_cur]]),
          (grouped_df_list_random[[group_cur]]$C31/.x + grouped_df_list_random[[group_cur]]$C26)/
            (1+1/.x),
          ylim=ylims,
          xlab="Fourier Amplitude",
          ylab="Approx C0",
          main=paste0("random, group ", group_cur, " with B=", round(.x, 2))))

plot_approxC0_vs_fourier_amp = function(dat, B, fourier_amp) {
  approxC0 = (dat$C31/B + dat$C26)/(1+1/B)
}

corr_C0_fourier_amp = function(dat, B, fourier_amp) {
  approxC0 = (dat$C31/B + dat$C26)/(1+1/B)
  return(cor(approxC0, fourier_amp))
}










#FIXME
# Try brute-forcing our way to optimal phi_k (for a given group). Include B values:
find_optimal_phi_AT_include_B = function(data, B_values, phi_values=seq(from=-pi,to=pi,by=pi/40)) {
  
  c0A_list = map(phi_values, ~find_c0A(data %>% select(C26, C29, C31), .x, B_values))
  
  k = 2*pi/10.4
  loss_by_component_list = map2(c0A_list, phi_values, ~c(sum((data %>% select(C26) - (.x[,1] + .x[,2]*sin(k*26+.y)))^2),
                                                         sum((data %>% select(C29) - (.x[,1] + .x[,2]*sin(k*29+.y)))^2),
                                                         sum((data %>% select(C31) - (.x[,1] + .x[,2]*sin(k*31+.y)))^2)))
  
  total_loss_list = map(loss_by_component_list, ~mean(.x)/dim(data)[1]) %>% unlist()
  
  to_return = list()
  
  optimal_index = which(total_loss_list == min(total_loss_list))
  
  if(length(optimal_index) > 1) {
    optimal_index=optimal_index[1]
  }
  
  to_return$optimal_phi = phi_values[optimal_index]
  to_return$c0A = c0A_list[[optimal_index]]
  to_return$phi_values = phi_values
  to_return$total_loss_list = total_loss_list
  to_return$loss_by_componenet_list = loss_by_component_list
  
  return(to_return)
}

# formatted_optimal_B_values_random = map(optimal_B_values_random, ~c(1, (1+.x)/2, .x))
# formatted_optimal_B_values_tiling = map(optimal_B_values_tiling, ~c(1, (1+.x)/2, .x))
# formatted_optimal_B_values_chrv = map(optimal_B_values_chrv, ~c(1, (1+.x)/2, .x))

# formatted_optimal_B_values_random = map(1:length(grouped_df_list_random), 
#                                         ~c(1, (1+optimal_B_val_random)/2, optimal_B_val_random))
# formatted_optimal_B_values_tiling = map(1:length(grouped_df_list_tiling), 
#                                         ~c(1, (1+optimal_B_val_tiling)/2, optimal_B_val_tiling))
# formatted_optimal_B_values_chrv = map(1:length(grouped_df_list_chrv), 
#                                         ~c(1, (1+optimal_B_val_chrv)/2, optimal_B_val_chrv))

formatted_optimal_B_values_random = map(1:length(grouped_df_list_random), 
                                        ~c(1, 0, optimal_B_val_random))
formatted_optimal_B_values_tiling = map(1:length(grouped_df_list_tiling), 
                                        ~c(1, 0, optimal_B_val_tiling))
formatted_optimal_B_values_chrv = map(1:length(grouped_df_list_chrv), 
                                      ~c(1, 0, optimal_B_val_chrv))

# B value optimized for library as a whole - not grouped
optimize_phi_list_groupedB_random = map2(grouped_df_list_random, formatted_optimal_B_values_random, 
                                         ~find_optimal_phi_AT_include_B(.x, .y))
optimize_phi_list_groupedB_tiling = map2(grouped_df_list_tiling, formatted_optimal_B_values_tiling, 
                                         ~find_optimal_phi_AT_include_B(.x, .y))
optimize_phi_list_groupedB_chrv = map2(grouped_df_list_chrv, formatted_optimal_B_values_chrv, 
                                       ~find_optimal_phi_AT_include_B(.x, .y))


# Optimal phi_k for each group:
plot(map(optimize_phi_list_groupedB_random, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, random",
     ylim=c(-pi,pi))
plot(map(optimize_phi_list_groupedB_tiling, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, tiling",
     ylim=c(-pi,pi))
plot(map(optimize_phi_list_groupedB_chrv, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, chrv",
     ylim=c(-pi,pi))


ylims=c(0,0.15)
# ylims=c(0,0.5)
# Plot loss vs phi for each group
map(1:length(optimize_phi_list_groupedB_random), ~plot(optimize_phi_list_groupedB_random[[.x]]$phi_values, 
                                                       optimize_phi_list_groupedB_random[[.x]]$total_loss_list,
                                                       ylab=paste(c("Loss for group ", .x)),
                                                       xlab="Phi Values, random",
                                                       ylim=ylims))
map(1:length(optimize_phi_list_groupedB_tiling), ~plot(optimize_phi_list_groupedB_tiling[[.x]]$phi_values, 
                                                       optimize_phi_list_groupedB_tiling[[.x]]$total_loss_list,
                                                       ylab=paste(c("Loss for group ", .x)),
                                                       xlab="Phi Values, tiling",
                                                       ylim=ylims))
map(1:length(optimize_phi_list_groupedB_chrv), ~plot(optimize_phi_list_groupedB_chrv[[.x]]$phi_values, 
                                                     optimize_phi_list_groupedB_chrv[[.x]]$total_loss_list,
                                                     ylab=paste(c("Loss for group ", .x)),
                                                     xlab="Phi Values, chrv",
                                                     ylim=ylims))


# Plot loss by component vs phi for each group:
plot_loss_by_component = function(optimize_phi_obj, group_number, library_name, ylims) {
  plot(optimize_phi_obj$phi_values,
       matrix(optimize_phi_obj$loss_by_componenet_list %>% unlist(), nrow=3)[1,]/length(optimize_phi_obj$phi_values),
       col="blue",
       ylab=paste(c("Loss for group ", group_number)),
       xlab=paste(c("Phi Values"), library_name),
       ylim=ylims)
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_componenet_list %>% unlist(), nrow=3)[2,]/length(optimize_phi_obj$phi_values),
         col="red")
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_componenet_list %>% unlist(), nrow=3)[3,]/length(optimize_phi_obj$phi_values),
         col="green")
  legend("topright", c("26", "29", "31"), fill=c("blue", "red", "green"))
}
map(1:length(optimize_phi_list_groupedB_random), ~plot_loss_by_component(optimize_phi_list_groupedB_random[[.x]], 
                                                                         .x, "random", c(0,3)))
map(1:length(optimize_phi_list_groupedB_tiling), ~plot_loss_by_component(optimize_phi_list_groupedB_tiling[[.x]], 
                                                                         .x, "tiling", c(0,30)))
map(1:length(optimize_phi_list_groupedB_chrv), ~plot_loss_by_component(optimize_phi_list_groupedB_chrv[[.x]], 
                                                                       .x, "chrv", c(0, 40)))
grouped_df_list_random = map(grouped_df_list_random, ~.x %>%
                               mutate(C0hat1 = mean(C26 + C31),
                                      C0hat2 = C31/optimal_B_val_random + 
                                        C26/(1 + 1/optimal_B_val_random)))
grouped_df_list_tiling = map(grouped_df_list_tiling, ~.x %>%
                               mutate(C0hat1 = mean(C26 + C31),
                                      C0hat2 = C31/optimal_B_val_tiling + 
                                        C26/(1 + 1/optimal_B_val_tiling)))
grouped_df_list_chrv = map(grouped_df_list_chrv, ~.x %>%
                             mutate(C0hat1 = mean(C26 + C31),
                                    C0hat2 = C31/optimal_B_val_chrv + 
                                      C26/(1 + 1/optimal_B_val_chrv)))


plot(1:length(grouped_df_list_random), map(grouped_df_list_random, ~mean(.x$C0hat1)),
     xlab="Group index",
     ylab="C0 hat 1",
     main="random",
     ylim = c(-0.5, 0.25))
plot(1:length(grouped_df_list_random), map(grouped_df_list_random, ~mean(.x$C0hat2)),
     xlab="Group index",
     ylab="C0 hat 2",
     main="random",
     ylim = c(-0.5, 0.25))

plot(1:length(grouped_df_list_tiling), map(grouped_df_list_tiling, ~mean(.x$C0hat1)),
     xlab="Group index",
     ylab="C0 hat 1",
     main="tiling",
     ylim = c(-0.5, 0.25))
plot(1:length(grouped_df_list_tiling), map(grouped_df_list_tiling, ~mean(.x$C0hat2)),
     xlab="Group index",
     ylab="C0 hat 2",
     main="tiling",
     ylim = c(-0.5, 0.25))

plot(1:length(grouped_df_list_chrv), map(grouped_df_list_chrv, ~mean(.x$C0hat1)),
     xlab="Group index",
     ylab="C0 hat 1",
     main="chrv",
     ylim = c(-0.5, 0.25))
plot(1:length(grouped_df_list_chrv), map(grouped_df_list_chrv, ~mean(.x$C0hat2)),
     xlab="Group index",
     ylab="C0 hat 2",
     main="chrv",
     ylim = c(-0.5, 0.25))



plot(map(optimize_phi_list_groupedB_random, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("random new C0")

plot(map(optimize_phi_list_groupedB_tiling, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling new C0")

plot(map(optimize_phi_list_groupedB_chrv, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv new C0")
