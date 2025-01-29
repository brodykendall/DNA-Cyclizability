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

cycle6_inc_fourier = chrv %>%
  select(x50mer, C26, C29, C31, C0) %>%
  rename(Sequence=x50mer, "n=26"=C26, "n=29"=C29, "n=31"=C31)
cycle6_inc_fourier$ps1_fourier_amp = Mod(all_ps1_fourier_chrv)
cycle6_inc_fourier$ps1_fourier_phase = Arg(all_ps1_fourier_chrv)

write.csv(cycle6_inc_fourier, "data/Created/cycle6_inc_fourier.csv", row.names=FALSE)

hist(Mod(all_ps1_fourier_random),xlim=c(0,13),ylim=c(0,0.3),freq=FALSE)
hist(Mod(all_ps1_fourier_tiling),xlim=c(0,13),ylim=c(0,0.3),freq=FALSE)
hist(Mod(all_ps1_fourier_chrv),xlim=c(0,13),ylim=c(0,0.3),freq=FALSE)

quantile(Mod(all_ps1_fourier_random), 0.75)
quantile(Mod(all_ps1_fourier_tiling), 0.75)
quantile(Mod(all_ps1_fourier_chrv), 0.75)

hist(Mod(all_ps1_fourier_random)[which(Mod(all_ps1_fourier_random) > 4)],xlim=c(3,13),ylim=c(0,0.7),freq=FALSE)
hist(Mod(all_ps1_fourier_tiling)[which(Mod(all_ps1_fourier_tiling) > 4)],xlim=c(3,13),ylim=c(0,0.7),freq=FALSE)
hist(Mod(all_ps1_fourier_chrv)[which(Mod(all_ps1_fourier_chrv) > 4)],xlim=c(3,13),ylim=c(0,0.7),freq=FALSE)

mod_cutoffs = seq(0.1, 5, by=0.1)
differences = vector()
for(mod_cutoff in mod_cutoffs) {
  avg_above_cutoff = random[which(Mod(all_ps1_fourier_random) > mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  avg_below_cutoff = random[which(Mod(all_ps1_fourier_random) < mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  differences = append(differences, (avg_above_cutoff - avg_below_cutoff)[[1]])
}
plot(mod_cutoffs,differences)

mod_cutoffs = seq(0.1, 5, by=0.1)
differences = vector()
for(mod_cutoff in mod_cutoffs) {
  avg_above_cutoff = tiling[which(Mod(all_ps1_fourier_tiling) > mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  avg_below_cutoff = tiling[which(Mod(all_ps1_fourier_tiling) < mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  differences = append(differences, (avg_above_cutoff - avg_below_cutoff)[[1]])
}
plot(mod_cutoffs,differences)

mod_cutoffs = seq(0.1, 5, by=0.1)
differences = vector()
for(mod_cutoff in mod_cutoffs) {
  avg_above_cutoff = chrv[which(Mod(all_ps1_fourier_chrv) > mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  avg_below_cutoff = chrv[which(Mod(all_ps1_fourier_chrv) < mod_cutoff), ] %>% 
    mutate(cn_range = abs(C26-C31)) %>% 
    summarise(mean_diff = mean(cn_range))
  differences = append(differences, (avg_above_cutoff - avg_below_cutoff)[[1]])
}
plot(mod_cutoffs,differences)

# mod_cutoff = 1
mod_cutoff = 2

random$below_cutoff = Mod(all_ps1_fourier_random) < mod_cutoff
tiling$below_cutoff = Mod(all_ps1_fourier_tiling) < mod_cutoff
chrv$below_cutoff = Mod(all_ps1_fourier_chrv) < mod_cutoff

num_groups = 20
# num_groups = 10
phase_lower_bounds = seq(from=-pi, to=pi, length.out=num_groups+1)[-(num_groups+1)]
phase_upper_bounds = seq(from=-pi, to=pi, length.out=num_groups+1)[-1]
# phase_lower_bounds = seq(from=-pi, to=pi, length.out=num_groups+1)[-(num_groups+1)] + pi/(2*num_groups)
# phase_upper_bounds = seq(from=-pi, to=pi, length.out=num_groups+1)[-1] + pi/(2*num_groups)

phase_centers = map2(phase_lower_bounds, phase_upper_bounds, ~mean(c(.x, .y))) %>%
  unlist()

phase_lower_bounds[1] = -Inf
phase_upper_bounds[num_groups] = Inf


group_assignment_list_random = map2(phase_lower_bounds, phase_upper_bounds, ~(Arg(all_ps1_fourier_random) >= .x) & 
                                      (Arg(all_ps1_fourier_random) < .y) & random$below_cutoff == FALSE)
group_assignment_list_tiling = map2(phase_lower_bounds, phase_upper_bounds, ~(Arg(all_ps1_fourier_tiling) >= .x) & 
                                      (Arg(all_ps1_fourier_tiling) < .y) & tiling$below_cutoff == FALSE)
group_assignment_list_chrv = map2(phase_lower_bounds, phase_upper_bounds, ~(Arg(all_ps1_fourier_chrv) >= .x) & 
                                    (Arg(all_ps1_fourier_chrv) < .y) & chrv$below_cutoff == FALSE)

get_group_assignment_index = function(group_assignment_list, num_groups, j) {
  for(i in 1:num_groups) {
    if(group_assignment_list[[i]][j]==TRUE) {
      return(i)
    }
  }
  return(num_groups + 1)
}

group_assignment_index_random = map(1:nrow(random), 
                                    ~get_group_assignment_index(
                                      group_assignment_list_random, 
                                      num_groups, 
                                      .x)) %>% unlist()
group_assignment_index_tiling = map(1:nrow(tiling), 
                                    ~get_group_assignment_index(
                                      group_assignment_list_tiling, 
                                      num_groups, 
                                      .x)) %>% unlist()
group_assignment_index_chrv = map(1:nrow(chrv), 
                                  ~get_group_assignment_index(
                                    group_assignment_list_chrv, 
                                    num_groups, 
                                    .x)) %>% unlist()

random$group_assignment_index = group_assignment_index_random
tiling$group_assignment_index = group_assignment_index_tiling
chrv$group_assignment_index = group_assignment_index_chrv

# Size of each group:
plot(sapply(group_assignment_list_random, sum))
plot(sapply(group_assignment_list_tiling, sum))
plot(sapply(group_assignment_list_chrv, sum))

hist(random$group_assignment_index)
hist(tiling$group_assignment_index)
hist(chrv$group_assignment_index)

# grouped_df_list_random = map(group_assignment_list_random, ~random[.x,])
# grouped_df_list_tiling = map(group_assignment_list_tiling, ~tiling[.x,])
# grouped_df_list_chrv = map(group_assignment_list_chrv, ~chrv[.x,])

# grouped_all_ps1_fourier_random = map(1:(num_groups+1), ~find_all_ps1_fourier(random %>% filter(group_assignment_index == .x)))
# grouped_all_ps1_fourier_tiling = map(1:(num_groups+1), ~find_all_ps1_fourier(tiling %>% filter(group_assignment_index == .x)))
# grouped_all_ps1_fourier_chrv = map(1:(num_groups+1), ~find_all_ps1_fourier(chrv %>% filter(group_assignment_index == .x)))
grouped_all_ps1_fourier_random = map(1:(num_groups), ~find_all_ps1_fourier(random %>% filter(group_assignment_index == .x)))
grouped_all_ps1_fourier_tiling = map(1:(num_groups), ~find_all_ps1_fourier(tiling %>% filter(group_assignment_index == .x)))
grouped_all_ps1_fourier_chrv = map(1:(num_groups), ~find_all_ps1_fourier(chrv %>% filter(group_assignment_index == .x)))

# grouped_all_ps1_fourier_random = sapply(grouped_df_list_random, find_all_ps1_fourier)
# grouped_all_ps1_fourier_tiling = sapply(grouped_df_list_tiling, find_all_ps1_fourier)
# grouped_all_ps1_fourier_chrv = sapply(grouped_df_list_chrv, find_all_ps1_fourier)

# Plot mean magnitude of each sequence in group:
plot(map(grouped_all_ps1_fourier_random, ~mean(Mod(.x)))%>%unlist())
plot(map(grouped_all_ps1_fourier_tiling, ~mean(Mod(.x)))%>%unlist())
plot(map(grouped_all_ps1_fourier_chrv, ~mean(Mod(.x)))%>%unlist())

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

plot(map(1:(num_groups+1), ~random %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("random C26")
plot(map(1:(num_groups+1), ~random %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("random C29")
plot(map(1:(num_groups+1), ~random %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("random C31")
plot(map(1:(num_groups+1), ~random %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("random C0")

plot(map(1:(num_groups+1), ~tiling %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C26")
plot(map(1:(num_groups+1), ~tiling %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C29")
plot(map(1:(num_groups+1), ~tiling %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C31")
plot(map(1:(num_groups+1), ~tiling %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling C0")

plot(map(1:(num_groups+1), ~chrv %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C26 = mean(C26)))%>%unlist(),
     ylab = "Mean C26 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C26")
plot(map(1:(num_groups+1), ~chrv %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C29 = mean(C29)))%>%unlist(),
     ylab = "Mean C29 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C29")
plot(map(1:(num_groups+1), ~chrv %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C31 = mean(C31)))%>%unlist(),
     ylab = "Mean C31 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C31")
plot(map(1:(num_groups+1), ~chrv %>% filter(group_assignment_index == .x) %>% 
           summarise(mean_C0 = mean(C0)))%>%unlist(),
     ylab = "Mean C0 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv C0")

# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
#      ylab = "Mean C26 Value by group",
#      ylim = c(-0.5, 0.25))
# title("random C26")
# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
#      ylab = "Mean C29 Value by group",
#      ylim = c(-0.5, 0.25))
# title("random C29")
# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
#      ylab = "Mean C31 Value by group",
#      ylim = c(-0.5, 0.25))
# title("random C31")
# plot(map(grouped_df_list_random, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
#      ylab = "Mean C0 Value by group",
#      ylim = c(-0.5, 0.25))
# title("random C0")

# plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
#      ylab = "Mean C26 Value by group",
#      ylim = c(-0.5, 0.25))
# title("tiling C26")
# plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
#      ylab = "Mean C29 Value by group",
#      ylim = c(-0.5, 0.25))
# title("tiling C29")
# plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
#      ylab = "Mean C31 Value by group",
#      ylim = c(-0.5, 0.25))
# title("tiling C31")
# plot(map(grouped_df_list_tiling, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
#      ylab = "Mean C0 Value by group",
#      ylim = c(-0.5, 0.25))
# title("tiling C0")

# plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C26 = mean(C26)))%>%unlist(),
#      ylab = "Mean C26 Value by group",
#      ylim = c(-0.5, 0.25))
# title("chrv C26")
# plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C29 = mean(C29)))%>%unlist(),
#      ylab = "Mean C29 Value by group",
#      ylim = c(-0.5, 0.25))
# title("chrv C29")
# plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C31 = mean(C31)))%>%unlist(),
#      ylab = "Mean C31 Value by group",
#      ylim = c(-0.5, 0.25))
# title("chrv C31")
# plot(map(grouped_df_list_chrv, ~.x%>%summarise(mean_C0 = mean(C0)))%>%unlist(),
#      ylab = "Mean C0 Value by group",
#      ylim = c(-0.5, 0.25))
# title("chrv C0")

# grouped_ps1_fourier_random = sapply(grouped_df_list_random, find_ps1_fourier)
# grouped_ps1_fourier_tiling = sapply(grouped_df_list_tiling, find_ps1_fourier)
# grouped_ps1_fourier_chrv = sapply(grouped_df_list_chrv, find_ps1_fourier)

# grouped_ps1_fourier_random = sapply(map(1:(num_groups+1), ~random%>%filter(group_assignment_index==.x)), 
#                                     find_ps1_fourier)
# grouped_ps1_fourier_tiling = sapply(map(1:(num_groups+1), ~tiling%>%filter(group_assignment_index==.x)), 
#                                     find_ps1_fourier)
# grouped_ps1_fourier_chrv = sapply(map(1:(num_groups+1), ~chrv%>%filter(group_assignment_index==.x)), 
#                                   find_ps1_fourier)
grouped_ps1_fourier_random = sapply(map(1:(num_groups), ~random%>%filter(group_assignment_index==.x)), 
                                    find_ps1_fourier)
grouped_ps1_fourier_tiling = sapply(map(1:(num_groups), ~tiling%>%filter(group_assignment_index==.x)), 
                                    find_ps1_fourier)
grouped_ps1_fourier_chrv = sapply(map(1:(num_groups), ~chrv%>%filter(group_assignment_index==.x)), 
                                  find_ps1_fourier)
# Plot magnitude of collapsed sequences by group:
plot(Mod(grouped_ps1_fourier_random))
plot(Mod(grouped_ps1_fourier_tiling))
plot(Mod(grouped_ps1_fourier_chrv))






# Given a phi_k, find C0 and A for each sequence:
find_c0A <- function(dat, phi, B=c(1,1,1)){
  # Find C0 and Amplitude for each sequence
  # dat: Each row corresponds to a sequence (only for one group).
  #   First three columns need to be C26, C29, and C31, in that order
  # 
  # Returns: columns C0, A
  
  # mat <- matrix(0, nrow = 3, ncol = 2)
  k <- 2*pi/10.4
  n <- c(26, 29, 31)
  
  # LEAST-SQUARES SOLUTION:
  mat[1:3,1] <- 1
  mat[1:3,2] <- B*sin(n*k+phi)
  
  mat2 <- solve(t(mat)%*%mat)%*%t(mat)
  
  c0A <-as.matrix(dat[,1:3]) %*% t(mat2)
  
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
  negA_rowidx = c0A[,2] < 0
  c0A[negA_rowidx,1] <- apply(dat[negA_rowidx, 1:3], 1, mean)
  c0A[negA_rowidx,2] <- 0
  return(c0A)
}

# Given a phi_k, find C0, m, d for each sequence:
find_c0A1A2 <- function(dat, phi){
  # Find C0, and 2 Amplitude parameters for each sequence
  # (Cn = C0 + (m*phi+d)*sin(b*n+phi))
  # dat: Each row corresponds to a sequence.
  #   First three columns need to be C26, C29, and C31, in that order
  #
  # Returns: columns C0, m, d
  mat <- matrix(0, nrow = 3, ncol = 3)
  k <- 2*pi/10.4
  mat[1:3,1] <- 1
  mat[1:3,2] <- sin(n*k + phi)
  mat[1:3,3] <- sin(n*k + phi)^2
  # mat[1:3,3] <- sin(phi)
  # mat[1,2:3] <- mat[1,2:3]*aa[1]
  # mat[2,2:3] <- mat[2,2:3]*aa[2]
  # mat[3,2:3] <- mat[3,2:3]*aa[3]
  inv_mat <- solve(mat)
  c0A1A2 <- as.matrix(dat%>%select(C26,C29,C31)) %*% t(inv_mat)
  return(c0A1A2)
}

min_loss_phi_from_B = function(data, phi_values, B_values) {
  k=2*pi/10.4
  c0A_list = map(phi_values, ~find_c0A(data %>% select(C26, C29, C31), .x, B=B_values))
  
  loss_by_component_list = map2(c0A_list, phi_values, ~c(sum((data %>% select(C26) - (.x[,1] + B_values[1]*.x[,2]*sin(k*26+.y)))^2),
                                                         sum((data %>% select(C29) - (.x[,1] + B_values[2]*.x[,2]*sin(k*29+.y)))^2),
                                                         sum((data %>% select(C31) - (.x[,1] + B_values[3]*.x[,2]*sin(k*31+.y)))^2)))
  
  total_loss_list = map(loss_by_component_list, ~mean(.x)/dim(data)[1]) %>% unlist()
  
  return(min(total_loss_list))
}

# Try brute-forcing our way to optimal phi_k (for a given group):
# find_optimal_phi_AT = function(data, phi_values=seq(from=-pi,to=pi,by=pi/40), B_values_grid) {
#   
#   min_total_loss_list = apply(B_values_grid, 1, function(row) {min_loss_phi_from_B(data,phi_values,row)})
#   
#   optimal_B = B_values_grid[which(min_total_loss_list == min(min_total_loss_list)),]
#   
#   c0A_list = map(phi_values, ~find_c0A(data %>% select(C26, C29, C31), .x, B=optimal_B))
#   
#   k = 2*pi/10.4
#   loss_by_component_list = map2(c0A_list, phi_values, ~c(sum((data %>% select(C26) - (.x[,1] + optimal_B[1]*.x[,2]*sin(k*26+.y)))^2),
#                                                          sum((data %>% select(C29) - (.x[,1] + optimal_B[2]*.x[,2]*sin(k*29+.y)))^2),
#                                                          sum((data %>% select(C31) - (.x[,1] + optimal_B[3]*.x[,2]*sin(k*31+.y)))^2)))
#   
#   total_loss_list = map(loss_by_component_list, ~mean(.x)/dim(data)[1]) %>% unlist()
#   
#   to_return = list()
#   
#   optimal_index = which(total_loss_list == min(total_loss_list))
#   
#   if(length(optimal_index) > 1) {
#     optimal_index=optimal_index[1]
#   }
#   
#   to_return$optimal_phi = phi_values[optimal_index]
#   to_return$optimal_B = optimal_B
#   to_return$c0A = c0A_list[[optimal_index]]
#   to_return$phi_values = phi_values
#   to_return$total_loss_list = total_loss_list
#   to_return$loss_by_component_list = loss_by_component_list
#   
#   return(to_return)
# }

find_optimal_phi_AT_dependentA = function(data, phi_values=seq(from=-pi,to=pi,by=pi/40)) {
  
  c0A1A2_list = map(phi_values, ~find_c0A1A2(data %>% select(C26, C29, C31), .x))
  # c0A0A1_list = map(phi_values, ~find_c0A0A1(data %>% select(C26, C29, C31), .x))
  
  k = 2*pi/10.4
  loss_by_component_list = map2(c0A1A2_list, phi_values, ~c(sum((data %>% select(C26) - (.x[,1] + .x[,2]*sin(k*26+.y) + .x[,3]*sin(k*26+.y)^2))^2),
                                                            sum((data %>% select(C29) - (.x[,1] + .x[,2]*sin(k*29+.y) + .x[,3]*sin(k*29+.y)^2))^2),
                                                            sum((data %>% select(C31) - (.x[,1] + .x[,2]*sin(k*31+.y) + .x[,3]*sin(k*31+.y)^2))^2)))
  
  total_loss_list = map(loss_by_component_list, ~mean(.x)/dim(data)[1]) %>% unlist()
  
  to_return = list()
  
  optimal_index = which(total_loss_list == min(total_loss_list))
  
  if(length(optimal_index) > 1) {
    optimal_index=optimal_index[1]
  }
  
  to_return$optimal_phi = phi_values[optimal_index]
  to_return$c0A1A2 = c0A1A2_list[[optimal_index]]
  to_return$phi_values = phi_values
  to_return$total_loss_list = total_loss_list
  to_return$loss_by_component_list = loss_by_component_list
  
  return(to_return)
}

optimize_phi_list_random = lapply(grouped_df_list_random, find_optimal_phi_AT_dependentA)

# optimize_phi_list_random = lapply(grouped_df_list_random, find_optimal_phi_AT)
# optimize_phi_list_tiling = lapply(grouped_df_list_tiling, find_optimal_phi_AT)
# optimize_phi_list_chrv = lapply(grouped_df_list_chrv, find_optimal_phi_AT)

# optimize_phi_list_random = lapply(map(1:(num_groups+1), ~random%>%filter(group_assignment_index==.x)), find_optimal_phi_AT)
# optimize_phi_list_tiling = lapply(map(1:(num_groups+1), ~tiling%>%filter(group_assignment_index==.x)), find_optimal_phi_AT)
# optimize_phi_list_chrv = lapply(map(1:(num_groups+1), ~chrv%>%filter(group_assignment_index==.x)), find_optimal_phi_AT)

B_values_grid_size = 7
B_values_grid = matrix(nrow=B_values_grid_size^2, ncol=3)
B_values_grid[,1] = 1
B_values_grid[,2] = rep(seq(from=1-(B_values_grid_size-1)/2*.1, 
                            to=1+(B_values_grid_size-1)/2*.1, 
                            by=0.1), 
                        each=B_values_grid_size)
B_values_grid[,3] = rep(seq(from=1-(B_values_grid_size-1)/2*.1, 
                            to=1+(B_values_grid_size-1)/2*.1, 
                            by=0.1), 
                        length.out=B_values_grid_size^2)

optimize_phi_list_random = map(map(1:(num_groups+1), ~random%>%filter(group_assignment_index==.x)), ~find_optimal_phi_AT(data=.x, B_values_grid = B_values_grid))
optimize_phi_list_tiling = map(map(1:(num_groups+1), ~tiling%>%filter(group_assignment_index==.x)), ~find_optimal_phi_AT(data=.x, B_values_grid = B_values_grid))
optimize_phi_list_chrv = map(map(1:(num_groups+1), ~chrv%>%filter(group_assignment_index==.x)), ~find_optimal_phi_AT(data=.x, B_values_grid = B_values_grid))

for(i in 1:10) {
  print(c(i,abs(pi-(optimize_phi_list_random[[i]]$optimal_phi - optimize_phi_list_random[[i+10]]$optimal_phi))))
}
for(i in 1:10) {
  print(c(i,abs(pi-(optimize_phi_list_tiling[[i]]$optimal_phi - optimize_phi_list_tiling[[i+10]]$optimal_phi))))
}


breaks=seq(from=-4.5, to=4.5, by=0.25)

# Group 7: C26 should be analogous to Group 17: C31, also,
# Group 7: C31 should be analogous to Group 17: C26
hist((random%>%filter(group_assignment_index==7)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==17)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==7)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==17)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
# 7 -> 17 shift left

# Group 8: C26 should be analogous to Group 18: C31, also,
# Group 8: C31 should be analogous to Group 18: C26
hist((random%>%filter(group_assignment_index==8)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==18)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==8)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((random%>%filter(group_assignment_index==18)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
# 8 -> 18 shift left


# Group 6: C26 should be analogous to Group 16: C31, also,
# Group 6: C31 should be analogous to Group 16: C26
hist((tiling%>%filter(group_assignment_index==6)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==16)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==6)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==16)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)

c(mean((tiling%>%filter(group_assignment_index==6)%>%select(C26))[,1]), var((tiling%>%filter(group_assignment_index==6)%>%select(C26))[,1]))
c(mean((tiling%>%filter(group_assignment_index==16)%>%select(C31))[,1]), var((tiling%>%filter(group_assignment_index==16)%>%select(C31))[,1]))
c(mean((tiling%>%filter(group_assignment_index==6)%>%select(C31))[,1]), var((tiling%>%filter(group_assignment_index==6)%>%select(C31))[,1]))
c(mean((tiling%>%filter(group_assignment_index==16)%>%select(C26))[,1]), var((tiling%>%filter(group_assignment_index==16)%>%select(C26))[,1]))

c(mean((tiling%>%filter(group_assignment_index==6)%>%select(C29))[,1]), var((tiling%>%filter(group_assignment_index==6)%>%select(C29))[,1]))
c(mean((tiling%>%filter(group_assignment_index==16)%>%select(C29))[,1]), var((tiling%>%filter(group_assignment_index==16)%>%select(C29))[,1]))

# 6 -> 16 shift left

# Group 7: C26 should be analogous to Group 17: C31, also,
# Group 7: C31 should be analogous to Group 17: C26
hist((tiling%>%filter(group_assignment_index==7)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==17)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==7)%>%select(C31))[,1],xlim=c(-3,3), breaks=breaks)
hist((tiling%>%filter(group_assignment_index==17)%>%select(C26))[,1],xlim=c(-3,3), breaks=breaks)

c(mean((tiling%>%filter(group_assignment_index==7)%>%select(C26))[,1]), var((tiling%>%filter(group_assignment_index==7)%>%select(C26))[,1]))
c(mean((tiling%>%filter(group_assignment_index==17)%>%select(C31))[,1]), var((tiling%>%filter(group_assignment_index==17)%>%select(C31))[,1]))
c(mean((tiling%>%filter(group_assignment_index==7)%>%select(C31))[,1]), var((tiling%>%filter(group_assignment_index==7)%>%select(C31))[,1]))
c(mean((tiling%>%filter(group_assignment_index==17)%>%select(C26))[,1]), var((tiling%>%filter(group_assignment_index==17)%>%select(C26))[,1]))
# 7 -> 17 shift left




# Optimal phi_k for each group:
plot(map(optimize_phi_list_random, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, random",
     ylim=c(-pi,pi))
plot(map(optimize_phi_list_tiling, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, tiling",
     ylim=c(-pi,pi))
plot(map(optimize_phi_list_chrv, ~.x$optimal_phi) %>% unlist(),
     ylab="optimal phi by group, chrv",
     ylim=c(-pi,pi))

# Optimal B[2] for each group:
plot(map(optimize_phi_list_random, ~.x$optimal_B[2]) %>% unlist(),
     ylab="optimal B[2] by group, random",
     ylim=c(.7,1.3))
plot(map(optimize_phi_list_tiling, ~.x$optimal_B[2]) %>% unlist(),
     ylab="optimal B[2] by group, tiling",
     ylim=c(.7,1.3))
plot(map(optimize_phi_list_chrv, ~.x$optimal_B[2]) %>% unlist(),
     ylab="optimal B[2] by group, chrv",
     ylim=c(.7,1.3))

# Optimal B[3] for each group:
plot(map(optimize_phi_list_random, ~.x$optimal_B[3]) %>% unlist(),
     ylab="optimal B[3] by group, random",
     ylim=c(.7,1.3))
plot(map(optimize_phi_list_tiling, ~.x$optimal_B[3]) %>% unlist(),
     ylab="optimal B[3] by group, tiling",
     ylim=c(.7,1.3))
plot(map(optimize_phi_list_chrv, ~.x$optimal_B[3]) %>% unlist(),
     ylab="optimal B[3] by group, chrv",
     ylim=c(.7,1.3))


ylims=c(0,0.15)
# ylims=c(0,0.5)
# Plot loss vs phi for each group
map(1:length(optimize_phi_list_random), ~plot(optimize_phi_list_random[[.x]]$phi_values, 
                                              optimize_phi_list_random[[.x]]$total_loss_list,
                                              ylab=paste(c("Loss for group ", .x)),
                                              xlab="Phi Values, random",
                                              ylim=ylims))
map(1:length(optimize_phi_list_tiling), ~plot(optimize_phi_list_tiling[[.x]]$phi_values, 
                                              optimize_phi_list_tiling[[.x]]$total_loss_list,
                                              ylab=paste(c("Loss for group ", .x)),
                                              xlab="Phi Values, tiling",
                                              ylim=ylims))
map(1:length(optimize_phi_list_chrv), ~plot(optimize_phi_list_chrv[[.x]]$phi_values, 
                                            optimize_phi_list_chrv[[.x]]$total_loss_list,
                                            ylab=paste(c("Loss for group ", .x)),
                                            xlab="Phi Values, chrv",
                                            ylim=ylims))


# Plot loss by component vs phi for each group:
plot_loss_by_component = function(optimize_phi_obj, group_number, library_name, ylims) {
  plot(optimize_phi_obj$phi_values,
       matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[1,]/length(optimize_phi_obj$phi_values),
       col="blue",
       ylab=paste(c("Loss for group ", group_number)),
       xlab=paste(c("Phi Values"), library_name),
       ylim=ylims)
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[2,]/length(optimize_phi_obj$phi_values),
         col="red")
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[3,]/length(optimize_phi_obj$phi_values),
         col="green")
  legend("topright", c("26", "29", "31"), fill=c("blue", "red", "green"))
}
map(1:length(optimize_phi_list_random), ~plot_loss_by_component(optimize_phi_list_random[[.x]], 
                                                                .x, "random", c(0,3)))
map(1:length(optimize_phi_list_tiling), ~plot_loss_by_component(optimize_phi_list_tiling[[.x]], 
                                                                .x, "tiling", c(0,30)))
map(1:length(optimize_phi_list_chrv), ~plot_loss_by_component(optimize_phi_list_chrv[[.x]], 
                                                              .x, "chrv", c(0, 40)))



plot(map(optimize_phi_list_random, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("random new C0")

plot(map(optimize_phi_list_tiling, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("tiling new C0")

plot(map(optimize_phi_list_chrv, ~mean(.x$c0A[,1]))%>%unlist(),
     ylab = "Mean New C0 Value by group",
     ylim = c(-0.5, 0.25))
title("chrv new C0")



for(i in 1:length(grouped_df_list_random)) {
  grouped_df_list_random[[i]]$C0_grouped_phi = optimize_phi_list_random[[i]]$c0A[,1]
  grouped_df_list_random[[i]]$A_grouped_phi = optimize_phi_list_random[[i]]$c0A[,2]
}
grouped_phi_df_random = bind_rows(grouped_df_list_random)

for(i in 1:length(grouped_df_list_tiling)) {
  grouped_df_list_tiling[[i]]$C0_grouped_phi = optimize_phi_list_tiling[[i]]$c0A[,1]
  grouped_df_list_tiling[[i]]$A_grouped_phi = optimize_phi_list_tiling[[i]]$c0A[,2]
}
grouped_phi_df_tiling = bind_rows(grouped_df_list_tiling)

for(i in 1:length(grouped_df_list_chrv)) {
  grouped_df_list_chrv[[i]]$C0_grouped_phi = optimize_phi_list_chrv[[i]]$c0A[,1]
  grouped_df_list_chrv[[i]]$A_grouped_phi = optimize_phi_list_chrv[[i]]$c0A[,2]
}
grouped_phi_df_chrv = bind_rows(grouped_df_list_chrv)


random_group1 = as.data.frame(random %>% filter(group_assignment_index==1))
random_group6 = as.data.frame(random %>% filter(group_assignment_index==6))
random_group11 = as.data.frame(random %>% filter(group_assignment_index==11))

construct_periodicity_plots_AT(random_group1, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group1, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group1, C31, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group1, C0, ylims_quartiles = c(0.3, 0.9))[[1]]

construct_periodicity_plots_AT(random_group6, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group6, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group6, C31, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group6, C0, ylims_quartiles = c(0.3, 0.9))[[1]]

construct_periodicity_plots_AT(random_group11, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group11, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(random_group11, C31, ylims_quartiles = c(0.3, 0.9))[[1]]



tiling_group1 = as.data.frame(tiling %>% filter(group_assignment_index==1))
tiling_group6 = as.data.frame(tiling %>% filter(group_assignment_index==6))
tiling_group11 = as.data.frame(tiling %>% filter(group_assignment_index==11))

construct_periodicity_plots_AT(tiling_group1, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group1, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group1, C31, ylims_quartiles = c(0.3, 0.9))[[1]]

construct_periodicity_plots_AT(tiling_group6, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group6, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group6, C31, ylims_quartiles = c(0.3, 0.9))[[1]]

construct_periodicity_plots_AT(tiling_group11, C26, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group11, C29, ylims_quartiles = c(0.3, 0.9))[[1]]
construct_periodicity_plots_AT(tiling_group11, C31, ylims_quartiles = c(0.3, 0.9))[[1]]





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
  to_return$loss_by_component_list = loss_by_component_list
  
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
       matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[1,]/length(optimize_phi_obj$phi_values),
       col="blue",
       ylab=paste(c("Loss for group ", group_number)),
       xlab=paste(c("Phi Values"), library_name),
       ylim=ylims)
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[2,]/length(optimize_phi_obj$phi_values),
         col="red")
  points(optimize_phi_obj$phi_values,
         matrix(optimize_phi_obj$loss_by_component_list %>% unlist(), nrow=3)[3,]/length(optimize_phi_obj$phi_values),
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
