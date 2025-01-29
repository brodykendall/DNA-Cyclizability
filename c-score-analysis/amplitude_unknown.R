library(gridExtra)

cycle3 = read_csv("cycle3.txt")
cycle3 = cycle3 %>%
  rename(C26 = "n=26", C29 = "n=29", C31 = "n=31")
# cycle3_seq1 = cycle3[1,]
ylims=c(-0.75, 0.75)

# Max C0:
cycle3_seq1 = cycle3[572,]
ylims=c(0, 5.5)


library(reshape2)
cycle3_long = melt(cycle3, id.vars=c("Sequence", "C0", "Amplitude", "Phase"))

#then plot
ggplot(cycle3_long, aes(x=value))+
  geom_histogram() + 
  facet_wrap(~variable)

ggplot(data=cycle3)+
  geom_boxplot(aes(y=C26))+
  geom_boxplot(aes(y=C29))+
  geom_boxplot(aes(y=C31))

point_size = 5
C0_linewidth = 1

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

dat_random$actual_C0 = find_c0Aphi(dat_random%>%select(C26,C29,C31), aa=c(1,0.82,0.7))

find_cn1 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  
  return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

cycle3_seq1_c0Aphi1 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1,1))
cycle3_seq1_cn1 = find_cn1(cycle3_seq1_c0Aphi1)
cycle3_seq1_df1 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn1)

constant_amplitude = ggplot() +
  geom_line(data=cycle3_seq1_df1, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi1[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)


find_cn2 <- function(c0Aphi, modulo, a=0.7) {
  k <- 2*pi/10.4
  n_list = 1:50
  # n_list = seq(1,50,by=0.1)
  a_ratio_list = 1 + sign(n_list-26)*(abs(n_list - 26)%%modulo)*(a-1)/(5)
  # a_ratio_list = 1 + (n_list - 26)*(a-1)/(5)
  
  print(a_ratio_list)
  
  a_list = c0Aphi[,2]/a_ratio_list
  
  a_list[a_ratio_list<0.05] = NA
  
  
  
  print(a_list)
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
  # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

cycle3_seq1_c0Aphi2 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/.82,1/.7))
cycle3_seq1_cn2 = find_cn2(cycle3_seq1_c0Aphi2, 5.2)
print(cycle3_seq1_cn2)
print(cycle3_seq1)
cycle3_seq1_df2 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn2)

ggplot() +
  geom_line(data=cycle3_seq1_df2, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi2[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)


cycle3_seq1_c0Aphi2_1 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/.844,1/.74))
cycle3_seq1_cn2_1 = find_cn2(cycle3_seq1_c0Aphi2_1, 5.2, a=0.74)
print(cycle3_seq1_cn2_1)
print(cycle3_seq1)
cycle3_seq1_df2_1 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn2_1)

ggplot() +
  geom_line(data=cycle3_seq1_df2_1, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi2_1[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)


cycle3_seq1_cn2.1 = find_cn2(cycle3_seq1_c0Aphi2, 10.4)
print(cycle3_seq1_cn2.1)
print(cycle3_seq1)
cycle3_seq1_df2.1 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn2.1)

ggplot() +
  geom_line(data=cycle3_seq1_df2.1, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi2[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)


cycle3_seq1_cn2.2 = find_cn2(cycle3_seq1_c0Aphi2, 50)
print(cycle3_seq1_cn2.2)
print(cycle3_seq1)
cycle3_seq1_df2.2 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn2.2)
# cycle3_seq1_df2.2 = data.frame("n" = seq(1,50,by=0.1), "Cn" = cycle3_seq1_cn2.2)

interpolated = ggplot() +
  geom_line(data=cycle3_seq1_df2.2, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi2[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)




# An/A26 follows Cauchy Distribution

find_cn3 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  a_ratio_list = 1/(1+((n_list-26)/(5*sqrt(7/3)))^2)
  
  print(a_ratio_list)
  
  a_list = c0Aphi[,2]/a_ratio_list
  
  # a_list[26] = c0Aphi[,2]
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
  # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

# b=0.1
cycle3_seq1_c0Aphi3 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/0.866,1/.7))
cycle3_seq1_cn3 = find_cn3(cycle3_seq1_c0Aphi3)
print(cycle3_seq1_cn3)
print(cycle3_seq1)
cycle3_seq1_df3 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn3)

cauchy = ggplot() +
  geom_line(data=cycle3_seq1_df3, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi3[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)

cauchy





find_cn4 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  # a_ratio_list = exp(-(n_list-26)/14.018)
  a_ratio_list = exp(-(n_list-26)/19.1304)
  
  print(a_ratio_list)
  
  a_list = c0Aphi[,2]/a_ratio_list
  # a_list = rev(c0Aphi[,2]/a_ratio_list)
  
  # a_list[26] = c0Aphi[,2]
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
  # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

# cycle3_seq1_c0Aphi4 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/.807,1/.7))
cycle3_seq1_c0Aphi4 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/0.8548592,1/.77))
cycle3_seq1_cn4 = find_cn4(cycle3_seq1_c0Aphi4)
print(cycle3_seq1_cn4)
print(cycle3_seq1)
cycle3_seq1_df4 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn4)

exp_decay = ggplot() +
  geom_line(data=cycle3_seq1_df4, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi4[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)

exp_decay



# NEED TO FIX THIS, BUT IT's CLOSE
# find_cn5 <- function(c0Aphi) {
#   k <- 2*pi/10.4
#   n_list = 1:50
#   a_ratio_list = exp(-(n_list-26)/-14.018)
#   
#   print(a_ratio_list)
#   
#   # a_list = c0Aphi[,2]/a_ratio_list
#   a_list = rev(c0Aphi[,2]*a_ratio_list)
#   
#   # a_list[26] = c0Aphi[,2]
#   
#   return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
#   # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
# }

# # b=0.1
# cycle3_seq1_c0Aphi5 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/.807,1/.7))
# cycle3_seq1_cn5 = find_cn5(cycle3_seq1_c0Aphi5)
# print(cycle3_seq1_cn5)
# print(cycle3_seq1)
# cycle3_seq1_df5 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn5)
# 
# ggplot() +
#   geom_line(data=cycle3_seq1_df5, mapping=aes(x=n, y=Cn)) + 
#   geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
#   geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
#   geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
#   geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi5[,1], color="C0"), linewidth=C0_linewidth) +
#   ylim(ylims)


find_cn6 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  # a_ratio_list = 1-0.3005483*sin(2*pi*(n_list-26)/20.8)
  a_ratio_list = 1-0.2304204*sin(2*pi*(n_list-26)/20.8)
  
  # print(a_ratio_list)
  
  a_list = c0Aphi[,2]/a_ratio_list
  
  # a_list[26] = c0Aphi[,2]
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
  # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

# b=0.1
# cycle3_seq1_c0Aphi6 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/0.7634133,1/.7))
cycle3_seq1_c0Aphi6 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/0.8186329,1/.77))
cycle3_seq1_cn6 = find_cn6(cycle3_seq1_c0Aphi6)
print(cycle3_seq1_cn6)
print(cycle3_seq1)
cycle3_seq1_df6 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn6)

ggplot() +
  geom_line(data=cycle3_seq1_df6, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi6[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)



find_cn7 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  a_ratio_list = 1-9.617547/n_list*sin(2*pi*(n_list-26)/20.8)
  
  # print(a_ratio_list)
  
  a_list = c0Aphi[,2]/a_ratio_list
  
  # a_list[26] = c0Aphi[,2]
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
  # return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
}

# b=0.1
cycle3_seq1_c0Aphi7 = find_c0Aphi((cycle3_seq1 %>% select(C26, C29, C31)), aa=c(1,1/0.7195878,1/.7))
cycle3_seq1_cn7 = find_cn7(cycle3_seq1_c0Aphi7)
print(cycle3_seq1_cn7)
print(cycle3_seq1)
cycle3_seq1_df7 = data.frame("n" = 1:50, "Cn" = cycle3_seq1_cn7)

ggplot() +
  geom_line(data=cycle3_seq1_df7, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=cycle3_seq1, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=cycle3_seq1, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  # geom_hline(mapping=aes(yintercept=cycle3_seq1_c0Aphi6[,1], color="C0"), linewidth=C0_linewidth) +
  ylim(ylims)




