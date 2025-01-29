library(tidyverse)
library(gridExtra)
library(reshape2)

source("scripts/functions/C0-functions.R")

dat_random = read.csv("cycle3.txt")
dat_tiling = read.csv("cycle5.txt")

colnames(dat_random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

dat_random_1seq = dat_random[572,]
ylims=c(1, 3.75)

dat_random_1seq = dat_random[6569,]
ylims=c(-3,3.5)
 
# dat_random_1seq = dat_random[10533,]
# ylims=c(-3,3)


dat_random_1seq = dat_tiling[45054,]
# ylims=c(-5,5)
ylims=c(-5,6)


point_size = 3
C0_linewidth = 1

# find_cn1 <- function(c0Aphi) {
#   k <- 2*pi/10.4
#   n_list = 1:50
#   
#   return(c0Aphi[,1] + c0Aphi[,2]*sin(k*n_list + c0Aphi[,3]))
# }
# 
# dat_random_1seq_c0Aphi1 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1,1))
# dat_random_1seq_cn1 = find_cn1(dat_random_1seq_c0Aphi1)
# dat_random_1seq_df1 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn1)
# 
# constant_amplitude = ggplot() +
#   geom_line(data=dat_random_1seq_df1, mapping=aes(x=n, y=Cn)) + 
#   geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
#   # geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi1[,1], color="C0"), linewidth=C0_linewidth) +
#   ylim(ylims) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# 
# constant_amplitude

find_cn <- function(c0Aphi, modulo, a=0.7) {
  k <- 2*pi/10.4
  n_list = 1:50
  a_ratio_list = 1 + sign(n_list-26)*(abs(n_list - 26)%%modulo)*(a-1)/(5)

  a_list = c0Aphi[,2]/a_ratio_list
  
  a_list[a_ratio_list<0.05] = NA
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
}

dat_random_1seq_c0Aphi = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.82,1/.7))
dat_random_1seq_cn = find_cn(dat_random_1seq_c0Aphi, 50)
dat_random_1seq_df = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn)

interpolated = ggplot() +
  geom_line(data=dat_random_1seq_df, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

interpolated

dat_random_1seq_c0Aphi[,1]


dat_random_1seq_c0Aphi1.2 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.94,1/.9))
dat_random_1seq_cn1.2 = find_cn(dat_random_1seq_c0Aphi1.2, 50, a=0.9)
dat_random_1seq_df1.2 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn1.2)

interpolated1.2 = ggplot() +
  geom_line(data=dat_random_1seq_df1.2, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi1.2[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

interpolated1.2

dat_random_1seq_c0Aphi1.2[,1]

# 
# 
# find_cn2 <- function(c0Aphi, modulo, a=0.7) {
#   k <- 2*pi/10.4
#   n_list = 1:50
#   a_ratio_list = 1 + sign(n_list-26)*(abs(n_list - 26)%%modulo)*(a-1)/(5)
#   
#   a_list = c0Aphi[,2]/a_ratio_list
#   
#   a_list[a_ratio_list<0.05] = NA
#   
#   return(c0Aphi[,1] + c0Aphi[,1]*a_list*sin(k*n_list + c0Aphi[,3]))
# }
# 
# 
# 
# 
# 
# dat_random_1seq_c0Aphi2 = find_c0Aphi2((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.82,1/.7))
# dat_random_1seq_cn2 = find_cn2(dat_random_1seq_c0Aphi2, 50)
# dat_random_1seq_df2 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn2)
# 
# interpolated2 = ggplot() +
#   geom_line(data=dat_random_1seq_df2, mapping=aes(x=n, y=Cn)) + 
#   geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
#   geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi2[,1]), color="red", linewidth=C0_linewidth) +
#   ylim(ylims) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# 
# interpolated2
# 
# dat_random_1seq_c0Aphi[,1]
# 
# 
# 
# dat_random_1seq_c0Aphi2.2 = find_c0Aphi2((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.82,1/.9))
# dat_random_1seq_cn2.2 = find_cn2(dat_random_1seq_c0Aphi2.2, 50)
# dat_random_1seq_df2.2 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn2.2)
# 
# interpolated2.2 = ggplot() +
#   geom_line(data=dat_random_1seq_df2.2, mapping=aes(x=n, y=Cn)) + 
#   geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
#   geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
#   geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi2.2[,1]), color="red", linewidth=C0_linewidth) +
#   ylim(ylims) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# 
# interpolated2.2


# A26/An = 1/(1+|n-26|/gamma)
# Or more simply, An=(1+|n-26|/gamma)*A26

find_cn3 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  # a_ratio_list = 1/(1+((n_list-26)/(5*sqrt(7/3)))^2)
  # a_ratio_list = 1/(1+abs((n_list-26)/16.6667))
  a_ratio_list = 1/(1+abs((n_list-26)/11.667))
  
  a_list = c0Aphi[,2]/a_ratio_list

  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
}

# dat_random_1seq_c0Aphi3 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.866,1/.7))
dat_random_1seq_c0Aphi3 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.8475,1/.7))
dat_random_1seq_cn3 = find_cn3(dat_random_1seq_c0Aphi3)
dat_random_1seq_df3 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn3)

cauchy = ggplot() +
  geom_line(data=dat_random_1seq_df3, mapping=aes(x=n, y=Cn)) +
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi3[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

cauchy

dat_random_1seq_c0Aphi3[,1]




# A31/An = 1/(1-|n-31|/gamma)
# Or more simply, An = (1-|n-31|/gamma)*A31

find_cn4 <- function(c0Aphi) {
  k <- 2*pi/10.4
  n_list = 1:50
  # a_ratio_list = 1/(1+((n_list-26)/(5*sqrt(7/3)))^2)
  # a_ratio_list = 1/(1-abs((n_list-31)/16.667))
  a_ratio_list = 1/(1-sqrt(abs(n_list-31))/7.454)
  
  a31 = c0Aphi[,2]/0.7
  
  a_list = a31/a_ratio_list
  
  # a_list[1:16] = 0
  
  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
}

# dat_random_1seq_c0Aphi3 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.866,1/.7))
# dat_random_1seq_c0Aphi4 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.7955,1/.7))
dat_random_1seq_c0Aphi4 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.8639,1/.7))
# dat_random_1seq_c0Aphi4 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/0.866,1/.7))
dat_random_1seq_cn4 = find_cn4(dat_random_1seq_c0Aphi4)
dat_random_1seq_df4 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn4)

cauchy2 = ggplot() +
  geom_line(data=dat_random_1seq_df4, mapping=aes(x=n, y=Cn)) +
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi4[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

cauchy2

dat_random_1seq_c0Aphi3[,1]



find_cn4 <- function(c0Aphi, C=14.0184) {
  k <- 2*pi/10.4
  n_list = 1:50
  a_ratio_list = exp(-(n_list-26)/C)
  
  a_list = c0Aphi[,2]/a_ratio_list

  return(c0Aphi[,1] + a_list*sin(k*n_list + c0Aphi[,3]))
}

dat_random_1seq_c0Aphi4 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.807,1/.7))
dat_random_1seq_cn4 = find_cn4(dat_random_1seq_c0Aphi4)
dat_random_1seq_df4 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn4)

exp_decay = ggplot() +
  geom_line(data=dat_random_1seq_df4, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi4[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

exp_decay

dat_random_1seq_c0Aphi4[,1]


dat_random_1seq_c0Aphi4.2 = find_c0Aphi((dat_random_1seq %>% select(C26, C29, C31)), aa=c(1,1/.936,1/.9))
dat_random_1seq_cn4.2 = find_cn4(dat_random_1seq_c0Aphi4.2, C=47.46)
dat_random_1seq_df4.2 = data.frame("n" = 1:50, "Cn" = dat_random_1seq_cn4.2)

exp_decay2 = ggplot() +
  geom_line(data=dat_random_1seq_df4.2, mapping=aes(x=n, y=Cn)) + 
  geom_point(data=dat_random_1seq, mapping=aes(x=26, y=C26), color="blue", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=29, y=C29), color="orange", size=point_size) +
  geom_point(data=dat_random_1seq, mapping=aes(x=31, y=C31), color="green", size=point_size) +
  geom_hline(mapping=aes(yintercept=dat_random_1seq_c0Aphi4.2[,1]), color="red", linewidth=C0_linewidth) +
  ylim(ylims) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

exp_decay2

dat_random_1seq_c0Aphi4.2[,1]


names <- c('C26', 'C29', 'C31', 'C0')
clrs <- c('blue', 'orange', 'green', 'red')
ltype <- c(3, 4, 5, 1)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend = names, bty='n', fill=clrs, ncol=4, x.intersp=0.2)



# grid.arrange(constant_amplitude, interpolated, cauchy, exp_decay)
grid.arrange(interpolated, interpolated2, interpolated1.2, interpolated2.2)

