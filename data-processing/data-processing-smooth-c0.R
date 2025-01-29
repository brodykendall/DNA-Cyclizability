library(tidyverse)

find_c0Aphi <- function(dat, n=c(26,29,31), aa=c(1,1/0.82,1/0.7)){
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

find_c0sAphi_lag <- function(dat, aa=c(1,1/0.82,1/0.7), lag_weight=1){
  # Find C0s, Amplitude, and Phase for each sequence
  # dat: Each row corresponds to a sequence.
  # 
  # Returns: columns C0 (of interest), C0 estimate (previous sequence), 
  #   C0 estimate (next sequence), A, phi
  
  # mat <- matrix(0, nrow = 9, ncol = 5)
  mat <- matrix(0, nrow = 9, ncol = 9)
  k <- 2*pi/10.4
  n <- c(26, 29, 31, 19, 22, 24, 33, 36, 38)
  
  # LEAST-SQUARES SOLUTION:
  mat[1:3,1] <- 1
  mat[4:6,2] <- 1
  mat[7:9,3] <- 1
  # mat[1:9,4] <- sin(n*k)
  # mat[1:9,5] <- cos(n*k)
  mat[1:3,4] <- sin(n[1:3]*k)
  mat[1:3,5] <- cos(n[1:3]*k)  
  mat[4:6,6] <- sin(n[4:6]*k)
  mat[4:6,7] <- cos(n[4:6]*k)
  mat[7:9,8] <- sin(n[7:9]*k)
  mat[7:9,9] <- cos(n[7:9]*k)
  # mat[1,4:5] <- mat[1,4:5]*aa[1]
  # mat[2,4:5] <- mat[2,4:5]*aa[2]
  # mat[3,4:5] <- mat[3,4:5]*aa[3]
  # mat[4,4:5] <- mat[4,4:5]*aa[1]
  # mat[5,4:5] <- mat[5,4:5]*aa[2]
  # mat[6,4:5] <- mat[6,4:5]*aa[3]
  # mat[7,4:5] <- mat[7,4:5]*aa[1]
  # mat[8,4:5] <- mat[8,4:5]*aa[2]
  # mat[9,4:5] <- mat[9,4:5]*aa[3]
  mat[1,4:9] <- mat[1,4:9]*aa[1]
  mat[2,4:9] <- mat[2,4:9]*aa[2]
  mat[3,4:9] <- mat[3,4:9]*aa[3]
  mat[4,4:9] <- mat[4,4:9]*aa[1]
  mat[5,4:9] <- mat[5,4:9]*aa[2]
  mat[6,4:9] <- mat[6,4:9]*aa[3]
  mat[7,4:9] <- mat[7,4:9]*aa[1]
  mat[8,4:9] <- mat[8,4:9]*aa[2]
  mat[9,4:9] <- mat[9,4:9]*aa[3]
  
  # mat2 <- solve(t(mat)%*%mat)%*%t(mat)
  weights <- diag(1, nrow = 9, ncol = 9)
  weights[4:9,4:9] = diag(lag_weight, 6, 6)
  mat2 <- solve(t(mat)%*%weights%*%mat)%*%t(mat)%*%weights
  
  c0sA1A2 <-as.matrix(dat %>% select(C26, C29, C31, C19, C22, C24, C33, C36, C38)) %*% t(mat2)
  
  c0sAphi <- c0sA1A2
  c0sAphi[,1:3] <- c0sA1A2[,1:3]
  # c0sAphi[,4] <- sqrt(c0sA1A2[,4]^2 + c0sA1A2[,5]^2)
  # c0sAphi[,5] <- sign(c0sA1A2[,5]) * acos(c0sA1A2[,4]/c0sAphi[,4]) - pi
  # Amplitude for main sequence:
  c0sAphi[,4] <- sqrt(c0sA1A2[,4]^2 + c0sA1A2[,5]^2)
  # Amplitude for previous sequence:
  c0sAphi[,5] <- sqrt(c0sA1A2[,6]^2 + c0sA1A2[,7]^2)
  # Amplitude for next sequence:
  c0sAphi[,6] <- sqrt(c0sA1A2[,8]^2 + c0sA1A2[,9]^2)
  # Phase
  c0sAphi[,7] <- atan2(c0sA1A2[,5] + c0sA1A2[,7] + c0sA1A2[,9], 
                       c0sA1A2[,4] + c0sA1A2[,6] + c0sA1A2[,8])
  c0sAphi = c0sAphi[,1:7]

  # If A < 0, then A <- 0 and C0 <- mean(Cn)
  # negA_rowidx = c0sAphi[,4] < 0
  # c0sAphi[negA_rowidx,1] <- apply(dat[negA_rowidx,] %>% select(C26, C29, C31), 1, mean)
  # c0sAphi[negA_rowidx,2] <- apply(dat[negA_rowidx,] %>% select(C19, C22, C24), 1, mean)
  # c0sAphi[negA_rowidx,3] <- apply(dat[negA_rowidx,] %>% select(C33, C36, C38), 1, mean)
  # c0sAphi[negA_rowidx,4] <- 0
  return(c0sAphi)
}

# TILING:
cycle5 <- read_csv("cycle5.txt")

dat_tiling <- cycle5

colnames(dat_tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

dat_tiling = dat_tiling %>%
  mutate(C19 = lag(C26),
         C22 = lag(C29),
         C24 = lag(C31),
         C33 = lead(C26),
         C36 = lead(C29),
         C38 = lead(C31))

# Remove first sequence from each gene:
dat_tiling = dat_tiling[-seq(1, nrow(dat_tiling), by=143),]
# Remove final sequence from each gene:
dat_tiling = dat_tiling[-seq(142, nrow(dat_tiling), by=142),]

tiling_C0sAphi_lag1 = find_c0sAphi_lag(dat_tiling, aa=c(1,1/0.82,1/0.7), lag_weight=1)

dat_tiling$C0_lag1 = tiling_C0sAphi_lag1[,1]
dat_tiling$Amplitude_lag1 = tiling_C0sAphi_lag1[,4]
dat_tiling$Phase_lag1 = tiling_C0sAphi_lag1[,5]

dat_tiling = dat_tiling %>%
  mutate(smooth_C26_v1 = (C19+C26+C33)/3,
         smooth_C29_v1 = (C22+C29+C36)/3,
         smooth_C31_v1 = (C24+C31+C38)/3)

dat_tiling = dat_tiling %>%
  mutate(smooth_C26_v2 = (43*C19+50*C26+43*C33)/(43+50+43),
         smooth_C29_v2 = (43*C22+50*C29+43*C36)/(43+50+43),
         smooth_C31_v2 = (43*C24+50*C31+43*C38)/(43+50+43))

dat_tiling = dat_tiling %>%
  mutate(smooth_C26_v3 = (43*C19+100*C26+43*C33)/(43+100+43),
         smooth_C29_v3 = (43*C22+100*C29+43*C36)/(43+100+43),
         smooth_C31_v3 = (43*C24+100*C31+43*C38)/(43+100+43))

tiling_C0Aphi_original = find_c0Aphi(dat_tiling %>% select(C26, C29, C31), aa=c(1,1/.82,1/.7))
tiling_C0Aphi_smooth_v1 = find_c0Aphi(dat_tiling %>% select(smooth_C26_v1, smooth_C29_v1, smooth_C31_v1), aa=c(1,1,1))
tiling_C0Aphi_smooth_v2 = find_c0Aphi(dat_tiling %>% select(smooth_C26_v2, smooth_C29_v2, smooth_C31_v2), aa=c(1,1,1))
tiling_C0Aphi_smooth_v3 = find_c0Aphi(dat_tiling %>% select(smooth_C26_v3, smooth_C29_v3, smooth_C31_v3), aa=c(1,1,1))

dat_tiling$C0_original = tiling_C0Aphi_original[,1]
dat_tiling$Amplitude_original = tiling_C0Aphi_original[,2]
dat_tiling$Phase_original = tiling_C0Aphi_original[,3]
dat_tiling$smooth_C0_v1 = tiling_C0Aphi_smooth_v1[,1]
dat_tiling$smooth_Amplitude_v1 = tiling_C0Aphi_smooth_v1[,2]
dat_tiling$smooth_Phase_v1 = tiling_C0Aphi_smooth_v1[,3]
dat_tiling$smooth_C0_v2 = tiling_C0Aphi_smooth_v2[,1]
dat_tiling$smooth_Amplitude_v2 = tiling_C0Aphi_smooth_v2[,2]
dat_tiling$smooth_Phase_v2 = tiling_C0Aphi_smooth_v2[,3]
dat_tiling$smooth_C0_v3 = tiling_C0Aphi_smooth_v3[,1]
dat_tiling$smooth_Amplitude_v3 = tiling_C0Aphi_smooth_v3[,2]
dat_tiling$smooth_Phase_v3 = tiling_C0Aphi_smooth_v3[,3]

write.csv(dat_tiling, "data/Created/tiling_smoothC0.csv")






# ChrV:
cycle6 <- read_csv("cycle6.txt")

dat_chrv <- cycle6

colnames(dat_chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

dat_chrv = dat_chrv %>%
  mutate(C19 = lag(C26),
         C22 = lag(C29),
         C24 = lag(C31),
         C33 = lead(C26),
         C36 = lead(C29),
         C38 = lead(C31))

# Remove first sequence:
dat_chrv = dat_chrv[-1,]
# Remove final sequence:
dat_chrv = dat_chrv[-nrow(dat_chrv),]

chrv_C0sAphi_lag1 = find_c0sAphi_lag(dat_chrv, aa=c(1,1,1), lag_weight=0.43/2)

dat_chrv$C0_lag1 = chrv_C0sAphi_lag1[,1]
dat_chrv$Amplitude_lag1 = chrv_C0sAphi_lag1[,4]
dat_chrv$Phase_lag1 = chrv_C0sAphi_lag1[,5]

dat_chrv = dat_chrv %>%
  mutate(smooth_C26_v1 = (C19+C26+C33)/3,
         smooth_C29_v1 = (C22+C29+C36)/3,
         smooth_C31_v1 = (C24+C31+C38)/3)

dat_chrv = dat_chrv %>%
  mutate(smooth_C26_v2 = (43*C19+50*C26+43*C33)/(43+50+43),
         smooth_C29_v2 = (43*C22+50*C29+43*C36)/(43+50+43),
         smooth_C31_v2 = (43*C24+50*C31+43*C38)/(43+50+43))

dat_chrv = dat_chrv %>%
  mutate(smooth_C26_v3 = (43*C19+100*C26+43*C33)/(43+100+43),
         smooth_C29_v3 = (43*C22+100*C29+43*C36)/(43+100+43),
         smooth_C31_v3 = (43*C24+100*C31+43*C38)/(43+100+43))

chrv_C0Aphi_original = find_c0Aphi(dat_chrv %>% select(C26, C29, C31), aa=c(1,1/.82,1/.7))
chrv_C0Aphi_smooth_v1 = find_c0Aphi(dat_chrv %>% select(smooth_C26_v1, smooth_C29_v1, smooth_C31_v1), aa=c(1,1,1))
chrv_C0Aphi_smooth_v2 = find_c0Aphi(dat_chrv %>% select(smooth_C26_v2, smooth_C29_v2, smooth_C31_v2), aa=c(1,1,1))
chrv_C0Aphi_smooth_v3 = find_c0Aphi(dat_chrv %>% select(smooth_C26_v3, smooth_C29_v3, smooth_C31_v3), aa=c(1,1,1))

dat_chrv$C0_original = chrv_C0Aphi_original[,1]
dat_chrv$Amplitude_original = chrv_C0Aphi_original[,2]
dat_chrv$Phase_original = chrv_C0Aphi_original[,3]
dat_chrv$smooth_C0_v1 = chrv_C0Aphi_smooth_v1[,1]
dat_chrv$smooth_Amplitude_v1 = chrv_C0Aphi_smooth_v1[,2]
dat_chrv$smooth_Phase_v1 = chrv_C0Aphi_smooth_v1[,3]
dat_chrv$smooth_C0_v2 = chrv_C0Aphi_smooth_v2[,1]
dat_chrv$smooth_Amplitude_v2 = chrv_C0Aphi_smooth_v2[,2]
dat_chrv$smooth_Phase_v2 = chrv_C0Aphi_smooth_v2[,3]
dat_chrv$smooth_C0_v3 = chrv_C0Aphi_smooth_v3[,1]
dat_chrv$smooth_Amplitude_v3 = chrv_C0Aphi_smooth_v3[,2]
dat_chrv$smooth_Phase_v3 = chrv_C0Aphi_smooth_v3[,3]

write.csv(dat_chrv, "data/Created/chrv_smoothC0.csv")






# dat = dat %>%
#   # mutate(smooth_C0 = (smooth_C26 + smooth_C29 + smooth_C31)/3)
#   mutate(smooth_C0 = (smooth_C26 + smooth_C31)/2)

plot_smoothC_by_library_position = function(data, starting_point=1, num_sequences=10, 
                                            include_C0=FALSE, include_smoothC0=FALSE) {
  row_sequence = seq(from=starting_point, by=1, length.out=num_sequences)
  position_sequence = (row_sequence-1)*7

  if (include_C0) { 
    alpha1=0.3
    alpha2=1
    
    ggplot(mapping=aes(x=position_sequence)) +
      geom_line(aes(y=data[row_sequence, "smooth_C26"] %>% unlist(), color="Smooth C26"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence, "smooth_C29"] %>% unlist(), color="Smooth C29"), alpha=alpha1) + 
      geom_line(aes(y=data[row_sequence, "smooth_C31"] %>% unlist(), color="Smooth C31"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence,"C0"] %>% unlist(), color="C0"), alpha=alpha2) +
      xlab("Position") +
      ylab("Value") +
      scale_color_manual(values=c("Smooth C26"="red", "Smooth C29"="blue", "Smooth C31"="green", "C0"="brown"))
  }
  else if(include_smoothC0) {
    alpha1=0.3
    alpha2=1
    
    ggplot(mapping=aes(x=position_sequence)) +
      geom_line(aes(y=data[row_sequence, "C26"] %>% unlist(), color="C26"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence, "C29"] %>% unlist(), color="C29"), alpha=alpha1) + 
      geom_line(aes(y=data[row_sequence, "C31"] %>% unlist(), color="C31"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence,"smooth_C0"] %>% unlist(), color="Smooth C0"), alpha=alpha2) +
      xlab("Position") +
      ylab("Value") +
      scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green", "Smooth C0"="brown"))
  }
  else {
    alpha1=0.3
    alpha2=1
    
    ggplot(mapping=aes(x=position_sequence)) +
      geom_line(aes(y=data[row_sequence, "C26"] %>% unlist(), color="C26"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence, "C29"] %>% unlist(), color="C29"), alpha=alpha1) + 
      geom_line(aes(y=data[row_sequence, "C31"] %>% unlist(), color="C31"), alpha=alpha1) +
      geom_line(aes(y=data[row_sequence, "smooth_C26"] %>% unlist(), color="C26"), alpha=alpha2) +
      geom_line(aes(y=data[row_sequence, "smooth_C29"] %>% unlist(), color="C29"), alpha=alpha2) + 
      geom_line(aes(y=data[row_sequence, "smooth_C31"] %>% unlist(), color="C31"), alpha=alpha2) +
      xlab("Position") +
      ylab("Value") +
      scale_color_manual(values=c("C26"="red", "C29"="blue", "C31"="green"))
  }
}

plot1_tiling_smooth_noC0 = plot_smoothC_by_library_position(dat_tiling, starting_point = 1, num_sequences=141)
plot1_tiling_smooth_noC0
plot1_tiling_smooth_C0 = plot_smoothC_by_library_position(dat_tiling, starting_point = 1, num_sequences=141, include_C0=TRUE)
plot1_tiling_smooth_C0
grid.arrange(plot1_tiling_smooth_noC0, plot1_tiling_smooth_C0, nrow=2)

plot1_tiling_smooth_smoothC0 = plot_smoothC_by_library_position(dat_tiling, starting_point = 1, num_sequences=141, include_smoothC0=TRUE)
plot1_tiling_smooth_smoothC0
grid.arrange(plot1_tiling_C0, plot1_tiling_smooth_smoothC0, nrow=2)



plot1_chrv_starting_point = 213500/7
plot1_chrv_smooth_noC0 = plot_smoothC_by_library_position(dat_chrv, starting_point = plot1_chrv_starting_point, num_sequences=50)
plot1_chrv_smooth_noC0
plot1_chrv_smooth_C0 = plot_smoothC_by_library_position(dat_chrv, starting_point = plot1_chrv_starting_point, num_sequences=50, include_C0=TRUE)
plot1_chrv_smooth_C0
grid.arrange(plot1_chrv_smooth_noC0, plot1_chrv_smooth_C0, nrow=2)

plot1_chrv_smooth_smoothC0 = plot_smoothC_by_library_position(dat_chrv, starting_point = plot1_chrv_starting_point, num_sequences=50, include_smoothC0=TRUE)
plot1_chrv_smooth_smoothC0
grid.arrange(plot1_chrv_C0, plot1_chrv_smooth_smoothC0, nrow=2)




ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(Amplitude_original,0.25)), aes(x=Phase_original, y=C0_original))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v3,0.25)), aes(x=smooth_Phase_v3, y=smooth_C0_v3))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v2,0.25)), aes(x=smooth_Phase_v2, y=smooth_C0_v2))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v1,0.25)), aes(x=smooth_Phase_v1, y=smooth_C0_v1))+
  geom_point(alpha=0.05)+
  geom_smooth()

ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(Amplitude_original,0.25)), aes(x=cos(Phase_original), y=C0_original))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v3,0.25)), aes(x=cos(smooth_Phase_v3), y=smooth_C0_v3))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v2,0.25)), aes(x=cos(smooth_Phase_v2), y=smooth_C0_v2))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_tiling %>% filter(Amplitude_original > quantile(smooth_Amplitude_v1,0.25)), aes(x=cos(smooth_Phase_v1), y=smooth_C0_v1))+
  geom_point(alpha=0.05)+
  geom_smooth()

library(energy)

dcor(dat_tiling$C0_original, dat_tiling$Phase_original)

cor(dat_tiling$C0_original, cos(dat_tiling$Phase_original))
# 0.1143306
cor(dat_tiling$smooth_C0_v3, cos(dat_tiling$smooth_Phase_v3))
# 0.2675596
cor(dat_tiling$smooth_C0_v2, cos(dat_tiling$smooth_Phase_v2))
# 0.3455052
cor(dat_tiling$smooth_C0_v1, cos(dat_tiling$smooth_Phase_v1))
# 0.3486609

ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(Amplitude_original,0.25)), aes(x=Phase_original, y=C0_original))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v3,0.25)), aes(x=smooth_Phase_v3, y=smooth_C0_v3))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v2,0.25)), aes(x=smooth_Phase_v2, y=smooth_C0_v2))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v1,0.25)), aes(x=smooth_Phase_v1, y=smooth_C0_v1))+
  geom_point(alpha=0.05)+
  geom_smooth()

ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(Amplitude_original,0.25)), aes(x=cos(Phase_original), y=C0_original))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v3,0.25)), aes(x=cos(smooth_Phase_v3), y=smooth_C0_v3))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v2,0.25)), aes(x=cos(smooth_Phase_v2), y=smooth_C0_v2))+
  geom_point(alpha=0.05)+
  geom_smooth()
ggplot(data=dat_chrv %>% filter(Amplitude_original > quantile(smooth_Amplitude_v1,0.25)), aes(x=cos(smooth_Phase_v1), y=smooth_C0_v1))+
  geom_point(alpha=0.05)+
  geom_smooth()

cor(dat_chrv$C0_original, cos(dat_chrv$Phase_original))
# 0.1877908
cor(dat_chrv$smooth_C0_v3, cos(dat_chrv$smooth_Phase_v3))
# 0.2117358
cor(dat_chrv$smooth_C0_v2, cos(dat_chrv$smooth_Phase_v2))
# 0.2386066
cor(dat_chrv$smooth_C0_v1, cos(dat_chrv$smooth_Phase_v1))
# 0.2407397

