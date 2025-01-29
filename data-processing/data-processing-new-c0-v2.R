n <- c(26, 29, 31)
aa <- c(1, 0.82, 0.7)
k <- 2*pi/10.4

find_c0Aphi <- function(dat){
  mat <- matrix(0, nrow = 3, ncol = 3)
  mat[1:3,1] <- 1
  mat[1:3,2] <- sin(n*k)
  mat[1:3,3] <- cos(n*k)
  mat[1,2:3] <- mat[1,2:3]*aa[1]
  mat[2,2:3] <- mat[2,2:3]*aa[2]
  mat[3,2:3] <- mat[3,2:3]*aa[3]
  inv_mat <- solve(mat)
  c0A1A2 <- as.matrix(dat[,c("C26", "C29", "C31")]) %*% t(inv_mat)
  c0Aphi <- c0A1A2
  c0Aphi[,1] <- c0A1A2[,1]
  c0Aphi[,2] <- sqrt(c0A1A2[,2]^2 + c0A1A2[,3]^2)
  c0Aphi[,3] <- sign(c0A1A2[,3]) * acos(c0A1A2[,2]/c0Aphi[,2])
  return(c0Aphi)
}





# Find the correct C0, Amplitude, Phase for the tiling data:

# Tiling library
cycle5 <- read_csv("cycle5.txt")
colnames(cycle5)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle5[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(cycle5)

set.seed(50)

train_indices_cycle5 = sample(1:nrow(cycle5), nrow(cycle5)*.9, replace=FALSE)

cycle5_train = cycle5[train_indices_cycle5,] %>%
  select(-Sequence, -C0)
cycle5_test = cycle5[-train_indices_cycle5,] %>%
  select(-Sequence, -C0)

dat_tiling = readRDS("data/Created/processed_tiling_newC0.rds")
dat_tiling_test = readRDS("data/Created/processed_tiling_test_newC0.rds")

dat_tiling[, colnames(cycle5_train)] = cycle5_train
dat_tiling_test[, colnames(cycle5_test)] = cycle5_test

saveRDS(dat_tiling, "data/Created/processed_tiling_newC0_v2.rds")
saveRDS(dat_tiling_test, "data/Created/processed_tiling_test_newC0_v2.rds")


# test the correctness of the algorithm: print the predicted value and compare it with the original data 
test_tiling <- function(ind){
  print(c(dat_tiling$C26[ind],dat_tiling$C29[ind],dat_tiling$C31[ind]))
  print(dat_tiling$C0_new[ind]+ aa * dat_tiling$Amplitude_new[ind] * sin(n*k+ dat_tiling$Phase_new[ind]))
}
ind <- sample(1:nrow(dat_tiling),1)
test_tiling(ind)







# Random library
cycle3 <- read_csv("cycle3.txt")
colnames(cycle3)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle3[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(cycle3)

set.seed(50)

train_indices_cycle3 = sample(1:nrow(cycle3),
                              nrow(cycle3)*.9, replace=FALSE)

cycle3_train = cycle3[train_indices_cycle3,] %>%
  select(-Sequence, -C0)
cycle3_test = cycle3[-train_indices_cycle3,] %>%
  select(-Sequence, -C0)

dat_random = readRDS("data/Created/processed_random_newC0.rds")
dat_random_test = readRDS("data/Created/processed_random_test_newC0.rds")

dat_random[, colnames(cycle3_train)] = cycle3_train
dat_random_test[, colnames(cycle3_test)] = cycle3_test

saveRDS(dat_random, "data/Created/processed_random_newC0_v2.rds")
saveRDS(dat_random_test, "data/Created/processed_random_test_newC0_v2.rds")

test_random <- function(ind){
  print(c(dat_random$C26[ind],dat_random$C29[ind],dat_random$C31[ind]))
  print(dat_random$C0_new[ind]+ aa * dat_random$Amplitude_new[ind] * sin(n*k+ dat_random$Phase_new[ind]))
}
ind <- sample(1:nrow(dat_random),1)
test_random(ind)










# ChrV library
cycle6 <- read_csv("cycle6.txt")
colnames(cycle6)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle6[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(cycle6)

set.seed(50)

train_indices_cycle6 = sample(1:nrow(cycle6),
                              nrow(cycle6)*.9, replace=FALSE)

cycle6_train = cycle6[train_indices_cycle6,] %>%
  select(-Sequence, -C0)
cycle6_test = cycle6[-train_indices_cycle6,] %>%
  select(-Sequence, -C0)

dat_chrV = readRDS("data/Created/processed_chrV_newC0.rds")
dat_chrV_test = readRDS("data/Created/processed_chrV_test_newC0.rds")

dat_chrV[, colnames(cycle6_train)] = cycle6_train
dat_chrV_test[, colnames(cycle6_test)] = cycle6_test

saveRDS(dat_chrV, "data/Created/processed_chrV_newC0_v2.rds")
saveRDS(dat_chrV_test, "data/Created/processed_chrV_test_newC0_v2.rds")

test_chrV <- function(ind){
  print(c(dat_chrV$C26[ind],dat_chrV$C29[ind],dat_chrV$C31[ind]))
  print(dat_chrV$C0_new[ind]+ aa * dat_chrV$Amplitude_new[ind] * sin(n*k+ dat_chrV$Phase_new[ind]))
}
ind <- sample(1:nrow(dat_chrV),1)
test_chrV(ind)








# Yeast library
cycle1 <- read_csv("cycle1.txt")
colnames(cycle1)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle1[,c("C0_new", "Amplitude_new", "Phase_new")] = find_c0Aphi(cycle1)

set.seed(50)

train_indices_cycle1 = sample(1:nrow(cycle1),
                              nrow(cycle1)*.9, replace=FALSE)

cycle1_train = cycle1[train_indices_cycle1,] %>%
  select(-Sequence, -C0)
cycle1_test = cycle1[-train_indices_cycle1,] %>%
  select(-Sequence, -C0)

dat_yeast = readRDS("data/Created/processed_yeast_newC0.rds")
dat_yeast_test = readRDS("data/Created/processed_yeast_test_newC0.rds")

dat_yeast[, colnames(cycle1_train)] = cycle1_train
dat_yeast_test[, colnames(cycle1_test)] = cycle1_test

saveRDS(dat_yeast, "data/Created/processed_yeast_newC0_v2.rds")
saveRDS(dat_yeast_test, "data/Created/processed_yeast_test_newC0_v2.rds")

test_yeast <- function(ind){
  print(c(dat_yeast$C26[ind],dat_yeast$C29[ind],dat_yeast$C31[ind]))
  print(dat_yeast$C0_new[ind]+ aa * dat_yeast$Amplitude_new[ind] * sin(n*k+ dat_yeast$Phase_new[ind]))
}
ind <- sample(1:nrow(dat_yeast),1)
test_yeast(ind)





