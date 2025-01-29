dat = readRDS("data/Created/processed_tiling_newC0.rds")
cycle1 <- read.delim("cycle1_may2023.txt")
colnames(cycle1)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle3 <- read_csv("cycle3.txt")
colnames(cycle3)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle5 <- read_csv("cycle5.txt")
colnames(cycle5)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")
cycle6 <- read_csv("cycle6.txt")
colnames(cycle6)=c("Sequence","C26","C29","C31", "C0", "Amplitude","Phase")

k=2*pi/10.4

print_Cs = function(ind, dat) {
  print(c(dat$C26[ind], dat$C0[ind] - dat$Amplitude[ind]*sin(26*k + dat$Phase[ind])))
  print(c(dat$C29[ind], dat$C0[ind] - dat$Amplitude[ind]/.82*sin(29*k + dat$Phase[ind])))
  print(c(dat$C31[ind], dat$C0[ind] - dat$Amplitude[ind]/.7*sin(31*k + dat$Phase[ind])))
}

test_ind = sample(1:nrow(dat), 1)
print_Cs(test_ind, dat)
# Inconsistent

test_ind = sample(1:nrow(cycle1), 1)
print_Cs(test_ind, cycle1)
# Inconsistent

test_ind = sample(1:nrow(cycle3), 1)
print_Cs(test_ind, cycle3)
# Inconsistent

test_ind = sample(1:nrow(cycle5), 1)
print_Cs(test_ind, cycle5)
# Inconsistent

test_ind = sample(1:nrow(cycle6), 1)
print_Cs(test_ind, cycle6)
# Inconsistent


dat = readRDS("data/Created/processed_tiling_newC0_v2.rds")

print_Cnews = function(ind, dat) {
  print(c(dat$C26[ind], dat$C0_new[ind] + dat$Amplitude_new[ind]*sin(26*k + dat$Phase_new[ind])))
  print(c(dat$C29[ind], dat$C0_new[ind] + dat$Amplitude_new[ind]*.82*sin(29*k + dat$Phase_new[ind])))
  print(c(dat$C31[ind], dat$C0_new[ind] + dat$Amplitude_new[ind]*.7*sin(31*k + dat$Phase_new[ind])))
}

test_ind = sample(1:nrow(dat), 1)
print_Cnews(1, dat)

