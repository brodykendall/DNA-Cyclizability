library(tidyverse)
library(infotheo)
library(ggplot2)

# Load data:
dat = readRDS("data/Created/processed_ratio.rds")

# Find cutoffs of c-scores for each quartile
cutoffs_0_25 = quantile(dat$C0, c(0, 0.25)) 
cutoffs_25_50 = quantile(dat$C0, c(0.25, 0.5))
cutoffs_50_75 = quantile(dat$C0, c(0.5, 0.75))
cutoffs_75_100 = quantile(dat$C0, c(0.75, 1))

# Position specific variable names 
# (e.g. X1di represents the dinucleotide in the 1st position of the sequence):
ps1 <- paste0("X", 1:50, "mono")
ps2 <- paste0("X", 1:49, "di")

# Divide data into quartiles:
dat_q1 = dat %>%
  filter(C0 >= cutoffs_0_25[1] & C0 <= cutoffs_0_25[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
dat_q2 = dat %>%
  filter(C0 >= cutoffs_25_50[1] & C0 <= cutoffs_25_50[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
dat_q3 = dat %>%
  filter(C0 >= cutoffs_50_75[1] & C0 <= cutoffs_50_75[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)
dat_q4 = dat %>%
  filter(C0 >= cutoffs_75_100[1] & C0 <= cutoffs_75_100[2]) %>%
  select(all_of(ps1), all_of(ps2), C0)

# Construct outline of mutual information dataframes:
d_q1 = as.data.frame(matrix(nrow=1225, ncol = 3))
colnames(d_q1) = c("Var1", "Var2", "MI")
d_q1[,c("Var1", "Var2")] = matrix(cross2(1:50, 1:50, .filter = ~.x >= .y) %>% unlist(),
                                  nrow=1225, ncol=2, byrow=TRUE)
d_q2 = as.data.frame(matrix(nrow=1225, ncol = 3))
colnames(d_q2) = c("Var1", "Var2", "MI")
d_q2[,c("Var1", "Var2")] = matrix(cross2(1:50, 1:50, .filter = ~.x >= .y) %>% unlist(),
                                  nrow=1225, ncol=2, byrow=TRUE)
d_q3 = as.data.frame(matrix(nrow=1225, ncol = 3))
colnames(d_q3) = c("Var1", "Var2", "MI")
d_q3[,c("Var1", "Var2")] = matrix(cross2(1:50, 1:50, .filter = ~.x >= .y) %>% unlist(),
                                  nrow=1225, ncol=2, byrow=TRUE)
d_q4 = as.data.frame(matrix(nrow=1225, ncol = 3))
colnames(d_q4) = c("Var1", "Var2", "MI")
d_q4[,c("Var1", "Var2")] = matrix(cross2(1:50, 1:50, .filter = ~.x >= .y) %>% unlist(),
                                  nrow=1225, ncol=2, byrow=TRUE)

# Calculate pairwise mutual information for each position:
d_q1[,"MI"] = apply(d_q1, 1, function(row) {
  mutinformation(dat_q1[,paste0("X",row[1],"mono")], dat_q1[,paste0("X",row[2],"mono")])
})
d_q2[,"MI"] = apply(d_q2, 1, function(row) {
  mutinformation(dat_q2[,paste0("X",row[1],"mono")], dat_q2[,paste0("X",row[2],"mono")])
})
d_q3[,"MI"] = apply(d_q3, 1, function(row) {
  mutinformation(dat_q3[,paste0("X",row[1],"mono")], dat_q3[,paste0("X",row[2],"mono")])
})
d_q4[,"MI"] = apply(d_q4, 1, function(row) {
  mutinformation(dat_q4[,paste0("X",row[1],"mono")], dat_q4[,paste0("X",row[2],"mono")])
})

# Filter out codon regions:
#   Note: could exclude this filtration to include codon regions if desired
d_q1.nc = d_q1 %>%
  filter(Var2 - Var1 >= 3)
d_q2.nc = d_q2 %>%
  filter(Var2 - Var1 >= 3)
d_q3.nc = d_q3 %>%
  filter(Var2 - Var1 >= 3)
d_q4.nc = d_q4 %>%
  filter(Var2 - Var1 >= 3)

# Find maximum pairwise mutual information across all four quartiles for scaling:
ps1_maxMI = max(d_q1.nc$MI, d_q2.nc$MI, d_q3.nc$MI, d_q4.nc$MI)

# Construct plots:
map_q1.nc <- d_q1.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps1_maxMI)) + 
  ylab("Nucleotide Position 2") +
  xlab("Nucleotide Position 1") +
  ggtitle("Mutual Information of Nucloetides by Position, First Quartile")
map_q2.nc <- d_q2.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps1_maxMI)) + 
  ylab("Nucleotide Position 2") +
  xlab("Nucleotide Position 1") +
  ggtitle("Mutual Information of Nucloetides by Position, Second Quartile")
map_q3.nc <- d_q3.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps1_maxMI)) + 
  ylab("Nucleotide Position 2") +
  xlab("Nucleotide Position 1") +
  ggtitle("Mutual Information of Nucloetides by Position, Third Quartile")
map_q4.nc <- d_q4.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps1_maxMI)) + 
  ylab("Nucleotide Position 2") +
  xlab("Nucleotide Position 1") +
  ggtitle("Mutual Information of Nucloetides by Position, Fourth Quartile")

# Display plots:
map_q1.nc
map_q2.nc
map_q3.nc
map_q4.nc







##########################################################################
# Dinucleotides:
##########################################################################

# Construct outline of mutual information dataframes:
d_q1.ps2 = as.data.frame(matrix(nrow=1176, ncol = 3))
colnames(d_q1.ps2) = c("Var1", "Var2", "MI")
d_q1.ps2[,c("Var1", "Var2")] = matrix(cross2(1:49, 1:49, .filter = ~.x >= .y) %>% unlist(),
                                  nrow=1176, ncol=2, byrow=TRUE)
d_q2.ps2 = as.data.frame(matrix(nrow=1176, ncol = 3))
colnames(d_q2.ps2) = c("Var1", "Var2", "MI")
d_q2.ps2[,c("Var1", "Var2")] = matrix(cross2(1:49, 1:49, .filter = ~.x >= .y) %>% unlist(),
                                      nrow=1176, ncol=2, byrow=TRUE)
d_q3.ps2 = as.data.frame(matrix(nrow=1176, ncol = 3))
colnames(d_q3.ps2) = c("Var1", "Var2", "MI")
d_q3.ps2[,c("Var1", "Var2")] = matrix(cross2(1:49, 1:49, .filter = ~.x >= .y) %>% unlist(),
                                      nrow=1176, ncol=2, byrow=TRUE)
d_q4.ps2 = as.data.frame(matrix(nrow=1176, ncol = 3))
colnames(d_q4.ps2) = c("Var1", "Var2", "MI")
d_q4.ps2[,c("Var1", "Var2")] = matrix(cross2(1:49, 1:49, .filter = ~.x >= .y) %>% unlist(),
                                      nrow=1176, ncol=2, byrow=TRUE)

# Calculate pairwise mutual information for each position:
d_q1.ps2[,"MI"] = apply(d_q1.ps2, 1, function(row) {
  mutinformation(dat_q1[,paste0("X",row[1],"di")], dat_q1[,paste0("X",row[2],"di")])
})
d_q2.ps2[,"MI"] = apply(d_q2.ps2, 1, function(row) {
  mutinformation(dat_q2[,paste0("X",row[1],"di")], dat_q2[,paste0("X",row[2],"di")])
})
d_q3.ps2[,"MI"] = apply(d_q3.ps2, 1, function(row) {
  mutinformation(dat_q3[,paste0("X",row[1],"di")], dat_q3[,paste0("X",row[2],"di")])
})
d_q4.ps2[,"MI"] = apply(d_q4.ps2, 1, function(row) {
  mutinformation(dat_q4[,paste0("X",row[1],"di")], dat_q4[,paste0("X",row[2],"di")])
})

# Filter out codon regions:
#   Note: could exclude this filtration to include codon regions if desired
d_q1.ps2.nc = d_q1.ps2 %>%
  filter(Var2 - Var1 >= 3)
d_q2.ps2.nc = d_q2.ps2 %>%
  filter(Var2 - Var1 >= 3)
d_q3.ps2.nc = d_q3.ps2 %>%
  filter(Var2 - Var1 >= 3)
d_q4.ps2.nc = d_q4.ps2 %>%
  filter(Var2 - Var1 >= 3)

# Find maximum pairwise mutual information across all four quartiles for scaling:
ps2_maxMI = max(d_q1.ps2.nc$MI, d_q2.ps2.nc$MI, d_q3.ps2.nc$MI, d_q4.ps2.nc$MI)

# Construct plots:
map_q1.ps2.nc <- d_q1.ps2.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps2_maxMI)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") +
  ggtitle("Mutual Information of Diucloetides by Position, First Quartile")
map_q2.ps2.nc <- d_q2.ps2.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps2_maxMI)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") +
  ggtitle("Mutual Information of Diucloetides by Position, Second Quartile")
map_q3.ps2.nc <- d_q3.ps2.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps2_maxMI)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") +
  ggtitle("Mutual Information of Diucloetides by Position, Third Quartile")
map_q4.ps2.nc <- d_q4.ps2.nc %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, ps2_maxMI)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") +
  ggtitle("Mutual Information of Diucloetides by Position, Fourth Quartile")

# Display plots:
map_q1.ps2.nc
map_q2.ps2.nc
map_q3.ps2.nc
map_q4.ps2.nc
