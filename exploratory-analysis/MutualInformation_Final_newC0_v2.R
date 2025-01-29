library(igraph)
library(infotheo)
library(ggplot2)

dat_tiling = readRDS("data/Created/processed_tiling_newC0_v2.rds")
dat_random = readRDS("data/Created/processed_random_newC0_v2.rds")


### TILING LIBRARY:


cutoffs_0_25_tiling = quantile(dat_tiling$C0_new, c(0, 0.25)) 
cutoffs_25_50_tiling = quantile(dat_tiling$C0_new, c(0.25, 0.5))
cutoffs_50_75_tiling = quantile(dat_tiling$C0_new, c(0.5, 0.75))
cutoffs_75_100_tiling = quantile(dat_tiling$C0_new, c(0.75, 1))

# First Quartile:
dat_q1_tiling = dat_tiling %>%
  filter(C0_new >= cutoffs_0_25_tiling[1] & C0_new <= cutoffs_0_25_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q1_tiling = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q1_tiling) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q1_tiling[z, "Var1"] = i
    d_q1_tiling[z, "Var2"] = j
    d_q1_tiling[z, "MI"] = mutinformation(dat_q1_tiling[,paste0("X", i, "di")], dat_q1_tiling[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q1_tiling = d_q1_tiling %>%
  arrange(desc(MI))

map_q1_tiling <- d_q1_tiling %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.0535)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") + 
  ggtitle("Bottom 25% of Tiling Library")



# Second Quartile:
dat_q2_tiling = dat_tiling %>%
  filter(C0_new >= cutoffs_25_50_tiling[1] & C0_new <= cutoffs_25_50_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q2_tiling = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q2_tiling) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q2_tiling[z, "Var1"] = i
    d_q2_tiling[z, "Var2"] = j
    d_q2_tiling[z, "MI"] = mutinformation(dat_q2_tiling[,paste0("X", i, "di")], dat_q2_tiling[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q2_tiling = d_q2_tiling %>%
  arrange(desc(MI))

map_q2_tiling <- d_q2_tiling %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.0535)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1")





# Third Quartile:
dat_q3_tiling = dat_tiling %>%
  filter(C0_new >= cutoffs_50_75_tiling[1] & C0_new <= cutoffs_50_75_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q3_tiling = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q3_tiling) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q3_tiling[z, "Var1"] = i
    d_q3_tiling[z, "Var2"] = j
    d_q3_tiling[z, "MI"] = mutinformation(dat_q3_tiling[,paste0("X", i, "di")], dat_q3_tiling[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q3_tiling = d_q3_tiling %>%
  arrange(desc(MI))

map_q3_tiling <- d_q3_tiling %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.0535)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1")





# Fourth Quartile:
dat_q4_tiling = dat_tiling %>%
  filter(C0_new >= cutoffs_75_100_tiling[1] & C0_new <= cutoffs_75_100_tiling[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q4_tiling = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q4_tiling) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q4_tiling[z, "Var1"] = i
    d_q4_tiling[z, "Var2"] = j
    d_q4_tiling[z, "MI"] = mutinformation(dat_q4_tiling[,paste0("X", i, "di")], dat_q4_tiling[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q4_tiling = d_q4_tiling %>%
  arrange(desc(MI))

map_q4_tiling <- d_q4_tiling %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.0535)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") + 
  ggtitle("Top 25% of Tiling Library")




# Print heatmaps:
map_q1_tiling
map_q2_tiling
map_q3_tiling
map_q4_tiling


ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/tiling_ps2_q1_heatmap_newC0_v2.png", plot = map_q1_tiling)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/tiling_ps2_q2_heatmap_newC0_v2.png", plot = map_q2_tiling)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/tiling_ps2_q3_heatmap_newC0_v2.png", plot = map_q3_tiling)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/tiling_ps2_q4_heatmap_newC0_v2.png", plot = map_q4_tiling)









### RANDOM LIBRARY:

cutoffs_0_25_random = quantile(dat_random$C0_new, c(0, 0.25)) 
cutoffs_25_50_random = quantile(dat_random$C0_new, c(0.25, 0.5))
cutoffs_50_75_random = quantile(dat_random$C0_new, c(0.5, 0.75))
cutoffs_75_100_random = quantile(dat_random$C0_new, c(0.75, 1))

# First Quartile:
dat_q1_random = dat_random %>%
  filter(C0_new >= cutoffs_0_25_random[1] & C0_new <= cutoffs_0_25_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q1_random = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q1_random) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q1_random[z, "Var1"] = i
    d_q1_random[z, "Var2"] = j
    d_q1_random[z, "MI"] = mutinformation(dat_q1_random[,paste0("X", i, "di")], dat_q1_random[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q1_random = d_q1_random %>%
  arrange(desc(MI))

map_q1_random <- d_q1_random %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.061)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") + 
  ggtitle("Bottom 25% of Random Library")



# Second Quartile:
dat_q2_random = dat_random %>%
  filter(C0_new >= cutoffs_25_50_random[1] & C0_new <= cutoffs_25_50_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q2_random = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q2_random) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q2_random[z, "Var1"] = i
    d_q2_random[z, "Var2"] = j
    d_q2_random[z, "MI"] = mutinformation(dat_q2_random[,paste0("X", i, "di")], dat_q2_random[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q2_random = d_q2_random %>%
  arrange(desc(MI))

map_q2_random <- d_q2_random %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.061)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1")





# Third Quartile:
dat_q3_random = dat_random %>%
  filter(C0_new >= cutoffs_50_75_random[1] & C0_new <= cutoffs_50_75_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q3_random = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q3_random) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q3_random[z, "Var1"] = i
    d_q3_random[z, "Var2"] = j
    d_q3_random[z, "MI"] = mutinformation(dat_q3_random[,paste0("X", i, "di")], dat_q3_random[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q3_random = d_q3_random %>%
  arrange(desc(MI))

map_q3_random <- d_q3_random %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.061)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1")





# Fourth Quartile:
dat_q4_random = dat_random %>%
  filter(C0_new >= cutoffs_75_100_random[1] & C0_new <= cutoffs_75_100_random[2]) %>%
  select(all_of(ps1), all_of(ps2), C0_new)

d_q4_random = as.data.frame(matrix(nrow = 1128, ncol = 3))
colnames(d_q4_random) = c("Var1", "Var2", "MI")
z=1
for(i in 1:47) {
  for(j in (i+2):49) {
    d_q4_random[z, "Var1"] = i
    d_q4_random[z, "Var2"] = j
    d_q4_random[z, "MI"] = mutinformation(dat_q4_random[,paste0("X", i, "di")], dat_q4_random[,paste0("X", j, "di")])
    z=z+1
  }
}

d_q4_random = d_q4_random %>%
  arrange(desc(MI))

map_q4_random <- d_q4_random %>%
  ggplot(aes(x=Var1, y=Var2, fill=MI)) + 
  geom_raster() + 
  scale_fill_gradient2(limits = c(0, 0.061)) + 
  ylab("Dinucleotide Position 2") +
  xlab("Dinucleotide Position 1") + 
  ggtitle("Top 25% of Random Library")




# Print heatmaps:
map_q1_random
map_q2_random
map_q3_random
map_q4_random


ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/random_ps2_q1_heatmap_newC0_v2.png", plot = map_q1_random)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/random_ps2_q2_heatmap_newC0_v2.png", plot = map_q2_random)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/random_ps2_q3_heatmap_newC0_v2.png", plot = map_q3_random)
ggsave("~/Documents/Northwestern/DNA_Cyclizability/figures/mutual_information/random_ps2_q4_heatmap_newC0_v2.png", plot = map_q4_random)

