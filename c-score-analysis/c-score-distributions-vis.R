nuc <- read_csv("cycle1.txt")
random <- read_csv("cycle3.txt")
tiling <- read_csv("cycle5.txt")
chrv <- read_csv("cycle6.txt")

colnames(nuc) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(random) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(tiling) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")
colnames(chrv) = c("x50mer", "C26", "C29", "C31", "C0", "Amplitude", "Phase")

nuc = nuc %>% 
  select(C26, C29, C31) %>%
  pivot_longer(cols=everything())

random = random %>% 
  select(C26, C29, C31) %>%
  pivot_longer(cols=everything())

tiling = tiling %>% 
  select(C26, C29, C31) %>%
  pivot_longer(cols=everything())

chrv = chrv %>% 
  select(C26, C29, C31) %>%
  pivot_longer(cols=everything())

# xmin_nuc <- min(nuc$C26, nuc$C29, nuc$C31)
# xmax_nuc <- max(nuc$C26, nuc$C29, nuc$C31)
# # par(mfrow=c(3,1))
# hist(nuc$C26, main = "C26", xlab = "Value", ylab = "Frequency",
#      col = adjustcolor( "dodgerblue", alpha.f = 0.2), xlim = c(xmin_nuc, xmax_nuc))
# hist(nuc$C29, main = "C29", xlab = "Value", ylab = "Frequency",
#      col = adjustcolor( "red", alpha.f = 0.2), add = TRUE, xlim = c(xmin_nuc, xmax_nuc))
# hist(nuc$C31, main = "C31", xlab = "Value", ylab = "Frequency",
#      col = adjustcolor( "green", alpha.f = 0.2), add = TRUE, xlim = c(xmin_nuc, xmax_nuc))
# legend()

library("ggplot2")
ggplot(nuc, aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  ggtitle("Nucleosome Library")

ggplot(random, aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  ggtitle("Random Library")

ggplot(tiling, aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  ggtitle("Tiling Library")

ggplot(chrv, aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  ggtitle("ChrV Library")




