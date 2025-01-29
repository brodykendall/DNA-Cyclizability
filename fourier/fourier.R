
# Correlation between new C0 and single nucleotides

# Throughout this document, we are only using the information of whether a nucleotide
# is A or T vs. C or G at the $i$th position.

# If the nucleotide is A or T at the $i$th position, then the feature $x_i$ is 1. Otherwise, it is 0.

#```{r load_package, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
library(tidyverse)
library(glmnet)
library(gtools)
#```

#```{r load_data, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
dat = readRDS("../DNA_cyclizability/data/processed_tiling_newC0.rds")
dat_test = readRDS("../DNA_cyclizability/data/processed_tiling_test_newC0.rds")
#```

#```{r get_mono_nucleotide_AT, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
ps1 = paste0("X", 1:50, "mono")
Xone = dat %>% select(all_of(ps1))
y_all_int = dat$C0_new
Y = y_all_int
# The histogram of $Y$ seems to be normally distributed. 
N1 = Y > quantile(Y)[4]
N2 = Y < quantile(Y)[2]
#```

# The Fourier coefficients are defined using:
# \[
# a_n = \sum_{k = 1}^N \cos (2\pi n(k-1)/N) x_k, 
# \]
# \[
# b_n = \sum_{k = 1}^N \sin (2\pi n(k-1)/N) x_k.
# \]

# ```{r construct_fourier_series, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
# A/T vs. C/G
Xone_AT = matrix(nrow=nrow(Xone), ncol=ncol(Xone))
colnames(Xone_AT) = colnames(Xone)
Xone_AT[] = ((Xone == "A") | (Xone == "T")) %>% as.matrix() %>% as.numeric()
N = 50
# rows: n (frequency), cols: k (location)
# period: N/n
cos_matrix = matrix(nrow = N, ncol = N)
sin_matrix = matrix(nrow = N, ncol = N)
for (n in 1:N) {
  for (k in 1:N) {
    cos_matrix[n, k] = cos(2*pi*n*(k-1)/N)
    sin_matrix[n, k] = sin(2*pi*n*(k-1)/N)
  }
}
Xone_AT_fourier_cos = Xone_AT %*% t(cos_matrix)
Xone_AT_fourier_sin = Xone_AT %*% t(sin_matrix)
Y_binary = rep(NA, length(Y))
Y_binary[N1] = 1
Y_binary[N2] = 0
# cos_matrix is symmetric around 25; sin_matrix is symmetric around 25 by a negative sign
opar = par(mfcol=c(2,1))
plot(cos_matrix[5,])
plot(sin_matrix[5,])
par(opar)
# ```

# ```{r correlations, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
# Check the correlation between features and the outcome
plot(drop(cor(Y, Xone_AT)), main="correlation between nucleotide A/T vs. C/G and new C0")
# We are particularly looking at 5th fourier coefficients, hence n = 5
cor_fourier_without_abs = function(theta, n = 5) {
  cor(Y, (Xone_AT %*% matrix(cos(2*pi*n*(0:(N-1))/N + theta), ncol = 1)))
}
cor_fourier_with_abs = function(theta, n = 5) {
  cor(Y, abs(Xone_AT %*% matrix(cos(2*pi*n*(0:(N-1))/N + theta), ncol = 1)))
}
thetas_len = 100
thetas = seq(0, pi, length.out=thetas_len)
result_without_abs = result_with_abs = rep(0, thetas_len)
for (j in 1:thetas_len) {
  result_with_abs[j] = cor_fourier_with_abs(thetas[j])
  result_without_abs[j] = cor_fourier_without_abs(thetas[j])
}
plot(result_with_abs, main="correlation betw. Fourier coefficients with different phases and new C0, w/ abs")
plot(result_without_abs, main="correlation betw. Fourier coefficients with different phases and new C0, w/o abs")
# cor(Y, rowSums(Xone_AT))
# ```


# 2. Regression models

# Adding the 5th Fourier cofficients (cosine and sine) offers substantial improvement to predicting new C0 from single nucleotides using a main effects model.

# However, adding the 5th Fourier cofficients (cosine and sine) offers negligible improvement to predicting new C0 from single nucleotides when we add the interaction terms.

# ```{r glmnet_fourier}
# Sparse logistic regression using glmnet
# Sparse logistic regression: Fourier
nonNAs = !is.na(Y_binary)
Xone_AT_fourier_logistic_glmnet = cv.glmnet(
  cbind(Xone_AT_fourier_cos[,1:25], Xone_AT_fourier_sin[,1:24])[nonNAs, ],
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_fourier_logistic_glmnet)
Xone_AT_fourier_logistic_glmnet = glmnet(
  cbind(Xone_AT_fourier_cos[,1:25], Xone_AT_fourier_sin[,1:24])[nonNAs, ],
  Y[nonNAs],
  family = "gaussian",
  lambda = 0.00326
  )
print(Xone_AT_fourier_logistic_glmnet)
coef(Xone_AT_fourier_logistic_glmnet)
fourier_logistic_glmnet_predict = 
  predict(Xone_AT_fourier_logistic_glmnet,
          newx = cbind(Xone_AT_fourier_cos[,1:25], Xone_AT_fourier_sin[,1:24])[nonNAs, ]) %>%
  as.numeric()
mean((fourier_logistic_glmnet_predict - Y[nonNAs])^2)
# ```

# ```{r glmnet_nucleotide}
# Sparse logistic regression: single nucleotide
Xone_AT_logistic_glmnet = cv.glmnet(
  Xone_AT[nonNAs,],
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet)
Xone_AT_logistic_glmnet = glmnet(
  Xone_AT[nonNAs,],
  Y[nonNAs],
  lambda = 0.00124,
  family = "gaussian")
print(Xone_AT_logistic_glmnet)
coef(Xone_AT_logistic_glmnet)
logistic_glmnet_predict = predict(Xone_AT_logistic_glmnet, type = "response", newx = Xone_AT[nonNAs,]) %>%
  as.numeric()
mean((logistic_glmnet_predict - Y[nonNAs])^2)
# View(data.frame(logistic_glmnet_predict, fourier_logistic_glmnet_predict, Y[nonNAs]))
# ```

# We start from a base model with only single nucleotides.

# ```{r improvements}
Xone_AT_logistic_glmnet_none = cv.glmnet(
  Xone_AT[nonNAs,],
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_none)
# ```

# ```{r improvement1}
Xone_AT_logistic_glmnet_freq5 = cv.glmnet(
  cbind(Xone_AT[nonNAs,], Xone_AT_fourier_cos[nonNAs,5], Xone_AT_fourier_sin[nonNAs,5]),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq5)
# ```

# ```{r improvement2}
Xone_AT_logistic_glmnet_freq5_squared = cv.glmnet(
  cbind(Xone_AT[nonNAs,], Xone_AT_fourier_cos[nonNAs,5]^2, Xone_AT_fourier_sin[nonNAs,5]^2),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq5_squared)
# ```

# ```{r improvement3}
Xone_AT_logistic_glmnet_freq5_abs = cv.glmnet(
  cbind(Xone_AT[nonNAs,], abs(Xone_AT_fourier_sin[nonNAs,5]), abs(Xone_AT_fourier_cos[nonNAs,5])),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq5_abs)
# ```

# ```{r improvements4}
Xone_AT_logistic_glmnet_freq_all = cv.glmnet(
  cbind(Xone_AT[nonNAs,], Xone_AT_fourier_cos[nonNAs, ]^2, Xone_AT_fourier_sin[nonNAs, ]^2),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq_all)
# ```

# ```{r improvements5}
Xone_AT_logistic_glmnet_freq_fourieronly = cv.glmnet(
  cbind(abs(Xone_AT_fourier_cos[nonNAs, ]), abs(Xone_AT_fourier_sin[nonNAs, ])),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq_fourieronly)

# ```

# ```{r improvements6}
Xone_AT_logistic_glmnet_freq_fourier5only = cv.glmnet(
  cbind(abs(Xone_AT_fourier_cos[nonNAs, 5]), abs(Xone_AT_fourier_sin[nonNAs, 5])),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_freq_fourier5only)

Y_nonNAs_predict = predict(Xone_AT_logistic_glmnet_freq_fourier5only,   cbind(abs(Xone_AT_fourier_cos[nonNAs, 5]), abs(Xone_AT_fourier_sin[nonNAs, 5])))

cor(Y_nonNAs_predict, Y[nonNAs])
# ```


# In this section, we are looking about whether adding the Fourier coefficients
# can improve the prediction when we have included the second-order interactions
# already. The conclusion is negative:

# ```{r interactions}
Xone_AT_interactions = matrix(nrow = nrow(Xone_AT), ncol = ncol(Xone_AT)*(ncol(Xone_AT)+1)*0.5)
k = 0
for (i in 1:N) {
  for (j in 1:i) {
    k = k + 1
    Xone_AT_interactions[, k] = Xone_AT[, i] * Xone_AT[, j]
  }
}
Xone_AT_logistic_glmnet_interactions = cv.glmnet(
  cbind(Xone_AT[nonNAs,], Xone_AT_interactions[nonNAs, ]),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_interactions)

Xone_AT_logistic_glmnet_interactions_5 = cv.glmnet(
  cbind(Xone_AT[nonNAs,], Xone_AT_interactions[nonNAs, ], abs(Xone_AT_fourier_cos[nonNAs, 5]), abs(Xone_AT_fourier_sin[nonNAs ,5])),
  Y[nonNAs],
  family = "gaussian")
print(Xone_AT_logistic_glmnet_interactions_5)
# ```

# ```{r do_not_use, echo=FALSE, results='hide', message=FALSE, warning=FALSE}
cor(Xone_AT_fourier_cos[,5]^2, Y)
cor(abs(Xone_AT_fourier_cos[,5]), Y)
plot(abs(Xone_AT_fourier_cos[,5]), Y, cex = 0.5, pty = 19, col = rgb(0,0,0,0.1))
plot(abs(Xone_AT_fourier_sin[,5]), Y)
abline(h=0, col="red")
cor(Xone_AT_fourier_sin[nonNAs,5]^2, Y[nonNAs])
cor(sign(Xone_AT_fourier_sin[nonNAs,5]), Y[nonNAs])
# ```

# ```{r maxcorr_do_not_use}
eigens = eigen(cov(Xone_AT))
U = eigens$vectors %*% diag(1/sqrt(eigens$values)) %*% t(eigens$vectors) %*% t(cov(Y, Xone_AT))
sqrt(sum(U * U)) / sqrt(cov(matrix(Y,ncol=1)))
# ```

# ```{r orthogonal_do_not_use}
Z = rbind(rep(1, 50), cos_matrix[1:25,], sin_matrix[1:24,])
W = Z %*% t(Z)
which(abs(W) >= 0.0001, arr.ind = TRUE)
# ```


# ![Bendability](https://3k8pb633fck715pnk94793lq-wpengine.netdna-ssl.com/wp-content/uploads/Red-bent-spring-spiral-leader.jpg)


# # References
# 
# 1. Marini, J. C., LEVENEt, S. D., CROTHERSt, D. M. & Englund, P. T. Bent helical structure in kinetoplast DNA. Proc. Natl. Acad. Sci. USA 5 (1982).
# 
# 2. Lipfert, J. et al. Double-stranded RNA under force and torque: Similarities to and striking differences from double-stranded DNA. Proc. Natl. Acad. Sci. U.S.A. 111, 15408–15413 (2014)
# 
# 3. Zhou, H. et al. m1A and m1G disrupt A-RNA structure through the intrinsic instability of Hoogsteen base pairs. Nat Struct Mol Biol 23, 803–810 (2016).
# 
# 4. Coll, M., Frederick, C. A., Wang, A. H. & Rich, A. A bifurcated hydrogen-bonded conformation in the d(A.T) base pairs of the DNA dodecamer d(CGCAAATTTGCG) and its complex with distamycin. Proceedings of the National Academy of Sciences 84, 8385–8389 (1987).
# 
# 5. Zhang, Y., Basu, A., Ha, T. & Bialek, W. Searching for sequence features that control DNA flexibility. Preprint at http://arxiv.org/abs/2012.06127 (2020).
# 
# 6. Khan, S. R., Sakib, S., Rahman, M. S. & Samee, M. A. H. DeepBend: An Interpretable Model of DNA Bendability. 2022.07.06.499067 Preprint at https://doi.org/10.1101/2022.07.06.499067 (2022).
# 
# 7. Trifonov, E. N. & Nibhani, R. Review fifteen years of search for strong nucleosomes: Search for Strong Nucleosomes. Biopolymers 103, 432–437 (2015).
# 
# 8. Brogaard, K., Xi, L., Wang, J.-P. & Widom, J. A map of nucleosome positions in yeast at base-pair resolution. Nature 486, 496–501 (2012).
# 
# 9. Jin, H., Rube, H. T. & Song, J. S. Categorical spectral analysis of periodicity in nucleosomal DNA. Nucleic Acids Res 44, 2047–2057 (2016).
# 
# 10. Basu, A. et al. Measuring DNA mechanics on the genome scale. Nature 589, 462–467 (2021).
# 
# 11. Tang, L. Sequencing DNA bendability. Nat Methods 18, 121–121 (2021).
# 
# 12. Li, K., Carroll, M., Vafabakhsh, R., Wang, X. A. & Wang, J.-P. DNAcycP: a deep learning tool for DNA cyclizability prediction. Nucleic Acids Research 50, 3142–3154 (2022).
# 
# 
# 
