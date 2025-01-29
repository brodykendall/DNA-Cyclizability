
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


# find_c0Aphi2 <- function(dat, n=c(26,29,31), aa=c(1,1,1)) {
#   # Cn = C0 + C0*A*sin(bn+phi)
#   mat <- matrix(0, nrow = 3, ncol = 3)
#   k <- 2*pi/10.4
#   mat[1:3,1] <- 1
#   mat[1:3,2] <- sin(n*k)
#   mat[1:3,3] <- cos(n*k)
#   mat[1,2:3] <- mat[1,2:3]*aa[1]
#   mat[2,2:3] <- mat[2,2:3]*aa[2]
#   mat[3,2:3] <- mat[3,2:3]*aa[3]
#   inv_mat <- solve(mat)
#   c0A1A2 <- as.matrix(dat[,1:3]) %*% t(inv_mat)
#   c0Aphi <- c0A1A2
#   c0Aphi[,1] <- c0A1A2[,1]
#   c0Aphi[,2] <- sqrt((c0A1A2[,2]/c0A1A2[,1])^2 + (c0A1A2[,3]/c0A1A2[,1])^2)
#   c0Aphi[,3] <- sign((c0A1A2[,3]/c0A1A2[,1])) * acos((c0A1A2[,2]/c0A1A2[,1])/c0Aphi[,2])
#   return(c0Aphi)
# }
