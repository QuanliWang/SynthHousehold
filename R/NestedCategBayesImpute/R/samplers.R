UpdatePhi <- function(IndividualData_all, z_Individual_all, K, L, p, d, maxd) {
  phi <- array(NA,dim = c(maxd,p, K*L))
  data = IndividualData_all[3:7,]
  groupIndex <- L*(z_Individual_all[1,]-1)+z_Individual_all[2,]
  for (j in 1:p) {
    phicount <- groupcount(groupIndex, data[j,], K*L, d[j])
    phi_j <- apply(phicount, c(1,2), function(x) rgamma(1,x+1,1))
    phi[1:d[j], j,] <- apply(phi_j, 1, function(x) x / sum(x))
  }
  dim(phi) <- c(maxd*p,K * L) #reshape to a 2D matrix
  return(phi)
}

UpdateW <- function(beta,z_Individual_all, K, L) {
  phicountcluster <- groupcount(z_Individual_all[1,],z_Individual_all[2,],K,L)

  cum <- t(apply(phicountcluster[,seq(L,1)],1, cumsum))[,seq(L,1)]
  v <- mapply(function(x,y) rbeta(1,x,y), 1 + phicountcluster[,1:(L-1)], beta + cum[,2:L])
  dim(v) <- c(K, L-1)
  v[v>1-1e-5] <- 1-1e-5
  v = cbind(v,1)
  w <- v * t(apply(cbind(1, 1-v[,1:(L-1)]), 1, cumprod))

  return(list(w = w, v = v))
}

UpdateLambda <- function(dHH,K,z_HH_all,HHdata_all) {
  lambda <- list()
  for (i in 1:length(dHH)) {
    lambdacount <- groupcount(z_HH_all,HHdata_all[i,],K, dHH[i])
    lam <- apply(lambdacount, c(1,2), function(x) rgamma(1,x+1,1))
    lam <- t(apply(lam, 1, function(x) x / sum(x)))
    lambda[[i]] = lam;
  }

  return(lambda)
}

UpdatePi <- function(alpha,z_HH_all,K) {

  kcount <- groupcount1D(z_HH_all, K)
  s <- seq(K,1)
  cum <- cumsum(kcount[s])[s]
  u <- mapply(function(x,y) rbeta(1,x,y), 1 + kcount[1:K-1], alpha + cum[2:K])
  u[u > 1-1e-5] <- 1-1e-5
  u <- c(u,1)
  u[K] <- 1

  pi  <- u* cumprod(c(1,1-u[1:K-1]))

  return(list(pi = pi, u = u))
}

UpdateAlpha <- function(aa,ab,u) {
  K <- length(u)
  alpha <- rgamma(1,aa + K - 1,scale = 1/(ab - sum(log(1-u[1:K-1]))))
  return(alpha)
}

UpdateBeta <- function(ba,bb,v) {
  K <- dim(v)[1]
  L <- dim(v)[2]
  beta <- rgamma(1,ba + K*(L-1), scale = 1/(bb - sum(log(1-v[,1:L-1])) ))
}



