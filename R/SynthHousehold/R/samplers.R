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
