initData <- function(household) {
  orig <- list()
  orig$origdata <- household
  orig$SS <- as.data.frame(table(household[,'Hhindex']))$Freq
  orig$n <- length(orig$SS)

  HHrowIndex <- c(1, cumsum(orig[["SS"]])+1)
  #column 8 and 9 are household level data
  orig$HHdataorigT <- t(household[HHrowIndex[1:orig$n],c(8,9)])
  orig$HHserial <- household[,1]

  orig$n_individuals <- dim(household)[1]
  #hard-coded data structure for now, as we don't know if it will apply to more general data set yet
  orig$p <- dim(household)[2] - 4

  #levels for each variable in the model. This might be different from the sun-sampled data in use
  orig$d <- c(2,9,5,94,12)
  orig$dataT <- t(household[,3:7])
  orig$maxd <- max(orig$d)

  counts <- as.data.frame(table(orig$SS))
  orig$ACS_count <- counts[order(counts[,1]),2]
  return(orig)
}


initParameters <- function(data,hyper) {
  para <- list()
  para$alpha <- 1 #hyperparameters for stick-breaking weights
  para$beta <- 1

  #intilize phi
  para$phi <- matrix(0, nrow  = data$maxd*data$p, ncol = hyper$K*hyper$L) #cell probabilities
  phi_1 <- matrix(0, nrow = data$maxd, ncol = data$p)
  for (i in 1:data$p) {
    for (j in 1:data$d[i]) {
      phi_1[j,i] <-sum(data$dataT[i,]==j)/data$n_individuals; #should we optimize this later?
    }
  }
  phi_1 <- as.vector(phi_1)
  for (i in 1:dim(para$phi)[2]) {
    para$phi[,i] <- phi_1
  }

  para$HHdata_all <- data$HHdataorigT;
  para$HHdata_all[2,] <- para$HHdata_all[2,] - 1

  #initialize lambda
  para$lambda <- list()
  for (i in 1:length(hyper$dHH)) {
    lambda <- matrix(0, nrow = hyper$K,ncol = hyper$dHH[i])
    for (j in 1:hyper$dHH[i]) {
      lambda[,j] <- sum(para$HHdata_all[i,]==j) / data$n
    }
    para$lambda[[i]] <- lambda
  }

  para$u <- c(rbeta(hyper$K-1, 1,para$alpha),1)
  para$pi <- para$u * cumprod(c(1,1.0-para$u[1:hyper$K-1]))

  ones <- matrix(1.0, hyper$K,1)
  para$v <- c(rbeta(hyper$K * (hyper$L-1), 1, para$beta),ones)
  dim(para$v) <- c(hyper$K, hyper$L)

  para$w <- matrix(0, nrow = hyper$K,ncol = hyper$L)
  for (i in 1:hyper$K) {
    v1 <- para$v[i,]
    para$w[i,]  <- v1 * cumprod(c(1, 1-v1[1:hyper$L-1]))
  }

  return(para)
}

initOutput <- function(data,hyper,mc) {
  output <- list()
  output$alphaout <- matrix(0,nrow = mc$nrun,ncol = 1)
  output$betaout <- matrix(0, nrow = mc$nrun,ncol = 1)
  output$piout <- matrix(0, mc$eff.sam,hyper$K)
  output$wout <- array(0, dim=c(mc$eff.sam,hyper$K,hyper$L))
  output$nout <- matrix(0,nrow = mc$nrun,ncol = 1)
  output$extrasize <- matrix(0,nrow = mc$nrun,ncol = 3)
  output$z_HH_save <- matrix(0, nrow = mc$nrun, ncol = data$n_individuals)
  output$z_member_save <- matrix(0, nrow = mc$nrun, ncol = data$n_individuals)
  output$elapsed_time <-  matrix(0,nrow = mc$nrun,ncol = 1)
  output$newphiout <- array(0, dim=c(mc$eff.sam,data$maxd*data$p,hyper$K*hyper$L))

  output$lambda1out = array(0, dim=c(mc$eff.sam,hyper$K, hyper$dHH[1]))
  output$lambda2out = array(0, dim=c(mc$eff.sam,hyper$K, hyper$dHH[2]))
  return(output)
}
