

for (i in 1:mc$nrun) {
  cat(paste("iteration ", i,"\n", sep = ""))
  t <- proc.time()

  z_household <- samplezHH(para$phi,orig$dataT,para$w,para$pi,orig$SS,t(para$HHdata_all[,1:orig$n]),para$lambda)

  z_Individuals <- samplezmember(para$phi,orig$dataT,para$w,z_household$z_HH,orig$HHserial)
  
  if(weight_option){
    data.extra <- GetImpossibleHouseholds(orig$d,ceiling(orig$ACS_count*struc_weight[-1]),para$lambda,para$w,para$phi,
                                          para$pi,hyper$blocksize,orig$n,is.element(i,synindex),format2)
    para$hh_size_new <- as.vector(data.extra$hh_size_new)
    DIM <- dim(data.extra$IndividualData_extra)[1]
    if (is.element(i,synindex)) {
      forsynData <- GetImpossibleHouseholds(orig$d,orig$ACS_count,para$lambda,para$w,para$phi,
                                            para$pi,hyper$blocksize,orig$n,is.element(i,synindex),format2) #synthetic data
      synData[[which(synindex ==i)]] <- forsynData$synIndividuals_all[1:DIM,]
    }
    
    #combine data and indicators -- use lists for weighting
    n_i <- as.data.frame(table(data.extra$IndividualData_extra[1,]))$Freq
    n_i_index <- rep(n_i,n_i)
    weighted_z_HH_all[[1]] <- z_household$z_HH
    temp_temp <- orig$HHdataorigT
    if (!format2) {
      temp_temp[2,] <- temp_temp[2,] -1
    }
    weighted_HHdata_all[[1]] <- temp_temp
    weighted_IndividualData_all[[1]] <- t(orig$origdata[,1:DIM])
    weighted_z_Individual_all[[1]] <- rbind(z_household$z_HH_Individuals,z_Individuals)
    
    for(w_i in 2:length(struc_weight)){
      weighted_z_HH_all[[w_i]] <- data.extra$z_HH_extra[which(n_i == (w_i-1))]
      weighted_HHdata_all[[w_i]] <- data.extra$HHdata_extra[,which(n_i == (w_i-1))]
      weighted_IndividualData_all[[w_i]] <- data.extra$IndividualData_extra[,which(n_i_index == (w_i-1))]
      weighted_z_Individual_all[[w_i]] <- data.extra$z_HHdata_individual_extra[,which(n_i_index == (w_i-1))]
    }
    
    # update phi
    para$phi <- UpdatePhiWeighted(weighted_IndividualData_all,weighted_z_Individual_all,
                                  hyper$K,hyper$L,orig$p,orig$d,orig$maxd,individual_varible_index,struc_weight)
    #update W
    W <- UpdateWWeighted(para$beta,weighted_z_Individual_all, hyper$K, hyper$L,struc_weight)
    para$w <- W$w
    para$v <- W$v
    
    # update lambda
    para$lambda <- UpdateLambdaWeighted(hyper$dHH,hyper$K,weighted_z_HH_all,weighted_HHdata_all,struc_weight)
    
    # update pi
    Pi <- UpdatePiWeighted(para$alpha,weighted_z_HH_all,hyper$K,struc_weight)
    para$pi <- Pi$pi
    para$u <- Pi$u
    
  } else {
    data.extra <- GetImpossibleHouseholds(orig$d,orig$ACS_count,para$lambda,para$w,para$phi,
                                          para$pi,hyper$blocksize,orig$n,is.element(i,synindex),format2)
    para$hh_size_new <- as.vector(data.extra$hh_size_new)
    DIM <- dim(data.extra$IndividualData_extra)[1]
    if (is.element(i,synindex)) {
      synData[[which(synindex ==i)]] <- data.extra$synIndividuals_all[1:DIM,]
    }
    
    #combine data and indicators
    para$z_HH_all <- c(z_household$z_HH, data.extra$z_HH_extra)
    para$HHdata_all <- orig$HHdataorigT
    if (!format2) {
      para$HHdata_all[2,] <- para$HHdata_all[2,] -1
    }
    para$HHdata_all <- cbind(para$HHdata_all,data.extra$HHdata_extra)
    para$IndividualData_all <- cbind(t(orig$origdata[,1:DIM]),data.extra$IndividualData_extra)
    
    #row 1 for K groups and row 2 for L groups
    temp <- rbind(z_household$z_HH_Individuals,z_Individuals)
    para$z_Individual_all  <- cbind(temp,data.extra$z_HHdata_individual_extra)
    
    # update phi
    para$phi <- UpdatePhi(para$IndividualData_all,para$z_Individual_all,
                          hyper$K,hyper$L,orig$p,orig$d,orig$maxd,individual_varible_index)
    
    #update W
    W <- UpdateW(para$beta,para$z_Individual_all, hyper$K, hyper$L)
    para$w <- W$w
    para$v <- W$v
    
    # update lambda
    para$lambda <- UpdateLambda(hyper$dHH,hyper$K,para$z_HH_all,para$HHdata_all)
    
    # update pi
    Pi <- UpdatePi(para$alpha,para$z_HH_all,hyper$K)
    para$pi <- Pi$pi
    para$u <- Pi$u
  }
 
  #update alpha
  para$alpha <- UpdateAlpha(hyper$aa,hyper$ab,para$u)

  #update beta
  para$beta <- UpdateBeta(hyper$ba,hyper$bb,para$v)

  #post save
  if (i %% mc$thin == 0 && i > mc$burn)  {
    index <- (i-mc$burn)/mc$thin
    output$piout[index,] <- para$pi
    output$wout[index,,] <- para$w
    output$newphiout[index,,] <- para$phi

    output$lambdaout[[1]][index,,] <- para$lambda[[1]]
    output$lambdaout[[2]][index,,] <- para$lambda[[2]]
  }

  total_household <- sum(c(orig$n,para$hh_size_new))
  cat(paste("sampled number of households is ", total_household, "\n", sep = ''))
  if(weight_option){
    est_total_household <- sum(c(orig$n,para$hh_size_new)/struc_weight)
    cat(paste("estimated (true) number of households is ", est_total_household, "\n", sep = ''))
  }
  

  output$elapsed_time[i] <- (proc.time() - t)[["elapsed"]]
  cat(paste("elapsed time = ", output$elapsed_time[i], "\n\n", sep = ' '))
  output$nout[i] <- total_household
  output$extrasize[i,] <- para$hh_size_new
  output$alphaout[i] <- para$alpha
  output$betaout[i] <- para$beta
}
