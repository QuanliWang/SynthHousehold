library(SynthHousehold)

#set up data set
data(household)
orig <- initData(household)

#mcmc parameters
mc <- list(nrun = 25, burn = 20, thin = 1)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin

#hyper parameters
#aa ab: gamma hyperparameters for alpha
hyper <- list(K=40 , L=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25, dHH = c(2,3), blocksize = 10000)

para <- initParameters(orig,hyper)
output <- initOutput(orig,hyper,mc)

synindex <- c(9910,9920,9930,9940,9950,9960,9970,9980,9990,10000)
synData <- list()

for (i in 1:mc$nrun) {
  # update zHH
  z_household <- samplezHH(para$phi,orig$dataT,para$w,para$pi,orig$SS,t(para$HHdata_all[,1:orig$n]),
            para$lambda[[1]],para$lambda[[2]])

  # update zIndividual
  z_Individuals <- samplezmember(para$phi,orig$dataT,para$w,z_household$z_HH,orig$HHserial)
  data.extra <- GetImpossibleHouseholds(orig$d,orig$ACS_count,para$lambda,para$w,para$phi,
                  para$pi,hyper$blocksize,orig$n,is.element(i,synindex))
  if (is.element(i,synindex)) {
    synData[[which(synindex ==i)]] <- data.extra$synIndividuals_all[1:8,]
  }
    #combine data and indicators
    para$z_HH_all <- c(z_household$z_HH, data.extra$z_HH_extra)
    para$HHdata_all <- orig$HHdataorigT
    para$HHdata_all[2,] <- para$HHdata_all[2,] -1
    para$HHdata_all <- rbind(z_household$z_HH_Individuals,z_Individuals)
    para$IndividualData_all <- cbind(t(orig$origdata[,1:8]),data.extra$IndividualData_extra)

    #row 1 for K groups and row 2 for L groups
    temp <- rbind(z_household$z_HH_Individuals,z_Individuals)
    para$z_Individual_all  <- cbind(temp,data.extra$HHdata_individual_extra)

    # update phi
    para$phi <- UpdatePhi(para$IndividualData_all,para$z_Individual_all, hyper$K,hyper$L,orig$p,orig$d,orig$maxd);

    W <- UpdateW(para$beta,para$z_Individual_all, hyper$K, hyper$L)
    para$w <- W$w
    para$v <- W$v


  #postsave
}
