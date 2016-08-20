rm(list = ls())
library(NestedCategBayesImpute)
library(dplyr)

orig.file <- system.file("extdata", "origdata_newFormat.txt", package = "NestedCategBayesImpute")


orig.data <- read.table(orig.file,header = TRUE, sep = " ")
names(orig.data) <- c("Hhindex","pernum", "sex", "race", "sthn","age","relate","ownership",
                      "headsex", "headrace", "headsthn","headage")
household.size <- as.data.frame(table(orig.data[,1]))
household.size[,1] <- as.numeric(household.size[,1])
names(household.size) <- c("Hhindex", 'householdsize')
household <- orig.data %>% inner_join(household.size)

individual_varible_index = c(3:7)

#column 8 to 13 are household level data
household_variable_index = c(8:13) #make sure the last one is household size

orig <- initData(household, individual_varible_index,household_variable_index)

#mcmc parameters
mc <- list(nrun = 10000, burn = 5000, thin = 50)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin
#hyper parameters
#aa ab: gamma hyperparameters for alpha

dHH <- rep(0,length(household_variable_index))
for (i in 1:length(dHH)) {
  dHH[i] <- max(household[,household_variable_index[i]])
}
hyper <- list(K=40 , L=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25,
              dHH = dHH, blocksize = 10000)
para <- initParameters(orig,hyper)
output <- initOutput(orig,hyper,mc)

synindex <- c(100,9920,9930,9940,9950,9960,9970,9980,9990,10000)
synData <- list()

#Rprof("R1.out",interval = 0.001)
for (i in 1:mc$nrun) {
  cat(paste("iteration ", i,"\n", sep = ""))
  t <- proc.time()
  # update zHH
  print("update zHH")

  stopifnot(!any(is.na(para$phi)))
  stopifnot(sum(para$phi) == 3000)
  stopifnot(length(para$w) == 600)
  stopifnot(sum(para$w) == 40)
  stopifnot(!any(is.na(para$w)))
  stopifnot(!any(is.na(para$pi)))
  stopifnot(length(unlist(para$lambda)) == 4680)
  z_household <- samplezHH(para$phi,orig$dataT,para$w,para$pi,orig$SS,t(para$HHdata_all[,1:orig$n]),para$lambda)

  stopifnot(max(z_household$z_HH)<=hyper$K && min(z_household$z_HH)>=1)
  # update zIndividual
  print("update zIndividual")
  z_Individuals <- samplezmember(para$phi,orig$dataT,para$w,z_household$z_HH,orig$HHserial)
  print("get new household")
  #save.image("debug.RData")
  data.extra <- GetImpossibleHouseholds(orig$d,orig$ACS_count,para$lambda,para$w,para$phi,
                  para$pi,hyper$blocksize,orig$n,is.element(i,synindex))

  para$hh_size_new <- as.vector(data.extra$hh_size_new)
  DIM <- dim(data.extra$IndividualData_extra)[1]
  if (is.element(i,synindex)) {
    synData[[which(synindex ==i)]] <- data.extra$synIndividuals_all[1:DIM,]
  }
    #combine data and indicators
   print("combine")
    para$z_HH_all <- c(z_household$z_HH, data.extra$z_HH_extra)
    para$HHdata_all <- orig$HHdataorigT
    para$HHdata_all <- cbind(para$HHdata_all,data.extra$HHdata_extra)
    para$IndividualData_all <- cbind(t(orig$origdata[,1:DIM]),data.extra$IndividualData_extra)

    #row 1 for K groups and row 2 for L groups
    temp <- rbind(z_household$z_HH_Individuals,z_Individuals)
    para$z_Individual_all  <- cbind(temp,data.extra$z_HHdata_individual_extra)

    # update phi
    print("update phi")
    para$phi <- UpdatePhi(para$IndividualData_all,para$z_Individual_all,
                    hyper$K,hyper$L,orig$p,orig$d,orig$maxd,individual_varible_index)

    #update W
    print("update W")
    stopifnot(all(para$z_Individual_all[1,]>0))
    stopifnot(all(para$z_Individual_all[1,]<=hyper$K))
    stopifnot(all(para$z_Individual_all[2,]>0))
    stopifnot(all(para$z_Individual_all[2,]<=hyper$L))
    W <- UpdateW(para$beta,para$z_Individual_all, hyper$K, hyper$L)
    para$w <- W$w
    para$v <- W$v
    stopifnot(length(para$w) == 600)
    stopifnot(sum(para$w) == 40)


  # update lambda
    print("update lambda")
  para$lambda <- UpdateLambda(hyper$dHH,hyper$K,para$z_HH_all,para$HHdata_all)

  # update pi
  print("update pi")
  Pi <- UpdatePi(para$alpha,para$z_HH_all,hyper$K)
  para$pi <- Pi$pi
  para$u <- Pi$u

  #update alpha
  print("update alpha")
  para$alpha <- UpdateAlpha(hyper$aa,hyper$ab,para$u)

  #update beta
  print("update beta")
  para$beta <- UpdateBeta(hyper$ba,hyper$bb,para$v)

  print("save")
  #post save
  if (i %% mc$thin == 0 && i > mc$burn)  {
    index <- (i-mc$burn)/mc$thin
    output$piout[index,] <- para$pi
    output$wout[index,,] <- para$w
    output$newphiout[index,,] <- para$phi

    output$lambda1out[index,,] <- para$lambda[[1]]
    output$lambda2out[index,,] <- para$lambda[[2]]
  }

  total_household <- dim(para$HHdata_all)[2]
  cat(paste("number of household is ", total_household, "\n", sep = ''))

  output$elapsed_time[i] <- (proc.time() - t)[["elapsed"]]
  cat(paste("elapsed time = ", output$elapsed_time[i], "\n\n", sep = ' '))
  output$nout[i] <- total_household
  output$extrasize[i,] <- para$hh_size_new
  #output$z_HH_save[i,] <- para$z_Individual_all[1,1:orig$n_individuals]
  #output$z_member_save[i,]  <- para$z_Individual_all[2,1:orig$n_individuals]
  output$alphaout[i] <- para$alpha
  output$betaout[i] <- para$beta

}
#Rprof(NULL)
#summaryRprof("R1.out")
