library(SynthHousehold)

#set up data set
data(household)
orig <- initData(household)

#mcmc parameters
mc <- list(nrun = 25, burn = 20, thin = 1)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin

#hyper parameters
#aa ab: gamma hyperparameters for alpha
hyper <- list(K=40 , L=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25, dHH = c(2,3), howmany = 10000)

para <- initParameters(orig,hyper)
output <- initOutput(orig,hyper,mc)
