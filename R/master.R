

rm(list = ls())
library(NestedCategBayesImpute)
library(dplyr)

format2 = TRUE #set format2 to TRUE for the new format

if (format2) {
  orig.file <- system.file("extdata", "origdata_newFormat.txt", package = "NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  names(orig.data) <- c("Hhindex","pernum", "sex", "race", "sthn","age","relate","ownership",
                        "headsex", "headrace", "headsthn","headage")
  household.size <- as.data.frame(table(orig.data[,1]))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)

  individual_varible_index = c(3:7)
  household_variable_index = c(8:13) #column 8 to 13 are household level data #make sure the last one is household size
} else {
  #set up data set
  #data(household)
  orig.file <- system.file("extdata", "origdata_oldFormat.txt", package = "NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  names(orig.data) <- c("Hhindex","pernum", "sex", "race", "sthn","age","relate","ownership")
  household.size <- as.data.frame(table(orig.data[,1]))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)
  
  individual_varible_index = c(3:7)
  household_variable_index = c(8,9) #make sure the last one is household size
}

orig <- initData(household,individual_varible_index,household_variable_index)

#mcmc parameters
mc <- list(nrun = 10000, burn = 5000, thin = 50)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin

#hyper parameters
#aa ab: gamma hyperparameters for alpha
dHH <- rep(0,length(household_variable_index))
for (i in 1:length(dHH)) {
  dHH[i] <- max(household[,household_variable_index[i]])
}

if (!format2) {
  dHH[length(dHH)] <- dHH[length(dHH)] - 1 # the old dat format has household size starting from 2
}

hyper <- list(K=40 , L=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25,
              dHH = dHH, blocksize = 10000)
para <- initParameters(orig,hyper,format2)
output <- initOutput(orig,hyper,mc)

#synindex <- c(9910,9920,9930,9940,9950,9960,9970,9980,9990,10000)
synindex <- sort(sample(seq((mc$nrun*0.8 +1),mc$nrun,mc$eff.sam),10,replace=F))
synData <- list()

#weighting
weight_option <- TRUE #set to true for weighting/capping option
struc_weight <- c(1,1,1/2,1/3,1/3) #set weights: must be ordered & no household size must be excluded

if(weight_option){
  struc_weight <- c(1,struc_weight) #add 1 for the weight of the observed data
  weighted_z_HH_all <- vector("list",length(struc_weight))
  weighted_HHdata_all <- vector("list",length(struc_weight))
  weighted_IndividualData_all <- vector("list",length(struc_weight))
  weighted_z_Individual_all <- vector("list",length(struc_weight))
}

#source("../mcmc.R")
source("mcmc.R")


#save synthetic data
writeFun <- function(LL){
  for(i in 1:length(LL)){
    if(format2){
      if(weight_option){
        write.table(t(LL[[i]]),paste0("synData_newFormat_weighted",i,".txt"),row.names = F,col.names = F)
      } else {
        write.table(t(LL[[i]]),paste0("synData_newFormat",i,".txt"),row.names = F,col.names = F)
      }
    } else {
      if(weight_option){
        write.table(t(LL[[i]]),paste0("synData_oldFormat_weighted",i,".txt"),row.names = F,col.names = F)
      } else {
        write.table(t(LL[[i]]),paste0("synData_oldFormat",i,".txt"),row.names = F,col.names = F)
      }
    }
  }
}
writeFun(synData)




