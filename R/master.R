

rm(list = ls())
library(NestedCategBayesImpute)
library(dplyr)

HHhead_at_group_level <- TRUE #set to TRUE to move household head to hthe group level
weight_option <- TRUE #set to true for weighting/capping option


#set up data set
if (HHhead_at_group_level) {
  orig.file <- system.file("extdata", "origdata_newFormat.txt", package = "NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  names(orig.data) <- c("Hhindex","pernum", "sex", "race", "hisp","age","relate","ownership","headsex", "headrace", "headhisp","headage")
  household.size <- as.data.frame(table(orig.data[,1]))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)

  individual_variable_index = c(3:7)
  household_variable_index = c(8:13) #column 8 to 13 are household level data #make sure the last one is household size
} else {
  orig.file <- system.file("extdata", "origdata_oldFormat.txt", package = "NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  names(orig.data) <- c("Hhindex","pernum", "sex", "race", "hisp","age","relate","ownership")
  household.size <- as.data.frame(table(orig.data[,1]))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)

  individual_variable_index = c(3:7)
  household_variable_index = c(8,9) #make sure the last one is household size
}

orig <- initData(household,individual_variable_index,household_variable_index)

#mcmc parameters
mc <- list(nrun = 20000, burn = 10000, thin = 50)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin

#hyper parameters
#aa ab: gamma hyperparameters for alpha
dHH <- rep(0,length(household_variable_index))
for (i in 1:length(dHH)) {
  dHH[i] <- max(household[,household_variable_index[i]])
}

if (!HHhead_at_group_level) {
  dHH[length(dHH)] <- dHH[length(dHH)] - 1 # the old dat format has household size starting from 2
}

hyper <- list(FF=40 , SS=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25,dHH = dHH, blocksize = 10000)
para <- initParameters(orig,hyper,HHhead_at_group_level)
output <- initOutput(orig,hyper,mc)

mm <- 50
synindex <- sort(sample(seq((mc$burn +1),mc$nrun,by=mc$thin),mm,replace=F))
synData <- list()

#weighting
struc_weight <- c(1/2,1/2,1/3,1/3,1/3) #set weights: must be ordered & no household size must be excluded
if(weight_option){
  struc_weight <- c(1,struc_weight) #add 1 for the weight of the observed data
  G_all_weighted <- vector("list",length(struc_weight))
  HHdata_all_weighted <- vector("list",length(struc_weight))
  IndividualData_all_weighted <- vector("list",length(struc_weight))
  M_all_weighted <- vector("list",length(struc_weight))
}

proc_t <- proc.time()
source("../mcmc.R")
#source("mcmc.R")
total_time <- (proc.time() - proc_t)[["elapsed"]]


#save synthetic data
writeFun <- function(LL){
  for(i in 1:length(LL)){
    if(HHhead_at_group_level){
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

#save computational time
if(HHhead_at_group_level){
  if(weight_option){
    write.table(total_time,"total_time_newFormat_weighted.txt",row.names = F,col.names = F)
  } else {
    write.table(total_time,"total_time_newFormat.txt",row.names = F,col.names = F)
  }
} else {
  if(weight_option){
    write.table(total_time,"total_time_oldFormat_weighted.txt",row.names = F,col.names = F)
  } else {
    write.table(total_time,"total_time_oldFormat.txt",row.names = F,col.names = F)
  }
}










