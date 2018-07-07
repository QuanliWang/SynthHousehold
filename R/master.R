### Empty environment and load required libraries
rm(list = ls())
library(NestedCategBayesImpute)
library(dplyr)

### Set indicator for whether of not to move the household head
### Also set indicator for the weighting/capping option
options <- list()
options$HHhead_at_group_level <- TRUE #set to TRUE to move household head to the group level
options$weight_option <- FALSE #set to TRUE for weighting/capping option. If TRUE, must supply weights

source("../GetExampleData.R")

ExampleData <- PrepareData(options)

### Initialize and set parameters for missing data
MissData <- initMissing(ExampleData,
                        struc_zero_variables=c("sex","age","relate","headsex","headage"),
                        miss_batch=10)
MissData$prop_batch <- 1.2


### Initialize the input data structure
orig <- initData(MissData)


### Supply weights; one for each household size
if(options$weight_option){
  struc_weight <- c(1/2,1/2,1/3,1/3,1/3) #must be ordered & no household size must be excluded
} else {
  struc_weight <- rep(1,length(orig$n_star_h)) #just a dummy column of ones if struc_weight=FALSE
}


### Set mcmc parameters

#mc <- list(nrun = 200, burn = 50, thin = 5)
mc <- list(nrun = 1000, burn = 500, thin = 5)
#mc <- list(nrun = 10000, burn = 5000, thin = 5)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin

### Set number of categories for each household level variable
dHH <- rep(0,length(ExampleData$household_variable_index))
for (i in 1:length(dHH)) {
  dHH[i] <- max(na.omit(ExampleData$household[,ExampleData$household_variable_index[i]]))
  if (i == length(dHH) & !options$HHhead_at_group_level) {
    dHH[length(dHH)] <- dHH[length(dHH)] - 1 #Household head within household assumes that the household size starts from 2
  }
}

### Set hyper parameters
#aa & ab are gamma hyperparameters for alpha while ba & bb are gamma hyperparameters for beta
#blocksize is the number of impossible households to sample at once (we use batch sampling to speed up mcmc)
#FF is the max number of group-level latent classes
#SS is the max number of individual-level classes
hyper <- list(FF=20 , SS=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25,dHH = dHH, blocksize = 10000)


### Initialize parameters and output
para <- initParameters(orig,hyper,options$HHhead_at_group_level)
output <- initOutput(orig,hyper,mc)


### Set number of synthetic data and the mcmc indexes for them
mm <- 50
synindex <- NULL
MissData$miss_index <- sort(sample(seq((mc$burn +1),mc$nrun,by=mc$thin),mm,replace=F))
#round(seq((mc$burn +1),mc$burn$nrun,length.out=mm))


### Run model
proc_t <- proc.time()
mc <- list(nrun = 200, burn = 50, thin = 5)
#Rprof()
library(profvis)
profvis({
ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,
                         ExampleData$individual_variable_index,
                         ExampleData$household_variable_index,
                         options$HHhead_at_group_level,options$weight_option,struc_weight,MissData)
})

#Rprof(NULL)
total_time <- (proc.time() - proc_t)[["elapsed"]]









### View first few lines of the first synthetic data.
head((ModelResults$synData)[[1]]) # Remember that the relate variable has been recoded to 11 levels


### Some posterior summaries and plots
library(coda)
names(ModelResults$output)
dim(ModelResults$output$alphaout)
alpha_output <- mcmc(ModelResults$output$alphaout)
plot(alpha_output)
summary(alpha_output)

dim(ModelResults$output$betaout)
beta_output <- mcmc(ModelResults$output$betaout)
plot(beta_output)
summary(beta_output)

dim(ModelResults$output$nout)
total_households <-mcmc(ModelResults$output$nout)
plot(total_households)
summary(total_households)

dim(ModelResults$output$extrasize)
impossible_households <-mcmc(ModelResults$output$extrasize)
plot(impossible_households)
summary(impossible_households)

dim(ModelResults$output$elapsed_time)
time_per_iteration <-mcmc(ModelResults$output$elapsed_time)
plot(time_per_iteration)
summary(time_per_iteration)

dim(ModelResults$output$F_occupied)
F_occupied <-mcmc(ModelResults$output$F_occupied)
plot(F_occupied)
summary(F_occupied)

dim(ModelResults$output$S_occupied_max)
S_occupied_max <-mcmc(ModelResults$output$S_occupied_max)
plot(S_occupied_max)
summary(S_occupied_max)

###################################################################################
####################################### END #######################################
###################################################################################


#save synthetic data
writeFun <- function(LL){
  for(i in 1:length(LL)){
    if(options$HHhead_at_group_level){
      if(options$weight_option){
        write.table(LL[[i]],paste0("synData_newFormat_weighted",i,".txt"),row.names = F,col.names = F)
      } else {
        write.table(LL[[i]],paste0("synData_newFormat",i,".txt"),row.names = F,col.names = F)
      }
    } else {
      if(options$weight_option){
        write.table(LL[[i]],paste0("synData_oldFormat_weighted",i,".txt"),row.names = F,col.names = F)
      } else {
        write.table(LL[[i]],paste0("synData_oldFormat",i,".txt"),row.names = F,col.names = F)
      }
    }
  }
}
writeFun(ModelResults$synData)

#save computational time
if(options$HHhead_at_group_level){
  if(options$weight_option){
    write.table(total_time,"total_time_newFormat_weighted.txt",row.names = F,col.names = F)
  } else {
    write.table(total_time,"total_time_newFormat.txt",row.names = F,col.names = F)
  }
} else {
  if(options$weight_option){
    write.table(total_time,"total_time_oldFormat_weighted.txt",row.names = F,col.names = F)
  } else {
    write.table(total_time,"total_time_oldFormat.txt",row.names = F,col.names = F)
  }
}

