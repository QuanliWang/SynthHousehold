### Run model
proc_t <- proc.time()
#mc <- list(nrun = 200, burn = 50, thin = 5)
#Rprof()
#library(profvis)
#profvis({
ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,
                         ExampleData$individual_variable_index,
                         ExampleData$household_variable_index,
                         options$HHhead_at_group_level,options$weight_option,struc_weight,MissData, Parallel = TRUE)
#})

#Rprof(NULL)
total_time <- (proc.time() - proc_t)[["elapsed"]]
total_time



#library(profvis)
#profvis({
#ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,
#                         ExampleData$individual_variable_index,
#                         ExampleData$household_variable_index,
#                         options$HHhead_at_group_level,options$weight_option,struc_weight,MissData, Parallel = TRUE)
#})

#library(rbenchmark)
#benchmark("serial" = {ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,
#                                               ExampleData$individual_variable_index,
#                                               ExampleData$household_variable_index,
#                                               options$HHhead_at_group_level,options$weight_option,struc_weight,MissData, Parallel = FALSE)},
#          "parallel" = {ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,
#                                                 ExampleData$individual_variable_index,
#                                                 ExampleData$household_variable_index,
#                                                 options$HHhead_at_group_level,options$weight_option,struc_weight,MissData, Parallel = TRUE)},
#          replications = 2,
#          columns = c("test", "replications", "elapsed",
#                      "relative", "user.self", "sys.self"))
#summaryRprof()

