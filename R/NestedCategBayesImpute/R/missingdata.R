SampleMissing <- function(MissData,para,orig,G_household,M,hyper) {
  MissData$n_batch_imp_sum <- MissData$n_batch_imp_sum +
    ceiling(MissData$n_0_reject*MissData$prop_batch)
  MissData$n_batch_imp <- ceiling(MissData$n_batch_imp_sum/i) + 1 #no. of batches of imputations to sample
  MissData$n_0_reject[] <- 0
  MissData$household <- as.matrix(MissData$household)
  #sample non structural zeros variables for everyone at once
  storage.mode(MissData$household) <- "integer" #very important if used to do in place update
  MissData <- SampleMissing_impC(MissData,para,orig,G_household,M,hyper)
  MissData$household <- as.data.frame(MissData$household)
  return(MissData)
}

initMissing <- function(data,struc_zero_variables,miss_batch){
  md <- list()
  #keep a copy of raw data
  md$household_with_miss <- data$household
  md$individual_variable_index <- data$individual_variable_index
  md$household_variable_index <- data$household_variable_index

  #find the index of structure zero varaibles
  struc_zero_variables_index <- which(is.element(colnames(data$household),struc_zero_variables))
  #find the index of structure zero varaibles for house level
  struc_zero_variables_house <- intersect(struc_zero_variables_index, data$household_variable_index)
  #find the index of structure zero varaibles for individual level
  struc_zero_variables_indiv <- intersect(struc_zero_variables_index,data$individual_variable_index)

  #indexing relative to house_variable index and invidual varible index
  md$house_szv_index <- match(struc_zero_variables_house,md$household_variable_index)
  md$indiv_szv_index <- match(struc_zero_variables_indiv,md$individual_variable_index)

  #precompute NA (misisng status)
  md$NA_indiv_missing_status <- is.na(as.matrix(md$household_with_miss[,md$individual_variable_index]))
  md$NA_house_missing_status <- is.na(as.matrix(md$household_with_miss[,md$household_variable_index]))

  all_variables_index <- c(md$individual_variable_index,md$household_variable_index)
  nonstruc_zero_variables_index <-
    all_variables_index[!is.element(all_variables_index,struc_zero_variables_index)]
  md$house_non_szv_index_raw <- intersect(nonstruc_zero_variables_index, data$household_variable_index)
  md$indiv_non_szv_index_raw <- intersect(nonstruc_zero_variables_index,data$individual_variable_index)
  md$house_non_szv_index <- match(md$house_non_szv_index_raw,md$household_variable_index)
  md$indiv_non_szv_index <- match(md$indiv_non_szv_index_raw,md$individual_variable_index)

  n <- length(unique(data$household$Hhindex)) #number of households
  md$n_batch_imp_sum <- rep(miss_batch,n)
  md$n_0_reject <- rep(1,n)
  md$miss_Hhindex <- sort(unique(data$household$Hhindex[!complete.cases(data$household[,struc_zero_variables_index])]))
  md$miss_Hh_invidual_index <- list()
  for(s in md$miss_Hhindex){
    md$miss_Hh_invidual_index[[s]] <- which(md$household_with_miss$Hhindex==s)
  }
  md$n_miss <- length(md$miss_Hhindex)

  if (sum(is.na(data$household)) == 0) {
    stop("No missing entries in data or missing data isn't coded as 'NA'")
  }
  if (sum(is.na(data$household[,c("Hhindex","pernum")])) > 0) {
    stop("Hhindex and pernum cannot contain missing entries")
  }
  for(j in md$individual_variable_index){
    levels_j <- sort(unique(data$household[,j]),na.last = NA)
    data$household[is.na(data$household[,j]),j] <-
      sample(levels_j,sum(is.na(data$household[,j])),replace=T,prob=summary(as.factor(na.omit(data$household[,j]))))
  }
  for(jj in md$household_variable_index){
    levels_j <- sort(unique(data$household[,jj]),na.last = NA)
    for(i in 1:n){
      which_indiv <- which(data$household$Hhindex==i)
      if(is.na(data$household[which_indiv[1],jj])==TRUE){
        data$household[which_indiv,jj] <-
          sample(levels_j,1,replace=T,prob=summary(as.factor(na.omit(data$household[,jj]))))
      }
    }
  }

  md$household <- data$household

  print("missing data must be coded as 'NA'")
  return(md)
}


