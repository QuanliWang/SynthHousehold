SampleMissingForOneHousehold_batch <- function(another_index,X_house_s_prop, X_indiv_s_prop,
                                         house_szv_index, NA_house_missing_status,
                                         indiv_szv_index, NA_indiv_missing_status,
                                         lambda,phi,G_household_G_s, index,
                                         orig_d,orig_maxd,
                                         batch) {
  n_0_reject = 0;
  first_valid <- 0
  while(first_valid == 0){
    for(real_k in house_szv_index){
      if(NA_house_missing_status[another_index[1],real_k]) {
        X_house_s_prop[,real_k] <- sampleW_multi(lambda[[real_k]][G_household_G_s,],runif(batch))
      }
    }
    for(real_k in indiv_szv_index){
      NA_current_indiv <- NA_indiv_missing_status[another_index,real_k]
      if (any(NA_current_indiv)) {
        pr_X_indiv_k <- phi[((1:orig_d[real_k]) + (real_k-1)*orig_maxd), index[NA_current_indiv]]
        if (is.matrix(pr_X_indiv_k)) {
          temp <- SampleMatrixByColumnC(pr_X_indiv_k,runif(batch * ncol(pr_X_indiv_k)), batch)
        } else {
          temp <- sampleW_multi(pr_X_indiv_k,runif(batch))
        }
        X_indiv_s_prop[rep(NA_current_indiv,batch),real_k] <- temp
      }
    }

    #Check edit rules; Need to make this part more general, very specific for this data and assumes head is
    #at the household level
    first_valid <- CheckSZ_batch(X_house_s_prop, X_indiv_s_prop)
    if (first_valid > 0) {
      n_0_reject <- n_0_reject + first_valid - 1
    } else {
      n_0_reject <- n_0_reject + batch
    }
  }
  return(list(n_0_reject = n_0_reject, first_valid = first_valid,
              X_house_s_prop = X_house_s_prop, X_indiv_s_prop = X_indiv_s_prop))
}

SampleMissing_imp <- function(MissData,para,orig,G_household,M,hyper){
  MissData$household <- as.matrix(MissData$household)

  #sample non structural zeros variables for everyone at once
  MissData$household <- SampleNonStructureZerosIndivC(MissData$household, MissData$NA_indiv_missing_status,
                                                    MissData$indiv_non_szv_index_raw,
                                                    (M + (G_household$G_Individuals-1)*hyper$SS),
                                                    MissData$indiv_non_szv_index,para$phi,orig$d,orig$maxd)

  MissData$household <- SampleNonStructureZerosHouseC(MissData$household, MissData$NA_house_missing_status,
                                                    MissData$house_non_szv_index_raw,
                                                    MissData$house_non_szv_index,para$lambda, G_household$G,
                                                    orig$n_i)


  relate_index <- which(colnames(MissData$household)[MissData$individual_variable_index] == "relate")
  stopifnot(relate_index == 5) #relate_index should always be 5

  for(s in MissData$miss_Hhindex){
    another_index <- MissData$miss_Hh_invidual_index[[s]] #the row index for all family members
    n_indiv <- length(another_index)
    X_house_s_prop <- MissData$household[rep(another_index[1],MissData$n_batch_imp[s]), MissData$household_variable_index]
    X_indiv_s_prop <- MissData$household[rep(another_index,   MissData$n_batch_imp[s]), MissData$individual_variable_index]

    index <- M[another_index] + (G_household$G[s]-1)*hyper$SS
    OneHousehold <- SampleMissingForOneHousehold_batch(another_index,X_house_s_prop, X_indiv_s_prop,
                                                 MissData$house_szv_index, MissData$NA_house_missing_status,
                                                 MissData$indiv_szv_index, MissData$NA_indiv_missing_status,
                                                 para$lambda,para$phi,G_household$G[s], index,
                                                 orig$d,orig$maxd,
                                                 MissData$n_batch_imp[s])
    MissData$n_0_reject[s] <- MissData$n_0_reject[s] + OneHousehold$n_0_reject
    MissData$household[another_index,MissData$household_variable_index] <-
      rep(OneHousehold$X_house_s_prop[OneHousehold$first_valid,],each=n_indiv)
    MissData$household[another_index,MissData$individual_variable_index] <-
      OneHousehold$X_indiv_s_prop[((OneHousehold$first_valid-1)*n_indiv+1) : (OneHousehold$first_valid*n_indiv),]
  }

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

SampleMissing <- function(MissData,para,orig,G_household,M,hyper) {
  MissData$n_batch_imp_sum <- MissData$n_batch_imp_sum +
    ceiling(MissData$n_0_reject*MissData$prop_batch)
  MissData$n_batch_imp <- ceiling(MissData$n_batch_imp_sum/i) + 1 #no. of batches of imputations to sample
  MissData$n_0_reject[] <- 0
  MissData <- SampleMissing_imp(MissData,para,orig,G_household,M,hyper)
  return(MissData)
}

