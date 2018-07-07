## All missing data realated code should stay in this file
SampleNonStructureZerosIndiv <- function(household, household_with_miss,
                                         nonstruc_zero_variables_index, phi_m_g_index,
                                         individual_variable_index,para_phi,
                                         orig_d, orig_maxd) {
  nonstruc_zero_variables_indiv <- intersect(nonstruc_zero_variables_index,individual_variable_index)
  if(any(is.na(household_with_miss[,nonstruc_zero_variables_indiv]))){
    phi_m_g <- para_phi[,phi_m_g_index]
    for(k in nonstruc_zero_variables_indiv){
      is_na_hwm_k <- is.na(household_with_miss[,k])
      if(any(is_na_hwm_k)){
        real_k <- which(individual_variable_index==k)
        pr_X_miss_p <- phi_m_g[((1:orig_d[real_k]) + (real_k-1)* orig_maxd),is_na_hwm_k]
        household[is_na_hwm_k,k] <- SampleMatrixByColumnC(pr_X_miss_p,runif(ncol(pr_X_miss_p)),1)
      }
    }
  }
  return(household)
}

SampleNonStructureZerosHouse <- function(household, household_with_miss,
                                         nonstruc_zero_variables_index, household_variable_index,
                                         para_lambda, G_household_G,orig_n_i, orig_n) {
  nonstruc_zero_variables_house <- intersect(nonstruc_zero_variables_index,household_variable_index)
  #if(any(is.na(household_with_miss[,nonstruc_zero_variables_house]))){
    lambda_g <- lapply(para_lambda,function(x) x[G_household_G,])
    for(k in nonstruc_zero_variables_house){
      if(any(is.na(household_with_miss[,k]))){
        real_k <- which(household_variable_index==k)
        lambda_g_k <- lambda_g[[real_k]]
        hwm_index <- c(1,cumsum(orig_n_i[-orig_n])+1)
        is_na_hwm_k <- is.na(household_with_miss[hwm_index,k])
        pr_X_miss_p <- lambda_g_k[is_na_hwm_k,]

        sampled_values <- SampleMatrixByRowC(pr_X_miss_p,runif(nrow(pr_X_miss_p)))
        household[is.na(household_with_miss[,k]),k] <- rep(sampled_values,orig_n_i[is_na_hwm_k])
      }
    }
  #}
  return(household)
}

SampleMissing <- function(MissData,para,orig,G_household,M,hyper) {
  MissData$n_batch_imp_sum <- MissData$n_batch_imp_sum +
    ceiling(MissData$n_0_reject*MissData$prop_batch)
  MissData$n_batch_imp <- ceiling(MissData$n_batch_imp_sum/i) + 1 #no. of batches of imputations to sample
  MissData$n_0_reject[] <- 0
  MissData <- SampleMissing_imp(MissData,para,orig,G_household,M,hyper)
  return(MissData)
}


SampleMissing_imp <- function(MissData,para,orig,G_household,M,hyper){

  MissData$household <- as.matrix(MissData$household)

  #sample non structural zeros variables for everyone at once
  MissData$household <- SampleNonStructureZerosIndiv(MissData$household, MissData$household_with_miss,
                                                    MissData$nonstruc_zero_variables_index,
                                                    (M + (G_household$G_Individuals-1)*hyper$SS),
                                                    MissData$individual_variable_index,para$phi,orig$d,orig$maxd)
  MissData$household <- SampleNonStructureZerosHouse(MissData$household, MissData$household_with_miss,
                                                    MissData$nonstruc_zero_variables_index,
                                                    MissData$household_variable_index,para$lambda, G_household$G,
                                                    orig$n_i,orig$n)


  household_variables <- MissData$household[,MissData$household_variable_index]
  individual_variables  <- MissData$household[,MissData$individual_variable_index]

  household_variables_with_miss <- as.matrix(MissData$household_with_miss[,MissData$household_variable_index])
  relate_index <- which(colnames(individual_variables)=="relate")

  for(s in MissData$miss_Hhindex){
    another_index <- MissData$miss_Hh_invidual_index[[s]] #the row index for all family members
    X_house_s_prop <- household_variables[rep(another_index[1],MissData$n_batch_imp[s]),]
    NA_error_house_s <- household_variables_with_miss[rep(another_index[1],MissData$n_batch_imp[s]),]

    X_indiv_s_prop <- individual_variables[rep(another_index,MissData$n_batch_imp[s]),]
    index <- M[another_index] + (G_household$G[s]-1)*hyper$SS
    check_counter_s <- 0;
    while(check_counter_s < 1){
      for(real_k in MissData$house_szv_index){
        if(any(is.na(NA_error_house_s[,real_k]))){
          pr_X_house_k <- para$lambda[[real_k]][G_household$G[s],]
          X_house_s_prop[,real_k] <- sampleW_multi(para$lambda[[real_k]][G_household$G[s],],
                                                     runif(MissData$n_batch_imp[s]))
        }
      }
      for(real_k in MissData$indiv_szv_index){
        NA_current_indiv <- MissData$NA_indiv_missing_status[another_index,real_k]
        if (any(NA_current_indiv)) {
          pr_X_indiv_k <- para$phi[((1:orig$d[real_k]) + (real_k-1)*orig$maxd), index[NA_current_indiv]]
          if (is.matrix(pr_X_indiv_k)) {
            temp <- SampleMatrixByColumnC(pr_X_indiv_k,
                                          runif(MissData$n_batch_imp[s] * ncol(pr_X_indiv_k)),MissData$n_batch_imp[s])
          } else {
            temp <- sampleW_multi(pr_X_indiv_k,runif(MissData$n_batch_imp[s]))
          }

          X_indiv_s_prop[rep(NA_current_indiv,MissData$n_batch_imp[s]),real_k] <- temp
        }
      }
      #Check edit rules; Need to make this part more general, very specific for this data and assumes head is
      #at the household level
      X_indiv_s_prop_orig <- X_indiv_s_prop
      X_indiv_s_prop_orig[,relate_index] <- X_indiv_s_prop_orig[,relate_index] + 1 #recode relate
      comb_to_check <- X_house_s_prop[,-1]
      comb_to_check[,relate_index] <- 1 #Set relate to 1
      comb_to_check <- cbind(comb_to_check,matrix(t(X_indiv_s_prop_orig),nrow=MissData$n_batch_imp[s],byrow=TRUE))
      check_counter <- checkSZ(comb_to_check,(length(another_index) + 1))
      check_counter_s <- check_counter_s + sum(check_counter)
      if(length(which(check_counter==1))>0){
        MissData$n_0_reject[s] <- MissData$n_0_reject[s] +
          length(which(check_counter[1:which(check_counter==1)[1]]==0))
      } else{
        MissData$n_0_reject[s] <- MissData$n_0_reject[s] + MissData$n_batch_imp[s]
      }
    }
    X_house <- X_house_s_prop[which(check_counter==1)[1],]
    X_indiv <- matrix(comb_to_check[which(check_counter==1)[1],-c(1:length(MissData$individual_variable_index))],
                      byrow=TRUE,nrow=length(another_index)) #remove household head
    X_indiv[,relate_index] <- X_indiv[,relate_index] - 1 #recode relate back
    MissData$household[another_index,MissData$household_variable_index] <- rep(X_house,each=length(another_index))
    MissData$household[another_index,MissData$individual_variable_index] <- X_indiv
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

  all_variables_index <- c(md$individual_variable_index,md$household_variable_index)
  nonstruc_zero_variables_index <-
    all_variables_index[!is.element(all_variables_index,struc_zero_variables_index)]
  nonstruc_zero_variables_house <- intersect(nonstruc_zero_variables_index, data$household_variable_index)
  nonstruc_zero_variables_indiv <- intersect(nonstruc_zero_variables_index,data$individual_variable_index)
  md$house_non_szv_index <- match(nonstruc_zero_variables_house,md$household_variable_index)
  md$indiv_non_szv_index <- match(nonstruc_zero_variables_indiv,md$individual_variable_index)

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



