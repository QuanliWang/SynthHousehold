SampleMissing <- function(MissData,para,orig, household_variable_index,individual_variable_index,
                          G_household,M,hyper){
  nonstruc_zero_variables <- MissData$nonstruc_zero_variables
  nonstruc_zero_variables_house <-
    nonstruc_zero_variables[is.element(nonstruc_zero_variables,household_variable_index)]
  nonstruc_zero_variables_indiv <-
    nonstruc_zero_variables[is.element(nonstruc_zero_variables,individual_variable_index)]
  struc_zero_variables <- MissData$struc_zero_variables
  struc_zero_variables_house <-
    struc_zero_variables[is.element(struc_zero_variables,household_variable_index)]
  struc_zero_variables_indiv <-
    struc_zero_variables[is.element(struc_zero_variables,individual_variable_index)]

  #sample non structural zeros variables for everyone at once
  if(sum(is.na(MissData$household_with_miss[,nonstruc_zero_variables_indiv])) > 0){
    phi_m_g <- t(para$phi[,(M + (G_household$G_Individuals-1)*hyper$SS)])
    for(k in nonstruc_zero_variables_indiv){
      if(sum(is.na(MissData$household_with_miss[,k]))>0){
        real_k <- which(individual_variable_index==k)
        pr_X_miss_p <- phi_m_g[which(is.na(MissData$household_with_miss[,k])==TRUE),
                               ((1:orig$d[real_k]) + (real_k-1)*orig$maxd)]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        #cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        cumul_miss_p <- t(apply(pr_X_miss_p,1,cumsum))
        MissData$household[is.na(MissData$household_with_miss[,k]),k] <- rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L
      }
    }
  }
  if(sum(is.na(MissData$household_with_miss[,nonstruc_zero_variables_house])) > 0){
    lambda_g <- lapply(para$lambda,function(x) x[G_household$G,])
    for(kk in nonstruc_zero_variables_house){
      if(sum(is.na(MissData$household_with_miss[,kk]))>0){
        real_kk <- which(household_variable_index==kk)
        lambda_g_kk <- lambda_g[[real_kk]]
        pr_X_miss_p <-
          lambda_g_kk[which(is.na(MissData$household_with_miss[c(1,cumsum(orig$n_i[-orig$n])+1),kk])==TRUE),]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        #cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        cumul_miss_p <- t(apply(pr_X_miss_p,1,cumsum))
        sampled_values <- rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L
        MissData$household[is.na(MissData$household_with_miss[,kk]),kk] <- rep(sampled_values,orig$n_i[
          which(is.na(MissData$household_with_miss[c(1,cumsum(orig$n_i[-orig$n])+1),kk])==TRUE)])
      }
    }
  }

  #SSS <- as.data.frame(MissData$miss_Hhindex)
  #names(SSS) <- "Hhindex"
  #df <- as.data.frame(MissData$household_with_miss$Hhindex)
  #names(df) <- "Hhindex"
  #df <- df %>% mutate(index = 1:length(MissData$household_with_miss$Hhindex))
  #df <- SSS %>% inner_join(df)

  MissData$household <- as.matrix(MissData$household)
  household_variables <- as.matrix(MissData$household[,household_variable_index])
  individual_variables  <- as.matrix(MissData$household[,individual_variable_index])

  household_variables_with_miss <- as.matrix(MissData$household_with_miss[,household_variable_index])
  individual_variables_with_miss  <- as.matrix(MissData$household_with_miss[,individual_variable_index])

  #sample structural zeros variables one household at a time
  for(sss in MissData$miss_Hhindex){
    another_index <- which(MissData$household_with_miss$Hhindex==sss)
    X_house_sss_prop <- household_variables[rep(another_index[1],MissData$n_batch_imp[sss]),]
    X_indiv_sss_prop <- individual_variables[rep(another_index,MissData$n_batch_imp[sss]),]
    relate_index <- which(colnames(X_indiv_sss_prop)=="relate")
    NA_error_house_sss <- household_variables_with_miss[rep(another_index[1],MissData$n_batch_imp[sss]),]
    NA_error_indiv_sss <- individual_variables_with_miss[rep(another_index,MissData$n_batch_imp[sss]),]

    G_prop <- rep(G_household$G[sss],MissData$n_batch_imp[sss])
    M_prop <- rep(M[another_index],MissData$n_batch_imp[sss])
    lambda_g <- lapply(para$lambda,function(x) x[G_prop,])
    phi_m_g <- t(para$phi[,(M_prop + (G_household$G[sss]-1)*hyper$SS)])
    check_counter_sss <- 0;
    while(check_counter_sss < 1){
      for(kkk in struc_zero_variables_house){
        real_kkk <- which(household_variable_index==kkk)
        if(sum(is.na(NA_error_house_sss[,real_kkk]))>0){
          pr_X_house_k <- lambda_g[[real_kkk]]
          Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
          #cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
          cumul_X_house_k <- t(apply(pr_X_house_k,1,cumsum))
          X_house_sss_prop[,real_kkk] <- rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L
        }
      }
      for(kkkk in struc_zero_variables_indiv){
        real_kkkk <- which(individual_variable_index==kkkk)
        if(sum(is.na(NA_error_indiv_sss[,real_kkkk]))>0){
          pr_X_indiv_k <- phi_m_g[is.na(NA_error_indiv_sss[,real_kkkk]),
                                  ((1:orig$d[real_kkkk]) + (real_kkkk-1)*orig$maxd)]

          #cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
          cumul_X_indiv_k <- t(apply(pr_X_indiv_k,1,cumsum))
          Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
          X_indiv_sss_prop[is.na(NA_error_indiv_sss[,real_kkkk]),real_kkkk] <-
            rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L
        }
      }
      #Check edit rules; Need to make this part more general, very specific for this data and assumes head is
      #at the household level
      X_indiv_sss_prop_orig <- X_indiv_sss_prop
      X_indiv_sss_prop_orig[,relate_index] <- X_indiv_sss_prop_orig[,relate_index] + 1 #recode relate
      comb_to_check <- X_house_sss_prop[,-1]
      comb_to_check[,relate_index] <- 1 #Set relate to 1
      comb_to_check <- cbind(comb_to_check,matrix(t(X_indiv_sss_prop_orig),nrow=MissData$n_batch_imp[sss],byrow=TRUE))
      check_counter <- checkSZ(comb_to_check,(length(another_index) + 1))
      check_counter_sss <- check_counter_sss + sum(check_counter)
      if(length(which(check_counter==1))>0){
        MissData$n_0_reject[sss] <- MissData$n_0_reject[sss] +
          length(which(check_counter[1:which(check_counter==1)[1]]==0))
      } else{
        MissData$n_0_reject[sss] <- MissData$n_0_reject[sss] + MissData$n_batch_imp[sss]
      }
    }
    X_house <- X_house_sss_prop[which(check_counter==1)[1],]
    X_indiv <- matrix(comb_to_check[which(check_counter==1)[1],-c(1:length(individual_variable_index))],
                      byrow=TRUE,nrow=length(another_index)) #remove household head
    X_indiv[,relate_index] <- X_indiv[,relate_index] - 1 #recode relate back
    MissData$household[another_index,household_variable_index] <- rep(X_house,each=length(another_index))
    MissData$household[another_index,individual_variable_index] <- X_indiv
  }
  MissData$household <- as.data.frame(MissData$household)
  return(MissData)
}




