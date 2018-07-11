#include "sampleW.h"
// [[Rcpp::export]]
NumericVector CheckSZ_batch(NumericMatrix X_house, NumericMatrix X_indiv) {
  if (X_house.ncol() != 6) {
    Rcout << "Household matrix must have 6 columns" << std::endl;
    return R_NilValue;
  }
  if (X_indiv.ncol() != 5) {
    Rcout << "Household matrix must have 6 columns" << std::endl;
    return R_NilValue;
  }
  int batch = X_house.nrow();
  int household_size = X_indiv.nrow() / batch + 1;
  int nvaraibles = X_indiv.ncol();
  NumericMatrix comb_to_check(batch, household_size * nvaraibles);
  for (int i = 0; i < batch; i++) {
    int count = 0;
    for (int j = 1; j < nvaraibles; j++) { //skip first household variables
      comb_to_check(i,count++) = X_house(i,j);
    }
    comb_to_check(i,count++) = 1; //Set relate to 1
    int indiv_start_row = i * (household_size -1); //household_size -1 is the number of other household members
    for (int j = indiv_start_row; j < indiv_start_row + (household_size -1); j++) {
      for (int k = 0; k < nvaraibles - 1; k++) {
        comb_to_check(i,count++) = X_indiv(j,k);
      }
      comb_to_check(i,count++) = X_indiv(j,nvaraibles -1) + 1; //recode relate
    }
  }
  return checkSZ2(comb_to_check,household_size);
}

// [[Rcpp::export]]
NumericMatrix SampleNonStructureZerosHouseC(NumericMatrix household,
                                            NumericMatrix NA_house_missing_status,
                                            NumericVector house_non_szv_index_raw,
                                            NumericVector house_non_szv_index,
                                            List para_lambda,
                                            NumericVector G_household_G,
                                            NumericVector orig_n_i
                                            ) {
  for (int count = 0; count < house_non_szv_index_raw.length(); count++) {
    int k = house_non_szv_index_raw[count] - 1; //index into raw data
    int real_k = house_non_szv_index[count] - 1; //index into invividual level variables
    NumericMatrix lambda_k = para_lambda[real_k];
    int cumcount = 0;

    for (int i = 0; i < orig_n_i.length(); i++) {
      if (NA_house_missing_status(cumcount,real_k)) {
        NumericVector pr_X_miss_p(lambda_k.row(G_household_G[count] - 1));
        int sample = samplew(pr_X_miss_p.begin(), pr_X_miss_p.length(), runif(1)[0]);
        for (int pos = cumcount; pos < cumcount + orig_n_i[i]; pos++) {
          household(pos,k) = sample;
        }
      }
      cumcount += orig_n_i[i];
    }

  }
  return(household);
}

// [[Rcpp::export]]
NumericMatrix SampleNonStructureZerosIndivC(NumericMatrix household,
                                           NumericMatrix NA_indiv_missing_status,
                                           NumericVector indiv_non_szv_index_raw,
                                           NumericVector phi_m_g_index,
                                           NumericVector indiv_non_szv_index,
                                           NumericMatrix para_phi,
                                           NumericVector orig_d,
                                           NumericVector orig_maxd) {
  for (int count = 0; count < indiv_non_szv_index_raw.length(); count++) {
    int k = indiv_non_szv_index_raw[count] - 1; //index into raw data
    int real_k = indiv_non_szv_index[count] - 1; //index into invividual level variables
    NumericVector r = runif(NA_indiv_missing_status.cols());
    for (int i = 0; i < NA_indiv_missing_status.cols(); i++) {
      if (NA_indiv_missing_status(i,real_k)) {
        int col = phi_m_g_index[i];
        int offset = (real_k-1)* orig_maxd[0];
        household(i,k) = samplew(para_phi.column(col -1).begin() + offset, orig_d[real_k], r[i]);
      }
    }
  }

  return household;
}

/*
SampleNonStructureZerosIndiv <- function(household, NA_indiv_missing_status,
                                         indiv_non_szv_index_raw, phi_m_g_index,
                                         indiv_non_szv_index,para_phi,
                                         orig_d, orig_maxd) {
  for (count in 1: length(indiv_non_szv_index_raw)) {
    k <- indiv_non_szv_index_raw[count] #index into raw data
    real_k <- indiv_non_szv_index[count] #index into invividual level variables
    is_na_hwm_k <- NA_indiv_missing_status[,real_k]
    if(any(is_na_hwm_k)){
      pr_X_miss_p <- para_phi[((1:orig_d[real_k]) + (real_k-1)* orig_maxd), phi_m_g_index[is_na_hwm_k]]
      household[is_na_hwm_k,k] <- SampleMatrixByColumnC(pr_X_miss_p,runif(ncol(pr_X_miss_p)),1)
    }
  }
  return(household)
}
*/

//hwm_index <- c(1,cumsum(orig$n_i[-orig$n])+1)
/*
 SampleNonStructureZerosHouse <- function(household, NA_house_missing_status,
 house_non_szv_index_raw, house_non_szv_index,
 para_lambda, G_household_G,hwm_index,orig_n_i) {
 for (count in 1: length(house_non_szv_index_raw)) {
 k <- house_non_szv_index_raw[count]
 real_k <- house_non_szv_index[count]
 if (any(NA_house_missing_status[,real_k])){

 is_na_hwm_k <- NA_house_missing_status[hwm_index,real_k]
 pr_X_miss_p <- para_lambda[[real_k]][G_household_G[is_na_hwm_k],]

 sampled_values <- SampleMatrixByRowC(pr_X_miss_p,runif(nrow(pr_X_miss_p)))
 household[NA_house_missing_status[,real_k],k] <- rep(sampled_values,orig_n_i[is_na_hwm_k])
 }
 }
 return(household)
 }
 */
/*
## All missing data realated code should stay in this file
 CheckSZ_HV_batch <- function(X_house, X_indiv, relate_index, household_size) {
X_indiv_orig <- X_indiv
X_indiv_orig[,relate_index] <- X_indiv_orig[,relate_index] + 1 #recode relate
comb_to_check <- X_house[,-1]
comb_to_check[,relate_index] <- 1 #Set relate to 1
comb_to_check <- cbind(comb_to_check,matrix(t(X_indiv_orig),nrow=dim(X_house)[1],byrow=TRUE))
check_counter <- checkSZ2(comb_to_check,household_size)
return(check_counter)
}
*/
