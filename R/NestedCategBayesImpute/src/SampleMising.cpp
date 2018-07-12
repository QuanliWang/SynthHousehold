#include "sampleW.h"

// [[Rcpp::export]]
IntegerVector CheckSZ_batch(IntegerMatrix X_house, IntegerMatrix X_indiv) {
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
  IntegerMatrix comb_to_check(batch, household_size * nvaraibles);
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
List SampleMissingForOneHousehold_batch(IntegerVector another_index,
                                        IntegerMatrix X_house_s_prop, IntegerMatrix X_indiv_s_prop,
                                        IntegerVector house_szv_index, LogicalMatrix NA_house_missing_status,
                                        IntegerVector indiv_szv_index, LogicalMatrix NA_indiv_missing_status,
                                        List lambda,NumericMatrix phi,int G_household_G_s, IntegerVector index,
                                        IntegerVector orig_d, int orig_maxd,
                                        int batch) {
  int n_0_reject = 0;
  int first_valid = 0;
  int nvidiv = another_index.length();
  while(first_valid == 0){
    for(int i = 0; i < house_szv_index.length(); i++) {
      int real_k =  house_szv_index[i] - 1; //to zero-based index
      if(NA_house_missing_status(another_index[0] - 1,real_k) != 0) {
        NumericMatrix lambda_real_k = lambda[real_k];
        NumericVector lambda_row = lambda_real_k.row(G_household_G_s - 1);
        X_house_s_prop.column(real_k) = sampleW_multi(lambda_row,runif(batch));
      }
    }
    for(int i = 0; i < indiv_szv_index.length(); i++) {
      int real_k =  indiv_szv_index[i] - 1; //to zero-based index
      int base = real_k*orig_maxd;
      for (int j = 0; j < nvidiv;j++) {
        int index_j = index[j] - 1;
        if (NA_indiv_missing_status(another_index[j]-1, real_k) != 0) {
          int n = orig_d[real_k];
          NumericVector p(n);
          for (int s = 0; s <  n; s++) {
            p[s] = phi(base+s,index_j);
          }
          IntegerVector samples = sampleW_multi(p,runif(batch));
          for (int s = 0; s < batch; s++) {
            X_indiv_s_prop(j + s * nvidiv, real_k) = samples[s];
          }
        }
      }
    }

    //Check edit rules; Need to make this part more general, very specific for
    //this data and assumes head is at the household level
    IntegerVector nv_first_valid = CheckSZ_batch(X_house_s_prop, X_indiv_s_prop);
    first_valid = nv_first_valid[0];
    if (first_valid > 0) {
      n_0_reject = n_0_reject + first_valid - 1;
    } else {
      n_0_reject = n_0_reject + batch;
    }
  }

  return List::create(Named("n_0_reject", n_0_reject),
                      Named("first_valid", first_valid),
                      Named("X_house_s_prop", X_house_s_prop),
                      Named("X_indiv_s_prop", X_indiv_s_prop));
}


// [[Rcpp::export]]
IntegerMatrix SampleNonStructureZerosHouseC(IntegerMatrix household,
                                            LogicalMatrix NA_house_missing_status,
                                            IntegerVector house_non_szv_index_raw,
                                            IntegerVector house_non_szv_index,
                                            List para_lambda,
                                            IntegerVector G_household_G,
                                            IntegerVector orig_n_i
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
IntegerMatrix SampleNonStructureZerosIndivC(IntegerMatrix household,
                                           LogicalMatrix NA_indiv_missing_status,
                                           IntegerVector indiv_non_szv_index_raw,
                                           NumericVector phi_m_g_index,
                                           IntegerVector indiv_non_szv_index,
                                           NumericMatrix para_phi,
                                           IntegerVector orig_d,
                                           IntegerVector orig_maxd) {
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

/*
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
*/
