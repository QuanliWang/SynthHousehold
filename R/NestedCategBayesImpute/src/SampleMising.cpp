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
    int nhouseholds = G_household_G.length();
    NumericVector r = runif(nhouseholds);
    for (int i = 0; i < nhouseholds; i++) {
      if (NA_house_missing_status(cumcount,real_k)) {
        NumericVector pr_X_miss_p(lambda_k.row(G_household_G[i] - 1));
        int sample = samplew(pr_X_miss_p.begin(), pr_X_miss_p.length(), r[i]);
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
                                           IntegerVector phi_m_g_index,
                                           IntegerVector indiv_non_szv_index,
                                           NumericMatrix para_phi,
                                           IntegerVector orig_d,
                                           IntegerVector orig_maxd) {
  for (int count = 0; count < indiv_non_szv_index_raw.length(); count++) {
    int k = indiv_non_szv_index_raw[count] - 1; //index into raw data
    int real_k = indiv_non_szv_index[count] - 1; //index into invividual level variables
    NumericVector r = runif(NA_indiv_missing_status.rows());
    for (int i = 0; i < NA_indiv_missing_status.rows(); i++) {
      if (NA_indiv_missing_status(i,real_k)) {
        int col = phi_m_g_index[i];
        int offset = (real_k-1)* orig_maxd[0];
        household(i,k) = samplew(para_phi.column(col -1).begin() + offset, orig_d[real_k], r[i]);
      }
    }
  }

  return(household);
}

// [[Rcpp::export]]
List SampleMissing_impC(List MissData, List para, List orig,List G_household, IntegerVector M, List hyper) {
  IntegerVector G_Individuals =  G_household["G_Individuals"];
  IntegerVector G = G_household["G"];
  int hyper_SS = hyper["SS"];
  List lambda = para["lambda"];
  NumericMatrix phi = para["phi"];
  IntegerVector d = orig["d"];
  IntegerVector maxd = orig["maxd"];
  IntegerMatrix household = MissData["household"];
  LogicalMatrix NA_house_missing_status = MissData["NA_house_missing_status"];
  LogicalMatrix NA_indiv_missing_status = MissData["NA_indiv_missing_status"];
  IntegerVector phi_m_g_index(M.length());
  for (int i = 0; i < M.length(); i++) {
    phi_m_g_index[i] = M[i] + (G_Individuals[i] -1) * hyper_SS;
  }

  household = SampleNonStructureZerosIndivC(household,
                                            NA_indiv_missing_status,
                                            as<IntegerVector>(MissData["indiv_non_szv_index_raw"]),
                                            phi_m_g_index,
                                            as<IntegerVector>(MissData["indiv_non_szv_index"]),
                                            phi,d,maxd);
  household = SampleNonStructureZerosHouseC(household,
                                            NA_house_missing_status,
                                            as<IntegerVector>(MissData["house_non_szv_index_raw"]),
                                            as<IntegerVector>(MissData["house_non_szv_index"]),
                                            lambda,G,
                                            as<IntegerVector>(orig["n_i"]));

  IntegerVector miss_Hhindex = MissData["miss_Hhindex"];
  List miss_Hh_invidual_index = MissData["miss_Hh_invidual_index"];
  IntegerVector batches = MissData["n_batch_imp"];
  IntegerVector household_variable_index = MissData["household_variable_index"];
  IntegerVector individual_variable_index = MissData["individual_variable_index"];
  IntegerVector house_szv_index = MissData["house_szv_index"];
  IntegerVector indiv_szv_index = MissData["indiv_szv_index"];

  IntegerVector n_0_reject = MissData["n_0_reject"];

  for (int i = 0; i < miss_Hhindex.length(); i++) {
    int s = miss_Hhindex[i] - 1;
    IntegerVector another_index = miss_Hh_invidual_index[s]; //the row index for all other family members
    int n_indiv = another_index.length();
    IntegerMatrix X_house(batches[s], household_variable_index.length());
    for (int b = 0; b < batches[s]; b++) {
        int house_row =  another_index[0] -1;
        for (int v = 0; v < household_variable_index.length();v++) {
          X_house(b,v) = household(house_row,household_variable_index[v] - 1);
        }
    }
    IntegerMatrix X_indiv(batches[s] * n_indiv, individual_variable_index.length());
    for (int b = 0; b < batches[s]; b++) {
      for (int ind = 0; ind < another_index.length(); ind++ ) {
        int indiv_row =  another_index[ind] -1;
        for (int v = 0; v < individual_variable_index.length();v++) {
          X_indiv( b*n_indiv + ind,v) = household(indiv_row,individual_variable_index[v] - 1);
        }
      }
    }

    IntegerVector index(another_index.length());
    for (int j = 0; j < index.length(); j++) {
      index[j] = M[another_index[j] - 1] + (G[s] -1) * hyper_SS;
    }

    List OneHousehold = SampleMissingForOneHousehold_batch(another_index, X_house, X_indiv,house_szv_index,
                                                           NA_house_missing_status,
                                                           indiv_szv_index, NA_indiv_missing_status,
                                                           lambda,phi, G[s], index,
                                                           d,maxd[0],batches[s]);

    n_0_reject[s] = n_0_reject[s] + as<int>(OneHousehold["n_0_reject"]);
    int first_valid = as<int>(OneHousehold["first_valid"]) - 1;
    IntegerMatrix House = OneHousehold["X_house_s_prop"];
    IntegerMatrix Indivs = OneHousehold["X_indiv_s_prop"];

    for (int ind = 0; ind < n_indiv; ind++) {
      int indiv_row =  another_index[ind] -1;
      for (int v = 0; v < household_variable_index.length();v++) {
        household(indiv_row,household_variable_index[v] - 1) = House(first_valid,v);
      }

      for (int v = 0; v < individual_variable_index.length();v++) {
        household(indiv_row,individual_variable_index[v] - 1) = Indivs( first_valid*n_indiv + ind,v);
      }
    }
  }
  MissData["household"] = household; //maybe not need as R/C++ all IntegerMatrix
  MissData["n_0_reject"] = n_0_reject; //numeric vector in R, integer vectort in Rcpp
  return(MissData);
}

/*
 * SampleMissing_imp <- function(MissData,para,orig,G_household,M,hyper){

 MissData$household <- SampleNonStructureZerosIndivC(MissData$household, MissData$NA_indiv_missing_status,
                                                     MissData$indiv_non_szv_index_raw,
(M + (G_household$G_Individuals-1)*hyper$SS),
MissData$indiv_non_szv_index,para$phi,orig$d,orig$maxd)
MissData$household <- SampleNonStructureZerosHouseC(MissData$household, MissData$NA_house_missing_status,
                                                    MissData$house_non_szv_index_raw,
MissData$house_non_szv_index,para$lambda, G_household$G,
orig$n_i)

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
return(MissData)
}
 */
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
}*/

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
