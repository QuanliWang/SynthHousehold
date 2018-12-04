#include <Rcpp.h>
using namespace Rcpp;
#include "SpecialFunctions.h"

// [[Rcpp::export]]
NumericVector gammarand(int n, double shape, double rate) {
  MTRand mt;
  mt.seed();
  vector<double> result;
  SpecialFunctions::gammarand(shape,1.0 /rate,n,mt,result);
  NumericVector r(result.begin(),result.end());
  return r;
}

NumericMatrix samplePhi(IntegerMatrix counts) {
  NumericMatrix result(counts.rows(),counts.cols());
  MTRand mt;
  mt.seed();
  for (int i = 0; i < counts.length(); i++) {
    result[i] = SpecialFunctions::gammarand(1 + counts[i], 1, mt);
  }

  return result;
}

// [[Rcpp::export]]
NumericMatrix UpdatePhi(IntegerMatrix data, IntegerMatrix M_all, int FF, int SS, IntegerVector d, int maxd) {
  MTRand mt;
  mt.seed();
  int p = d.length();
  int groups =  FF * SS;
  int phi_rows = maxd * p;
  NumericMatrix phi(phi_rows , groups);
  int n = M_all.ncol();
  IntegerVector groupIndex(n);
  for (int i = 0; i <n; i++) {
    groupIndex[i] = SS*(M_all(0,i)-1) + M_all(1,i) - 1;
  }
  for (int j = 0; j < p; j++) {

    NumericMatrix counts(groups,d[j]);
    for (int i = 0; i < n; i++) { //group count
      counts[groupIndex[i] + (data(j,i)-1) * groups]++;
    }
    for (int i = 0; i < counts.length();i++) { //gammarand sampling
        counts[i] = SpecialFunctions::gammarand(1 + counts[i], 1, mt);
    }

    for (int k =0; k < groups; k++) { //normalization
      int base = k * phi_rows + j * maxd;
      double dsum = 0;
      for (int i = 0; i < d[j];i++) {
        dsum+=counts(k,i);
      }
      if (dsum <=0 ) {dsum =1;}
      for (int i = 0; i < d[j];i++) {
        phi[base + i] = counts(k,i) / dsum;
        //counts(k,i) /= dsum;
      }
    }
  }
  return phi;
}

// [[Rcpp::export]]
NumericMatrix UpdatePhiWeighted(List data, List M_all, int FF, int SS, IntegerVector d, int maxd, NumericVector struc_weight) {
  MTRand mt;
  mt.seed();
  int p = d.length();
  int groups =  FF * SS;
  int phi_rows = maxd * p;
  NumericMatrix phi(phi_rows , groups);

  std::vector<IntegerVector> groupIndexes;
  for (int i = 0; i < struc_weight.length(); i++) {
    IntegerMatrix M = M_all[i];
    int n = M.ncol();
    IntegerVector groupIndex(n);
    for (int g = 0; g <n; g++) {
      groupIndex[g] = SS*(M(0,g)-1) + M(1,g) - 1;
    }
    groupIndexes.push_back(groupIndex);
  }

  for (int j = 0; j < p; j++) {
    NumericMatrix counts(groups,d[j]);
    for (int s = 0; s < struc_weight.length(); s++) {
      double weight = 1.0 / struc_weight[s];
      IntegerMatrix current_data = data[s];
      for (int i = 0; i < current_data.ncol(); i++) { //group count
        counts[groupIndexes[s][j] + (current_data(j,i)-1) * groups] += weight;
      }
    }

    for (int i = 0; i < counts.length();i++) { //gammarand sampling
      counts[i] = SpecialFunctions::gammarand(1 + counts[i], 1, mt);
    }

    for (int k =0; k < groups; k++) { //normalization
      int base = k * phi_rows + j * maxd;
      double dsum = 0;
      for (int i = 0; i < d[j];i++) {
        dsum+=counts(k,i);
      }
      if (dsum <=0 ) {dsum =1;}
      for (int i = 0; i < d[j];i++) {
        phi[base + i] = counts(k,i) / dsum;
        //counts(k,i) /= dsum;
      }
    }
  }
  return phi;
}


/*
UpdatePhiWeighted <- function(data, M_all, FF, SS, p, d, maxd, struc_weight) {
  phi <- array(0,dim = c(maxd,p, FF*SS))

  groupIndex <- lapply(M_all,function(x) SS*(x[1,]-1)+x[2,])
  for (j in 1:p) {

    phicount <- 0
    for(w_i in 1:length(struc_weight)){
      data_w_i <- data[[w_i]]
      phicount <- phicount + (groupcount(groupIndex[[w_i]], data_w_i[j,], FF*SS, d[j]) / struc_weight[w_i])
    }

    phi_j <- apply(phicount, c(1,2), function(x) rgamma(1,x+1,1))
      phi[1:d[j], j,] <- apply(phi_j, 1, function(x) x / sum(x))
  }
  dim(phi) <- c(maxd*p,FF * SS) #reshape to a 2D matrix
    return(phi)
}
*/

/*
UpdatePhi <- function(data, M_all, FF, SS, d, maxd) {
  p <- length(d)
  phi <- array(0,dim = c(maxd,p, FF*SS))
  groupIndex <- SS*(M_all[1,]-1)+M_all[2,]
  for (j in 1:p) {
    phicount <- groupcount(groupIndex, data[j,], FF*SS, d[j])
    phi_j <- samplePhi(phicount);
    phi[1:d[j], j,] <- apply(phi_j, 1, function(x) x / sum(x))
  }
  dim(phi) <- c(maxd*p,FF * SS) #reshape to a 2D matrix
    return(phi)
}
 */


