#include <Rcpp.h>
using namespace Rcpp;
#include "MersenneTwister.h"
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

