#include <Rcpp.h>
using namespace Rcpp;
#include "samplehouseholds.h"

// [[Rcpp::export]]
NumericMatrix samplehouseholds(NumericMatrix phi, NumericMatrix w, NumericVector pi,
               NumericVector d, NumericMatrix lambda1,NumericMatrix lambda2,
               int currrentbatch, int nHouseholds,  int householdsize) {

  int K = w.nrow();
  int L = w.ncol();
  int p = d.length();
  int lambda1_columns = lambda1.ncol();

  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  int ncol = householdsize * DIM + 1 + householdsize;
  NumericMatrix data(nHouseholds, ncol);

  NumericVector rand = runif(nHouseholds * (householdsize *(1+p) + 2));
  sampleHouseholds_imp(data.begin(), rand.begin(), lambda1.begin(), lambda2.begin(), w.begin(),
                   phi.begin(), pi.begin(),d.begin(),
                   nHouseholds, householdsize, K, L,maxdd,p,lambda1_columns, currrentbatch);

  return data;
}
