#include <Rcpp.h>
using namespace Rcpp;
#include "samplehouseholds.h"

// [[Rcpp::export]]
NumericMatrix samplehouseholds(NumericMatrix phi, NumericMatrix w, NumericVector pi,
               NumericVector d, List lambda,
               int currrentbatch, int nHouseholds,  int householdsize) {

  int K = w.nrow();
  int L = w.ncol();
  int p = d.length();
  int n_lambdas = lambda.length();
  int *lambda_columns = new int[n_lambdas];
  double **lambdas = new double*[n_lambdas];
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  int ncol = householdsize * DIM + 1 + householdsize;
  NumericMatrix data(nHouseholds, ncol);

  //copy data from list of matrices to C++
  for (int i = 0; i < n_lambdas; i++) {
    NumericMatrix l = lambda[i];
    lambda_columns[i] = l.ncol();
    lambdas[i] = new double[l.length()];
    std::copy(l.begin(), l.end(), lambdas[i]);
  }

  NumericVector rand = runif(nHouseholds * (householdsize *(1+p) + 2));
  //n_lambdas is not used for now, but might be useful when there are more than two household level variables
  sampleHouseholds_imp(data.begin(), rand.begin(), lambdas, lambda_columns, w.begin(),
                   phi.begin(), pi.begin(),d.begin(),
                   nHouseholds, householdsize, K, L,maxdd,p, currrentbatch,n_lambdas);

  //clean up
  delete [] lambda_columns;
  for (int i = 0; i < n_lambdas; i++) {
    delete [] lambdas[i];
  }
  delete [] lambdas;
  return data;
}
