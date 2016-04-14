#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"

// [[Rcpp::export]]
NumericVector samplezmember(NumericMatrix phi, NumericMatrix data,
               NumericMatrix w, NumericVector zHH, NumericVector serial) {

  int p = data.nrow();
  int nIndividuals = data.ncol();
  int K = w.nrow();
  int L = w.ncol();
  int maxdd = phi.nrow() / p;

  NumericVector indi(nIndividuals);
  NumericVector rand = runif(nIndividuals);

  int maxDDtp = maxdd*p;
  double *zupdateprob2= new double[L];
  for (int m = 0; m < nIndividuals; m++) {
    int zHHasg = zHH[int(serial[m])-1];
    int base = m*p;
    double updatesum = 0.0;
    for (int l = 0; l < L; l++) {
      double phiprod = 1.0;
      for (int j = 0; j < p; j++) {
        int u = (int)data[base+j]-1;
        phiprod *= phi[maxDDtp*((zHHasg-1)*L+l)+j*maxdd+u];
      }
      zupdateprob2[l] = w[K*l+zHHasg-1]*phiprod;
    }
    indi[m] = samplew(zupdateprob2, L, rand[m]);
  }
  delete [] zupdateprob2;
  return indi;
}

