#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"

// [[Rcpp::export]]
NumericVector samplezmember(NumericMatrix phi, NumericMatrix data,
               NumericMatrix w, NumericVector zHH, NumericVector serial) {

  printf("in samplezmember\n");
  int p = data.nrow();
  int nIndividuals = data.ncol();
  int K = w.nrow();
  int L = w.ncol();
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;
  //printf("maxdd = %d\n", maxdd);
  NumericVector indi(nIndividuals);
  NumericVector rand = runif(nIndividuals);

  double *zupdateprob2= new double[L];
  for (int m = 0; m < nIndividuals; m++) {

    int zHHasg = zHH[int(serial[m])-1];
    int base = m*p;
    for (int l = 0; l < L; l++) {
      try {
        double phiprod = 1.0;
        for (int j = 0; j < p; j++) {
          int u = (int)data[base+j]-1;
          phiprod *= phi[maxDDtp*((zHHasg-1)*L+l)+j*maxdd+u];
        }
        zupdateprob2[l] = w[K*l+zHHasg-1]*phiprod;
      } catch(...) {
        zupdateprob2[l] = 0;
      }
    }
    indi[m] = samplew(zupdateprob2, L, rand[m]);
    printf("indi = (%d,%d)\n",m, (int)indi[m]);
  }
  delete [] zupdateprob2;
  printf("done samplezmember\n");
  return indi;
}

