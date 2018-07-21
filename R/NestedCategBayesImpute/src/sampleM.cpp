#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"


// [[Rcpp::export]]
NumericVector sampleM(NumericMatrix phi, IntegerMatrix data,
               NumericMatrix omega, IntegerVector G, IntegerVector serial) {
  int p = data.nrow();
  int nIndividuals = data.ncol();
  int FF = omega.nrow();
  int SS = omega.ncol();
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;
  NumericVector indi(nIndividuals);
  NumericVector rand = runif(nIndividuals);
  double *Gupdateprob2= new double[SS];
  for (int m = 0; m < nIndividuals; m++) {

    int G_asg = G[int(serial[m])-1];
    int base = m*p;
    int base2 = maxDDtp*((G_asg-1)*SS);
    for (int l = 0; l < SS; l++) {
      try {
        double phiprod = 1.0;
        for (int j = 0; j < p; j++) {
          int u = data[base+j]-1;
          phiprod *= phi[base2+j*maxdd+u];
        }
        Gupdateprob2[l] = omega[FF*l+G_asg-1]*phiprod; //can use logrithm to speed up here and also for accuracy
      } catch(...) {
        Gupdateprob2[l] = 0;
      }
    }
    indi[m] = samplew(Gupdateprob2, SS, rand[m]);
  }
  delete [] Gupdateprob2;
  return indi;
}

