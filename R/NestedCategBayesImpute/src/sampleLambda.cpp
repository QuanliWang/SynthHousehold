#include <Rcpp.h>
using namespace Rcpp;

#include "MersenneTwister.h"
#include "SpecialFunctions.h"

// [[Rcpp::export]]
List UpdateLambda(IntegerMatrix HHdata_all, IntegerVector G_all, IntegerVector dHH, int FF) {
  MTRand mt;
  mt.seed();
  int p = dHH.length();
  List lambda(p);
  for (int  j = 0; j < p; j++) {
    NumericMatrix counts(FF,dHH[j]);
    int n = G_all.length();
    for (int i = 0; i < n; i++) { //group count
      counts[G_all[i] - 1 + (HHdata_all(j,i)-1) * FF]++;
    }
    for (int i = 0; i < counts.length();i++) { //gammarand sampling
      counts[i] = SpecialFunctions::gammarand(1 + counts[i], 1, mt);
    }
    for (int k =0; k < FF; k++) { //normalization
      double dsum = 0;
      for (int i = 0; i < dHH[j];i++) {
        dsum+=counts(k,i);
      }
      if (dsum <=0 ) {dsum =1;}
      for (int i = 0; i < dHH[j];i++) {
        counts(k,i) /= dsum;
      }
    }
    lambda[j] = counts;
  }
  return lambda;
}

/*
 UpdateLambda <- function(dHH,FF,G_all,HHdata_all) {
lambda <- list()
for (i in 1:length(dHH)) {
lambdacount <- groupcount(G_all,HHdata_all[i,],FF, dHH[i])
lam <- apply(lambdacount, c(1,2), function(x) rgamma(1,x+1,1))
lam <- t(apply(lam, 1, function(x) x / sum(x)))
lambda[[i]] = lam;
}
return(lambda)
}
 */
