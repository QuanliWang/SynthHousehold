#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"

// [[Rcpp::export]]
List samplezHH(NumericMatrix phi, NumericMatrix data,
                          NumericMatrix w, NumericVector pi, NumericVector S,
                          NumericMatrix HHdata, List lambda) {
  int p = data.nrow();
  int nIndividuals = data.ncol();
  int K = w.nrow();
  int L = w.ncol();
  int maxdd = phi.nrow() / p;
  int n = S.length();

  std::vector<NumericMatrix> Lambdas;
  for (int i = 0; i < lambda.length(); i++) {
    Lambdas.push_back(lambda[i]);
  }

  //printf("K = %d, L = %d, p = %d, maxd = %d, nIndividuals = %d, n=%d\n", K, L, p, maxdd, nIndividuals,n);

  NumericVector group(n);
  NumericVector indi(nIndividuals);
  int count = 0;

  NumericVector rand = runif(n);
  //use one-d indexing here to be consistant with Matalb code
  //might need to abandon this if we are going to abondon the Matlab version
  double *zupdateprob1 = new double[K];
  int *cumS = new int[n];
  cumS[0] = 0;
  for (int i = 1; i < n; i++) {
    cumS[i] = cumS[i-1] + (int)S[i-1];
  }
  int maxDDtP = maxdd*p;
  for (int h = 0; h < n; h++) {
    for (int k=0; k < K; k++) {
      double zupdateprod = 1.0;
      for (int memberindex=0; memberindex < S[h]; memberindex++){
        int base = (cumS[h]+memberindex)*p; //base for data
        double add = 0.0;
        for (int l=0; l < L; l++) {
          double phiprod = 1.0;
          int phi_base = (int)(maxDDtP*(k*L+l));
          for (int j=0; j < p; j++) {
            int u = (int)data[base+j]-1;
            phiprod *= phi[phi_base+j*maxdd+u];
          }
          add += w[K*l+k]*phiprod;
        } // closing l++
        zupdateprod *= add;
      } // closing member++

      for (int hv = 0; hv < Lambdas.size(); hv++) {
        zupdateprod *= Lambdas[hv][(HHdata[h+hv*n]-1)*K+k];
      }
      zupdateprob1[k] = pi[k]*zupdateprod;
    } // closing k++
    group[h] = samplew(zupdateprob1, K, rand[h]);
    for (int m=0; m < S[h];m++) {
      indi[count++] = group[h];
    }
  }
  delete [] cumS;
  delete [] zupdateprob1;
  return List::create(Named("z_HH", group), Named("z_HH_Individuals", indi));
}

