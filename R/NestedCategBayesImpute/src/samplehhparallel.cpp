#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

#include "samplehouseholds.h"
#include "utils.h"

struct HeadAtGroupLevelHHSampler : public Worker
{
  // source matrix
  RMatrix<double> phi;
  RMatrix<double> omega;
  RVector<double> pi;
  RVector<int> d;
  List lambda;
  int currrentbatchbase;
  int householdsize;

  int FF;
  int SS;
  int p;
  int n_lambdas;
  int *lambda_columns = NULL;
  double **lambdas = NULL;
  int maxDDtp;
  int maxdd;

  int DIM ; //not output the household size
  int ncol;

  NumericMatrix omegaT;

  // destination matrix
  RMatrix<int> data;
  NumericVector r; //at most this many

  // initialize with source and destination
  HeadAtGroupLevelHHSampler(NumericMatrix phi, NumericMatrix omega,
                            NumericVector pi,IntegerVector d,
                            List lambda,
                            int currrentbatchbase, int householdsize,
                            IntegerMatrix data)
    : phi(phi), omega(omega), pi(pi), d(d), lambda(lambda), currrentbatchbase(currrentbatchbase), householdsize(householdsize),data(data) {
    FF = omega.nrow();
    SS = omega.ncol();
    p = d.length();
    n_lambdas = lambda.length();
    lambda_columns = new int[n_lambdas];
    lambdas = new double*[n_lambdas];
    maxDDtp = phi.nrow();
    maxdd = maxDDtp / p;

    DIM = 2 + p + n_lambdas - 1; //not output the household size
    ncol = DIM * householdsize + householdsize  + 1;

    //copy data from list of matrices to C++
    for (int i = 0; i < n_lambdas; i++) {
      NumericMatrix l = lambda[i];
      lambda_columns[i] = l.ncol();
      lambdas[i] = new double[l.length()];
      std::copy(l.begin(), l.end(), lambdas[i]);
    }
    //prepare omega for group sampling, first need to transpose omega
    omegaT = transposeAndNormalize(omega);

    r = runif(data.nrow()  * ncol);
  }

  void cleanup() {
    delete [] lambda_columns;
    for (int i = 0; i < n_lambdas; i++) {
      delete [] lambdas[i];
    }
    delete [] lambdas;
  }

  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    int nHouseholds = end - begin;
    ::sampleHouseholds_imp_HHhead_at_group_level(data.begin() + begin * ncol, r.begin() +  begin * ncol,
                                                 lambdas, lambda_columns, omegaT.begin(),
                                                 phi.begin(), pi.begin(),d.begin(),
                                                 nHouseholds, householdsize, FF, SS,maxdd,p, currrentbatchbase + begin,n_lambdas);

  }
};

// [[Rcpp::export]]
IntegerMatrix sampleHH_HHhead_at_group_level(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
                                             IntegerVector d, List lambda,
                                             int currrentbatch, int nHouseholds,  int householdsize) {
  int DIM = 2 + d.length() + lambda.length() - 1;
  IntegerMatrix data(nHouseholds, DIM * householdsize + householdsize  + 1);
  HeadAtGroupLevelHHSampler worker(phi, omega, pi, d, lambda, currrentbatch*nHouseholds, householdsize, data);
  parallelFor(0, data.nrow(), worker);
  //worker.cleanup();
  return data;
}
