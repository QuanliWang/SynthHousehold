#include <Rcpp.h>
using namespace Rcpp;

//#include "MersenneTwister.h"
//#include "SpecialFunctions.h"

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

#include "samplehouseholds.h"

// [[Rcpp::export]]
IntegerMatrix samplehouseholds(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
               IntegerVector d, List lambda,
               int currrentbatch, int nHouseholds,  int householdsize) {

  int FF = omega.nrow();
  int SS = omega.ncol();
  int p = d.length();
  int n_lambdas = lambda.length();
  int *lambda_columns = new int[n_lambdas];
  double **lambdas = new double*[n_lambdas];
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  //int ncol = householdsize * DIM + 1 + householdsize;
  //output data: zero-based
  //column 0: household unique index
  //column 1: member index within the household (pernum:person number?)
  //column 2 to 2+p-1: individual data
  //column 2+p: 2+p+n_lambdas-2: household level data
  //column 2+p+n_lambdas-1: household group indicator
  //column last hh_size: individual group indicator
  int DIM = 2 + p + n_lambdas - 1; //not output the household size
  int ncol = DIM * householdsize + householdsize  + 1;
  IntegerMatrix data(nHouseholds, ncol);

  //copy data from list of matrices to C++
  for (int i = 0; i < n_lambdas; i++) {
    NumericMatrix l = lambda[i];
    lambda_columns[i] = l.ncol();
    lambdas[i] = new double[l.length()];
    std::copy(l.begin(), l.end(), lambdas[i]);
  }
  //printf("in samplehouseholds\n");
  NumericVector rand = runif(nHouseholds * ncol); //at most this many
  sampleHouseholds_imp(data.begin(), rand.begin(), lambdas, lambda_columns, omega.begin(),
                   phi.begin(), pi.begin(),d.begin(),
                   nHouseholds, householdsize, FF, SS,maxdd,p, currrentbatch,n_lambdas);


  //clean up
  delete [] lambda_columns;
  for (int i = 0; i < n_lambdas; i++) {
    delete [] lambdas[i];
  }
  delete [] lambdas;
  //printf("done samplehouseholds\n");
  return data;
}

NumericMatrix prepareOmegaT(NumericMatrix omega) {
  int FF = omega.nrow();
  int SS = omega.ncol();
  NumericMatrix omegaT(SS,FF);
  //prepare omega for group sampling, first need to transpose omega
  double *currentrow = omegaT.begin();
  for (int k =0; k < FF; k++) {
    double dsum = 0.0;
    for (int l = 0; l <SS; l++) {
      //transpose first
      currentrow[l] = omega[l*FF+k];  //omegat[k * SS + l] = omega[l*FF+k];
      dsum += currentrow[l];
    }
    currentrow[0] /= dsum;
    for (int l = 1; l <SS; l++) {
      currentrow[l] = currentrow[l]/dsum + currentrow[l-1]; //normilized cum_sum
    }
    currentrow += SS;
  }
  return omegaT;
}

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
    omegaT = prepareOmegaT(omega);

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

// [[Rcpp::export]]
IntegerMatrix samplehouseholds_HHhead_at_group_level(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
                               IntegerVector d, List lambda,
                               int currrentbatch, int nHouseholds,  int householdsize) {

  int FF = omega.nrow();
  int SS = omega.ncol();
  int p = d.length();
  int n_lambdas = lambda.length();
  int *lambda_columns = new int[n_lambdas];
  double **lambdas = new double*[n_lambdas];
  int maxDDtp = phi.nrow();
  int maxdd = maxDDtp / p;

  //int ncol = householdsize * DIM + 1 + householdsize;
  //output data: zero-based
  //column 0: household unique index
  //column 1: member index within the household (pernum:person number?)
  //column 2 to 2+p-1: individual data
  //column 2+p: 2+p+n_lambdas-2: household level data
  //column 2+p+n_lambdas-1: household group indicator
  //column last hh_size: individual group indicator
  int DIM = 2 + p + n_lambdas - 1; //not output the household size
  int ncol = DIM * householdsize + householdsize  + 1;
  IntegerMatrix data(nHouseholds, ncol);

  //copy data from list of matrices to C++
  for (int i = 0; i < n_lambdas; i++) {
    NumericMatrix l = lambda[i];
    lambda_columns[i] = l.ncol();
    lambdas[i] = new double[l.length()];
    std::copy(l.begin(), l.end(), lambdas[i]);
  }
  //printf("in samplehouseholds\n");
  NumericMatrix omegaT = prepareOmegaT(omega);
  NumericVector rand = runif(nHouseholds * ncol); //at most this many
  sampleHouseholds_imp_HHhead_at_group_level(data.begin(), rand.begin(), lambdas, lambda_columns, omegaT.begin(),
                       phi.begin(), pi.begin(),d.begin(),
                       nHouseholds, householdsize, FF, SS,maxdd,p, currrentbatch*nHouseholds,n_lambdas);

  //clean up
  delete [] lambda_columns;
  for (int i = 0; i < n_lambdas; i++) {
    delete [] lambdas[i];
  }
  delete [] lambdas;
  //printf("done samplehouseholds\n");
  return data;
}
// [[Rcpp::export]]
IntegerMatrix households2individuals(IntegerMatrix data, int hh_size){

  int nHouseholds = data.ncol();

  //use the raw data instead, which has hh_size * DIM + 1 + hh_size columns (in C)
  int columns = data.nrow();
  //int hh_size = (columns - 1) / (DIM+1);
  int DIM = (columns -1) / hh_size -1;
  IntegerMatrix Individuals(DIM + 2, nHouseholds*hh_size);

  int c9 = hh_size * DIM;
  int count = 0;
  for (int i = 0; i < nHouseholds; i++) {
    int base = i * columns;
    for (int j = 0; j < hh_size; j++) {
      for (int k = 0; k < DIM;k++) {
        Individuals[count++] = data[base + j*DIM+k];
      }
      Individuals[count++] = data[base + c9];

      Individuals[count++] = data[base + c9 + 1 + j];
    }
  }
  return(Individuals);
}
