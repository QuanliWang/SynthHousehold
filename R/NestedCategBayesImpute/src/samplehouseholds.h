//  sampleHouseholds.h
#include <Rcpp.h>
using namespace Rcpp;
void sampleHouseholds_imp(int* data, double* rand,  double** lambda, int* lambda_columns, double* omega, double* phi,
                      double *pi, int* d,int nHouseholds, int householdsize, int FF,int SS,
                      int maxdd, int p, int currrentbatch, int n_lambdas);
void sampleHouseholds_imp_HHhead_at_group_level(int* data, double* rand,  double** lambda, int* lambda_columns, double* omega, double* phi,
                          double *pi, int* d,int nHouseholds, int householdsize, int FF,int SS,
                          int maxdd, int p, int currrentbatch, int n_lambdas);
IntegerMatrix samplehouseholds_HHhead_at_group_level(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
                                                     IntegerVector d, List lambda,
                                                     int currrentbatch, int nHouseholds,  int householdsize);
IntegerMatrix samplehouseholds(NumericMatrix phi, NumericMatrix omega, NumericVector pi,
                               IntegerVector d, List lambda,
                               int currrentbatch, int nHouseholds,  int householdsize);
IntegerMatrix households2individuals(IntegerMatrix data, int hh_size);
