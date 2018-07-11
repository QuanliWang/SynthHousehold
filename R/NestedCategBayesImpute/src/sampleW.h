#include <Rcpp.h>
using namespace Rcpp;
//  sampleW.h
int samplew(double *p, int n, double d);
void samplew_multi(double *p, int n, double *d,int howmany);
//this version put results into a different place
void samplew_multi2(double *p, int n, double *d, double* result,int howmany);

NumericVector checkSZ2(NumericMatrix Data_to_check, int h);
