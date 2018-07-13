#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"
#include "MersenneTwister.h"
#include "SpecialFunctions.h"

//  sampleW.cpp
//  Created by Quanli Wang on 2/20/16.
int samplew(double *p, int n, double d) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }

    for(k=0;k < n && d>myw[k];k++)
        ;
    delete [] myw;
    if (k == n) {k = n-1;}
     return k+1;
}

void samplew_multi(double *p, int n, double *d,int howmany) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    for (int h=0; h < howmany; h++) {
        for(k=0;k < n && d[h]>myw[k];k++)
            ;
        if (k == n) {k = n-1;}
        d[h] = k+1;

    }
    delete [] myw;
}

//this version put results into a different place
void samplew_multi2(double *p, int n, double *d, int* result,int howmany) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    if (dsum <=0 ) {dsum =1;}
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    for (int h=0; h < howmany; h++) {
        for(k=0;k < n && d[h]>myw[k];k++)
            ;
        if (k == n) {k = n-1;}
        result[h] = k+1;
    }
    delete [] myw;
}

// [[Rcpp::export]]
IntegerVector sampleW_multi(NumericVector p, NumericVector d) {
  int howmany = d.length();
  IntegerVector samples(howmany);
  int n = p.length();
  samplew_multi2(p.begin(), n, d.begin(), samples.begin(), howmany);
  return samples;
}

// [[Rcpp::export]]
NumericVector gammarand(int n, double shape, double rate) {
  MTRand mt;
  mt.seed();
  vector<double> result;
  SpecialFunctions::gammarand(shape,1.0 /rate,n,mt,result);
  NumericVector r(result.begin(),result.end());
  return r;
}

// [[Rcpp::export]]
NumericMatrix samplePhi(IntegerMatrix counts) {
  NumericMatrix result(counts.rows(),counts.cols());
  MTRand mt;
  mt.seed();
  for (int i = 0; i < counts.length(); i++) {
    result[i] = SpecialFunctions::gammarand(1 + counts[i], 1, mt);
  }

  return result;
}
