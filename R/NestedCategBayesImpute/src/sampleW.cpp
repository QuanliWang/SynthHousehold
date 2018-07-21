#include <Rcpp.h>
using namespace Rcpp;
#include "sampleW.h"
#include <algorithm>


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
    delete [] myw;
    return std::distance(myw, std::lower_bound(myw, myw+n,d)) + 1;
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
        d[h] = std::distance(myw, std::lower_bound(myw, myw+n,d[h])) + 1;
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
        result[h] = std::distance(myw, std::lower_bound(myw, myw+n,d[h])) + 1;
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

