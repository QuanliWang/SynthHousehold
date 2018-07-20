#include <Rcpp.h>
using namespace Rcpp;

#include "sampleW.h"
#include "samplehouseholds.h"
#include <cstdio>
#include "utils.h"

void sampleHHindex(double** lambda, int n_lambdas, int householdsize, double* pi, int FF, double* nextrand, int* hhindexh, int nHouseholds) {
  double* currentlambdacolumn = lambda[n_lambdas-1] + (householdsize - 1) * FF; //column hh_size-1, addjusted to zero based
  double* pi_lambda_last = new double[FF];
  //note that now household size start from 1, instead of 2
  for (int i = 0; i < FF; i++) {
    pi_lambda_last[i] = pi[i] * currentlambdacolumn[i];
  }
  samplew_multi2(pi_lambda_last, FF, nextrand,hhindexh,nHouseholds);
  delete [] pi_lambda_last;
}

//output: householdsize of data columns, starting at column base (0-based)
void sampleIndivMemberIndex(int* data,int* hhindexh, int nHouseholds, int base, int householdsize,
                            double* omegat, int SS, double* nextrand) {
  int** columns = new int*[householdsize];
  double* currentrand = nextrand;
  for (int j = 0; j < householdsize; j++) {
    columns[j] = data + (base + j) * nHouseholds; //zero-based column
  }
  for (int i = 0; i < nHouseholds; i++) {
    int group = hhindexh[i]-1;
    double* currentp = omegat + group * SS;

    for (int j = 0; j < householdsize; j++) {
      double rn = *currentrand++;
      int k;
      for(k=0;k < SS && rn>currentp[k];k++) //see sampleW for algorithm, currentp allreay cumulitive
        ;
      columns[j][i] = k + 1;
    }
  }
}

void sampleHH_level_data(int* data, int* hhindexh, double* nextrand, int nHouseholds, int DIM,  double* lambda, int n_lambda,
                         int FF, int householdsize,  int p, int g) {
  double* currentrand = nextrand;
  int** columns = new int*[householdsize];
  for (int j = 0; j < householdsize; j++) {
    columns[j] = data + (j * DIM + 2 + p + g) * nHouseholds; //zero-based column
  }
  //prepare lambdas for for group sampling, first need to transpose lambda
  //the code  here duplicate the lines above for omega
  double* lambda_t = new double[FF * n_lambda];
  transposeAndNormalize(lambda, FF, n_lambda,  lambda_t);

  for (int i = 0; i < nHouseholds; i++) {
    int group = hhindexh[i]-1;
    double* currentp = lambda_t + group * n_lambda;
    double rn = *currentrand++;
    int k;
    for(k=0;k < n_lambda && rn>currentp[k];k++) //see sampleW for algorithm
      ;
    columns[0][i] = k + 1;
  }
  for (int j = 1; j < householdsize; j++) {
    std::copy(columns[0], columns[0] + nHouseholds, columns[j]);
  }
  delete [] lambda_t;
  delete [] columns;
}

//output ps, which needs to be deallocated from calling function
void preparePhis(double** ps, double* phi, int* d, int maxdd, int p, int FF, int SS) {
  int maxDDtp = maxdd * p;
  int groups = FF * SS;
  for (int i = 0; i < p; i++) {
    int currentd = d[i];
    ps[i] = new double[currentd*groups];
    double* currentp = ps[i];
    for (int j = 0; j <groups; j++) { //for each group/cluster
      double dsum = 0.0; //copy phi's and normalize them
      int base = maxDDtp*j+i*maxdd;
      for (int k = 0; k <currentd; k++) {
        currentp[k] =  phi[base+k];
        dsum += currentp[k];
      }
      currentp[0] /= dsum; //convert p-values to cum-p-values
      for (int k = 1; k <currentd; k++) {
        currentp[k] = currentp[k]/dsum + currentp[k-1];
      }
      currentp += currentd; //advance pointer to next group
    }
  }
}
void sampleIndiv_level_data(int* data, int* hhindexh, double* nextrand, int nHouseholds,
                            double** ps, int* d, int p, int SS,int householdsize, int DIM,
                            int currrentbatchbase) {
  double* currentrand = nextrand;
  int** datacolumns = new int*[2+p]; // 2+p variables to genetate data
  int* groupindex = new int[nHouseholds];
  //now generate individual data (column 2 through 2+p-1, zero-based )
  for (int hh =0; hh < householdsize; hh++) {
    //set columns
    for (int i =0; i < 2+p ; i++) {
      datacolumns[i] = data + (hh * DIM + i) * nHouseholds; //zero-based column
    }
    int houseIndex = currrentbatchbase + 1; //one based houseIndex
    for (int j = 0; j < nHouseholds; j++) {
      datacolumns[0][j] = houseIndex + j;
    }
    std::fill(datacolumns[1], datacolumns[1] + nHouseholds, hh + 1);

    //set groupindex
    int* hh_column = data + ((householdsize * DIM + 1) + hh) * nHouseholds;
    for (int i = 0; i < nHouseholds; i++) {
      groupindex[i] = (hhindexh[i]-1)*SS + hh_column[i];
    }

    for (int i = 0; i < p; i++) {
      int n = d[i];
      for (int j = 0; j < nHouseholds; j++) {
        int group = int(groupindex[j])-1;
        double* cum_curentphi_j = ps[i] + group * n;
        double rn = *currentrand++;
        int k;
        for(k=0;k < n && rn>cum_curentphi_j[k];k++) //see sampleW for algorithm
          ;
        datacolumns[i+2][j] = k+1; //start at column 2, zero-based
      }
    }

  }
  delete [] datacolumns;
  delete [] groupindex;
}

void sampleHouseholds_imp_HHhead_at_group_level(int* data, double* rand,
                                                double** lambda, int* lambda_columns,
                                                double* omegat, double* phi,
                                                double *pi, int* d,int nHouseholds,
                                                int householdsize, int FF,int SS,
                                                int maxdd, int p,
                                                int currrentbatchbase,int n_lambdas) {

  //number of columns in the final output
  int DIM = 2 + p + n_lambdas - 1;  //total number of variables
  double* nextrand = rand; //traverse through random numbers

  //sampling hhindexh, column: householdsize * DIM + 1 (one-based)
  int column = (householdsize * DIM + 1) - 1; //zero-based column
  int* hhindexh = data + column * nHouseholds;
  sampleHHindex(lambda, n_lambdas, householdsize, pi, FF, nextrand, hhindexh, nHouseholds);
  nextrand += nHouseholds; //advance nHouseholds random numbers

  //now sampling from each group for each individual memberindexhh
  //do random samples for the same probs at the same time
  int base = (householdsize * DIM + 1);
  sampleIndivMemberIndex(data, hhindexh, nHouseholds, base, householdsize, omegat, SS, nextrand);
  nextrand += householdsize * nHouseholds;

  //generate household level data
  for (int g = 0; g < n_lambdas-1; g++) {
    sampleHH_level_data(data, hhindexh, nextrand, nHouseholds, DIM, lambda[g], lambda_columns[g], FF, householdsize, p, g);
    nextrand += nHouseholds;
  }

  //extract p values for each individual variable
  double** ps = new double*[p];
  preparePhis(ps, phi, d, maxdd, p, FF, SS);
  sampleIndiv_level_data(data, hhindexh, nextrand, nHouseholds, ps, d, p, SS, householdsize, DIM, currrentbatchbase);
  nextrand += p * nHouseholds; //for record keeping only

  //clearn up the memory
  for (int i = 0; i < p; i++) {
    delete [] ps[i];
  }
  delete [] ps;
}
