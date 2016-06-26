#include "sampleW.h"
#include "samplehouseholds.h"
#include <cstdio>

void sampleHouseholds_imp(double* data, double* rand,  double** lambda, int* lambda_columns, double* w, double* phi,
                          double *pi, double* d,int nHouseholds, int householdsize, int K,int L,
                          int maxdd, int p, int currrentbatch,int n_lambdas) {

    //number of columns in the final output
    int groups = K * L;
    int maxDDtp = maxdd * p;
    double* nextrand = rand; //traverse through random numbers

    int column;
    int DIM = 2 + p + n_lambdas - 1;
    //sampling hhindexh, column: householdsize * DIM + 1
    column = (householdsize * DIM + 1) - 1; //zero-based column
    double* hhindexh = data + column * nHouseholds;

    double* pi_lambda_last = new double[K];
    double* currentlambdacolumn = lambda[n_lambdas-1] + (householdsize - 1 -1) * K; //column hh_size-1, addjusted to zero based
    for (int i = 0; i < K; i++) {
      pi_lambda_last[i] = pi[i] * currentlambdacolumn[i];
    }
    samplew_multi2(pi_lambda_last, K, nextrand,hhindexh,nHouseholds);
    nextrand += nHouseholds; //advance nHouseholds random numbers
    delete [] pi_lambda_last;

    //prepare w for group sampling, first need to transpose w
    double* wt = new double[K * L];
    double *currentrow = wt;
    for (int k =0; k < K; k++) {
        double dsum = 0.0;
        for (int l = 0; l <L; l++) {
            //transpose first
            currentrow[l] = w[l*K+k];  //wt[k * L + l] = w[l*K+k];
            dsum += currentrow[l];
        }
        currentrow[0] /= dsum;
        for (int l = 1; l <L; l++) {
            currentrow[l] = currentrow[l]/dsum + currentrow[l-1]; //normilized cum_sum
        }
        currentrow += L;
    }

    //now sampling from each group for each individual
    //memberindexhh
    //do random samples for the same probs at the same time
    double** columns = new double*[householdsize];
    for (int i = 0; i < householdsize; i++) {
        columns[i] = data + ((householdsize * DIM + 1) + i) * nHouseholds; //zero-based column
    }
    for (int j = 0; j < householdsize; j++) {
        for (int i = 0; i < nHouseholds; i++) {
            int group = (int)hhindexh[i]-1;
            double* currentp = wt + group * L;
            double rn = *nextrand++;
            int k;
            for(k=0;k < L && rn>currentp[k];k++) //see sampleW for algorithm
                ;
            columns[j][i] = k + 1;
        }
    }

    //generate household level data
    for (int g = 0; g < n_lambdas-1; g++) {
      //printf("g = %d\n", g);
      //prepare lambdas for for group sampling, first need to transpose lambda
      //the code  here duplicate the lines above for w
      double* lambda_t = new double[K * lambda_columns[g]];
      currentrow = lambda_t;
      for (int k =0; k < K; k++) {
        double dsum = 0.0;
        for (int l = 0; l <lambda_columns[g]; l++) {
          //transpose first
          currentrow[l] = lambda[g][l*K+k];
          dsum += currentrow[l];
        }
        currentrow[0] /= dsum;
        for (int l = 1; l <lambda_columns[g]; l++) {
          currentrow[l] = currentrow[l]/dsum + currentrow[l-1]; //normilized cum_sum
        }
        currentrow += lambda_columns[g];
      }

      for (int i = 0; i < householdsize; i++) {
          columns[i] = data + (i * DIM + 2 + p + g) * nHouseholds; //zero-based column
      }

      for (int i = 0; i < nHouseholds; i++) {
          int group = (int)hhindexh[i]-1;
          double* currentp = lambda_t + group * lambda_columns[g];
          double rn = *nextrand++;
          int k;
          for(k=0;k < lambda_columns[g] && rn>currentp[k];k++) //see sampleW for algorithm
              ;
          for (int j = 0; j < householdsize; j++) {
              columns[j][i] = k + 1;
          }
      }
      delete [] lambda_t;
    }
    delete [] wt;
    delete [] columns;

    //extract p values for each variable
    double** ps = new double*[p];
    for (int i = 0; i < p; i++) {
        int currentd = (int)d[i];
        ps[i] = new double[currentd*groups];

        double* currentp = ps[i];
        for (int j = 0; j <groups; j++) { //for each group/cluster
            double dsum = 0.0; //copy phi's and normalize them
            for (int k = 0; k <currentd; k++) {
                double currentphi = phi[maxDDtp*j+i*maxdd+k];
                currentp[k] = currentphi;
                dsum += currentphi;
            }
            currentp[0] /= dsum; //convert p-values to cum-p-values
            for (int k = 1; k <currentd; k++) {
                currentp[k] = currentp[k]/dsum + currentp[k-1];
            }
            currentp += currentd; //advance pointer to next group
        }
    }

    double** datacolumns = new double*[2+p]; // 2+p variables to genetate data
    int* groupindex = new int[nHouseholds];
    //now generate individual data (column 2 through 2+p-1, zero-based )
    for (int hh =0; hh < householdsize; hh++) {
        //set columns
        for (int i =0; i < 2+p ; i++) {
            datacolumns[i] = data + (hh * DIM + i) * nHouseholds; //zero-based column
        }
        //set groupindex
        double* hh_column = data + ((householdsize * DIM + 1) + hh) * nHouseholds;
        for (int i = 0; i < nHouseholds; i++) {
            groupindex[i] = (hhindexh[i]-1)*L + hh_column[i];
        }

        int houseIndex = currrentbatch *nHouseholds;
        for (int j = 0; j < nHouseholds; j++) {
          datacolumns[0][j] = houseIndex + j + 1; //one based houseIndex
          datacolumns[1][j] = hh+1;
        }

        for (int i = 0; i < p; i++) {
            int n = (int)d[i];
            for (int j = 0; j < nHouseholds; j++) {
                int group = int(groupindex[j])-1;
                double* cum_curentphi_j = ps[i] + group * n;
                double rn = *nextrand++;
                int k;
                for(k=0;k < n && rn>cum_curentphi_j[k];k++) //see sampleW for algorithm
                    ;
                datacolumns[i+2][j] = k+1; //start at column 2, zero-based
            }
        }


    }
    delete [] datacolumns;
    delete [] groupindex;

    //clearn up the memory
    for (int i = 0; i < p; i++) {
        delete [] ps[i];
    }
    delete [] ps;

}
