#include "mex.h"
#include <vector>
#include "sampleW.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    
    const int DIM = 8;

    /* pointers to input matrices/vectors */
    const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
    const double* const w = mxGetPr(prhs[1]);
    const double* const pi = mxGetPr(prhs[2]);
    const double* const d = mxGetPr(prhs[3]); //pointer to levels
    double* lambda1 = mxGetPr(prhs[4]);
    double* lambda2 = mxGetPr(prhs[5]);
    
    const int currrentbatch = (int)mxGetScalar(prhs[6]);
    const int nHouseholds = (int)mxGetScalar(prhs[7]);
    const int householdsize = (int)mxGetScalar(prhs[8]);
    double* rand = mxGetPr(prhs[9]); // pointer to random number
    
    int K = (int)mxGetM(prhs[4]);
    int L = (int)mxGetN(prhs[1]);
    const int p = (int)mxGetNumberOfElements(prhs[3]);
    int lambda1_columns = (int)mxGetN(prhs[4]);
    
    //mexPrintf("p = %d, L = %d, nHouseholds = %d, currrentbatch = %d, householdsize = %d, lambda1_columns =%d\n", p, L, nHouseholds,currrentbatch,householdsize,lambda1_columns);

    int maxdd, maxDDtp, groups;
    const mwSize *dims;
    mwSize number_of_dimensions;
    dims = mxGetDimensions(prhs[0]);
    number_of_dimensions = mxGetNumberOfDimensions(prhs[0]);
    if (number_of_dimensions == 2) {
        maxDDtp = dims[0];
        groups = dims[1];
        maxdd = ((int)maxDDtp) / p;
    } else if (number_of_dimensions == 3) {
        if (dims[1] != p) {
            mexPrintf("phi's second dim has to match p");
            return;
        }
        maxdd = (int)dims[0];
        maxDDtp = (int)dims[0] * dims[1];
        groups = (int)dims[2];
    } else {
        mexPrintf("phi had to be a 2D or 3D array");
        return;
    }
    
    if (mxGetNumberOfElements(prhs[9]) != nHouseholds * (p*householdsize+householdsize+2)) {
        mexPrintf("Too few random numbers provided");
        return;
    }
    
    //mexPrintf("maxDDtp = %d, groups = %d, maxdd = %d\n", maxDDtp, groups, maxdd);
    
    int ncol = householdsize * DIM + 1 + householdsize;
    
    plhs[0] = mxCreateDoubleMatrix(nHouseholds, ncol, mxREAL);
    double *data = mxGetPr(plhs[0]);
    
    double* nextrand = rand; //traverse through random numbers
    
    int column;
    //sampling hhindexh, column: householdsize householdsize * DIM + 1
    column = (householdsize * DIM + 1) - 1; //zero-based column
    double* hhindexh = data + column * nHouseholds;
    
    double* pi_lambda2 = new double[K];
    double* currentlambdacolumn = lambda2 + (householdsize - 1 -1) * K; //column hh_size-1, addjusted to zero based
    for (int i = 0; i < K; i++) {
        pi_lambda2[i] = pi[i] * currentlambdacolumn[i];
    }
    samplew_multi2(pi_lambda2, K, nextrand,hhindexh,nHouseholds);
    nextrand += nHouseholds; //advance nHouseholds random numbers
    delete [] pi_lambda2;
    
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
    
    //prepare lambda1 for group sampling, first need to transpose lambda1
    //the code  here duplicate the lines above for w
    double* lambda1t = new double[K * lambda1_columns];
    currentrow = lambda1t;
    for (int k =0; k < K; k++) {
        double dsum = 0.0;
        for (int l = 0; l <lambda1_columns; l++) {
            //transpose first
            currentrow[l] = lambda1[l*K+k];
            dsum += currentrow[l];
        }
        currentrow[0] /= dsum;
        for (int l = 1; l <lambda1_columns; l++) {
            currentrow[l] = currentrow[l]/dsum + currentrow[l-1]; //normilized cum_sum
        }
        currentrow += lambda1_columns;
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
    
    for (int i = 0; i < householdsize; i++) {
        columns[i] = data + ((i+1) * DIM -1) * nHouseholds; //zero-based column
    }
    
    for (int i = 0; i < nHouseholds; i++) {
        int group = (int)hhindexh[i]-1;
        double* currentp = lambda1t + group * lambda1_columns;
        double rn = *nextrand++;
        int k;
        for(k=0;k < lambda1_columns && rn>currentp[k];k++) //see sampleW for algorithm
            ;
        for (int j = 0; j < householdsize; j++) {
            columns[j][i] = k + 1;
        }
    }

    delete [] wt;
    delete [] lambda1t;
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
                double currentphi = newphi[maxDDtp*j+i*maxdd+k];
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
    
    double** datacolumns = new double*[7]; //7 variables to genetate data
    int* groupindex = new int[nHouseholds];
    //now generate individual data (column 2 through 6, zero-based )
    for (int hh =0; hh < householdsize; hh++) {
        //set columns
        for (int i =0; i < 7; i++) {
            datacolumns[i] = data + (hh * DIM + i) * nHouseholds; //zero-based column
        }
        //set groupindex
        //groupindex = (hhindexh-1)*L+data_to_check(:,hh_size * 8 +1 + hh);
        double* hh_column = data + ((householdsize * DIM + 1) + hh) * nHouseholds;
        for (int i = 0; i < nHouseholds; i++) {
            groupindex[i] = (hhindexh[i]-1)*L + hh_column[i];
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
        int houseIndex = currrentbatch *nHouseholds;
        for (int j = 0; j < nHouseholds; j++) {
            datacolumns[0][j] = houseIndex + j + 1; //i based houseIndex
            datacolumns[1][j] = hh+1;
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