#include "mex.h"
#include "samplehouseholds_imp.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* pointers to input matrices/vectors */
    double* phi = mxGetPr(prhs[0]); // pointer to newphi
    double* w = mxGetPr(prhs[1]);
    double* pi = mxGetPr(prhs[2]);
    double* d = mxGetPr(prhs[3]); //pointer to levels
    double* lambda1 = mxGetPr(prhs[4]);
    double* lambda2 = mxGetPr(prhs[5]);
    
    int currrentbatch = (int)mxGetScalar(prhs[6]);
    int nHouseholds = (int)mxGetScalar(prhs[7]);
    int householdsize = (int)mxGetScalar(prhs[8]);
    double* rand = mxGetPr(prhs[9]); // pointer to random number
    
    int K = (int)mxGetM(prhs[4]);
    int L = (int)mxGetN(prhs[1]);
    int p = (int)mxGetNumberOfElements(prhs[3]);
    int lambda1_columns = (int)mxGetN(prhs[4]);
    
    //mexPrintf("p = %d, L = %d, nHouseholds = %d, currrentbatch = %d, householdsize = %d, lambda1_columns =%d\n", p, L, nHouseholds,currrentbatch,householdsize,lambda1_columns);

    int maxdd, maxDDtp, groups;
    const mwSize *dims;
    mwSize number_of_dimensions;
    dims = mxGetDimensions(prhs[0]);
    number_of_dimensions = mxGetNumberOfDimensions(prhs[0]);
    if (number_of_dimensions == 2) {
        maxDDtp = dims[0];
        maxdd = ((int)maxDDtp) / p;
    } else if (number_of_dimensions == 3) {
        if (dims[1] != p) {
            mexPrintf("phi's second dim has to match p");
            return;
        }
        maxdd = (int)dims[0];
        maxDDtp = (int)dims[0] * dims[1];
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
    //end of matlab related
    
    sampleHouseholds(data, rand, lambda1, lambda2, w, phi, pi,d,
                     nHouseholds, householdsize, K, L,maxdd,p,lambda1_columns, currrentbatch);
    }