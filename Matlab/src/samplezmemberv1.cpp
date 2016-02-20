#include "mex.h"
#include <vector>
#include "sampleW.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* pointers to input matrices/vectors */
    const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
    const double* const data = mxGetPr(prhs[1]); // pointer to data (transposed)
    const double* const w = mxGetPr(prhs[2]); // pointer to w

    const double* const zHH = mxGetPr(prhs[3]); // pointer to HHdata
    const double* const serial = mxGetPr(prhs[4]); // pointer to lambda
    const double* const rand = mxGetPr(prhs[5]); // pointer to random number
    const int nIndividuals = mxGetNumberOfElements(prhs[4]); // number of individuals
    const int K = (int)mxGetM(prhs[2]);
    const int L = (int)mxGetN(prhs[2]);
    const int p = (int)mxGetM(prhs[1]);
    const int maxdd =(int) mxGetM(prhs[0]) / p;
    //mexPrintf("K = %d, L = %d, p = %d, maxd = %d, nIndividuals = %d\n", K, L, p, maxdd, nIndividuals);
    plhs[0] = mxCreateDoubleMatrix(nIndividuals, 1, mxREAL);
    double *group = mxGetPr(plhs[0]);
    int maxDDtp = maxdd*p;
    double *zupdateprob2= new double[L];
    for (int m = 0; m < nIndividuals; m++) {
        int zHHasg = zHH[int(serial[m])-1];
        int base = m*p;
        double updatesum = 0.0;
        for (int l = 0; l < L; l++) {
            double phiprod = 1.0;
            for (int j = 0; j < p; j++) {
                int u = (int)data[base+j]-1;
                phiprod *= newphi[maxDDtp*((zHHasg-1)*L+l)+j*maxdd+u];
            }
            zupdateprob2[l] = w[K*l+zHHasg-1]*phiprod;
        }
        group[m] = samplew(zupdateprob2, L, rand[m]);
    }
    delete [] zupdateprob2;
}