#include "mex.h"
#include <vector>
#include "sampleW.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* pointers to input matrices/vectors */
    const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
    const double* const d = mxGetPr(prhs[1]); //pointer to levels
    const double* const groupindex = mxGetPr(prhs[2]); //one based
    
    const int nIndividuals = (int)mxGetNumberOfElements(prhs[2]); // number of individuals
    const int p = (int)mxGetScalar(prhs[3]);
     double* rand = mxGetPr(prhs[4]); // pointer to random number
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
    if (mxGetNumberOfElements(prhs[4]) != nIndividuals * p) {
        mexPrintf("The number of random numbers has to be ");
        return;
    }
    
    //mexPrintf("p = %d, maxd = %d, nIndividuals = %d, maxDDtp = %d\n", p, maxdd, nIndividuals,maxDDtp);
    
    plhs[0] = mxCreateDoubleMatrix(nIndividuals, p, mxREAL);
    double *syn = mxGetPr(plhs[0]);
    
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
    
    //now sampling from each group for each individual
    double* nextrand = rand; //traverse through random numbers
    double* nextsyn = syn;
    
    for (int i = 0; i < p; i++) {
        int n = (int)d[i];
        for (int j = 0; j < nIndividuals; j++) {
        int group = int(groupindex[j])-1;
            double* cum_curentphi_j = ps[i] + group * n;
            double rn = *nextrand++;
            int k;
            for(k=0;k < n && rn>cum_curentphi_j[k];k++) //see sampleW for algorithm
                ;
            *nextsyn++ = k+1;
        }
    }
    
    //clearn up the memory
    for (int i = 0; i < p; i++) {
        delete [] ps[i];
    }
    delete [] ps;
}