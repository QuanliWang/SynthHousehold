#include "mex.h"
#include <cstring>
#define DIM 8
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
	/* pointers to input matrices/vectors */
	double* data1 = mxGetPr(prhs[0]);
    double* data2 = mxGetPr(prhs[1]);
    
    int nRecordLength = DIM + 2;
    int nindividuals1 = (int)mxGetNumberOfElements(prhs[0]) / nRecordLength;
    int nindividuals2 = (int)mxGetNumberOfElements(prhs[1]) / nRecordLength;
    
    plhs[0] = mxCreateDoubleMatrix(nRecordLength,nindividuals1 + nindividuals2, mxREAL);
	double* Individuals = mxGetPr(plhs[0]);
    if (nindividuals1 > 0) {
        std::memcpy(Individuals, data1, nindividuals1*nRecordLength*sizeof(double));
        Individuals += nindividuals1*nRecordLength;
    }
    if (nindividuals2 > 0) {
        std::memcpy(Individuals, data2, nindividuals2*nRecordLength*sizeof(double));
    }
}