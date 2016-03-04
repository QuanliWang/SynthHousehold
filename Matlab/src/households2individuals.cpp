#include "mex.h"
#include <cstring>
#define DIM 8
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
	/* pointers to input matrices/vectors */
	double* data = mxGetPr(prhs[0]);
    int nHouseholds = (int)mxGetN(prhs[0]);
    //use the raw data instead, which has hh_size * 8 + 1 + hh_size
    int columns = (int)mxGetM(prhs[0]);
    int hh_size = (columns -1) / (DIM+1);
    
    plhs[0] = mxCreateDoubleMatrix(DIM + 2,nHouseholds*hh_size, mxREAL);
	double* Individuals = mxGetPr(plhs[0]);
    
    int bytes = DIM * sizeof(double);
    int c9 = hh_size * DIM;
    for (int i = 0; i < nHouseholds; i++) {
        double column9 = data[c9];
        for (int j = 0; j < hh_size; j++) {
            std::memcpy(Individuals, data + j * DIM, bytes); //the first DIM columns
            Individuals += DIM;
            *Individuals++ = data[c9];
            *Individuals++ = data[c9+1+j];
        }
        data += columns; //move data to the next household
    }
}