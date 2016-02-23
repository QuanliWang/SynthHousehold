#include "checkconstraints_imp.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
	/* pointers to input matrices/vectors */
	double* const data = mxGetPr(prhs[0]);
    int nHouseholds = (int)mxGetM(prhs[0]);
    //use the raw data instead, which has hh_size * 8 + 1 + hh_size
    int hh_size = (int) (mxGetN(prhs[0]) -1)  / (DIM+1);
     
	plhs[0] = mxCreateDoubleMatrix(nHouseholds, 1, mxREAL);
	double *isPossible = mxGetPr(plhs[0]);
    
    checkconstraints(data, isPossible,hh_size, nHouseholds);
    
}