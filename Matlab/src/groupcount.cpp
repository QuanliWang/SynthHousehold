
#include "mex.h"
#include <math.h>

//note the second parameter will be overwriten
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *g;  //group index
    double *d;  //data
    int n;      //number of data points
    int row, col;   //number of row and columns
	
	g = mxGetPr(prhs[0]);
    n = (int)(mxGetM( prhs[0] ) *  mxGetN( prhs[0]));
    d =mxGetPr(prhs[1]);
    row = (int)mxGetScalar(prhs[2]);
    col = (int)mxGetScalar(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(row, col, mxREAL);
    double *r = mxGetPr(plhs[0]);
    for (int i = 0; i < n; i++) {
        r[((int)g[i] -1) + ((int)d[i]-1) * row ]++;
    }
}
