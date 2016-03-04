
#include "mex.h"
#include <math.h>
#include "sampleW.cpp"

//note the second parameter will be overwriten
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *p;
    double *d;
	int n,howmany;
	if( nrhs != 2 )
		mexErrMsgTxt( "Usage: ..." );
		
	p = mxGetPr(prhs[0]);
    n = (int)(mxGetM( prhs[0] ) *  mxGetN( prhs[0])); 
    
    d = mxGetPr(prhs[1]);
    howmany = (int)(mxGetM( prhs[1] ) *  mxGetN( prhs[1]));
    
    samplew_multi(p,n,d,howmany);
    plhs[0] = mxCreateDoubleMatrix(howmany, 1, mxREAL);
    double *C = mxGetPr(plhs[0]);
    for (int k =0; k<howmany;k++) {
        C[k] = d[k];
    }
}
