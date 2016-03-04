
#include "mex.h"
#include <math.h>
#include "sampleW.cpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *p;
    double d;
	int n,k;
	if( nrhs != 2 )
		mexErrMsgTxt( "Usage: ..." );
		
	p = mxGetPr(prhs[0]);
	d = mxGetScalar(prhs[1]);
    n = (int)(mxGetM( prhs[0] ) *  mxGetN( prhs[0])); 
    k = samplew(p,n,d);
    plhs[0] = mxCreateDoubleScalar(k);
	
}
