
#include "mex.h"
#include <math.h>
int samplew(double *p, int n, double d) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    
    for(k=0;k < n && d>myw[k];k++)
        ;
    delete [] myw;
    return k+1;
}
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
