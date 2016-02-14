
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
void samplew_multi(double *p, int n, double *d,int howmany) {
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
    for (int h=0; h < howmany; h++) {
        for(k=0;k < n && d[h]>myw[k];k++)
            ;
        d[h] = k+1;
    }
    delete [] myw;
}

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
