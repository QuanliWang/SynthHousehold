#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* pointers to input matrices/vectors */
	const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
	const double* const data = mxGetPr(prhs[1]); // pointer to data (transposed)
	const double* const w = mxGetPr(prhs[2]); // pointer to w
	const int K = (int)mxGetScalar(prhs[3]);
	const int L = (int)mxGetScalar(prhs[4]);
	const int p = (int)mxGetScalar(prhs[5]);
	const int maxdd = (int)mxGetScalar(prhs[6]);
	const int n_s = (int)mxGetScalar(prhs[7]);
	const double* const zHH = mxGetPr(prhs[8]); // pointer to HHdata
	const double* const serial = mxGetPr(prhs[9]); // pointer to lambda
	
	plhs[0] = mxCreateDoubleMatrix(n_s*L, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);
	int maxDDtp = maxdd*p;
	double *zupdateprob2= new double[L];
	for (int m = 0; m < n_s; m++) {
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
				updatesum += zupdateprob2[l];
		}
		for (int l = 0; l < L; l++){
			*coef++ = zupdateprob2[l]/updatesum;
		}
	}
	delete [] zupdateprob2;
}