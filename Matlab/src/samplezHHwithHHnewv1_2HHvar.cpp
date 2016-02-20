#include "mex.h"
#include "sampleW.cpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* pointers to input matrices/vectors */
	const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
	const double* const data = mxGetPr(prhs[1]); // pointer to data (transposed)
	const double* const w = mxGetPr(prhs[2]); // pointer to w
	const double* const pi = mxGetPr(prhs[3]); // pointer to pi
	const double* const S = mxGetPr(prhs[4]); // pointer to S
    const double* const HHdata1 = mxGetPr(prhs[5]); // pointer to HHdata1
    const double* const lambda1 = mxGetPr(prhs[6]); // pointer to lambda1
    const double* const HHdata2 = mxGetPr(prhs[7]); // pointer to HHdata2
    const double* const lambda2 = mxGetPr(prhs[8]); // pointer to lambda2
    const double* const rand = mxGetPr(prhs[9]); // pointer to random number
    
    const int p = (int)mxGetM(prhs[1]);
    const int nIndividuals = (int)mxGetN(prhs[1]); // number of individuals
    const int K = (int)mxGetM(prhs[2]);
    const int L = (int)mxGetN(prhs[2]);
    const int maxdd =(int) mxGetM(prhs[0]) / p;
    const int n = (int)mxGetNumberOfElements(prhs[4]);
    
    //mexPrintf("K = %d, L = %d, p = %d, maxd = %d, nIndividuals = %d, n=%d\n", K, L, p, maxdd, nIndividuals,n);

	
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nIndividuals, 1, mxREAL);
    
    int count = 0;
	double *group = mxGetPr(plhs[0]);
    double *indi = mxGetPr(plhs[1]);

	double *zupdateprob1 = new double[K];
	int *cumS = new int[n];
	cumS[0] = 0;
	for (int i = 1; i < n; i++) {
		cumS[i] = cumS[i-1] + (int)S[i-1];
	}
	int maxDDtP = maxdd*p;
	for (int h = 0; h < n; h++) {
		int HHdata_base1 = (HHdata1[h]-1)*K;
		int HHdata_base2 = (HHdata2[h]-1)*K;
		for (int k=0; k < K; k++) {
			double zupdateprod = 1.0;
			for (int memberindex=0; memberindex < S[h]; memberindex++){
				int base = (cumS[h]+memberindex)*p; //base for data
				double add = 0.0;
				for (int l=0; l < L; l++) {
					double phiprod = 1.0;
                    int phi_base = (int)(maxDDtP*(k*L+l));
					for (int j=0; j < p; j++) {
						int u = (int)data[base+j]-1;
						phiprod *= newphi[phi_base+j*maxdd+u];
					}
					add += w[K*l+k]*phiprod;
				} // closing l++
				zupdateprod *= add;
			} // closing member++
			zupdateprod *= lambda1[HHdata_base1+k];
			zupdateprod *= lambda2[HHdata_base2+k];

			zupdateprob1[k] = pi[k]*zupdateprod;
		} // closing k++
        group[h] = samplew(zupdateprob1, K, rand[h]);
        for (int m=0; m < S[h];m++) {
            indi[count++] = group[h];
        }
	}
	delete [] cumS;
	delete [] zupdateprob1;
}