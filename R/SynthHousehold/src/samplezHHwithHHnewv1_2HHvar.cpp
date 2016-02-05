void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* pointers to input matrices/vectors */
	const double* const newphi = mxGetPr(prhs[0]); // pointer to newphi
	const double* const data = mxGetPr(prhs[1]); // pointer to data (transposed)
	const double* const w = mxGetPr(prhs[2]); // pointer to w
	const double* const pi = mxGetPr(prhs[3]); // pointer to pi
	const double* const S = mxGetPr(prhs[4]); // pointer to S
	const int K = (int)mxGetScalar(prhs[5]);
	const int L = (int)mxGetScalar(prhs[6]);
	const int p = (int)mxGetScalar(prhs[7]);
	const int maxdd = (int)mxGetScalar(prhs[8]);
	const int n = (int)mxGetScalar(prhs[9]);
	const double* const HHdata1 = mxGetPr(prhs[10]); // pointer to HHdata1
	const double* const lambda1 = mxGetPr(prhs[11]); // pointer to lambda1
	const double* const HHdata2 = mxGetPr(prhs[12]); // pointer to HHdata2
	const double* const lambda2 = mxGetPr(prhs[13]); // pointer to lambda2

	plhs[0] = mxCreateDoubleMatrix(n*K, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);

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
		double updatesum = 0.0;
		for (int k=0; k < K; k++) {
			double zupdateprod = 1.0;
			for (int memberindex=0; memberindex < S[h]; memberindex++){
				int base = (cumS[h]+memberindex)*p;
				double add = 0.0;
				for (int l=0; l < L; l++) {
					double phiprod = 1.0;
					for (int j=0; j < p; j++) {
						int u = (int)data[base+j]-1;
						phiprod *= newphi[maxDDtP*(k*L+l)+j*maxdd+u];
					}
					add += w[K*l+k]*phiprod;
				} // closing l++
				zupdateprod *= add;
			} // closing member++
			zupdateprod *= lambda1[HHdata_base1+k];
			zupdateprod *= lambda2[HHdata_base2+k];

			zupdateprob1[k] = pi[k]*zupdateprod;
			updatesum += zupdateprob1[k];
		} // closing k++

		for (int k=0; k < K; k++){
			*coef++ = zupdateprob1[k]/updatesum;
		}
	}
	delete [] cumS;
	delete [] zupdateprob1;
}
