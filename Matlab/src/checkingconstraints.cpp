#include "utils.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* pointers to input matrices/vectors */
	const double* const datah1 = mxGetPr(prhs[0]);
	const int number = (int)mxGetScalar(prhs[1]);
    const int hh_size = (int)mxGetScalar(prhs[2]);
    
	plhs[0] = mxCreateDoubleMatrix(number, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);

	double *datah = new double[hh_size * 3 + 1];
    int record_length = hh_size * 8;
    
    //column 3, 6, 7 = sex, age and relte
    
	for (int m = 1; m <= number; m++){
        for (int j = 1; j <= hh_size; j++) {
            datah[j] = datah1[record_length*(m-1)+3-1 + (j-1) * 8];
            datah[hh_size + j] = datah1[record_length*(m-1)+6-1 + (j-1) * 8];
            datah[2 * hh_size + j] = datah1[record_length*(m-1)+7-1 + (j-1) * 8];
        }
        /*
		datah[1] = datah1[32*(m-1)+3-1]; //sex1
		datah[2] = datah1[32*(m-1)+11-1]; //sex2
		datah[3] = datah1[32*(m-1)+19-1]; //sex3
		datah[4] = datah1[32*(m-1)+25-1]; //sex4
		datah[5] = datah1[32*(m-1)+6-1]; //age1
		datah[6] = datah1[32*(m-1)+14-1]; //age2
		datah[7] = datah1[32*(m-1)+22-1]; //age3
		datah[8] = datah1[32*(m-1)+30-1]; //age4
		datah[9] = datah1[32*(m-1)+7-1]; //relate1
		datah[10] = datah1[32*(m-1)+15-1]; //relate2
		datah[11] = datah1[32*(m-1)+23-1]; //relate3
		datah[12] = datah1[32*(m-1)+31-1]; //relate4
        */
		coef[m-1] = isValid(datah, hh_size);
	}
	delete [] datah;
}