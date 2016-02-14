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
            datah[j] = datah1[record_length*(m-1)+3-1 + (j-1) * 8]; //sex
            datah[hh_size + j] = datah1[record_length*(m-1)+6-1 + (j-1) * 8]; //age
            datah[2 * hh_size + j] = datah1[record_length*(m-1)+7-1 + (j-1) * 8];  //relate
        }
		coef[m-1] = isValid(datah, hh_size);
	}
	delete [] datah;
}