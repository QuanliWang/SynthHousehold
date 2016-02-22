#include "utils.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    const int DIM = 8;
    const int COL = 3;
	/* pointers to input matrices/vectors */
	const double* const datah1 = mxGetPr(prhs[0]);
    const int number = (int)mxGetM(prhs[0]);
    //use the raw data instead, which has hh_size * 8 + 1 + hh_size
    const int hh_size = (int) (mxGetN(prhs[0]) -1)  / (DIM+1);
     
	plhs[0] = mxCreateDoubleMatrix(number, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);

	double *datah = new double[hh_size * 3 + 1];
    
    //column 3, 6, 7 = sex, age and relte
    int column[COL]; column[0] = 3; column[1] = 6; column[2] = 7;
    
	for (int m = 1; m <= number; m++){
        for (int j = 1; j <= hh_size; j++) {
            for (int k = 0; k < COL; k++) {
                datah[k * hh_size + j] = datah1[((j-1) * 8 + column[k] -1) * number + (m-1)];
            }
        }
		coef[m-1] = isValid(datah, hh_size);
	}
    
	delete [] datah;
}