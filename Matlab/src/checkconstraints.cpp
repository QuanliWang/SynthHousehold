#include <cstring>
#include "checkconstraints_imp.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
	/* pointers to input matrices/vectors */
	double* const data = mxGetPr(prhs[0]);
    int nHouseholds = (int)mxGetM(prhs[0]);
    //use the raw data instead, which has hh_size * 8 + 1 + hh_size
    int columns = (int)mxGetN(prhs[0]);
    int hh_size = (columns -1) / (DIM+1);
    
    int TRANSPOSE = nrhs > 1? (int)mxGetScalar(prhs[1]) : 0;
    
	int *isPossible = new int[nHouseholds];
    int totalpossible = checkconstraints(data, isPossible,hh_size, nHouseholds);
    
    int rows = nHouseholds-totalpossible;
    if (TRANSPOSE) {
        plhs[0] = mxCreateDoubleMatrix(columns, rows, mxREAL);
    } else {
        plhs[0] = mxCreateDoubleMatrix(rows, columns, mxREAL);
    }
    
    plhs[1] = mxCreateDoubleMatrix(totalpossible, 1, mxREAL);
    double* newdata = mxGetPr(plhs[0]);
    double* impossible_counts = mxGetPr(plhs[1]);
    int count1 = 0;
    for (int i = 0; i < nHouseholds; i++) {
        if (isPossible[i] == 0) { //found an impossible household
            if (TRANSPOSE) {
                //std::memcpy(newdata + count1 * columns, source, sizeof dest); not working, need data to be transposed too. Do it later.
                for (int j = 0; j < columns; j++) { //put impossible one at the begning
                    newdata[count1*columns+j] = data[j*nHouseholds+i];
                }
            } else {
                for (int j = 0; j < columns; j++) { //put impossible one at the begning
                    newdata[j*rows+count1] = data[j*nHouseholds+i];
                }
            }
            
            count1++;
        } else {
            *impossible_counts++ = count1;
        }

    }
    delete [] isPossible;
}