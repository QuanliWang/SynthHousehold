#include <cstring>
#include "checkconstraints_imp.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
	/* pointers to input matrices/vectors */
	double* const data = mxGetPr(prhs[0]);
    int nHouseholds = (int)mxGetM(prhs[0]);
    //use the raw data instead, which has hh_size * 8 + 1 + hh_size
    int columns = (int)mxGetN(prhs[0]);
    int hh_size = (columns -1) / (DIM+1);
    
    int neededpossiblehh = (int)mxGetScalar(prhs[1]);
	int *isPossible = new int[nHouseholds];
    int totalpossible = checkconstraints(data, isPossible,hh_size, nHouseholds);
    
    int rows = nHouseholds-totalpossible;
    plhs[0] = mxCreateDoubleMatrix(columns, rows, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(totalpossible, 1, mxREAL);
    double* newdata = mxGetPr(plhs[0]);
    
    double* impossible_counts = mxGetPr(plhs[1]);
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < nHouseholds && count2 < neededpossiblehh; i++) {
        if (isPossible[i] == 0) { //found an impossible household
            //std::memcpy(newdata + count1 * columns, source, sizeof dest); not working, need data to be transposed too. Do it later.
            for (int j = 0; j < columns; j++) { //put impossible one at the begning
                newdata[count1*columns+j] = data[j*nHouseholds+i];
            }
            count1++;
        } else {
            *impossible_counts++ = count1;
            count2++;
        }
    }
    
    if (count1 < rows) { //rows in C are the columns in the return matrix
        //need to resize the output matrix
        void *newptr = mxRealloc(newdata, count1 * columns * sizeof(double));
        mxSetPr(plhs[0], (double *)newptr);
        mxSetN(plhs[0], count1);
    }
    
    delete [] isPossible;
    plhs[2] = mxCreateDoubleScalar(count2);
}