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
    
    plhs[0] = mxCreateDoubleMatrix(nHouseholds, 1, mxREAL);
	double* isPossible = mxGetPr(plhs[0]);
    
    int totalpossible = checkconstraints(data, isPossible,hh_size, nHouseholds);
    
    int rows = nHouseholds-totalpossible;
    plhs[1] = mxCreateDoubleMatrix(columns, rows, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(totalpossible, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(columns,totalpossible, mxREAL);
    double* newdata = mxGetPr(plhs[1]);
    double* syndata = mxGetPr(plhs[3]);
    
    double* impossible_counts = mxGetPr(plhs[2]);
    
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < nHouseholds && count2 < neededpossiblehh; i++) {
        if (isPossible[i] == 0) { //found an impossible household
            for (int j = 0; j < columns; j++) { //impossible ones
                newdata[count1*columns+j] = data[j*nHouseholds+i];
            }
            count1++;
        } else {
            *impossible_counts++ = count1;
            for (int j = 0; j < columns; j++) { //possible ones (syndata)
                syndata[count2*columns+j] = data[j*nHouseholds+i];
            }
            count2++;
        }
    }
    
    if (count1 < rows) { //rows in C are the columns in the return matrix
        //need to resize the output matrix
        void *newptr = mxRealloc(newdata, count1 * columns * sizeof(double));
        mxSetPr(plhs[1], (double *)newptr);
        mxSetN(plhs[1], count1);
    }
    
    if (count2 < totalpossible) { //truncate possible households if too many
        //need to resize the output matrix
        void *newptr = mxRealloc(syndata, count2 * columns * sizeof(double));
        mxSetPr(plhs[3], (double *)newptr);
        mxSetN(plhs[3], count2);
    }
    
    plhs[4] = mxCreateDoubleScalar(count2);
    
}