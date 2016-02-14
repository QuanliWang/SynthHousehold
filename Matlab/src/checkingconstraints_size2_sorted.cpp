#include "mex.h"
#include <vector>
#include <cmath>
inline bool IsHead(double relate, double age) {
    return (relate == 1 && age >=17);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* pointers to input matrices/vectors */
	const double* const datah1 = mxGetPr(prhs[0]); 
	const int number = (int)mxGetScalar(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(number, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);

	double *datah = new double[10];
    //1 = head/householder, 2 = spouse, 3 = child, 4 = child-in-law, 5 = parent, 6 = parent-in- law, 7 = sibling, 8 = sibling-in-law,
    //9 = grandchild, 
    //10 = other relatives,
    //11 = partner, friend, visitor,
    //12 = other non-relatives
	for (int m = 1; m <= number; m++){		

		datah[1] = datah1[16*(m-1)+3-1];
		datah[2] = datah1[16*(m-1)+11-1];
		datah[3] = datah1[16*(m-1)+6-1];
		datah[4] = datah1[16*(m-1)+14-1];
		datah[5] = datah1[16*(m-1)+7-1];
		datah[6] = datah1[16*(m-1)+15-1];

		coef[m-1] = 0;
        if (!IsHead(datah[5], datah[3])) { continue; }
        //head/spouse
		if (datah[6]==2 && (datah[1] != datah[2]) && (datah[4] >= 17) && (std::abs(datah[3]-datah[4])<=52))  // type 1
			coef[m-1] = 1;
        
        //head/child
		else if (datah[6]==3 && ((datah[3] - datah[4]) >= 12)) // type 2
			coef[m-1] = 1;
        
        //head/child-in-law
		else if (datah[6]==4 && ((datah[3] - datah[4]) >= 10)) // type 3
			coef[m-1] = 1;

        //head/grand parent
		else if (datah[6]==5 && ((datah[4] - datah[3]) >= 13)) // type 4
			coef[m-1] = 1;
        
        //head/parent-in-law
		else if (datah[6]==6 && ((datah[4] - datah[3]) >= 9))// type 5
			coef[m-1] = 1;
        
        //head/sibling
		else if (datah[6]==7 && (std::abs(datah[3] - datah[4]) <= 33)) // type 6
			coef[m-1] = 1;
        
        //head/sibling-in-law
		else if (datah[6]==8 && (std::abs(datah[3] - datah[4]) <= 33))// type 7
			coef[m-1] = 1;
        
        //head/grand child
		else if (datah[6]==9 && ((datah[3] - datah[4]) >= 30))// type 8
			coef[m-1] = 1;
        
        //head/other relatives
		else if (datah[6]==10)// type 9
			coef[m-1] = 1;
        
        //head/partner, friend, visitor
		else if (datah[6]==11)// type 10
			coef[m-1] = 1;
        
        //head/other non-relatives
		else if (datah[6]==12) // type 11
			coef[m-1] = 1;
		}

		delete [] datah;
}