#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	/* pointers to input matrices/vectors */
	const double* const datah1 = mxGetPr(prhs[0]); 
	const int number = (int)mxGetScalar(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(number, 1, mxREAL);
	double *coef = mxGetPr(plhs[0]);

	double *datah = new double[10];

	for (int m = 1; m <= number; m++){		

		datah[1] = datah1[24*(m-1)+3-1];
		datah[2] = datah1[24*(m-1)+11-1];
		datah[3] = datah1[24*(m-1)+19-1];
		datah[4] = datah1[24*(m-1)+6-1];
		datah[5] = datah1[24*(m-1)+14-1];
		datah[6] = datah1[24*(m-1)+22-1];
		datah[7] = datah1[24*(m-1)+7-1];
		datah[8] = datah1[24*(m-1)+15-1];
		datah[9] = datah1[24*(m-1)+23-1];

		coef[m-1] = 0;
	
		if(datah[7]==1 && datah[8]==2 && datah[9] ==3 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && ((datah[4] - datah[6])>=12) && (abs(datah[4]-datah[5])<=52)) // type 1
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==4 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && ((datah[4] - datah[6])>=10) && (abs(datah[4]-datah[5])<=52))// type 2
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==5 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && ((datah[6] - datah[4])>=13) && (abs(datah[4]-datah[5])<=52)) // type 3
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==6 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && ((datah[6] - datah[4])>=9) && (abs(datah[4]-datah[5])<=52))// type 4 
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==7 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && (abs(datah[4] - datah[6])<=33) && (abs(datah[4]-datah[5])<=52))// type 5
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==8 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && (abs(datah[4] - datah[6])<=33) && (abs(datah[4]-datah[5])<=52))// type 6
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==9 && (datah[1] != datah[2]) && (datah[4] >= 33) && (datah[5] >= 33) && ((datah[4] - datah[6])>=30) && (abs(datah[4]-datah[5])<=52)) // type 7
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==10 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && (abs(datah[4]-datah[5])<=52))// type 8
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==11 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && (abs(datah[4]-datah[5])<=52)) // type 9
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==2 && datah[9] ==12 && (datah[1] != datah[2]) && (datah[4] >= 17) && (datah[5] >= 17) && (abs(datah[4]-datah[5])<=52)) // type 10
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==3 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && ((datah[4] - datah[6])>=12)) // type 11
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==4 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && ((datah[4] - datah[6])>=10))// type 12
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==5 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && ((datah[6] - datah[4])>=13))// type 13
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==6 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && ((datah[6] - datah[4])>=9))// type 14
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==7 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && (abs(datah[6] - datah[4])<=33))// type 15
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==8 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12) && (abs(datah[6] - datah[4])<=33))// type 16
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==9 && (datah[4] >= 33) && ((datah[4] - datah[5])>=12) && ((datah[4] - datah[6])>=30))// type 17
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==10 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12)) // type 18
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==11 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12)) // type 19
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==3 && datah[9] ==12 && (datah[4] >= 17) && ((datah[4] - datah[5])>=12))// type 20
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==4 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10) && ((datah[4] - datah[6])>=10))// type 21
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==5 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10) && ((datah[6] - datah[4])>=13)) // type 22
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==6 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10) && ((datah[6] - datah[4])>=9)) // type 23
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==7 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10) && (abs(datah[6] - datah[4])<=33))// type 24
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==8 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10) && (abs(datah[6] - datah[4])<=33)) // type 25
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==9 && (datah[4] >= 33) && ((datah[4] - datah[5])>=10) && ((datah[4] - datah[6])>=30)) // type 26
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==10 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10))// type 27
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==11 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10))// type 28
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==4 && datah[9] ==12 && (datah[4] >= 17) && ((datah[4] - datah[5])>=10)) // type 29
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==5 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13) && ((datah[6] - datah[4])>=13)) // type 30
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==6 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13) && ((datah[6] - datah[4])>=9))// type 31
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==7 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13) && (abs(datah[6] - datah[4])<=33)) // type 32
			 coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==8 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13) && (abs(datah[6] - datah[4])<=33))// type 33
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==9 && (datah[4] >= 33) && ((datah[5] - datah[4])>=13) && ((datah[4] - datah[6])>=30))// type 34
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==10 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13))// type 35
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==11 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13)) // type 36
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==5 && datah[9] ==12 && (datah[4] >= 17) && ((datah[5] - datah[4])>=13))// type 37
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==6 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9) && ((datah[6] - datah[4])>=9))// type 38
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==7 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9) && (abs(datah[6] - datah[4])<=33))// type 39
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==8 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9) && (abs(datah[6] - datah[4])<=33)) // type 40
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==9 && (datah[4] >= 33) && ((datah[5] - datah[4])>=9) && ((datah[4] - datah[6])>=30)) // type 41
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==10 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9)) // type 42
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==11 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9))// type 43
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==6 && datah[9] ==12 && (datah[4] >= 17) && ((datah[5] - datah[4])>=9))// type 44
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==7 && (datah[4] >= 17) && (abs(datah[5] - datah[4])<=33) && (abs(datah[6] - datah[4])<=33)) // type 45
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==8 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33) && (abs(datah[4] - datah[6])<=33)) // type 46
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==9 && (datah[4] >= 8) && (abs(datah[4] - datah[5])<=33) && ((datah[4] - datah[6])>=30))// type 47
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==10 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33))// type 48
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==11 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33)) // type 49
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==7 && datah[9] ==12 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33)) // type 50
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==8 && datah[9] ==8 && (datah[4] >= 17) && (abs(datah[5] - datah[4])<=33) && (abs(datah[6] - datah[4])<=33))// type 51
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==8 && datah[9] ==9 && (datah[4] >= 33) && (abs(datah[4] - datah[5])<=33) && ((datah[4] - datah[6])>=30)) // type 52
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==8 && datah[9] ==10 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33)) // type 53
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==8 && datah[9] ==11 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33))// type 54
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==8 && datah[9] ==12 && (datah[4] >= 17) && (abs(datah[4] - datah[5])<=33)) // type 55
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==9 && datah[9] ==9 && (datah[4] >= 33) && ((datah[4] - datah[5])>=30) && ((datah[4] - datah[6])>=30))// type 56
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==9 && datah[9] ==10 && (datah[4] >= 33) && ((datah[4] - datah[5])>=30))// type 57
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==9 && datah[9] ==11 && (datah[4] >= 33) && ((datah[4] - datah[5])>=30))// type 58
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==9 && datah[9] ==12 && (datah[4] >= 33) && ((datah[4] - datah[5])>=30)) // type 59
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==10 && datah[9] ==10 && (datah[4] >= 17)) // type 60
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==10 && datah[9] ==11 && (datah[4] >= 17))// type 61
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==10 && datah[9] ==12 && (datah[4] >= 17)) // type 62
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==11 && datah[9] ==11 && (datah[4] >= 17))// type 63
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==11 && datah[9] ==12 && (datah[4] >= 17))// type 64
			coef[m-1] = 1;

		else if (datah[7]==1 && datah[8]==12 && datah[9] ==12 && (datah[4] >= 17))// type 65
			coef[m-1] = 1;
		}

		delete [] datah;
}