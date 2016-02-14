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

	double *datah = new double[13];
	for (int m = 1; m <= number; m++){		

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
       
		coef[m-1] = 0;
        if (!IsHead(datah[9], datah[5])) { continue; }
		//type1
		if(datah[10]==2 && datah[11] ==3 && datah[12]==3 && datah[1] != datah[2] && datah[6] >= 17 && (datah[5] - datah[7])>=12 && (datah[5] - datah[8])>=12  && std::abs(datah[5]-datah[6])<=52)
			coef[m-1] = 1;
		//type2
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==4 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && ((datah[5] - datah[8])>=10)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type3
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==5  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && ((datah[8] - datah[5])>=13) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type4
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==6  && (datah[1] != datah[2]) &&  (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && ((datah[8] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type5
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==7  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type6
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==8  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type7
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=12) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type8
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type9
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type10
		else if (datah[10]==2 && datah[11] ==3 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=12) && (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type11
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==4 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=10)&& (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type12
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==5 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=13)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type13
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==6 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=9)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type14
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==7 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type15
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==8 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type16
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==9 && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=30)&& (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type17
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==10 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type18
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==11 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type19
		else if (datah[10]==2 && datah[11] ==4 && datah[12]==12 && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[5] - datah[7])>=10) && (std::abs(datah[5]-datah[6])<=52))
			coef[m-1] = 1;
		//type20
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==5  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=13) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type21
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==6  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type22
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==7  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type23
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==8  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type24
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[7] - datah[5])>=13) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type25
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type26
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type27
		else if (datah[10]==2 && datah[11] ==5 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=13) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type28
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==6  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && ((datah[8] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type29
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==7  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type30
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==8  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type31
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[7] - datah[5])>=9) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type32
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type33
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type34
		else if (datah[10]==2 && datah[11] ==6 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && ((datah[7] - datah[5])>=9) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type35
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==7  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type36
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==8  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type37
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type38
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type39
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type40
		else if (datah[10]==2 && datah[11] ==7 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type41
		else if (datah[10]==2 && datah[11] ==8 && datah[12]==8  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type42
		else if (datah[10]==2 && datah[11] ==8 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type43
		else if (datah[10]==2 && datah[11] ==8 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type44
		else if (datah[10]==2 && datah[11] ==8 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type45
		else if (datah[10]==2 && datah[11] ==8 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type46
		else if (datah[10]==2 && datah[11] ==9 && datah[12]==9  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=30) && ((datah[5] - datah[8])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type47
		else if (datah[10]==2 && datah[11] ==9 && datah[12]==10  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type48
		else if (datah[10]==2 && datah[11] ==9 && datah[12]==11  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type49
		else if (datah[10]==2 && datah[11] ==9 && datah[12]==12  && (datah[1] != datah[2]) && (datah[5] >= 33) && (datah[6] >= 33) && ((datah[5] - datah[7])>=30) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type50
		else if (datah[10]==2 && datah[11] ==10 && datah[12]==10  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
	    //type51
		else if (datah[10]==2 && datah[11] ==10 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type52
		else if (datah[10]==2 && datah[11] ==10 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type53
		else if (datah[10]==2 && datah[11] ==11 && datah[12]==11  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type54
		else if (datah[10]==2 && datah[11] ==11 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type55
		else if (datah[10]==2 && datah[11] ==12 && datah[12]==12  && (datah[1] != datah[2]) && (datah[6] >= 17) && (std::abs(datah[5]-datah[6])<=52)) 
			coef[m-1] = 1;
		//type56
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==3 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && ((datah[5] - datah[8])>=12))
			coef[m-1] = 1;
		//type57
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==4 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && ((datah[5] - datah[8])>=10))
			coef[m-1] = 1;
		//type58
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==5 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type59
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==6 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type60
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==7 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type61
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type62
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type63
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12))
			coef[m-1] = 1;
		//type64
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12))
			coef[m-1] = 1;
		//type65
		else if(datah[10]==3 && datah[11] ==3 && datah[12]==12 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=12))
			coef[m-1] = 1;
		//type66
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==4 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=10))
			coef[m-1] = 1;
		//type67
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==5 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type68
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==6 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type69
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==7 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type70
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type71
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type72
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type73
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type74
		else if(datah[10]==3 && datah[11] ==4 && datah[12]==12 && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type75
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==5 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type76
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==6 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type77
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==7 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type78
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type79
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type80
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type81
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type82
		else if(datah[10]==3 && datah[11] ==5 && datah[12]==12 && (datah[5] >= 17) && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type83
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==6 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type84
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==7 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type85
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type86
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type87
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type88
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type89
		else if(datah[10]==3 && datah[11] ==6 && datah[12]==12 && ((datah[5] - datah[6])>=12)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type90
		else if(datah[10]==3 && datah[11] ==7 && datah[12]==7 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
	}
	delete [] datah;
}