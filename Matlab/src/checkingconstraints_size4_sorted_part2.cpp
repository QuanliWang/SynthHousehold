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
		//type91
		if(datah[10]==3 && datah[11] ==7 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type92
		else if(datah[10]==3 && datah[11] ==7 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type93
		else if(datah[10]==3 && datah[11] ==7 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type94
		else if(datah[10]==3 && datah[11] ==7 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type95
		else if(datah[10]==3 && datah[11] ==7 && datah[12]==12 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type96
		else if(datah[10]==3 && datah[11] ==8 && datah[12]==8 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type97
		else if(datah[10]==3 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type98
		else if(datah[10]==3 && datah[11] ==8 && datah[12]==10 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type99
		else if(datah[10]==3 && datah[11] ==8 && datah[12]==11 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type100
		else if(datah[10]==3 && datah[11] ==8 && datah[12]==12 && ((datah[5] - datah[6])>=12)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type101
		else if(datah[10]==3 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=30) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type102
		else if(datah[10]==3 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type103
		else if(datah[10]==3 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type104
		else if(datah[10]==3 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=12)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type105
		else if(datah[10]==3 && datah[11] ==10 && datah[12]==10 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type106
		else if(datah[10]==3 && datah[11] ==10 && datah[12]==11 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type107
		else if(datah[10]==3 && datah[11] ==10 && datah[12]==12 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type108
		else if(datah[10]==3 && datah[11] ==11 && datah[12]==11 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type109
		else if(datah[10]==3 && datah[11] ==11 && datah[12]==12 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type110
		else if(datah[10]==3 && datah[11] ==12 && datah[12]==12 && ((datah[5] - datah[6])>=12))
			coef[m-1] = 1;
		//type111
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==4 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=10))
			coef[m-1] = 1;
		//type112
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==5 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type113
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==6 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type114
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==7 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type115
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==8 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type116
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type117
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==10 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type118
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==11 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type119
		else if(datah[10]==4 && datah[11] ==4 && datah[12]==12 && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=10))
			coef[m-1] = 1;
		//type120
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==5 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=13))
			coef[m-1] = 1;
		//type121
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==6 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type122
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==7 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type123
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==8 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type124
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type125
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==10 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type126
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==11 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type127
		else if(datah[10]==4 && datah[11] ==5 && datah[12]==12 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type128
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==6 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type129
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==7 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type130
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==8 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type131
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type132
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==10 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type133
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==11 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type134
		else if(datah[10]==4 && datah[11] ==6 && datah[12]==12 && ((datah[5] - datah[6])>=10)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type135
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==7 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type136
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==8 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type137
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type138
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==10 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type139
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==11 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type140
		else if(datah[10]==4 && datah[11] ==7 && datah[12]==12 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type141
		else if(datah[10]==4 && datah[11] ==8 && datah[12]==8 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type142
		else if(datah[10]==4 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type143
		else if(datah[10]==4 && datah[11] ==8 && datah[12]==10 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type144
		else if(datah[10]==4 && datah[11] ==8 && datah[12]==11 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type145
		else if(datah[10]==4 && datah[11] ==8 && datah[12]==12 && ((datah[5] - datah[6])>=10)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type146
		else if(datah[10]==4 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=30) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type147
		else if(datah[10]==4 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type148
		else if(datah[10]==4 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type149
		else if(datah[10]==4 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=10)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type150
		else if(datah[10]==4 && datah[11] ==10 && datah[12]==10 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type151
		else if(datah[10]==4 && datah[11] ==10 && datah[12]==11 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type152
		else if(datah[10]==4 && datah[11] ==10 && datah[12]==12 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type153
		else if(datah[10]==4 && datah[11] ==11 && datah[12]==11 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type154
		else if(datah[10]==4 && datah[11] ==11 && datah[12]==12 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type155
		else if(datah[10]==4 && datah[11] ==12 && datah[12]==12 && ((datah[5] - datah[6])>=10))
			coef[m-1] = 1;
		//type156
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==6 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type157
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==7 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type158
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==8 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type159
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type160
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==10 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type161
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==11 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type162
		else if(datah[10]==5 && datah[11] ==5 && datah[12]==12 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=13))
			coef[m-1] = 1;
		//type163
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==6 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9) && ((datah[8] - datah[5])>=9))
			coef[m-1] = 1;
		//type164
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==7 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type165
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==8 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type166
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type167
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==10 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type168
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==11 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type169
		else if(datah[10]==5 && datah[11] ==6 && datah[12]==12 && ((datah[6] - datah[5])>=13)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type170
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==7 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type171
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==8 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type172
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type173
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==10 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type174
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==11 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type175
		else if(datah[10]==5 && datah[11] ==7 && datah[12]==12 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type176
		else if(datah[10]==5 && datah[11] ==8 && datah[12]==8 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33) && (std::abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type177
		else if(datah[10]==5 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type178
		else if(datah[10]==5 && datah[11] ==8 && datah[12]==10 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type179
		else if(datah[10]==5 && datah[11] ==8 && datah[12]==11 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type180
		else if(datah[10]==5 && datah[11] ==8 && datah[12]==12 && ((datah[6] - datah[5])>=13)  && (std::abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		}
		delete [] datah;
}