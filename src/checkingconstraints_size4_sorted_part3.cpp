#include "mex.h"
#include <vector>

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
	
		//type181
		if(datah[9]==1 && datah[10]==5 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[5] - datah[7])>=30) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type182
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type183
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type184
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 33) && ((datah[6] - datah[5])>=13)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type185
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type186
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type187
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type188
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type189
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type190
		else if(datah[9]==1 && datah[10]==5 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=13))
			coef[m-1] = 1;
		//type191
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==7 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type182
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==8 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type193
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type194
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==10 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type195
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type196
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==6 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && ((datah[7] - datah[5])>=9))
			coef[m-1] = 1;
		//type197
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==7 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type198
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==8 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type199
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type200
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==10 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type201
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type202
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==7 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type203
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==8 && datah[12]==8 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type204
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type205
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==8 && datah[12]==10 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type206
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==8 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type207
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==8 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type208
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && ((datah[5] - datah[7])>=30) && ((datah[5] - datah[8])>=30))
			coef[m-1] = 1;
		//type209
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type210
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type211
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 33) && ((datah[6] - datah[5])>=9)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type212
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type213
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type214
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type215
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type216
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type217
		else if(datah[9]==1 && datah[10]==6 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17) && ((datah[6] - datah[5])>=9))
			coef[m-1] = 1;
		//type218
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==7 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type219
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==8 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type220
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==9 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type221
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type222
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type223
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==7 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type224
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==8 && datah[12]==8 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type225
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type226
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==8 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type227
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==8 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type228
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==8 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type229
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type230
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type231
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type232
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type233
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type234
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type235
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type236
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type237
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type238
		else if(datah[9]==1 && datah[10]==7 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type239
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==8 && datah[12]==8 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (abs(datah[5] - datah[8])<=33))
			coef[m-1] = 1;
		//type240
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==8 && datah[12]==9 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type241
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==8 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type242
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==8 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type243
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==8 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && (abs(datah[5] - datah[7])<=33))
			coef[m-1] = 1;
		//type244
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type245
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type246
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type247
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type248
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type249
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type250
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type251
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type252
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type253
		else if(datah[9]==1 && datah[10]==8 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17) && (abs(datah[5] - datah[6])<=33))
			coef[m-1] = 1;
		//type254
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==9 && datah[12]==9 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30)  && ((datah[5] - datah[7])>=30) && (datah[5] - datah[8])>=30)
			coef[m-1] = 1;
		//type255
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==9 && datah[12]==10 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type256
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==9 && datah[12]==11 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type257
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==9 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30)  && ((datah[5] - datah[7])>=30))
			coef[m-1] = 1;
		//type258
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type259
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type260
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type261
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type262
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type263
		else if(datah[9]==1 && datah[10]==9 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 33) && ((datah[5] - datah[6])>=30))
			coef[m-1] = 1;
		//type264
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==10 && datah[12]==10 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type265
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==10 && datah[12]==11 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type266
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==10 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type267
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type268
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type269
		else if(datah[9]==1 && datah[10]==10 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type270
		else if(datah[9]==1 && datah[10]==11 && datah[11] ==11 && datah[12]==11 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type271
		else if(datah[9]==1 && datah[10]==11 && datah[11] ==11 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type272
		else if(datah[9]==1 && datah[10]==11 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		//type273
		else if(datah[9]==1 && datah[10]==12 && datah[11] ==12 && datah[12]==12 && (datah[5] >= 17))
			coef[m-1] = 1;
		}
		delete [] datah;
}