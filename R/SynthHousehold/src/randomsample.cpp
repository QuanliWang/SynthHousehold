#include <math.h>
int samplew(double *p, int n, double d) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }

    for(k=0;k < n && d>myw[k];k++)
        ;
    delete [] myw;
    return k+1;
}

