#include <Rcpp.h>
using namespace Rcpp;
#include "utils.h"

NumericMatrix transposeAndNormalize(NumericMatrix mat) {
  int row = mat.nrow();
  int col = mat.ncol();
  NumericMatrix omegaT(col,row);

  double *currentrow = omegaT.begin();
  for (int k =0; k < row; k++) {
    double dsum = 0.0;
    for (int l = 0; l <col; l++) {
      //transpose first
      currentrow[l] = mat[l*row+k];  //omegat[k * col + l] = mat[l*row+k];
      dsum += currentrow[l];
    }
    currentrow[0] /= dsum;
    for (int l = 1; l <col; l++) {
      currentrow[l] = currentrow[l]/dsum + currentrow[l-1]; //normilized cum_sum
    }
    currentrow += col;
  }
  return omegaT;
}


