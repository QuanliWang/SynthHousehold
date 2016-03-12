#include <cstring>
#include <Rcpp.h>
using namespace Rcpp;
#include "checkconstraints.h"

// [[Rcpp::export]]
List checkconstraints(NumericMatrix data,int neededpossiblehh) {
  int nHouseholds = data.nrow(); //data item in rows. !!

  //use the raw data instead, which has hh_size * 8 + 1 + hh_size
  int columns = data.ncol();
  int hh_size = (columns -1) / (DIM+1);

  NumericVector isPossible(nHouseholds);
  int totalpossible = checkconstraints_imp(data.begin(), isPossible.begin(), hh_size, nHouseholds);

  int rows = nHouseholds-totalpossible;
  NumericMatrix newdata(columns,rows);
  NumericMatrix syndata(columns,totalpossible);
  NumericVector impossible_counts(totalpossible);

  int count1 = 0;
  int count2 = 0;
  for (int i = 0; i < nHouseholds && count2 < neededpossiblehh; i++) {
    if (isPossible[i] == 0) { //found an impossible household
      for (int j = 0; j < columns; j++) { //impossible ones
        newdata[count1*columns+j] = data[j*nHouseholds+i];
      }
      count1++;
    } else {
      impossible_counts[count2] = count1;
      for (int j = 0; j < columns; j++) { //possible ones (syndata)
        syndata[count2*columns+j] = data[j*nHouseholds+i];
      }
      count2++;
    }
  }

  if (count1 < rows) { //rows in C are the columns in the return matrix
    //need to resize the output matrix
    newdata = newdata(Range(0,columns-1), Range(0,count1-1));
  }

  if (count2 < totalpossible) { //truncate possible households if too many
    //need to resize the output matrix
    syndata = syndata(Range(0,columns-1), Range(0,count2-1));
  }

  return List::create(Named("outcome", isPossible),
                      Named("Households", newdata),
                      Named("Index", impossible_counts),
                      Named("synHouseholds", syndata),
                      Named("possible", count2));
}

// [[Rcpp::export]]
NumericMatrix households2individuals(NumericMatrix data){

  int nHouseholds = data.ncol();

  //use the raw data instead, which has hh_size * 8 + 1 + hh_size columns (in C)
  int columns = data.nrow();
  int hh_size = (columns - 1) / (DIM+1);

  NumericMatrix Individuals(DIM + 2, nHouseholds*hh_size);

  int c9 = hh_size * DIM;
  int count = 0;
  for (int i = 0; i < nHouseholds; i++) {
    int base = i * columns;
    for (int j = 0; j < hh_size; j++) {
      for (int k = 0; k < DIM;k++) {
        Individuals[count++] = data[base + j*DIM+k];
      }
      Individuals[count++] = data[base + c9];
      Individuals[count++] = data[base + c9 + 1 + j];
    }
  }
  return(Individuals);
}
