#include <cstring>
#include <Rcpp.h>
using namespace Rcpp;
#include "checkconstraints.h"

// [[Rcpp::export]]
List checkconstraints(IntegerMatrix data,int neededpossiblehh, int hh_size) {
  int nHouseholds = data.nrow(); //data item in rows. !!

  //use the raw data instead, which has hh_size * DIM + 1 + hh_size
  int columns = data.ncol();
  int DIM = (columns -1) / hh_size -1;
  //int hh_size = (columns -1) / (DIM+1);

  IntegerVector isPossible(nHouseholds);
  int totalpossible = checkconstraints_imp(data.begin(), isPossible.begin(), hh_size, DIM, nHouseholds);

  int rows = nHouseholds-totalpossible;
  IntegerMatrix newdata(columns,rows);
  IntegerMatrix syndata(columns,totalpossible);
  IntegerVector impossible_counts(totalpossible);

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
    if (count1 >= 1) {
      newdata = newdata(Range(0,columns-1), Range(0,count1-1));
    } else {
      newdata = R_NilValue;
    }
  }

  if (count2 < totalpossible) { //truncate possible households if too many
    //need to resize the output matrix
    if (count2 >= 1) {
      syndata = syndata(Range(0,columns-1), Range(0,count2-1));
    } else {
      syndata = R_NilValue;
    }
  }

  return List::create(Named("outcome", isPossible),
                      Named("Households", newdata),
                      Named("Index", impossible_counts),
                      Named("synHouseholds", syndata),
                      Named("possible", count2));
}


// [[Rcpp::export]]
List checkconstraints_HHhead_at_group_level(IntegerMatrix data,int neededpossiblehh, int hh_size) {
  int nHouseholds = data.nrow(); //data item in rows. !!

  //use the raw data instead, which has hh_size * DIM + 1 + hh_size
  int columns = data.ncol();
  int DIM = (columns -1) / hh_size -1;

  IntegerVector isPossible(nHouseholds);
  int totalpossible = checkconstraints_imp_HHhead_at_group_level(data.begin(), isPossible.begin(), hh_size, DIM, nHouseholds);

  int rows = nHouseholds-totalpossible;
  IntegerMatrix newdata(columns,rows);
  IntegerMatrix syndata(columns,totalpossible);
  IntegerVector impossible_counts(totalpossible);

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
    if (count1 >= 1) {
      newdata = newdata(Range(0,columns-1), Range(0,count1-1));
    } else {
      newdata = R_NilValue;
    }
  }

  if (count2 < totalpossible) { //truncate possible households if too many
    //need to resize the output matrix
    if (count2 >= 1) {
      syndata = syndata(Range(0,columns-1), Range(0,count2-1));
    } else {
      syndata = R_NilValue;
    }
  }

  return List::create(Named("outcome", isPossible),
                      Named("Households", newdata),
                      Named("Index", impossible_counts),
                      Named("synHouseholds", syndata),
                      Named("possible", count2));
}

