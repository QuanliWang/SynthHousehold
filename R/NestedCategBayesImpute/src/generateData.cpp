#include <Rcpp.h>
using namespace Rcpp;
#include "samplehouseholds.h"
#include "checkconstraints.h"

IntegerMatrix Concatenate(List data) {
  int totalcolumns = 0;
  int rows;
  for (int i = 0; i < data.length(); i++) {
    IntegerMatrix current = data[i];
    totalcolumns += current.ncol();
    rows = current.nrow();
  }
  IntegerMatrix result(rows, totalcolumns);
  int offset = 0;
  for (int i = 0; i < data.length(); i++) {
    IntegerMatrix current = data[i];
    std::copy(current.begin(),current.end(),result.begin() + offset);
    offset += current.length();
  }
  return result;
}

// [[Rcpp::export]]
List GenerateData(int hh_size,List lambda, NumericMatrix omega, NumericMatrix phi,
                  NumericVector pi, IntegerVector d, int batches_done,
                         int valid_hh_needed, int blocksize, int synindex, bool HHhead_at_group_level) {
  List Individuals_extra;
  List G_extra;
  List HHData_extra;
  List synIndividuals;
  int batch_index = 0;
  int valid_hh_found = 0;
  int p = d.length();

  while (valid_hh_found < valid_hh_needed) {
    batch_index ++;
    //generate a batch of 10K household
    List checked_households;
    if (HHhead_at_group_level) {
      IntegerMatrix data_to_check = samplehouseholds_HHhead_at_group_level(phi,omega, pi, d, lambda,
                                                              batch_index+batches_done, blocksize,hh_size);
      checked_households = checkconstraints_HHhead_at_group_level(data_to_check,valid_hh_needed-valid_hh_found, hh_size);
    } else {
      IntegerMatrix data_to_check = samplehouseholds(phi,omega, pi, d,
                                                                lambda,batch_index+batches_done, blocksize,hh_size);
      checked_households = checkconstraints(data_to_check,valid_hh_needed-valid_hh_found, hh_size);
    }

    IntegerMatrix Households = checked_households["Households"];
    int possible = checked_households["possible"];
    IntegerMatrix synHouseholds = checked_households["synHouseholds"];
    if (Households != R_NilValue) {
      int DIM = p + lambda.length() + 1;
      IntegerVector hhrow(Households.ncol());
      G_extra.push_back(Households(Range(hh_size * DIM,hh_size * DIM), Range(0,Households.ncol()-1)));
      HHData_extra.push_back(Households(Range(p+2,DIM-1), Range(0,Households.ncol()-1)));
      Individuals_extra.push_back(households2individuals(Households, hh_size));
    }

    valid_hh_found += possible;
    if (synindex > 0) {
      if (synHouseholds != R_NilValue) {
        synIndividuals.push_back(households2individuals(synHouseholds,hh_size));
      }
    }
  }

  batch_index += batches_done;
  return(List::create(Named("Individuals_extra", Concatenate(Individuals_extra)),
                      Named("G_extra", Concatenate(G_extra)),
                      Named("HHData_extra", Concatenate(HHData_extra)),
                      Named("synIndividuals", Concatenate(synIndividuals)),
                      Named("batch.index", batch_index)
                        ));
}

/*
 GenerateData <- function(hh_size,lambda, omega, phi,pi, d, batches_done,
 valid_hh_needed,blocksize,synindex,HHhead_at_group_level) {
Individuals_extra <- list()
G_extra <- list()
HHData_extra <- list()
batch.index <- 0
valid_hh_found <- 0
p <- length(d)
synIndividuals <- list()

while (valid_hh_found< valid_hh_needed) {
batch.index <- batch.index + 1
#generate a batch of 10K household
if (HHhead_at_group_level) {
data_to_check <- samplehouseholds_HHhead_at_group_level(phi,omega, pi, d, lambda,
                                                        batch.index+batches_done, blocksize,hh_size)
checked_households <- checkconstraints_HHhead_at_group_level(data_to_check,valid_hh_needed-valid_hh_found, hh_size)
} else {
data_to_check <- samplehouseholds(phi,omega, pi, d, lambda,batch.index+batches_done, blocksize,hh_size)
checked_households <- checkconstraints(data_to_check,valid_hh_needed-valid_hh_found, hh_size)
}
if (length(checked_households$Households) > 0) {
DIM <- p + length(lambda) + 1
G_extra[[batch.index]] <- checked_households$Households[hh_size * DIM +1,]
HHData_extra[[batch.index]] <- checked_households$Households[(p+3): DIM,]
Individuals_extra[[batch.index]] <- households2individuals(checked_households$Households, hh_size)
}

valid_hh_found <- valid_hh_found + checked_households$possible
if (synindex > 0) {
if (length(checked_households$synHouseholds) > 0) {
synIndividuals[[batch.index]]  <- households2individuals(checked_households$synHouseholds,hh_size)
}
}
}

Individuals_extra <- do.call(cbind, Individuals_extra)
G_extra <- unlist(G_extra)
if (HHhead_at_group_level) {
HHData_extra <- do.call(cbind, HHData_extra)
} else {
HHData_extra <- unlist(HHData_extra)
}

if (synindex > 0) {
synIndividuals <- do.call(cbind, synIndividuals)
}
batch.index <- batch.index + batches_done
return(list(Individuals_extra = Individuals_extra,
            G_extra = G_extra,
HHData_extra = HHData_extra,
synIndividuals = synIndividuals,
batch.index = batch.index))
}

*/
