#include <Rcpp.h>
using namespace Rcpp;
#include "samplehouseholds.h"
#include "checkconstraints.h"

//input index is for debugging only, can be taken out later on
IntegerMatrix Concatenate(List data, int index) {
  //Rcout << index << std::endl;
  int totalcolumns = 0;
  int rows;
  for (int i = 0; i < data.length(); i++) {
    IntegerMatrix current = data[i];
    if (current != R_NilValue) {
      totalcolumns += current.ncol();
      rows = current.nrow();
    }
  }
  if (rows <= 0 || totalcolumns <=0) {
    IntegerMatrix Empty(0,0);
    return Empty;
  }
  IntegerMatrix result(rows, totalcolumns);
  int offset = 0;
  for (int i = 0; i < data.length(); i++) {
    IntegerMatrix current = data[i];
    if (current != R_NilValue) {
      std::copy(current.begin(),current.end(),result.begin() + offset);
    }
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
  return(List::create(Named("Individuals_extra", Concatenate(Individuals_extra,1)),
                      Named("G_extra", Concatenate(G_extra,2)),
                      Named("HHData_extra", Concatenate(HHData_extra,3)),
                      Named("synIndividuals", Concatenate(synIndividuals,4)),
                      Named("batch.index", batch_index)
                        ));
}

// [[Rcpp::export]]
List GetImpossibleHouseholds(IntegerVector d,IntegerVector n_star_h, List lambda,
                             NumericMatrix omega, NumericMatrix phi, NumericVector pi,
                             int blocksize,  int n, int synindex, bool HHhead_at_group_level) {
  int cumsize = 0;
  NumericMatrix hh_size_new(n_star_h.length(),1);
  List hh_index;
  List ImpossibleIndividuals;
  List G_extra;
  List HHdata_extra;
  List synIndividuals_all;

  int batches_done = 0;
  for (int hh_size  = 1; hh_size <=n_star_h.length(); hh_size++) {
    int hh_size_real;
    if (HHhead_at_group_level) {
      hh_size_real = hh_size;
    } else {
      hh_size_real = hh_size + 1;
    }
    List batch = GenerateData(hh_size_real,lambda, omega, phi,pi, d, batches_done,
                          n_star_h[hh_size - 1],blocksize,synindex,HHhead_at_group_level);

    IntegerMatrix G_extra_1 = batch["G_extra"];
    IntegerMatrix Individuals_extra = batch["Individuals_extra"];
    IntegerMatrix HHData_extra_1 = batch["HHData_extra"];
    if (G_extra_1 != R_NilValue) {
      hh_size_new[hh_size - 1] = G_extra_1.length();
    } else {
      hh_size_new[hh_size - 1] = 0;
    }
    IntegerMatrix hh_temp(hh_size_new[hh_size-1] * hh_size_real,1);
    int count = 0;
    for (int j = 0; j < hh_size_new[hh_size-1]; j++) {
      for (int k = 0; k < hh_size_real; k++) {
        hh_temp[count++] = cumsize + 1 + j;
      }
    }
    cumsize += hh_size_new[hh_size-1];
    ImpossibleIndividuals.push_back(Individuals_extra);
    G_extra.push_back(G_extra_1);

    if (HHData_extra_1 != R_NilValue) {
      IntegerMatrix HHData_extra_wth_size(HHData_extra_1.nrow() + 1, HHData_extra_1.ncol());
      for (int i = 0; i < HHData_extra_1.ncol(); i++) {
        for (int j = 0; j < HHData_extra_1.nrow(); j++) {
          HHData_extra_wth_size(j,i) = HHData_extra_1(j,i);
        }
        HHData_extra_wth_size(HHData_extra_1.nrow(),i) = hh_size;
      }
      HHdata_extra.push_back(HHData_extra_wth_size);
    }
    if (synindex > 0) {
      IntegerMatrix synIndividuals = batch["synIndividuals"];
      synIndividuals_all.push_back(synIndividuals);
    }
    int batch_index = batch["batch.index"];
    batches_done = batch_index;
  }

  IntegerMatrix hh_index_combined = Concatenate(hh_index,5);
  IntegerMatrix ImpossibleIndividuals_combined = Concatenate(ImpossibleIndividuals,6);
  for (int i = 0; i < hh_index.length(); i++) {
    ImpossibleIndividuals_combined(0,i) = hh_index(0,i);
    ImpossibleIndividuals_combined(0,i) += n;
  }

  int DIM = ImpossibleIndividuals_combined.nrow() - 2;
  IntegerMatrix IndividualData_extra(DIM,ImpossibleIndividuals_combined.ncol());
  IntegerMatrix G_Individuals_and_M_extra(2,ImpossibleIndividuals_combined.ncol());

  for (int i = 0; i < ImpossibleIndividuals_combined.ncol(); i++) {
    for (int j = 0; j < DIM; j++) {
      IndividualData_extra(j,i) = ImpossibleIndividuals_combined(j,i);
    }
    G_Individuals_and_M_extra(0,i) = ImpossibleIndividuals_combined(DIM,i);
    G_Individuals_and_M_extra(1,i) = ImpossibleIndividuals_combined(DIM+1,i);
  }
  return(List::create(Named("G_Individuals_and_M_extra", G_Individuals_and_M_extra),
                      Named("G_extra", Concatenate(G_extra,7)),
                      Named("IndividualData_extra", IndividualData_extra),
                      Named("HHdata_extra", Concatenate(HHdata_extra,8)),
                      Named("hh_size_new", hh_size_new),
                      Named("synIndividuals_all", Concatenate(synIndividuals_all,9))
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
/*
 GetImpossibleHouseholds <- function(d,n_star_h,lambda,omega,phi,pi,blocksize,n,synindex,HHhead_at_group_level) {
#save(d,n_star_h,lambda,omega,phi,pi,blocksize,n,synindex, file = "GetImpossibleHouseholds.RData")
 cumsize <- 0
hh_size_new <-  matrix(0,nrow = length(n_star_h), ncol = 1)
hh_index <- list()
ImpossibleIndividuals <- list()
G_extra <- list()
HHdata_extra <- list()
synIndividuals_all <- list()

##
batches.done <- 0
for (hh_size in  1:(length(n_star_h))) {
if (HHhead_at_group_level) {
hh_size_real <- hh_size
} else {
hh_size_real <- hh_size + 1
}
batch <- GenerateData(hh_size_real,lambda, omega, phi,pi, d, batches.done,
                      n_star_h[hh_size],blocksize,synindex,HHhead_at_group_level)

hh_size_new[hh_size] <- length(batch$G_extra)
hh_index[[hh_size]] <- cumsize + rep(1:hh_size_new[hh_size], each = hh_size_real)
cumsize <- cumsize + hh_size_new[hh_size]
ImpossibleIndividuals[[hh_size]] <- batch$Individuals_extra
G_extra[[hh_size]] <-  batch$G_extra
HHdata_extra[[hh_size]] <- rbind(batch$HHData_extra,rep(hh_size, times = hh_size_new[hh_size])) # length(lambda) by ...
if (synindex > 0) {
synIndividuals_all[[hh_size]] <- batch$synIndividuals
}
batches.done <- batch$batch.index
}

##
hh_index <- unlist(hh_index)
ImpossibleIndividuals <- do.call(cbind, ImpossibleIndividuals)
G_extra <- unlist(G_extra)
HHdata_extra <- do.call(cbind, HHdata_extra)
if (synindex > 0) {
synIndividuals_all <- do.call(cbind, synIndividuals_all)
}

ImpossibleIndividuals[1,] <- n + hh_index
DIM <- dim(ImpossibleIndividuals)[1] - 2
IndividualData_extra <- ImpossibleIndividuals[1:DIM,]
G_Individuals_and_M_extra <- ImpossibleIndividuals[(DIM+1):(DIM+2),]
return(list(G_Individuals_and_M_extra = G_Individuals_and_M_extra,
            G_extra = G_extra,
IndividualData_extra = IndividualData_extra,
HHdata_extra = HHdata_extra,
hh_size_new = hh_size_new,
synIndividuals_all = synIndividuals_all))
}
 */
