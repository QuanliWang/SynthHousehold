#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix groupcount(NumericVector g1, NumericVector g2, int row, int col) {

  NumericMatrix counts(row, col);
  for (int i = 0; i < g1.length(); i++) {
    counts[((int)g1[i] -1) + ((int)g2[i]-1) * row ]++;
  }
  return counts;
}

// [[Rcpp::export]]
NumericVector groupcount1D(NumericVector g, int groups) {
  NumericVector counts(groups);
  for (int i = 0; i < g.length(); i++) {
    counts[((int)g[i]-1)]++;
  }
  return counts;
}
