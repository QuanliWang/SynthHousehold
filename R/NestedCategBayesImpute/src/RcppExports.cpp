// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// checkconstraints
List checkconstraints(NumericMatrix data, int neededpossiblehh);
RcppExport SEXP NestedCategBayesImpute_checkconstraints(SEXP dataSEXP, SEXP neededpossiblehhSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type neededpossiblehh(neededpossiblehhSEXP);
    __result = Rcpp::wrap(checkconstraints(data, neededpossiblehh));
    return __result;
END_RCPP
}
// households2individuals
NumericMatrix households2individuals(NumericMatrix data);
RcppExport SEXP NestedCategBayesImpute_households2individuals(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    __result = Rcpp::wrap(households2individuals(data));
    return __result;
END_RCPP
}
// groupcount
NumericMatrix groupcount(NumericVector g1, NumericVector g2, int n1, int n2);
RcppExport SEXP NestedCategBayesImpute_groupcount(SEXP g1SEXP, SEXP g2SEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    __result = Rcpp::wrap(groupcount(g1, g2, n1, n2));
    return __result;
END_RCPP
}
// groupcount1D
NumericVector groupcount1D(NumericVector g, int n);
RcppExport SEXP NestedCategBayesImpute_groupcount1D(SEXP gSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    __result = Rcpp::wrap(groupcount1D(g, n));
    return __result;
END_RCPP
}
// samplehouseholds
NumericMatrix samplehouseholds(NumericMatrix phi, NumericMatrix w, NumericVector pi, NumericVector d, List lambda, int currrentbatch, int nHouseholds, int householdsize);
RcppExport SEXP NestedCategBayesImpute_samplehouseholds(SEXP phiSEXP, SEXP wSEXP, SEXP piSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP currrentbatchSEXP, SEXP nHouseholdsSEXP, SEXP householdsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type currrentbatch(currrentbatchSEXP);
    Rcpp::traits::input_parameter< int >::type nHouseholds(nHouseholdsSEXP);
    Rcpp::traits::input_parameter< int >::type householdsize(householdsizeSEXP);
    __result = Rcpp::wrap(samplehouseholds(phi, w, pi, d, lambda, currrentbatch, nHouseholds, householdsize));
    return __result;
END_RCPP
}
// samplezHH
List samplezHH(NumericMatrix phi, NumericMatrix data, NumericMatrix w, NumericVector pi, NumericVector S, NumericMatrix HHdata, NumericMatrix lambda1, NumericMatrix lambda2);
RcppExport SEXP NestedCategBayesImpute_samplezHH(SEXP phiSEXP, SEXP dataSEXP, SEXP wSEXP, SEXP piSEXP, SEXP SSEXP, SEXP HHdataSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type HHdata(HHdataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda2(lambda2SEXP);
    __result = Rcpp::wrap(samplezHH(phi, data, w, pi, S, HHdata, lambda1, lambda2));
    return __result;
END_RCPP
}
// samplezmember
NumericVector samplezmember(NumericMatrix phi, NumericMatrix data, NumericMatrix w, NumericVector zHH, NumericVector serial);
RcppExport SEXP NestedCategBayesImpute_samplezmember(SEXP phiSEXP, SEXP dataSEXP, SEXP wSEXP, SEXP zHHSEXP, SEXP serialSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zHH(zHHSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type serial(serialSEXP);
    __result = Rcpp::wrap(samplezmember(phi, data, w, zHH, serial));
    return __result;
END_RCPP
}