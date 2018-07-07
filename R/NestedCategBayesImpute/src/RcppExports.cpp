// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// checkconstraints
List checkconstraints(NumericMatrix data, int neededpossiblehh, int hh_size);
RcppExport SEXP _NestedCategBayesImpute_checkconstraints(SEXP dataSEXP, SEXP neededpossiblehhSEXP, SEXP hh_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type neededpossiblehh(neededpossiblehhSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(checkconstraints(data, neededpossiblehh, hh_size));
    return rcpp_result_gen;
END_RCPP
}
// checkconstraints_HHhead_at_group_level
List checkconstraints_HHhead_at_group_level(NumericMatrix data, int neededpossiblehh, int hh_size);
RcppExport SEXP _NestedCategBayesImpute_checkconstraints_HHhead_at_group_level(SEXP dataSEXP, SEXP neededpossiblehhSEXP, SEXP hh_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type neededpossiblehh(neededpossiblehhSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(checkconstraints_HHhead_at_group_level(data, neededpossiblehh, hh_size));
    return rcpp_result_gen;
END_RCPP
}
// households2individuals
NumericMatrix households2individuals(NumericMatrix data, int hh_size);
RcppExport SEXP _NestedCategBayesImpute_households2individuals(SEXP dataSEXP, SEXP hh_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(households2individuals(data, hh_size));
    return rcpp_result_gen;
END_RCPP
}
// checkSZ
NumericVector checkSZ(NumericMatrix Data_to_check, int h);
RcppExport SEXP _NestedCategBayesImpute_checkSZ(SEXP Data_to_checkSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Data_to_check(Data_to_checkSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(checkSZ(Data_to_check, h));
    return rcpp_result_gen;
END_RCPP
}
// groupcount
NumericMatrix groupcount(NumericVector g1, NumericVector g2, int n1, int n2);
RcppExport SEXP _NestedCategBayesImpute_groupcount(SEXP g1SEXP, SEXP g2SEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(groupcount(g1, g2, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// groupcount1D
NumericVector groupcount1D(NumericVector g, int n);
RcppExport SEXP _NestedCategBayesImpute_groupcount1D(SEXP gSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(groupcount1D(g, n));
    return rcpp_result_gen;
END_RCPP
}
// sampleG
List sampleG(NumericMatrix phi, NumericMatrix data, NumericMatrix omega, NumericVector pi, NumericVector ni, NumericMatrix HHdata, List lambda);
RcppExport SEXP _NestedCategBayesImpute_sampleG(SEXP phiSEXP, SEXP dataSEXP, SEXP omegaSEXP, SEXP piSEXP, SEXP niSEXP, SEXP HHdataSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ni(niSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type HHdata(HHdataSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleG(phi, data, omega, pi, ni, HHdata, lambda));
    return rcpp_result_gen;
END_RCPP
}
// samplehouseholds
NumericMatrix samplehouseholds(NumericMatrix phi, NumericMatrix omega, NumericVector pi, NumericVector d, List lambda, int currrentbatch, int nHouseholds, int householdsize);
RcppExport SEXP _NestedCategBayesImpute_samplehouseholds(SEXP phiSEXP, SEXP omegaSEXP, SEXP piSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP currrentbatchSEXP, SEXP nHouseholdsSEXP, SEXP householdsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type currrentbatch(currrentbatchSEXP);
    Rcpp::traits::input_parameter< int >::type nHouseholds(nHouseholdsSEXP);
    Rcpp::traits::input_parameter< int >::type householdsize(householdsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(samplehouseholds(phi, omega, pi, d, lambda, currrentbatch, nHouseholds, householdsize));
    return rcpp_result_gen;
END_RCPP
}
// samplehouseholds_HHhead_at_group_level
NumericMatrix samplehouseholds_HHhead_at_group_level(NumericMatrix phi, NumericMatrix omega, NumericVector pi, NumericVector d, List lambda, int currrentbatch, int nHouseholds, int householdsize);
RcppExport SEXP _NestedCategBayesImpute_samplehouseholds_HHhead_at_group_level(SEXP phiSEXP, SEXP omegaSEXP, SEXP piSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP currrentbatchSEXP, SEXP nHouseholdsSEXP, SEXP householdsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type currrentbatch(currrentbatchSEXP);
    Rcpp::traits::input_parameter< int >::type nHouseholds(nHouseholdsSEXP);
    Rcpp::traits::input_parameter< int >::type householdsize(householdsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(samplehouseholds_HHhead_at_group_level(phi, omega, pi, d, lambda, currrentbatch, nHouseholds, householdsize));
    return rcpp_result_gen;
END_RCPP
}
// sampleM
NumericVector sampleM(NumericMatrix phi, NumericMatrix data, NumericMatrix omega, NumericVector G, NumericVector serial);
RcppExport SEXP _NestedCategBayesImpute_sampleM(SEXP phiSEXP, SEXP dataSEXP, SEXP omegaSEXP, SEXP GSEXP, SEXP serialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type serial(serialSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleM(phi, data, omega, G, serial));
    return rcpp_result_gen;
END_RCPP
}
// SampleMatrixByColumnC
NumericVector SampleMatrixByColumnC(NumericMatrix data, NumericVector r, NumericVector dup);
RcppExport SEXP _NestedCategBayesImpute_SampleMatrixByColumnC(SEXP dataSEXP, SEXP rSEXP, SEXP dupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dup(dupSEXP);
    rcpp_result_gen = Rcpp::wrap(SampleMatrixByColumnC(data, r, dup));
    return rcpp_result_gen;
END_RCPP
}
// SampleMatrixByRowC
NumericVector SampleMatrixByRowC(NumericMatrix data, NumericVector r);
RcppExport SEXP _NestedCategBayesImpute_SampleMatrixByRowC(SEXP dataSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(SampleMatrixByRowC(data, r));
    return rcpp_result_gen;
END_RCPP
}
// sampleW_multi
NumericVector sampleW_multi(NumericVector p, NumericVector d);
RcppExport SEXP _NestedCategBayesImpute_sampleW_multi(SEXP pSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleW_multi(p, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NestedCategBayesImpute_checkconstraints", (DL_FUNC) &_NestedCategBayesImpute_checkconstraints, 3},
    {"_NestedCategBayesImpute_checkconstraints_HHhead_at_group_level", (DL_FUNC) &_NestedCategBayesImpute_checkconstraints_HHhead_at_group_level, 3},
    {"_NestedCategBayesImpute_households2individuals", (DL_FUNC) &_NestedCategBayesImpute_households2individuals, 2},
    {"_NestedCategBayesImpute_checkSZ", (DL_FUNC) &_NestedCategBayesImpute_checkSZ, 2},
    {"_NestedCategBayesImpute_groupcount", (DL_FUNC) &_NestedCategBayesImpute_groupcount, 4},
    {"_NestedCategBayesImpute_groupcount1D", (DL_FUNC) &_NestedCategBayesImpute_groupcount1D, 2},
    {"_NestedCategBayesImpute_sampleG", (DL_FUNC) &_NestedCategBayesImpute_sampleG, 7},
    {"_NestedCategBayesImpute_samplehouseholds", (DL_FUNC) &_NestedCategBayesImpute_samplehouseholds, 8},
    {"_NestedCategBayesImpute_samplehouseholds_HHhead_at_group_level", (DL_FUNC) &_NestedCategBayesImpute_samplehouseholds_HHhead_at_group_level, 8},
    {"_NestedCategBayesImpute_sampleM", (DL_FUNC) &_NestedCategBayesImpute_sampleM, 5},
    {"_NestedCategBayesImpute_SampleMatrixByColumnC", (DL_FUNC) &_NestedCategBayesImpute_SampleMatrixByColumnC, 3},
    {"_NestedCategBayesImpute_SampleMatrixByRowC", (DL_FUNC) &_NestedCategBayesImpute_SampleMatrixByRowC, 2},
    {"_NestedCategBayesImpute_sampleW_multi", (DL_FUNC) &_NestedCategBayesImpute_sampleW_multi, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_NestedCategBayesImpute(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
