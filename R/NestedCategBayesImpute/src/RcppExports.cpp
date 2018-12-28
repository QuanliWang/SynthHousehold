// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// checkSZ
IntegerVector checkSZ(IntegerMatrix Data_to_check, int h);
RcppExport SEXP _NestedCategBayesImpute_checkSZ(SEXP Data_to_checkSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Data_to_check(Data_to_checkSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(checkSZ(Data_to_check, h));
    return rcpp_result_gen;
END_RCPP
}
// checkSZ2
IntegerVector checkSZ2(IntegerMatrix Data_to_check, int h);
RcppExport SEXP _NestedCategBayesImpute_checkSZ2(SEXP Data_to_checkSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Data_to_check(Data_to_checkSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(checkSZ2(Data_to_check, h));
    return rcpp_result_gen;
END_RCPP
}
// checkconstraints_HHhead_at_group_level
List checkconstraints_HHhead_at_group_level(IntegerMatrix data, int neededpossiblehh, int hh_size, int parallel);
RcppExport SEXP _NestedCategBayesImpute_checkconstraints_HHhead_at_group_level(SEXP dataSEXP, SEXP neededpossiblehhSEXP, SEXP hh_sizeSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type neededpossiblehh(neededpossiblehhSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type parallel(parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(checkconstraints_HHhead_at_group_level(data, neededpossiblehh, hh_size, parallel));
    return rcpp_result_gen;
END_RCPP
}
// checkconstraints
List checkconstraints(IntegerMatrix data, int neededpossiblehh, int hh_size);
RcppExport SEXP _NestedCategBayesImpute_checkconstraints(SEXP dataSEXP, SEXP neededpossiblehhSEXP, SEXP hh_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type neededpossiblehh(neededpossiblehhSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(checkconstraints(data, neededpossiblehh, hh_size));
    return rcpp_result_gen;
END_RCPP
}
// GetImpossibleHouseholds
List GetImpossibleHouseholds(IntegerVector d, IntegerVector n_star_h, List lambda, NumericMatrix omega, NumericMatrix phi, NumericVector pi, int blocksize, int n, int synindex, bool HHhead_at_group_level, bool Parallel);
RcppExport SEXP _NestedCategBayesImpute_GetImpossibleHouseholds(SEXP dSEXP, SEXP n_star_hSEXP, SEXP lambdaSEXP, SEXP omegaSEXP, SEXP phiSEXP, SEXP piSEXP, SEXP blocksizeSEXP, SEXP nSEXP, SEXP synindexSEXP, SEXP HHhead_at_group_levelSEXP, SEXP ParallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_star_h(n_star_hSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< int >::type blocksize(blocksizeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type synindex(synindexSEXP);
    Rcpp::traits::input_parameter< bool >::type HHhead_at_group_level(HHhead_at_group_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type Parallel(ParallelSEXP);
    rcpp_result_gen = Rcpp::wrap(GetImpossibleHouseholds(d, n_star_h, lambda, omega, phi, pi, blocksize, n, synindex, HHhead_at_group_level, Parallel));
    return rcpp_result_gen;
END_RCPP
}
// groupcount
IntegerMatrix groupcount(IntegerVector g1, IntegerVector g2, int n1, int n2);
RcppExport SEXP _NestedCategBayesImpute_groupcount(SEXP g1SEXP, SEXP g2SEXP, SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(groupcount(g1, g2, n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// groupcount1D
IntegerVector groupcount1D(IntegerVector g, int n);
RcppExport SEXP _NestedCategBayesImpute_groupcount1D(SEXP gSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(groupcount1D(g, n));
    return rcpp_result_gen;
END_RCPP
}
// UpdateBeta
double UpdateBeta(double ba, double bb, NumericMatrix v);
RcppExport SEXP _NestedCategBayesImpute_UpdateBeta(SEXP baSEXP, SEXP bbSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type ba(baSEXP);
    Rcpp::traits::input_parameter< double >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateBeta(ba, bb, v));
    return rcpp_result_gen;
END_RCPP
}
// UpdateAlpha
double UpdateAlpha(double aa, double ab, NumericVector u);
RcppExport SEXP _NestedCategBayesImpute_UpdateAlpha(SEXP aaSEXP, SEXP abSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aa(aaSEXP);
    Rcpp::traits::input_parameter< double >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateAlpha(aa, ab, u));
    return rcpp_result_gen;
END_RCPP
}
// sampleG
List sampleG(NumericMatrix phi, IntegerMatrix data, NumericMatrix omega, NumericVector pi, IntegerVector ni, IntegerMatrix HHdata, List lambda, int Parallel);
RcppExport SEXP _NestedCategBayesImpute_sampleG(SEXP phiSEXP, SEXP dataSEXP, SEXP omegaSEXP, SEXP piSEXP, SEXP niSEXP, SEXP HHdataSEXP, SEXP lambdaSEXP, SEXP ParallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ni(niSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type HHdata(HHdataSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type Parallel(ParallelSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleG(phi, data, omega, pi, ni, HHdata, lambda, Parallel));
    return rcpp_result_gen;
END_RCPP
}
// households2individuals
IntegerMatrix households2individuals(IntegerMatrix data, int hh_size);
RcppExport SEXP _NestedCategBayesImpute_households2individuals(SEXP dataSEXP, SEXP hh_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type hh_size(hh_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(households2individuals(data, hh_size));
    return rcpp_result_gen;
END_RCPP
}
// samplehouseholds
IntegerMatrix samplehouseholds(NumericMatrix phi, NumericMatrix omega, NumericVector pi, IntegerVector d, List lambda, int currrentbatch, int nHouseholds, int householdsize, int HeadAtGroupLevel, int Parallel);
RcppExport SEXP _NestedCategBayesImpute_samplehouseholds(SEXP phiSEXP, SEXP omegaSEXP, SEXP piSEXP, SEXP dSEXP, SEXP lambdaSEXP, SEXP currrentbatchSEXP, SEXP nHouseholdsSEXP, SEXP householdsizeSEXP, SEXP HeadAtGroupLevelSEXP, SEXP ParallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< List >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type currrentbatch(currrentbatchSEXP);
    Rcpp::traits::input_parameter< int >::type nHouseholds(nHouseholdsSEXP);
    Rcpp::traits::input_parameter< int >::type householdsize(householdsizeSEXP);
    Rcpp::traits::input_parameter< int >::type HeadAtGroupLevel(HeadAtGroupLevelSEXP);
    Rcpp::traits::input_parameter< int >::type Parallel(ParallelSEXP);
    rcpp_result_gen = Rcpp::wrap(samplehouseholds(phi, omega, pi, d, lambda, currrentbatch, nHouseholds, householdsize, HeadAtGroupLevel, Parallel));
    return rcpp_result_gen;
END_RCPP
}
// UpdateLambda
List UpdateLambda(IntegerMatrix HHdata_all, IntegerVector G_all, IntegerVector dHH, int FF);
RcppExport SEXP _NestedCategBayesImpute_UpdateLambda(SEXP HHdata_allSEXP, SEXP G_allSEXP, SEXP dHHSEXP, SEXP FFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type HHdata_all(HHdata_allSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type G_all(G_allSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dHH(dHHSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateLambda(HHdata_all, G_all, dHH, FF));
    return rcpp_result_gen;
END_RCPP
}
// UpdateLambdaWeighted
List UpdateLambdaWeighted(List HHdata_all, List G_all, IntegerVector dHH, int FF, NumericVector struc_weight);
RcppExport SEXP _NestedCategBayesImpute_UpdateLambdaWeighted(SEXP HHdata_allSEXP, SEXP G_allSEXP, SEXP dHHSEXP, SEXP FFSEXP, SEXP struc_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type HHdata_all(HHdata_allSEXP);
    Rcpp::traits::input_parameter< List >::type G_all(G_allSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dHH(dHHSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type struc_weight(struc_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateLambdaWeighted(HHdata_all, G_all, dHH, FF, struc_weight));
    return rcpp_result_gen;
END_RCPP
}
// sampleM
IntegerVector sampleM(NumericMatrix phi, IntegerMatrix data, NumericMatrix omega, IntegerVector G, IntegerVector serial, int Parallel);
RcppExport SEXP _NestedCategBayesImpute_sampleM(SEXP phiSEXP, SEXP dataSEXP, SEXP omegaSEXP, SEXP GSEXP, SEXP serialSEXP, SEXP ParallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type serial(serialSEXP);
    Rcpp::traits::input_parameter< int >::type Parallel(ParallelSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleM(phi, data, omega, G, serial, Parallel));
    return rcpp_result_gen;
END_RCPP
}
// SampleMissing
List SampleMissing(List MissData, List para, List orig, List G_household, IntegerVector M, List hyper);
RcppExport SEXP _NestedCategBayesImpute_SampleMissing(SEXP MissDataSEXP, SEXP paraSEXP, SEXP origSEXP, SEXP G_householdSEXP, SEXP MSEXP, SEXP hyperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type MissData(MissDataSEXP);
    Rcpp::traits::input_parameter< List >::type para(paraSEXP);
    Rcpp::traits::input_parameter< List >::type orig(origSEXP);
    Rcpp::traits::input_parameter< List >::type G_household(G_householdSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type M(MSEXP);
    Rcpp::traits::input_parameter< List >::type hyper(hyperSEXP);
    rcpp_result_gen = Rcpp::wrap(SampleMissing(MissData, para, orig, G_household, M, hyper));
    return rcpp_result_gen;
END_RCPP
}
// UpdateOmega
List UpdateOmega(double beta, IntegerMatrix M_all, int FF, int SS);
RcppExport SEXP _NestedCategBayesImpute_UpdateOmega(SEXP betaSEXP, SEXP M_allSEXP, SEXP FFSEXP, SEXP SSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type M_all(M_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< int >::type SS(SSSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateOmega(beta, M_all, FF, SS));
    return rcpp_result_gen;
END_RCPP
}
// UpdateOmegaWeighted
List UpdateOmegaWeighted(double beta, List M_all, int FF, int SS, NumericVector struc_weight);
RcppExport SEXP _NestedCategBayesImpute_UpdateOmegaWeighted(SEXP betaSEXP, SEXP M_allSEXP, SEXP FFSEXP, SEXP SSSEXP, SEXP struc_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type M_all(M_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< int >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type struc_weight(struc_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateOmegaWeighted(beta, M_all, FF, SS, struc_weight));
    return rcpp_result_gen;
END_RCPP
}
// gammarand
NumericVector gammarand(int n, double shape, double rate);
RcppExport SEXP _NestedCategBayesImpute_gammarand(SEXP nSEXP, SEXP shapeSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(gammarand(n, shape, rate));
    return rcpp_result_gen;
END_RCPP
}
// UpdatePhi
NumericMatrix UpdatePhi(IntegerMatrix data, IntegerMatrix M_all, int FF, int SS, IntegerVector d, int maxd);
RcppExport SEXP _NestedCategBayesImpute_UpdatePhi(SEXP dataSEXP, SEXP M_allSEXP, SEXP FFSEXP, SEXP SSSEXP, SEXP dSEXP, SEXP maxdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type M_all(M_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< int >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type maxd(maxdSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdatePhi(data, M_all, FF, SS, d, maxd));
    return rcpp_result_gen;
END_RCPP
}
// UpdatePhiWeighted
NumericMatrix UpdatePhiWeighted(List data, List M_all, int FF, int SS, IntegerVector d, int maxd, NumericVector struc_weight);
RcppExport SEXP _NestedCategBayesImpute_UpdatePhiWeighted(SEXP dataSEXP, SEXP M_allSEXP, SEXP FFSEXP, SEXP SSSEXP, SEXP dSEXP, SEXP maxdSEXP, SEXP struc_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type M_all(M_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< int >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type maxd(maxdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type struc_weight(struc_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdatePhiWeighted(data, M_all, FF, SS, d, maxd, struc_weight));
    return rcpp_result_gen;
END_RCPP
}
// UpdatePi
List UpdatePi(double alpha, IntegerVector G_all, int FF);
RcppExport SEXP _NestedCategBayesImpute_UpdatePi(SEXP alphaSEXP, SEXP G_allSEXP, SEXP FFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type G_all(G_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdatePi(alpha, G_all, FF));
    return rcpp_result_gen;
END_RCPP
}
// UpdatePiWeighted
List UpdatePiWeighted(double alpha, List G_all, int FF, NumericVector struc_weight);
RcppExport SEXP _NestedCategBayesImpute_UpdatePiWeighted(SEXP alphaSEXP, SEXP G_allSEXP, SEXP FFSEXP, SEXP struc_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< List >::type G_all(G_allSEXP);
    Rcpp::traits::input_parameter< int >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type struc_weight(struc_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdatePiWeighted(alpha, G_all, FF, struc_weight));
    return rcpp_result_gen;
END_RCPP
}
// sampleW_multi
IntegerVector sampleW_multi(NumericVector p, NumericVector d);
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
    {"_NestedCategBayesImpute_checkSZ", (DL_FUNC) &_NestedCategBayesImpute_checkSZ, 2},
    {"_NestedCategBayesImpute_checkSZ2", (DL_FUNC) &_NestedCategBayesImpute_checkSZ2, 2},
    {"_NestedCategBayesImpute_checkconstraints_HHhead_at_group_level", (DL_FUNC) &_NestedCategBayesImpute_checkconstraints_HHhead_at_group_level, 4},
    {"_NestedCategBayesImpute_checkconstraints", (DL_FUNC) &_NestedCategBayesImpute_checkconstraints, 3},
    {"_NestedCategBayesImpute_GetImpossibleHouseholds", (DL_FUNC) &_NestedCategBayesImpute_GetImpossibleHouseholds, 11},
    {"_NestedCategBayesImpute_groupcount", (DL_FUNC) &_NestedCategBayesImpute_groupcount, 4},
    {"_NestedCategBayesImpute_groupcount1D", (DL_FUNC) &_NestedCategBayesImpute_groupcount1D, 2},
    {"_NestedCategBayesImpute_UpdateBeta", (DL_FUNC) &_NestedCategBayesImpute_UpdateBeta, 3},
    {"_NestedCategBayesImpute_UpdateAlpha", (DL_FUNC) &_NestedCategBayesImpute_UpdateAlpha, 3},
    {"_NestedCategBayesImpute_sampleG", (DL_FUNC) &_NestedCategBayesImpute_sampleG, 8},
    {"_NestedCategBayesImpute_households2individuals", (DL_FUNC) &_NestedCategBayesImpute_households2individuals, 2},
    {"_NestedCategBayesImpute_samplehouseholds", (DL_FUNC) &_NestedCategBayesImpute_samplehouseholds, 10},
    {"_NestedCategBayesImpute_UpdateLambda", (DL_FUNC) &_NestedCategBayesImpute_UpdateLambda, 4},
    {"_NestedCategBayesImpute_UpdateLambdaWeighted", (DL_FUNC) &_NestedCategBayesImpute_UpdateLambdaWeighted, 5},
    {"_NestedCategBayesImpute_sampleM", (DL_FUNC) &_NestedCategBayesImpute_sampleM, 6},
    {"_NestedCategBayesImpute_SampleMissing", (DL_FUNC) &_NestedCategBayesImpute_SampleMissing, 6},
    {"_NestedCategBayesImpute_UpdateOmega", (DL_FUNC) &_NestedCategBayesImpute_UpdateOmega, 4},
    {"_NestedCategBayesImpute_UpdateOmegaWeighted", (DL_FUNC) &_NestedCategBayesImpute_UpdateOmegaWeighted, 5},
    {"_NestedCategBayesImpute_gammarand", (DL_FUNC) &_NestedCategBayesImpute_gammarand, 3},
    {"_NestedCategBayesImpute_UpdatePhi", (DL_FUNC) &_NestedCategBayesImpute_UpdatePhi, 6},
    {"_NestedCategBayesImpute_UpdatePhiWeighted", (DL_FUNC) &_NestedCategBayesImpute_UpdatePhiWeighted, 7},
    {"_NestedCategBayesImpute_UpdatePi", (DL_FUNC) &_NestedCategBayesImpute_UpdatePi, 3},
    {"_NestedCategBayesImpute_UpdatePiWeighted", (DL_FUNC) &_NestedCategBayesImpute_UpdatePiWeighted, 4},
    {"_NestedCategBayesImpute_sampleW_multi", (DL_FUNC) &_NestedCategBayesImpute_sampleW_multi, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_NestedCategBayesImpute(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
