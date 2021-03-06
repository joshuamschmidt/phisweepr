// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Arma_colSums
arma::rowvec Arma_colSums(const arma::mat& x);
RcppExport SEXP _phisweepr_Arma_colSums(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Arma_colSums(x));
    return rcpp_result_gen;
END_RCPP
}
// Sugar_colSums
NumericVector Sugar_colSums(const NumericMatrix& x);
RcppExport SEXP _phisweepr_Sugar_colSums(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Sugar_colSums(x));
    return rcpp_result_gen;
END_RCPP
}
// Cpp_colSums
NumericVector Cpp_colSums(const NumericMatrix& x);
RcppExport SEXP _phisweepr_Cpp_colSums(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpp_colSums(x));
    return rcpp_result_gen;
END_RCPP
}
// Arma_colSums_spM
arma::sp_mat Arma_colSums_spM(const arma::sp_mat& x);
RcppExport SEXP _phisweepr_Arma_colSums_spM(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Arma_colSums_spM(x));
    return rcpp_result_gen;
END_RCPP
}
// theta_from_pi
double theta_from_pi(IntegerVector derivedcount, int nSam, int locusLength);
RcppExport SEXP _phisweepr_theta_from_pi(SEXP derivedcountSEXP, SEXP nSamSEXP, SEXP locusLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type derivedcount(derivedcountSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< int >::type locusLength(locusLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_from_pi(derivedcount, nSam, locusLength));
    return rcpp_result_gen;
END_RCPP
}
// CbTable
arma::sp_mat CbTable(int nSam);
RcppExport SEXP _phisweepr_CbTable(SEXP nSamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    rcpp_result_gen = Rcpp::wrap(CbTable(nSam));
    return rcpp_result_gen;
END_RCPP
}
// subPop_freqSpec
arma::sp_mat subPop_freqSpec(const arma::sp_mat& CbTable, int nSam);
RcppExport SEXP _phisweepr_subPop_freqSpec(SEXP CbTableSEXP, SEXP nSamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type CbTable(CbTableSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    rcpp_result_gen = Rcpp::wrap(subPop_freqSpec(CbTable, nSam));
    return rcpp_result_gen;
END_RCPP
}
// calc_pHomo
NumericVector calc_pHomo(const arma::sp_mat& CbTable, int nSam, int n1Max, int n1Min);
RcppExport SEXP _phisweepr_calc_pHomo(SEXP CbTableSEXP, SEXP nSamSEXP, SEXP n1MaxSEXP, SEXP n1MinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type CbTable(CbTableSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_pHomo(CbTable, nSam, n1Max, n1Min));
    return rcpp_result_gen;
END_RCPP
}
// get_equilibrium2dSfsEntry
NumericVector get_equilibrium2dSfsEntry(const arma::sp_mat& CbTable, int nSam, int n1, int k1, int k2);
RcppExport SEXP _phisweepr_get_equilibrium2dSfsEntry(SEXP CbTableSEXP, SEXP nSamSEXP, SEXP n1SEXP, SEXP k1SEXP, SEXP k2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type CbTable(CbTableSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    rcpp_result_gen = Rcpp::wrap(get_equilibrium2dSfsEntry(CbTable, nSam, n1, k1, k2));
    return rcpp_result_gen;
END_RCPP
}
// standardNeutralJointSFSArmaCube
arma::cube standardNeutralJointSFSArmaCube(const arma::sp_mat& CbTable, NumericVector pHomo, int n1Min, int n1Max, int nSam, double theta, bool monomorphic, bool fixedDerived);
RcppExport SEXP _phisweepr_standardNeutralJointSFSArmaCube(SEXP CbTableSEXP, SEXP pHomoSEXP, SEXP n1MinSEXP, SEXP n1MaxSEXP, SEXP nSamSEXP, SEXP thetaSEXP, SEXP monomorphicSEXP, SEXP fixedDerivedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type CbTable(CbTableSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pHomo(pHomoSEXP);
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< bool >::type monomorphic(monomorphicSEXP);
    Rcpp::traits::input_parameter< bool >::type fixedDerived(fixedDerivedSEXP);
    rcpp_result_gen = Rcpp::wrap(standardNeutralJointSFSArmaCube(CbTable, pHomo, n1Min, n1Max, nSam, theta, monomorphic, fixedDerived));
    return rcpp_result_gen;
END_RCPP
}
// listMatricesToArmaCube
arma::cube listMatricesToArmaCube(Rcpp::List sfslist, int n1Min, int n1Max);
RcppExport SEXP _phisweepr_listMatricesToArmaCube(SEXP sfslistSEXP, SEXP n1MinSEXP, SEXP n1MaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sfslist(sfslistSEXP);
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    rcpp_result_gen = Rcpp::wrap(listMatricesToArmaCube(sfslist, n1Min, n1Max));
    return rcpp_result_gen;
END_RCPP
}
// compute_phiS_field
arma::field<arma::cube> compute_phiS_field(int n1Min, int n1Max, int nSam, const arma::cube& neutralSFS, const NumericVector& alphaD);
RcppExport SEXP _phisweepr_compute_phiS_field(SEXP n1MinSEXP, SEXP n1MaxSEXP, SEXP nSamSEXP, SEXP neutralSFSSEXP, SEXP alphaDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type neutralSFS(neutralSFSSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alphaD(alphaDSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_phiS_field(n1Min, n1Max, nSam, neutralSFS, alphaD));
    return rcpp_result_gen;
END_RCPP
}
// returnCubeFromField
arma::cube returnCubeFromField(const arma::field<arma::cube>& inField, int index);
RcppExport SEXP _phisweepr_returnCubeFromField(SEXP inFieldSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::field<arma::cube>& >::type inField(inFieldSEXP);
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(returnCubeFromField(inField, index));
    return rcpp_result_gen;
END_RCPP
}
// vectorCubes
std::vector<arma::cube> vectorCubes(int n1Min, int n1Max, int nSam, const arma::cube& neutralSFS, const NumericVector& alphaD);
RcppExport SEXP _phisweepr_vectorCubes(SEXP n1MinSEXP, SEXP n1MaxSEXP, SEXP nSamSEXP, SEXP neutralSFSSEXP, SEXP alphaDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type neutralSFS(neutralSFSSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alphaD(alphaDSEXP);
    rcpp_result_gen = Rcpp::wrap(vectorCubes(n1Min, n1Max, nSam, neutralSFS, alphaD));
    return rcpp_result_gen;
END_RCPP
}
// mDimVector
Rcpp::NumericVector mDimVector(NumericVector dimensions);
RcppExport SEXP _phisweepr_mDimVector(SEXP dimensionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dimensions(dimensionsSEXP);
    rcpp_result_gen = Rcpp::wrap(mDimVector(dimensions));
    return rcpp_result_gen;
END_RCPP
}
// listMatrices3dVector
Rcpp::NumericVector listMatrices3dVector(Rcpp::List sfslist, int n1Min, int n1Max);
RcppExport SEXP _phisweepr_listMatrices3dVector(SEXP sfslistSEXP, SEXP n1MinSEXP, SEXP n1MaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sfslist(sfslistSEXP);
    Rcpp::traits::input_parameter< int >::type n1Min(n1MinSEXP);
    Rcpp::traits::input_parameter< int >::type n1Max(n1MaxSEXP);
    rcpp_result_gen = Rcpp::wrap(listMatrices3dVector(sfslist, n1Min, n1Max));
    return rcpp_result_gen;
END_RCPP
}
// groupInfo
arma::mat groupInfo(const arma::sp_mat& geno_Matrix, int coreIdx, int n1min, int n1max, int nSam);
RcppExport SEXP _phisweepr_groupInfo(SEXP geno_MatrixSEXP, SEXP coreIdxSEXP, SEXP n1minSEXP, SEXP n1maxSEXP, SEXP nSamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type geno_Matrix(geno_MatrixSEXP);
    Rcpp::traits::input_parameter< int >::type coreIdx(coreIdxSEXP);
    Rcpp::traits::input_parameter< int >::type n1min(n1minSEXP);
    Rcpp::traits::input_parameter< int >::type n1max(n1maxSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    rcpp_result_gen = Rcpp::wrap(groupInfo(geno_Matrix, coreIdx, n1min, n1max, nSam));
    return rcpp_result_gen;
END_RCPP
}
// probEscape_Sample_C
double probEscape_Sample_C(int n, int k, double alpha, double d, double beta);
RcppExport SEXP _phisweepr_probEscape_Sample_C(SEXP nSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP dSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(probEscape_Sample_C(n, k, alpha, d, beta));
    return rcpp_result_gen;
END_RCPP
}
// p_jH_C
double p_jH_C(NumericVector pvec, int j, int H, int n);
RcppExport SEXP _phisweepr_p_jH_C(SEXP pvecSEXP, SEXP jSEXP, SEXP HSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pvec(pvecSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(p_jH_C(pvec, j, H, n));
    return rcpp_result_gen;
END_RCPP
}
// p_Phi_Selection_second_term_inner_left
double p_Phi_Selection_second_term_inner_left(int n1, int k1, int n2, int k2, const arma::sp_mat& ptable, int i);
RcppExport SEXP _phisweepr_p_Phi_Selection_second_term_inner_left(SEXP n1SEXP, SEXP k1SEXP, SEXP n2SEXP, SEXP k2SEXP, SEXP ptableSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type ptable(ptableSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(p_Phi_Selection_second_term_inner_left(n1, k1, n2, k2, ptable, i));
    return rcpp_result_gen;
END_RCPP
}
// p_Phi_Selection_second_term_inner_right
double p_Phi_Selection_second_term_inner_right(int n1, int k1, int n2, int k2, const arma::sp_mat& ptable, int i);
RcppExport SEXP _phisweepr_p_Phi_Selection_second_term_inner_right(SEXP n1SEXP, SEXP k1SEXP, SEXP n2SEXP, SEXP k2SEXP, SEXP ptableSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type ptable(ptableSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(p_Phi_Selection_second_term_inner_right(n1, k1, n2, k2, ptable, i));
    return rcpp_result_gen;
END_RCPP
}
// phi_S_alphad_C
double phi_S_alphad_C(int n1, int k1, int n2, int k2, const arma::sp_mat& ptable, double alphad, double beta);
RcppExport SEXP _phisweepr_phi_S_alphad_C(SEXP n1SEXP, SEXP k1SEXP, SEXP n2SEXP, SEXP k2SEXP, SEXP ptableSEXP, SEXP alphadSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type ptable(ptableSEXP);
    Rcpp::traits::input_parameter< double >::type alphad(alphadSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_S_alphad_C(n1, k1, n2, k2, ptable, alphad, beta));
    return rcpp_result_gen;
END_RCPP
}
// phi_S_alphad_lookupGenerator_C
NumericVector phi_S_alphad_lookupGenerator_C(int n1, NumericVector k1, int n2, NumericVector k2, const arma::sp_mat& ptable, NumericVector alphad, NumericVector beta);
RcppExport SEXP _phisweepr_phi_S_alphad_lookupGenerator_C(SEXP n1SEXP, SEXP k1SEXP, SEXP n2SEXP, SEXP k2SEXP, SEXP ptableSEXP, SEXP alphadSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type ptable(ptableSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphad(alphadSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_S_alphad_lookupGenerator_C(n1, k1, n2, k2, ptable, alphad, beta));
    return rcpp_result_gen;
END_RCPP
}
// test4dArrayReturnFirstSlot
xt::rarray<int> test4dArrayReturnFirstSlot(xt::rarray<int>& t, NumericVector alphaD);
RcppExport SEXP _phisweepr_test4dArrayReturnFirstSlot(SEXP tSEXP, SEXP alphaDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< xt::rarray<int>& >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphaD(alphaDSEXP);
    rcpp_result_gen = Rcpp::wrap(test4dArrayReturnFirstSlot(t, alphaD));
    return rcpp_result_gen;
END_RCPP
}
// arraySubviewTest
xt::rarray<int> arraySubviewTest(xt::rarray<int>& t);
RcppExport SEXP _phisweepr_arraySubviewTest(SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< xt::rarray<int>& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(arraySubviewTest(t));
    return rcpp_result_gen;
END_RCPP
}
// rcppExpandGridFromZero
NumericMatrix rcppExpandGridFromZero(int n1, int n2);
RcppExport SEXP _phisweepr_rcppExpandGridFromZero(SEXP n1SEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(rcppExpandGridFromZero(n1, n2));
    return rcpp_result_gen;
END_RCPP
}
// makePhiSTable
NumericVector makePhiSTable(NumericVector testN1s, int nSam, Rcpp::List sfsTable, NumericVector alphaD);
RcppExport SEXP _phisweepr_makePhiSTable(SEXP testN1sSEXP, SEXP nSamSEXP, SEXP sfsTableSEXP, SEXP alphaDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type testN1s(testN1sSEXP);
    Rcpp::traits::input_parameter< int >::type nSam(nSamSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type sfsTable(sfsTableSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphaD(alphaDSEXP);
    rcpp_result_gen = Rcpp::wrap(makePhiSTable(testN1s, nSam, sfsTable, alphaD));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_phisweepr_Arma_colSums", (DL_FUNC) &_phisweepr_Arma_colSums, 1},
    {"_phisweepr_Sugar_colSums", (DL_FUNC) &_phisweepr_Sugar_colSums, 1},
    {"_phisweepr_Cpp_colSums", (DL_FUNC) &_phisweepr_Cpp_colSums, 1},
    {"_phisweepr_Arma_colSums_spM", (DL_FUNC) &_phisweepr_Arma_colSums_spM, 1},
    {"_phisweepr_theta_from_pi", (DL_FUNC) &_phisweepr_theta_from_pi, 3},
    {"_phisweepr_CbTable", (DL_FUNC) &_phisweepr_CbTable, 1},
    {"_phisweepr_subPop_freqSpec", (DL_FUNC) &_phisweepr_subPop_freqSpec, 2},
    {"_phisweepr_calc_pHomo", (DL_FUNC) &_phisweepr_calc_pHomo, 4},
    {"_phisweepr_get_equilibrium2dSfsEntry", (DL_FUNC) &_phisweepr_get_equilibrium2dSfsEntry, 5},
    {"_phisweepr_standardNeutralJointSFSArmaCube", (DL_FUNC) &_phisweepr_standardNeutralJointSFSArmaCube, 8},
    {"_phisweepr_listMatricesToArmaCube", (DL_FUNC) &_phisweepr_listMatricesToArmaCube, 3},
    {"_phisweepr_compute_phiS_field", (DL_FUNC) &_phisweepr_compute_phiS_field, 5},
    {"_phisweepr_returnCubeFromField", (DL_FUNC) &_phisweepr_returnCubeFromField, 2},
    {"_phisweepr_vectorCubes", (DL_FUNC) &_phisweepr_vectorCubes, 5},
    {"_phisweepr_mDimVector", (DL_FUNC) &_phisweepr_mDimVector, 1},
    {"_phisweepr_listMatrices3dVector", (DL_FUNC) &_phisweepr_listMatrices3dVector, 3},
    {"_phisweepr_groupInfo", (DL_FUNC) &_phisweepr_groupInfo, 5},
    {"_phisweepr_probEscape_Sample_C", (DL_FUNC) &_phisweepr_probEscape_Sample_C, 5},
    {"_phisweepr_p_jH_C", (DL_FUNC) &_phisweepr_p_jH_C, 4},
    {"_phisweepr_p_Phi_Selection_second_term_inner_left", (DL_FUNC) &_phisweepr_p_Phi_Selection_second_term_inner_left, 6},
    {"_phisweepr_p_Phi_Selection_second_term_inner_right", (DL_FUNC) &_phisweepr_p_Phi_Selection_second_term_inner_right, 6},
    {"_phisweepr_phi_S_alphad_C", (DL_FUNC) &_phisweepr_phi_S_alphad_C, 7},
    {"_phisweepr_phi_S_alphad_lookupGenerator_C", (DL_FUNC) &_phisweepr_phi_S_alphad_lookupGenerator_C, 7},
    {"_phisweepr_test4dArrayReturnFirstSlot", (DL_FUNC) &_phisweepr_test4dArrayReturnFirstSlot, 2},
    {"_phisweepr_arraySubviewTest", (DL_FUNC) &_phisweepr_arraySubviewTest, 1},
    {"_phisweepr_rcppExpandGridFromZero", (DL_FUNC) &_phisweepr_rcppExpandGridFromZero, 2},
    {"_phisweepr_makePhiSTable", (DL_FUNC) &_phisweepr_makePhiSTable, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_phisweepr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
