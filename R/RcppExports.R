# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Arma_colSums <- function(x) {
    .Call(`_phisweepr_Arma_colSums`, x)
}

Sugar_colSums <- function(x) {
    .Call(`_phisweepr_Sugar_colSums`, x)
}

Cpp_colSums <- function(x) {
    .Call(`_phisweepr_Cpp_colSums`, x)
}

Arma_colSums_spM <- function(x) {
    .Call(`_phisweepr_Arma_colSums_spM`, x)
}

theta_from_pi <- function(derivedcount, nSam, locusLength) {
    .Call(`_phisweepr_theta_from_pi`, derivedcount, nSam, locusLength)
}

CbTable <- function(nSam) {
    .Call(`_phisweepr_CbTable`, nSam)
}

subPop_freqSpec <- function(CbTable, nSam) {
    .Call(`_phisweepr_subPop_freqSpec`, CbTable, nSam)
}

calc_pHomo <- function(CbTable, nSam, n1Max, n1Min) {
    .Call(`_phisweepr_calc_pHomo`, CbTable, nSam, n1Max, n1Min)
}

get_equilibrium2dSfsEntry <- function(CbTable, nSam, n1, k1, k2) {
    .Call(`_phisweepr_get_equilibrium2dSfsEntry`, CbTable, nSam, n1, k1, k2)
}

standardNeutralJointSFSArmaCube <- function(CbTable, pHomo, n1Min, n1Max, nSam, theta, monomorphic, fixedDerived) {
    .Call(`_phisweepr_standardNeutralJointSFSArmaCube`, CbTable, pHomo, n1Min, n1Max, nSam, theta, monomorphic, fixedDerived)
}

listMatricesToArmaCube <- function(sfslist, n1Min, n1Max) {
    .Call(`_phisweepr_listMatricesToArmaCube`, sfslist, n1Min, n1Max)
}

compute_phiS_field <- function(n1Min, n1Max, nSam, neutralSFS, alphaD) {
    .Call(`_phisweepr_compute_phiS_field`, n1Min, n1Max, nSam, neutralSFS, alphaD)
}

returnCubeFromField <- function(inField, index) {
    .Call(`_phisweepr_returnCubeFromField`, inField, index)
}

vectorCubes <- function(n1Min, n1Max, nSam, neutralSFS, alphaD) {
    .Call(`_phisweepr_vectorCubes`, n1Min, n1Max, nSam, neutralSFS, alphaD)
}

mDimVector <- function(dimensions) {
    .Call(`_phisweepr_mDimVector`, dimensions)
}

listMatrices3dVector <- function(sfslist, n1Min, n1Max) {
    .Call(`_phisweepr_listMatrices3dVector`, sfslist, n1Min, n1Max)
}

groupInfo <- function(geno_Matrix, coreIdx, n1min, n1max, nSam) {
    .Call(`_phisweepr_groupInfo`, geno_Matrix, coreIdx, n1min, n1max, nSam)
}

probEscape_Sample_C <- function(n, k, alpha, d, beta) {
    .Call(`_phisweepr_probEscape_Sample_C`, n, k, alpha, d, beta)
}

p_jH_C <- function(pvec, j, H, n) {
    .Call(`_phisweepr_p_jH_C`, pvec, j, H, n)
}

p_Phi_Selection_second_term_inner_left <- function(n1, k1, n2, k2, ptable, i) {
    .Call(`_phisweepr_p_Phi_Selection_second_term_inner_left`, n1, k1, n2, k2, ptable, i)
}

p_Phi_Selection_second_term_inner_right <- function(n1, k1, n2, k2, ptable, i) {
    .Call(`_phisweepr_p_Phi_Selection_second_term_inner_right`, n1, k1, n2, k2, ptable, i)
}

phi_S_alphad_C <- function(n1, k1, n2, k2, ptable, alphad, beta) {
    .Call(`_phisweepr_phi_S_alphad_C`, n1, k1, n2, k2, ptable, alphad, beta)
}

phi_S_alphad_lookupGenerator_C <- function(n1, k1, n2, k2, ptable, alphad, beta) {
    .Call(`_phisweepr_phi_S_alphad_lookupGenerator_C`, n1, k1, n2, k2, ptable, alphad, beta)
}

test4dArrayReturnFirstSlot <- function(t, alphaD) {
    .Call(`_phisweepr_test4dArrayReturnFirstSlot`, t, alphaD)
}

arraySubviewTest <- function(t) {
    .Call(`_phisweepr_arraySubviewTest`, t)
}

rcppExpandGridFromZero <- function(n1, n2) {
    .Call(`_phisweepr_rcppExpandGridFromZero`, n1, n2)
}

makePhiSTable <- function(testN1s, nSam, sfsTable, alphaD) {
    .Call(`_phisweepr_makePhiSTable`, testN1s, nSam, sfsTable, alphaD)
}

