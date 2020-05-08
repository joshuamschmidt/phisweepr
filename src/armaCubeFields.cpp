#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::field<arma::cube> testCubeField(int n_rows, int n_cols, int nN1, NumericVector alphad) {
  int nAlphad = alphad.size();
  arma::field<arma::cube> phiStable(nAlphad);
  phiStable.fill(arma::cube(n_rows, n_cols, nN1, arma::fill::zeros));
  return phiStable;
}

// [[Rcpp::export]]
arma::mat PhiSMatFromPhiSCubeField(arma::field<arma::cube> phiStable, int n_slice, int n_matrix) {
  arma::cube phiSslice = phiStable[n_slice]; 
  arma::mat phiSMat = phiSslice.slice(n_matrix);
  return phiSMat;
}
