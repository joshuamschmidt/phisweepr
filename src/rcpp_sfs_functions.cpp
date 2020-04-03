// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::sp_mat row_slice(const arma::sp_mat& x, const int n) {
  return x.row(n - 1);
}