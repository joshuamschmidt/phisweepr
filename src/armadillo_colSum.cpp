#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

// [[Rcpp::export]]
NumericVector Sugar_colSums(const NumericMatrix& x) {
  return colSums(x);
}

// [[Rcpp::export]]
NumericVector Cpp_colSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += x(i, j);
    }
    ans[j] = sum;
  }
  return ans;
}

// [[Rcpp::export]]
arma::sp_mat Arma_colSums_spM(const arma::sp_mat& x) {
  return arma::sum(x, 0);
}

// // [[Rcpp::export]]
// NumericVector Sugar_colSums_spM(const NumericMatrix& x) {
//   return colSums(x);
// }
// 
// // [[Rcpp::export]]
// NumericVector Cpp_colSums(const NumericMatrix& x) {
//   int nr = x.nrow(), nc = x.ncol();
//   NumericVector ans(nc);
//   for (int j = 0; j < nc; j++) {
//     double sum = 0.0;
//     for (int i = 0; i < nr; i++) {
//       sum += x(i, j);
//     }
//     ans[j] = sum;
//   }
//   return ans;
// }
// 



// // [[Rcpp::export]]
// void sum_by_iterator(const arma::sp_mat& x, int c) {
//   for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
//     if (i.col()==c) {
//       Rcpp::Rcout << " row is " << i.row() << std::endl;
//     }
//   }
// }
// 
// // [[Rcpp::export]]
// arma::sp_mat group(const arma::sp_mat& x, int c) {
//   for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
//     if (i.col()==c) {
//       Rcpp::Rcout << " row is " << i.row() << std::endl;
//     }
//   }
// }


