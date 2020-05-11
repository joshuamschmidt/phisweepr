#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector getProbVectorPhiS(arma::field<arma::cube> phitable, int n1, int minN1, NumericVector k1, NumericVector k2, NumericVector sub_alphaD_idx){
  if(k1.size()!=k2.size() || k1.size()!=sub_alphaD_idx.size()){
    stop("vectors for k1, k2 and/or sub_alphaD_idx are not the same length!");
  }
  NumericVector outvec(sub_alphaD_idx.size(), 0);
  int n=n1-minN1;
  arma::cube nCube = phitable[n];
  for(int i=0; i<sub_alphaD_idx.size();i++){
    outvec[i] = nCube(k2[i],k1[i],sub_alphaD_idx[i]);
  }
  return (outvec);
}

// [[Rcpp::export]]
NumericVector getProbVectorPhiN(arma::cube phitable, int n1, int minN1, NumericVector k1, NumericVector k2){
  if(k1.size()!=k2.size()){
    stop("vectors k1 and k2 are not the same length!");
  }
  NumericVector outvec(k1.size(), 0);
  int n=n1-minN1;
  for(int i=0; i<k1.size();i++){
    outvec[i] = phitable(k2[i],k1[i],n);
  }
  return (outvec);
}