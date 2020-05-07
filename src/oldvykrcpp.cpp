#define STRICT_R_HEADERS
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(xtensor)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <xtensor/xmath.hpp>
#include <xtensor-r/rarray.hpp>
#include <xtensor-r/rtensor.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xslice.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xshape.hpp>
#include <Rcpp.h>

using namespace Rcpp;

// // [[Rcpp::export]]
// arma::field<arma::cube> fourDObject_as_armaFieldofCubes(int N1Length, int k1Length, int k2Length, int alphaDLength) {
//   arma::field<arma::cube> A(lengthNRange);
//   
//   A.fill(arma::cube(n, L, max(N), arma::fill::zeros));
//   
//   return A;
// }

;

// [[Rcpp::export]]
double probEscape_Sample_C(int n, int k, double alpha, double d, double beta){
  double returnValue = 0;
  double prob_escape = 1-beta*exp(-alpha*std::abs(d));
  returnValue = R::dbinom(k, n, prob_escape, false);
  return(returnValue);
}

// [[Rcpp::export]]
double p_Phi_Selection_second_term_inner_right(int n1, int k1, int n2, int k2,NumericMatrix ptable, int i){
  double returnValue = 0;
  if (1+i-k1 <= 0) {
    return(returnValue);
  }
  double first_term = R::choose(double(i)+1, k1) / R::choose(n2+i+1, k1+k2) * R::choose(n2, k2);
  double second_term = (i+1-k1)/(double(i)+1);
  returnValue = first_term * second_term * ptable(n2+i+1,k1+k2);
  return(returnValue);
}

// [[Rcpp::export]]
double p_Phi_Selection_second_term_inner_left(int n1, int k1, int n2, int k2, NumericMatrix ptable,int i){
  double returnValue = 0;
  if (k1+1-n1+i <= 0) {
    return(returnValue);
  }
  double first_term = R::choose(double(i)+1, k1+1-n1+i) / R::choose(n2+i+1, k1+1-n1+i+k2) * R::choose(n2, k2);
  double second_term = (k1+1-n1+i)/(1+double(i));
  returnValue = first_term * second_term * ptable(n2+i+1, k1+1-n1+i+k2);
  return(returnValue);
}

// [[Rcpp::export]]
double p_jH_C(NumericVector pvec, int j, int H, int n){
  int pvec_size = pvec.size();
  if (pvec_size != n+1) return(-9999);
  double returnValue = 0;
  if (j<0 || j>H) {
    return(returnValue);
  }
  for (int i=j; i<=n; ++i) {
    returnValue += pvec(i)* (R::dhyper(j, i, n-i, H, false));
  }
  return(returnValue);
}

// [[Rcpp::export]]
double phi_S_alphad_C(int n1, int k1, int n2, int k2, NumericMatrix ptable, double alphad, double beta){
  int ntt,ktt,i;
  double p1,p2,p_total,p_outer,p_inner;
  
  ntt = n1+n2;
  ktt = k1+k2;
  p_outer = probEscape_Sample_C(n1, n1, alphad, 1, beta) * ((R::choose(double(n1), k1) * R::choose(double(n2), k2)) / R::choose(double(ntt), ktt)) * double(ptable(ntt,ktt));
  p_inner = 0;
  for(i=0; i< n1; i++) {
    p1=0;
    p2=0;
    if((k1+1+i-n1) > 0) {
      p1 = p_Phi_Selection_second_term_inner_left(n1, k1, n2, k2,ptable, i);
    }
    if((1+i-k1) > 0) {
      p2 = p_Phi_Selection_second_term_inner_right(n1, k1, n2, k2,ptable, i);
    }
    p_inner += probEscape_Sample_C(n1, i, alphad, 1, beta)*(p1 + p2);
  }
  p_total = p_outer + p_inner;
  return(p_total);
}


// [[Rcpp::export]]
NumericVector phi_S_alphad_lookupGenerator_C(int n1, int nSam, NumericVector k1, NumericVector k2, NumericMatrix ptable, NumericVector alphad, NumericVector beta){
  int permSize = alphad.size();
  NumericVector outputVector(alphad.size(), 0);
  for (int i=0; i<permSize; i++) {
    outputVector(i) = phi_S_alphad_C(n1, nSam, k1(i), k2(i), ptable, alphad(i), beta(i));
  }
  return(outputVector);
}

// [[Rcpp::export]]
xt::rarray<int> generate_phiS_4darray(NumericMatrix ptable, int nSam, int minN1, int maxN1, NumericVector alphad, int beta){
  int permSize = alphad.size();
  
  rcppExpandGridFromZero;
  phi_S_alphad_lookupGenerator_C
}
  
  
  
  
  
  
  