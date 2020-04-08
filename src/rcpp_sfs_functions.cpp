#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double theta_from_pi(const IntegerVector& derivedcount, int nSam, int locusLength) {
  int i;
  double frequency = 0.0;
  double theta = 0.0; // theta is in this case estimated as pi
  for(i=0; i<derivedcount.size(); i++){	
    frequency = derivedcount(i)/static_cast<double>(nSam);
    theta += 2.0*frequency*(1-frequency);
  }
  // theta/seq_length normalises theta summed over sites to theta per site.	
  theta = theta/locusLength;
  return theta;
}

// [[Rcpp::export]]
arma::sp_mat CbTable(int nSam) {
  //calculate combination(n,k)
  // literally just a table of binomial coefficients of n choose k
  int i, k, m, n;
  float C;
  arma::sp_mat Cb(2 * nSam + 1, 2 * nSam + 1);
  for (n = 0; n <= 2 * nSam; n++)
    for (k = 0; k <= n / 2; k++) {
      m = n - k;
      if (m == 0 || k == 0) {
        C = 1.0;
      } else if (m < 0 || n == 0) {
        C = 0.0;
      } else {
        C = 1.0;
        for (i = m + 1; i <= n; i++) {
          C = C * i / (i - m);
        }
      }
      Cb(n, k) = C;
    }
    for (n = 0; n <= 2 * nSam; n++) {
      for (k = n / 2 + 1; k <= 2 * nSam; k++) {
        if (k > n) {
          Cb(n, k) = 0.0;
        } else Cb(n, k) = Cb(n, n - k);
      }
    }
    return (Cb);
}

// [[Rcpp::export]]
arma::sp_mat subPop_freqSpec(const arma::sp_mat& CbTable, int nSam){
  //calculate neutral sampling probability of having k derived alleles in sample of size n 
  //(n<nSam due to missing data)
  //subP is the neutral SFS
  // can replicate as 1/(1:nSam)
  // doesn't have fixed derived?
  int i, j, n;
  arma::sp_mat subP(nSam + 1, nSam + 1); 
  for (n = 1; n <= nSam; n++) {
    if (n == nSam){
      for (j = 1; j < n; j++){
        subP(n, j) = 1.0 / j;
      }
    }
    else
      for (j = 1; j < n; j++) {
        float pnj = 0.0;
        for (i = j; i < nSam; i++)
          if ((nSam - i - n + j) >= 0)
            pnj += CbTable(i, j) * CbTable(nSam - i, n - j) / CbTable(nSam, n) / i;
          subP(n, j) = pnj;
      }
  }
  return(subP);
}

