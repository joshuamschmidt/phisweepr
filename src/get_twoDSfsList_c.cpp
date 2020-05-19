#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector seq_c(int l){
  NumericVector x (1);
  for (int i=1;i<l;i++)
    x.push_back(i);
  return(x);
}


// [[Rcpp::export]]
// loop over sites to calculate the observed 2d-SFS
arma::sp_mat get_two_dimensionalSFSlist_C(List dataObject,int n1Min, int n1Max, char monomorphic, char fixedDerived){
  NumericVector dCounts = dataObject["derived_allele_counts"];
  //Rcout << "made it here " << std::endl;
  arma::sp_mat genotypes = dataObject["genotypes"];
  int nSam = dataObject["sample.size"];
  arma::sp_mat nMatrix(n1Max+1,n1Max+1);
  for(int n=n1Min; n <= n1Min; n++) {
     //arma::sp_mat nMatrix(n1Max+1,n1Max+1);
     //int n2 = nSam - n;
     for(int i=0; i<dCounts.size();i++){
       // i is index in dCounts and genotypes matrix
       if(dCounts(i)==n){
         // find indxs of rows of derived alleles
        NumericVector d_idxs;
        for(arma::sp_mat::const_iterator it = genotypes.begin_col(i); it != genotypes.end_col(i); ++it){
          d_idxs.push_back(it.row());
          }
        NumericVector a_idxs = setdiff(seq_c(nSam),d_idxs);
        std::sort(a_idxs.begin(), a_idxs.end());
        for(int j=0; j<genotypes.n_cols;j++){
          if(j!=i){
            int k1 = 0;
            int k2 = 0;
            for(int d_idx =0;d_idx< d_idxs.size();d_idx++){
              k1 += genotypes(d_idxs(d_idx),j);
            }
            for(int a_idx =0;a_idx< a_idxs.size();a_idx++){
              k2 += genotypes(a_idxs(a_idx),j);
            }
            // nMatrix is k2 rows, k1 cols
            nMatrix(k2,k1)+=1;
          }
         }
       }
     }
   }
  // return(dCounts);
  // int i = 0;
  // NumericVector d_idxs;
  // for(arma::sp_mat::const_iterator it = genotypes.begin_col(i); it != genotypes.end_col(i); ++it){
  //   d_idxs.push_back(it.row());
  // }
  // NumericVector a_idxs = setdiff(seq_c(nSam),d_idxs);
  // std::sort(a_idxs.begin(), a_idxs.end());
  // //arma::sp_mat sub = genotypes.col(i);
  //arma::uvec A = find(genotypes.col(i) < 1.0);
  //arma::sp_mat B = genotypes.rows(A);
  return(nMatrix);
}

// [[Rcpp::export]]
arma::sp_mat get_two_dimensionalSFS_C(const arma::sp_mat& genotypes,int n1, int nSam, NumericVector n_idxs, char monomorphic, char fixedDerived){
  int n2 = nSam - n1;
  arma::sp_mat nMatrix(n2+1,n1+1);
  for(int i=0; i<n_idxs.size();i++){
    // find indxs of rows of derived alleles
    int idx = n_idxs[i]-1;
    NumericVector d_idxs;
    for(arma::sp_mat::const_iterator it = genotypes.begin_col(idx); it != genotypes.end_col(idx); ++it){
      d_idxs.push_back(it.row());
    }
    NumericVector a_idxs = setdiff(seq_c(nSam),d_idxs);
    std::sort(a_idxs.begin(), a_idxs.end());
    for(int j=0; j<genotypes.n_cols;j++){
      if(j!=idx){
        int k1 = 0;
        int k2 = 0;
        for(int d_idx =0;d_idx< d_idxs.size();d_idx++){
          k1 += genotypes(d_idxs(d_idx),j);
        }
        for(int a_idx =0;a_idx< a_idxs.size();a_idx++){
          k2 += genotypes(a_idxs(a_idx),j);
        }
        // nMatrix is k2 rows, k1 cols
        nMatrix(k2,k1)+=1;
      }
    }
  }
  return(nMatrix);
}

