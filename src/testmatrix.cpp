#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// // [[Rcpp::export]]
// void testiter(const arma::sp_mat& x) {
//   // Make const iterator
//   //arma::sp_mat::const_iterator i = x.begin();
//   //arma::sp_mat::const_iterator j = x.end();
//   for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
//     if (i.col()==0) {
//       Rcpp::Rcout << "it's row " << i.row() << std::endl;
//     }
//   }
//   
// } 
// [[Rcpp::export]]
arma::mat groupInfo(const arma::sp_mat& geno_Matrix, int coreIdx, int n1min, int n1max, int nSam) {
  //int i,j,groups,nSites;
  // initialise matrix to hold 2-d SFS remember (row,col)
  int outrows = n1max - n1min +1;
  arma::mat outMatrix(outrows, nSam, arma::fill::zeros);
  outMatrix(0,0) = 1;
  outMatrix(0,1) += 5;
  return(outMatrix);
  //return(arma::sum( geno_Matrix, 0 ));
}
//   groups = 2;
//   // change idx from 1-based to 0-based array
//   coreIdx = coreIdx - 1;
//   // allocate space of group info 2-d array
//   nSites = geno_Matrix.n_cols;
//   NumericMatrix groupInfo(nSites, groups);
//   
//   nSam = geno_Matrix.n_rows;
//   // loop over the sample haplotypes            
//   for(i=0;i<nSam;i++){
//     if(geno_Matrix(i,coreIdx)==1)
//       for(j=0; j<nSites; j++) groupInfo(j,0) += geno_Matrix(i,j);	// get count of k derived alleles over all other sites.
//     else
//       for(j=0; j<nSites; j++) groupInfo(j,1) += geno_Matrix(i,j);
//   }
//   return(groupInfo);
// }


// [[Rcpp::export]]
NumericVector test_colSums(arma::sp_mat& geno_Matrix) { 
  size_t cols = geno_Matrix.n_cols; 
  NumericVector res(cols);
  for (size_t i=0; i<cols; i++) { 
    res[i] = sum(geno_Matrix.col(i));
  }
  return(res); 
}