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

// [[Rcpp::export]]
NumericMatrix groupInfo( NumericMatrix geno_Matrix, int coreIdx) {
  int i,j,groups,nSites,nSam;
  groups = 2;
  // change idx from 1-based to 0-based array
  coreIdx = coreIdx -1;
  // allocate space of group info 2-d array
  nSam = geno_Matrix.nrow();
  nSites = geno_Matrix.ncol();
  NumericMatrix groupInfo(nSites, groups);
  
  // loop over the sample haplotypes            
  for(i=0;i<nSam;i++){
    if(geno_Matrix(i,coreIdx)==1)
      for(j=0; j<nSites; j++) groupInfo(j,0) += geno_Matrix(i,j);	// get count of k derived alleles over all other sites.
    else
      for(j=0; j<nSites; j++) groupInfo(j,1) += geno_Matrix(i,j);
  }
  return(groupInfo);
}




// phiS porbability functions
// [[Rcpp::export]]
double probEscape_Sample_C(int n, int k, double alpha, double d, double beta){
  double returnValue = 0;
  double prob_escape = 1-beta*exp(-alpha*std::abs(d));
  returnValue = R::dbinom(k, n, prob_escape, false);
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
double p_Phi_Selection_second_term_inner_left(int n1, int k1, int n2, int k2, const arma::sp_mat& ptable,int i){
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
double p_Phi_Selection_second_term_inner_right(int n1, int k1, int n2, int k2,const arma::sp_mat& ptable, int i){
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
double phi_S_alphad_C(int n1, int k1, int n2, int k2, const arma::sp_mat& ptable, double alphad, double beta){
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
NumericVector phi_S_alphad_lookupGenerator_C(int n1,NumericVector k1,int n2,NumericVector k2, const arma::sp_mat& ptable, NumericVector alphad, NumericVector beta){
  int permSize = alphad.size();
  NumericVector outputVector(alphad.size(), 0);
  for (int i=0; i<permSize; i++) {
    outputVector(i) = phi_S_alphad_C(n1, k1(i), n2, k2(i), ptable, alphad(i), beta(i));
  }
  return(outputVector);
}


// // [[Rcpp::export]]
// xt::rarray<int> test4dArray(xt::rarray<int>& t) {
//   t(0,0,0,0) = 90;
//   return t;
// }
// 
// // [[Rcpp::export]]
// xt::rarray<int> test4dArrayIndex(xt::rarray<int>& t) {
//   auto select = t(0,0,0,0);
//   return select;
// }

// [[Rcpp::export]]
xt::rarray<int> test4dArrayReturnFirstSlot(xt::rarray<int>& t, NumericVector alphaD) {
  int dim = t.dimension();
  if(dim!=4){
    Rcpp::stop("check that that input array has four dimensions!");
  }
  //xt::rarray<int> subt = xt::view(t, xt::all(),xt::all(),xt::all(), 0);
  const int n_rows = t.shape(0);
  int n_cols = t.shape(1);
  int n_matrices = t.shape(2);
  int n_slices = t.shape(3);
  int nAlphD = alphaD.size();
  if(nAlphD != n_slices){
    Rcpp::stop("check that that alphD length array has four dimensions!");
  }
  for(int s=0;s<n_slices;s++){
    for(int m=0;m<n_matrices;m++){
      for(int c=0;c<n_cols;c++){
        for(int r=0;r<n_rows;r++){
          t(r,c,m,s) = t(r,c,m,s)+1;
        }
      }
    }
  }
  // Rcpp::Rcout << "t has " << dim << " dimensions" << std::endl;
  // Rcpp::Rcout << "t has " << n_rows << " rows" << std::endl;
  // Rcpp::Rcout << "t has " << n_cols << " cols" << std::endl;
  // Rcpp::Rcout << "t has " << n_matrices << " matrices" << std::endl;
  // Rcpp::Rcout << "t has " << n_slices << " slices" << std::endl;
  return t;
}

// [[Rcpp::export]]
xt::rarray<int> arraySubviewTest(xt::rarray<int>& t) {
  xt::rarray<int> subt = xt::view(t,  1,xt::all(),xt::all(), 0);
  return subt;
}

// [[Rcpp::export]]
NumericMatrix rcppExpandGridFromZero(int n1, int n2) {
  // creates every pairwise combination of 0:n1,0:n2
  int combinations = (n1+1) * (n2+1);
  NumericMatrix k1k2(combinations,2);
  int jstart = 0;
  for(int i=0;i<=n1;i++){
    if(i!=0) {
      jstart = i * (n2+1);
    }
    for(int j=0; j<=n2;j++) {
      k1k2(jstart+j,0) = i;
      k1k2(jstart+j,1) = j;
    }
  }
  return k1k2;
}

//xt::rarray<double>

// [[Rcpp::export]]
NumericVector makePhiSTable(NumericVector testN1s, int nSam, Rcpp::List sfsTable, NumericVector alphaD) {
  // grab expand.grid function from R
  Rcpp::Function expGrid("expand.grid");
  // get dimensions of objects
  int nN1 = testN1s.size();
  int maxN1 = max(testN1s);
  arma::sp_mat mat = sfsTable[0];
  unsigned long int n_rows = mat.n_rows;
  unsigned long int n_cols = mat.n_cols;
  unsigned long int n_matrices = sfsTable.size();
  unsigned long int nAlphaD = alphaD.size();
  // check dimensions and sizes of inputs match
  if(n_rows != (maxN1+1) || n_rows != n_cols){
     Rcpp::stop("check that input matrices are square\n and/or dimensions macth max testN1s");
  }
  if(n_matrices!= nN1){
    Rcpp::stop("check that n matrices in sfsTable (dim[3]) == length(testN1s)!");
  }
  //Rcout << " all dims ok!" << std::endl;
  
  // create xarray to store phiS
  const xt::rarray<double>::shape_type& shape = { n_rows, n_cols, n_matrices, nAlphaD };
  xt::rarray<double> phiSout(shape);
  // claculate phiS and
  NumericVector nPhiSVector(1);
  for(int n =0; n <1; n++ ) {
    int n1 = testN1s[n];
    int n2 = nSam - n1;
    NumericVector beta = {1}; // phiS this is hard coded = 1, to preserve matching beta to partialSweepFinder.
    NumericMatrix k1k2 = rcppExpandGridFromZero(n1, n2);
    NumericVector k1 = k1k2( _ , 0 );
    NumericVector k2 = k1k2( _ , 1 );
    arma::sp_mat n1mat = sfsTable[n];
    Rcout << n1 << " is the current n1 " << std::endl;
    nPhiSVector = phi_S_alphad_lookupGenerator_C(n1,k1,n2,k2,n1mat,alphaD, 1);
  }
  return nPhiSVector;
}


