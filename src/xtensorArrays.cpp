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
double probEscape_Single_C(double alpha, double beta, double d){
  return(1-beta*exp(-alpha*std::abs(d)));
}

// [[Rcpp::export]]
NumericVector probEscape_Single_C_plot(double alpha, double beta, NumericVector d){
  NumericVector outvec(d.size(), 0);
  for(int i=0;i<d.size();i++){
    outvec(i) = 1-beta*exp(-alpha*std::abs(d(i)));
  }
  return(outvec);
}



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
double p_Phi_Selection_second_term_inner_right(int n1, int k1, int n2, int k2, NumericMatrix ptable, int i){
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
NumericVector phi_S_alphad_lookupGenerator_C(int n1, NumericVector k1, int n2, NumericVector k2, NumericMatrix ptable, double alphad, double beta){
  int permSize = k1.size();
  NumericVector outputVector(k1.size(), 0);
  for (int i=0; i<permSize; i++) {
    outputVector(i) = phi_S_alphad_C(n1, k1(i), n2, k2(i), ptable, alphad, beta);
  }
  return(outputVector);
}


// [[Rcpp::export]]
xt::rarray<double> phi_S_alphad_lookupGenerator_Cx(int n1, NumericVector k1, int n2, NumericVector k2, NumericMatrix ptable, double alphad, double beta){
  int permSize = k1.size();
  xt::xarray<double> outputVector = xt::zeros<double>({k1.size()});
  for (int i=0; i<permSize; i++) {
    outputVector(i) = phi_S_alphad_C(n1, k1(i), n2, k2(i), ptable, alphad, beta);
  }
  // unsigned long int n_rows = n2+1;
  // unsigned long int n_cols = n1+1;
  // std::array<std::size_t, 2> sh1 = {n_rows, n_cols};;
  // outputVector.reshape(sh1);
  return(outputVector);
}



// [[Rcpp::export]]
xt::rarray<double> makePhiSTable(int nSam, IntegerVector testN1s, NumericMatrix ptable, NumericVector alphad, double beta) {
   // get dimensions of objects
  int nN1 = testN1s.size();
  //int maxN1 = max(testN1s);
  unsigned long int n_rows = nSam+1;
  unsigned long int n_cols = nSam+1;
  unsigned long int n_matrices = nN1;
  unsigned long int nAlphad = alphad.size();
  // check dimensions and sizes of inputs match
  if(n_rows != ptable.nrow() || ptable.nrow() != ptable.ncol()){
     Rcpp::stop("check that input matrices are square\n and/or dimensions match the sample size");
  }
  // create xarray to store phiS
  const xt::rarray<double>::shape_type& shape = { n_rows, n_cols, n_matrices, nAlphad };
  xt::rarray<double> phiSout = xt::zeros<double>(shape);
  // claculate phiS and
  // NumericVector nPhiSVector(1);
   for(int n =0; n < 6; n++ ) {
     int n1 = testN1s[n];
     int n2 = nSam - n1;
     NumericMatrix k1k2 = rcppExpandGridFromZero(n1, n2);
     NumericVector k1 = k1k2( _ , 0 );
     NumericVector k2 = k1k2( _ , 1 );
     for(int i=0; i <nAlphad;i++){
       // vector comes out col by col....
       NumericVector alphadPhiS = phi_S_alphad_lookupGenerator_C(n1,k1,n2,k2,ptable,alphad[i], beta);
       // loop j k over the vector to fill the
       // for(int j=0;j<=n2;j++){
       //   int offset = 0;
       //   if(j > 0){
       //    offset = j*n1+1;
       //   }
       //   for(int k=0;k<=n1;k++){
       //     phiSout(j,k,n,i) = alphadPhiS(offset+k);
       //   }
       // }
       int offset = 0;
       for(int k=0;k<=n1;k++){
         if(k > 0){
           offset = k*(n2+1);
         }
         for(int j=0;j<=n2;j++){
           phiSout(j,k,n,i) = alphadPhiS(offset+j);
         }
       }
       // //phiSout( all(), all(), n, i) = NumericMatrix(n2+1, n1+1, alphadPhiS.begin());
       // unsigned long int sub_n_rows = nSam+1-n1;
       // unsigned long int sub_n_cols = nSam+1-n2;
       // const xt::rarray<double>::shape_type& shape = { sub_n_rows, sub_n_cols };
       // auto submatrix = xt::view(phiSout, xt::all(), xt::all(), n ,i);
       // submatrix = xt::adapt(alphadPhiS, shape);
       // //xt::view(phiSout,  xt::all(),xt::all(),n, i) = NumericMatrix(n2+1, n1+1, alphadPhiS.begin());
     }
  }
  return(phiSout);
}


