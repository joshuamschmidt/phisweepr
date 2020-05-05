#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]
double theta_from_pi( IntegerVector derivedcount, int nSam, int locusLength) {
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

// [[Rcpp::export]]
NumericVector calc_pHomo(const arma::sp_mat& CbTable, int nSam, int n1Max, int n1Min){
  int n1size = n1Max - n1Min + 1;
  NumericVector pHomo(n1size);
  for(int n=0; n< n1size; n++) {
    int n1=n+n1Min;
    int n2=nSam-n1;
    pHomo[n] = 0.0;
    for(int i=0; i<=n1; i++){
      for(int j=0; j<=n2; j++) {
        if( (i+j)!=0 && (i+j)!=nSam ) {
          pHomo[n] += CbTable(n1,i)*CbTable(n2,j)/CbTable(nSam,i+j)/(i+j);
        }
      }
    }
  }
  return(pHomo);
}

// [[Rcpp::export]]
NumericVector get_equilibrium2dSfsEntry(const arma::sp_mat& CbTable, int nSam, int n1, int k1, int k2) {
  int ktt = k1 + k2;
  int n2 = nSam - n1;
  NumericVector sfs_entry(1);
  sfs_entry = CbTable(n1,k1)*CbTable(n2,k2)/CbTable(nSam,ktt)/ktt;
  return(sfs_entry);
}

// [[Rcpp::export]]
arma::cube standardNeutralJointSFSArmaCube(const arma::sp_mat& CbTable, NumericVector pHomo, int n1Min, int n1Max, int nSam, double theta, bool monomorphic, bool fixedDerived) {
  int n1Length = n1Max - n1Min + 1;
  int n_rows = n1Max + 1; // rows hold k2
  int n_cols = n_rows; // cols hold k1
  int n_slices = n1Length; // length range of n1 values to be tested
  arma::cube sfsCube(n_rows, n_cols, n_slices);
  for(int i=0; i< n_slices; i++) {
    arma::sp_mat nMatrix(n_rows,n_cols);
    int n1 = i + n1Min;
    int n2 = nSam - n1;
    for(int k1=0; k1<= n1; k1++) {
      for(int k2=0; k2<= n2; k2++) {
        int ktt = k1+k2;
        if( ktt==0 || ktt==nSam ){
          nMatrix(k2,k1) = 1.0-theta*pHomo[i];
        }
        else{
          nMatrix(k2,k1) = theta*CbTable(n1,k1)*CbTable(n2,k2)/CbTable(nSam,ktt)/ktt;
        }
      }
    }
    if(fixedDerived == FALSE) {
      nMatrix(n2,n1) = 0.0;
    }
    if(monomorphic == FALSE) {
      nMatrix(0,0) = 0.0;
    }
    nMatrix = nMatrix/accu(nMatrix);
    sfsCube.slice(i) = nMatrix;
  }
  return(sfsCube);
}

// [[Rcpp::export]]
arma::cube listMatricesToArmaCube(Rcpp::List sfslist, int n1Min, int n1Max) {
  // create a cube structre (3d array)
  // slices are n1 counts, scaled so that the first slice is n1Min
  // colums are k1, rows are k2 counts.
  // so n1, k1, k2 => cube(k2, k1, n1)
  int sfslistSize = sfslist.size();
  int n_rows = n1Max + 1; // rows hold k2
  int n_cols = n_rows; // cols hold k1
  int n_slices = n1Max - n1Min + 1; // length range of n1 values to be tested
  // error check? is slice dimension size correct
  if (n_slices != sfslistSize) {
    stop("check that n1 range and length of sfslist match");
  }
  arma::cube sfsCube(n_rows, n_cols, n_slices);
  for(int i=0; i< n_slices; i++) {
    arma::sp_mat nMatrix = sfslist[i];
    sfsCube.slice(i) = nMatrix;
  }
  //arma::sp_mat M = sfslist[15];
  return sfsCube;
}

// [[Rcpp::export]]
arma::field<arma::cube> compute_phiS_field( int n1Min, int n1Max, int nSam, const arma::cube& neutralSFS, const NumericVector& alphaD) {
  int nAlphaD = alphaD.size();
  int n_slices = neutralSFS.n_slices;
  int n_rows = n1Max + 1; // rows hold k2
  int n_cols = n_rows; // cols hold k1
  int n1Length = n1Max - n1Min + 1;
  if (n_slices != n1Length) {
    stop("check that n1 range defined by n1Min and n1Max\n and n.slices of neutralSFS match");
  }
 // define field of length nAlphaD.
  arma::field<arma::cube> outField(nAlphaD);
  outField.fill(arma::cube(n_rows, n_cols, n_slices, arma::fill::randu));
  // each index of the field is an alphaD/ outField[0] ~ alphaD 0.001 etc.
  // contains a cube of k2 rows, k1 cols, n1 slices.
  // loop over alphaD; loop over cubes and calculate the new probability!
  // put new cube into field.
  //outField(0) = neutralSFS;
  //outField(1) = sqrt(neutralSFS);
  // phi_S_alphad_lookupGenerator_C
  return(outField);
}


// [[Rcpp::export]]
arma::cube returnCubeFromField(const arma::field<arma::cube>& inField, int index){
  arma::cube outCube = inField(0);
  return outCube;
}

// [[Rcpp::export]]
std::vector<arma::cube> vectorCubes(int n1Min, int n1Max, int nSam, const arma::cube& neutralSFS, const NumericVector& alphaD) {
  int nAlphaD = alphaD.size();
  int n_slices = neutralSFS.n_slices;
  int n_rows = n1Max + 1; // rows hold k2
  int n_cols = n_rows; // cols hold k1
  int n1Length = n1Max - n1Min + 1;
  if (n_slices != n1Length) {
    stop("check that n1 range defined by n1Min and n1Max\n and n.slices of neutralSFS match");
  }
  // define field of length nAlphaD.
  std::vector<arma::cube> outVecCube(nAlphaD, arma::cube(n_rows, n_cols, n_slices, arma::fill::randu));
  return outVecCube;
}


// [[Rcpp::export]]
Rcpp::NumericVector mDimVector(NumericVector dimensions){
  Rcpp::NumericVector vec = Rcpp::NumericVector( Rcpp::Dimension(dimensions));
  //int dims = vec.attr("dim");
  return vec;
}
  

