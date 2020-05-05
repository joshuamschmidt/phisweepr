#include <Rcpp.h>

using namespace Rcpp;

class ThreeDIndex {
private:
  int nrows, ncols, nmatrices;
  
public:
  ThreeDIndex( int nrows_, int ncols_, int nmatrices_) : nrows(nrows_), 
  ncols(ncols_), nmatrices(nmatrices_){}
  
  int operator()( int i, int j, int k){
    return i + j * nrows + k * ( nrows * ncols ) ;
  }
  
} ;




/// [[Rcpp::export]]
Rcpp::NumericVector listMatrices3dVector(Rcpp::List sfslist, int n1Min, int n1Max) {
  // create a 3d vector
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
  Rcpp::NumericVector outVec = Rcpp::NumericVector( Rcpp::Dimension(n_rows,n_cols,n_slices));
  ThreeDIndex vecIndex( n_rows, n_cols, n_slices ) ;
  for(int s=0; s< n_slices; s++) {
    arma::sp_mat nMatrix = sfslist[s];
    for(int j=0; j < n_cols; j++) {
      for(int i=0; i < n_rows; i++) {
        outVec[ vecIndex(i,j,s) ] = nMatrix(i,j);
      }
    }
  }
  return outVec;
}
