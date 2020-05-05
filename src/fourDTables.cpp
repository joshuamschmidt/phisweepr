// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// //#include "FourDimension.h"
// using namespace Rcpp;
// 
// class FourDIndex {
// private:
//   int nrows, ncols, nmatrices, nslices;
// 
// public:
//   FourDIndex( IntegerVector dim ) : nrows(dim[0]), ncols(dim[1]),
//   nmatrices(dim[2]), nslices(dim[3]){}
// 
//   FourDIndex( int nrows_, int ncols_, int nmatrices_, int nslices_) :
//     nrows(nrows_), ncols(ncols_), nmatrices(nmatrices_), nslices(nslices_){}
// 
//   int operator()( int i, int j, int k, int l){
//     return i + j * nrows + k * ( nrows * ncols ) + l * ( nrows * ncols * nmatrices);
//   }
// 
// } ;
// 
// class Dimfour {
// public:
//   typedef std::vector<int>::reference reference ;
//   typedef std::vector<int>::const_reference const_reference ;
//   
//   Dimfour() : dims(){}
//   
//   Dimfour(SEXP dims) ;
//   
//   Dimfour( const Dimfour& other ) : dims(other.dims){}
//   Dimfour& operator=( const Dimfour& other ) {
//     if( *this != other )
//       dims = other.dims ;
//     return *this ;
//   }
//   Dimfour(const size_t& n1, const size_t& n2, const size_t& n3, const size_t& n4) : dims(4){
//     dims[0] = static_cast<int>(n1) ;
//     dims[1] = static_cast<int>(n2) ;
//     dims[2] = static_cast<int>(n3) ;
//     dims[3] = static_cast<int>(n4) ;
//   }
//   operator SEXP() const ;
//   
//   inline int size() const {
//     return (int) dims.size() ;
//   }
//   inline R_xlen_t prod() const {
//     return std::accumulate( dims.begin(), dims.end(), static_cast<R_xlen_t>(1), std::multiplies<R_xlen_t>() );
//   }
//   inline reference operator[](int i){
//     if( i < 0 || i>=static_cast<int>(dims.size()) ) throw std::range_error("index out of bounds") ;
//     return dims[i] ;
//   }
//   inline const_reference operator[](int i) const{
//     if( i < 0 || i>=static_cast<int>(dims.size()) ) throw std::range_error("index out of bounds") ;
//     return dims[i] ;
//   }
//   
// private:
//   std::vector<int> dims;
// };
// 
// 
// 
// 
// 
// // [[Rcpp::export]]
// Rcpp::NumericVector listMatrices4dVector(Rcpp::List sfslist, int n1Min, int n1Max, NumericVector alphaD) {
//   // create a 4d vector
//   // slices are n1 counts, scaled so that the first slice is n1Min
//   // colums are k1, rows are k2 counts. alphaD, indexed so need a way to convert from alphaD and its value.
//   // so n1, k1, k2 => cube(k2, k1, n1, alphaD)
//   int sfslistSize = sfslist.size();
//   int n_rows = n1Max + 1; // rows hold k2
//   int n_cols = n_rows; // cols hold k1
//   int n_matrices = n1Max - n1Min + 1; // length range of n1 values to be tested
//   int n_slices = alphaD.size();
//   // error check? is slice dimension size correct
//   if (n_matrices != sfslistSize) {
//     stop("check that n1 range and length of sfslist match");
//   }
//   Rcpp::NumericVector outVec = Rcpp::NumericVector( Dimfour::Dimfour( n_rows, n_cols, n_matrices, n_slices ));
//   FourDIndex vecIndex( n_rows, n_cols, n_matrices, n_slices ) ;
//   for(int l=0; l<n_slices; l++) {
//     for(int k=0; k< n_matrices; k++) {
//       arma::sp_mat nMatrix = sfslist[k];
//       for(int j=0; j < n_cols; j++) {
//         for(int i=0; i < n_rows; i++) {
//           outVec[ vecIndex(i,j,k,l) ] = nMatrix(i,j)*l;
//         }
//       }
//     }
//   }
//   return outVec;
// }
