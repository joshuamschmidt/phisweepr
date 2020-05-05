// #include <Rcpp.h>
// using namespace Rcpp ;
// 
// class Offset{
// private:
//   int nrows, ncols, nmats ;
//   
// public:
//   Offset( IntegerVector dim ) : nrows(dim[0]), ncols(dim[1]), 
//   nmats(dim[2]){}
//   
//   Offset( int nrows_, int ncols_, int nmats_) : nrows(nrows_), 
//   ncols(ncols_), nmats(nmats_){}
//   
//   int operator()( int i, int j, int k){
//     return i + j * nrows + k * ( nrows * ncols ) ;
//   }
//   
// } ;
// 
// class Array3 : public NumericVector {
// private:
//   Offset offset ;
//   
// public:
//   Array3( SEXP x) : NumericVector(x), offset( 
//       (IntegerVector)((RObject)x).attr("dim") ) { }
//   Array3( Dimension dim ): NumericVector( dim ), offset( dim[0], 
//           dim[1], dim[2] ) {}
//   Array3( int nrows_, int ncols_, int nmats_ ): NumericVector( 
//       Dimension(nrows_, ncols_, nmats_ ) ), offset( nrows_, ncols_, nmats_ ) {}
//   
//   inline double& operator()( int i, int j, int k){
//     return ( (NumericVector*)(this) )->operator[]( offset( i, j, k) 
//     ) ;
//   }
//   
// } ;
// // [[Rcpp::export]]
// Array3 foo(Array3 arr){
//   arr(0,0,0) = 1.0 ;
//   arr(1,1,1) = 2.0 ;
//   return arr ;
// }
// 
// /*** R
// a <- array( 0, dim = c(4,5,6 ) )
// foo(a) 
// */