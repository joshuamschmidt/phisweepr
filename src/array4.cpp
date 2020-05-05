// #include <Rcpp.h>
// using namespace Rcpp ;
// 
// // following example from http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-August/006309.html
// 
// class Offset{
// private:
//   int nrows, ncols, nmats , nslcs;
// 
// public:
//   Offset( IntegerVector dim ) : nrows(dim[0]), ncols(dim[1]),
//   nmats(dim[2]), nslcs(dim[3]){}
// 
//   Offset( int nrows_, int ncols_, int nmats_, int nslcs_) : nrows(nrows_),
//   ncols(ncols_), nmats(nmats_), nslcs(nslcs_){}
// 
//   int operator()( int i, int j, int k, int l){
//     return i + j * nrows + k * ( nrows * ncols ) + l * ( nrows * ncols * nmats );
//   }
// 
// } ;
// 
// class Array4 : public NumericVector {
// private:
//   Offset offset ;
// 
// public:
//   Array4( SEXP x) : NumericVector(x), offset((IntegerVector)((RObject)x).attr("dim") ) { }
//   Array4( Dimension dim ): NumericVector( dim ), offset( dim[0], dim[1], dim[2] , dim[3]) {}
//   Array4( int nrows_, int ncols_, int nmats_, int nslcs_ ): NumericVector(Dimension(nrows_, ncols_, nmats_, nslcs_ ) ), offset( nrows_, ncols_, nmats_, nslcs_ ) {}
// 
//   inline double& operator()( int i, int j, int k, int l){
//     return ( (NumericVector*)(this) )->operator[]( offset( i, j, k, l)
//     ) ;
//   }
// 
// } ;
// 
// // [[Rcpp::export]]
// Array4 foo(Array4 arr){
//   arr(0,0,0,0) = 1.0 ;
//   arr(1,1,1,1) = 2.0 ;
//   return arr ;
// }
// 
// /*** R
// a <- array( 0, dim = c(4,5,6, 7) )
// foo(a)
// */