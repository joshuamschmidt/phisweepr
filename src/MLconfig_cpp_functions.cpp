#include <Rcpp.h>
using namespace Rcpp;

// ML config cpp functions
// [[Rcpp::export]]
double distance_pESingle_C(double alpha, double beta, double pEscape){
  return(-(log(-(pEscape-1/beta))/alpha));
}
