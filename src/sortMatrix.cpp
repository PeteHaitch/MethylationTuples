// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h> 

using namespace Rcpp;

// Read http://thread.gmane.org/gmane.comp.lang.r.rcpp/5818 for a discussion 
// on passing R matrices to Armadillo by reference.

//' Sort rows or columns of a numeric matrix.
//' 
//' This simply uses the \code{sort} routine in \code{Armadillo} via 
//' \code{RcppArmadillo}.
//' 
//' @param x A numeric matrix.
//' @param sort_direction The sort_direction argument is optional; 
//' sort_direction is either "ascend" or "descend"; by default "ascend" is used.
//' @param dim The dim argument is optional; by default dim = 0 is used.
//' 
//' @return A matrix with the elements of the input matrix sorted in each column 
//' (dim = 0), or each row (dim = 1).
//' 
//' @keywords internal
// [[Rcpp::export(".sortMatrixCpp")]]
arma::mat sortMatrix(const arma::mat& x, 
                        const char* sort_direction = "ascend", 
                        const int dim = 0L) {
  return arma::sort(x, sort_direction, dim);
}
