#include <Rcpp.h>
using namespace Rcpp;

// TODO: Remove this function; it's just a placeholder so that 
// devtools::load_all(".") will work.

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}
