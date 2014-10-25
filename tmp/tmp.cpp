#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
std::map<std::string, int> test(DataFrame idDF) {
    int nid = idDF.nrows();
    std::map<std::string, int> idMap;
    CharacterVector key = idDF["key"];
    IntegerVector val = idDF["val"];
    for (int i = 0; i < nid; i++) {
      String tmp_key = key[i];
      idMap[tmp_key] = val[i];
    }
    return idMap;
}


// [[Rcpp::export]]
const bool test2(IntegerVector ipd, int ipd_) {
    //int max_ipd = max(ipd);
    std::set<int> ipd_set(ipd.begin(), ipd.end());
    if (ipd_set.find(ipd_) != ipd_set.end()) {
      return true;
    } else {
      return false;
    }
}

// [[Rcpp::export]]
const bool test3(IntegerVector ipd, int ipd_) {
    //int max_ipd = max(ipd);
    if (std::find(ipd.begin(), ipd.end(), ipd_) != ipd.end()) {
      return true;
    } else {
      return false;
    }
}

// [[Rcpp::export]]
std::vector<int> test4(const std::vector<int>& x) {
  return std::vector<int>(x.begin(), x.end() - 1);
}

// [[Rcpp::export]]
IntegerVector test5(IntegerVector x) {
  IntegerVector y(x.begin(), x.end() - 1);
  return y;
}

// [[Rcpp::export]]
IntegerVector test6(int n) {
  return seq_len(n);
}


