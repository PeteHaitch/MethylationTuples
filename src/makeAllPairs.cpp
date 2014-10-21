#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".makeAllPairsCpp")]]

List makeAllPairs(IntegerVector s, int max_ipd) {
  
  // Initiate vectors to store results.
  std::vector<int> xx;
  std::vector<int> yy;

  // Create variables used in the for-loop.
  int n = s.size();
  int j = 0;
  int k = 1;
  
  // Loop over start positions, s, and find pairs with ipd <= max_ipd
  for (int i = 0; i < (n - 1); i++) {
    j = i + 1;
    while (j <= (n - 1)) {
      /* 
      Ideally, would use (s[j] - s[i]) %in% ipd but there is no %in% sugar 
      operator. Instead, use this monstrosity. 
      Note need to wrap in is_true() (see 
      http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-May/005807.html)
      */
      if ( (s[j] - s[i]) <= max_ipd) {
        
        /* 
        (i, j) are a pair, so add to x and y.
        Record (i + 1, j + 1) because R vectors are 1-based and C++ vectors 
        are 0-based. 
        */
        xx.push_back(i + 1); 
        yy.push_back(j + 1);
        j += 1;
        k += 1;
      } else {
        /*
        (i, j) aren't a pair and (i, j + 1) can't be, so move onto pairs 
        with x = i + 1.
        */
        break;
      }
    }
  }
  
  return List::create(_["xx"] = xx, _["yy"] = yy);
  
}
