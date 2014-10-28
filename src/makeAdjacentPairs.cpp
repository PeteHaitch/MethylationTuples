#include <Rcpp.h>
using namespace Rcpp;

// TODO: Fail gracefully if output is likely to exceed R's vector limits.

// NOTE: Could require seqnames and strand to be IntegerVector rather than 
// CharacterVector to possibly save memory but with added minor complication of 
// back-converting integers to seqnames and strand.
// NOTE: Assumes that all loci are either in a feature (in_feature is TRUE or 
// FALSE) or all loci don't have feature information available (in_feature is 
// NA).

//' Get row numbers to create all pairs of adjacent methylation loci.
//' 
//' In the description below, \code{x} is a \code{\link{MethPat}} object 
//' containing 1-tuples. The argument descriptions are chosen for clarity, not 
//' because they are necessarily an efficient way to obtain said arguments.
//'
//' @param methpat_order \code{order(x)}.
//' @param seqnames \code{as.character(seqnames(sort(x)))}.
//' @param strand \code{as.character(sort(x))}.
//' @param pos \code{start(sort(x))}.
//' @param in_feature \code{overlapsAny(x, feature)} or \code{rep(NA, nrow(x))} 
//' if no feature is supplied to \code{methLevelCor}.
//' @param id_dt A \code{\link[data.table]{data.table}} mapping the
//' \code{seqnames-strand-IPD-in_feature} combination to an integer ID.
//' 
//' @keywords internal
//' 
//' @return A \code{list} with the \code{ID}, \code{i} and \code{j}
//' \code{j} of each pair, where \code{i} (resp. \code{j}) is the row number of 
//' the first (resp. second) loci in the pair with respect to \code{x}.
// [[Rcpp::export(".makeAdjacentPairsCpp")]]
List makeAdjacentPairs(IntegerVector methpat_order,
                       std::vector<std::string> seqnames,
                       std::vector<std::string> strand, 
                       IntegerVector pos,
                       LogicalVector in_feature,
                       DataFrame id_dt) {
  
  // Initialise vectors to store results.
  std::vector<int> id_out;
  std::vector<int> i_out;
  std::vector<int> j_out;
  
  // Reserve memory for output vectors.
  // n is an upper bound on the number of pairs created.
  double nr = seqnames.size();
  double n = nr;
  id_out.reserve(n);
  i_out.reserve(n);
  j_out.reserve(n);
  
  // Create id_map from id_dt
  int nid = id_dt.nrows();
  std::map<std::string, int> id_map;
  CharacterVector key = id_dt["KEY"];
  IntegerVector val = id_dt["ID"];
  for (int i = 0; i < nid; i++) {
    String tmp_key = key[i];
    id_map[tmp_key] = val[i];
  }
  
  // Initialise variables used in the for-loop.
  std::string pair_feature_status_string = "";
  
  // Loop over loci and make adjacent pairs.
  for (double i = 0; i < nr - 1; i++) {
    if (seqnames[i] == seqnames[i + 1] and strand[i] == strand[i + 1]) {
      int ipd_ = pos[i + 1] - pos[i];
      // (i, i + 1) are a pair, so add to id_out, i_out and j_out
      std::string ipd_string = Rcpp::toString(ipd_);
      // pair_feature_status_string: 
      // NA; out/out ("0"); in/out or out/in ("1"); in/in ("2")
      // NA occurs when no feature was supplied to methLevelCor
      if (LogicalVector::is_na(in_feature[i])) {
        pair_feature_status_string = "NA";
      } else {
        pair_feature_status_string = Rcpp::toString(in_feature[i] + 
        in_feature[i + 1]);
      }
      std::string id_key = ipd_string + strand[i] + pair_feature_status_string;
      // Look-up id_key in id_map to get the value and store in id_out
      int id_val = id_map[id_key];
      id_out.push_back(id_val);
      i_out.push_back(methpat_order[i]);
      j_out.push_back(methpat_order[i + 1]);
    }
  }
  
  return List::create(_["ID"] = id_out, _["i"] = i_out, _["j"] = j_out); 
}
