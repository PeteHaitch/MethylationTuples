#include <Rcpp.h>
using namespace Rcpp;

// TODO: Fail gracefully if output is likely to exceed R's vector limits.

// NOTE: Could require seqnames and strand to be IntegerVector rather than 
// CharacterVector to possibly save memory but with added minor complication of 
// back-converting integers to seqnames and strand
// NOTE: Assumes that all loci are either in a feature (in_feature is TRUE or 
// FALSE) or all loci don't have feature information available (in_feature is 
// NA).

//' Get row numbers to create all pairs of methylation loci with given IPDs.
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
//' if no feature is supplied to \code{betaCor}.
//' \code{rep(NA_integer_, nrow(x))} if no feature is supplied to 
//' \code{betaCor}.
//' @param ipd An integer vector of IPD, e.g., \code{ipd = 1:100}.
//' @param id_dt A \code{\link[data.table]{data.table}} mapping the
//' \code{seqnames-strand-IPD-in_feature} combination to an integer ID.
//' 
//' @keywords internal
//' 
//' @return A \code{list} with the \code{ID}, \code{i} and \code{j}
//' \code{j} of each pair, where \code{i} (resp. \code{j}) is the row number of 
//' the first (resp. second) loci in the pair with respect to \code{x}.
// [[Rcpp::export(".makeAllPairsCpp")]]
List makeAllPairs(IntegerVector methpat_order,
                  std::vector<std::string> seqnames,
                  std::vector<std::string> strand, 
                  IntegerVector pos,
                  LogicalVector in_feature,
                  IntegerVector ipd,
                  DataFrame id_dt) {
  
  // Initialise vectors to store results.
  std::vector<int> id_out;
  std::vector<double> i_out;
  std::vector<double> j_out;
  
  // Reserve memory for output vectors.
  // Hard to estimate this given each loci can be part of multiple pairs and 
  // the number of pairs depends on the genome and the ipd vector.
  // n is an initial guess that assumes each loci is involved in 50 pairs, 
  // which is simply based on my (limited) experience. 
  // TODO: There should be a better way to estimate "50".
  int nr = seqnames.size();
  int n = nr * 50;
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
  int j = 0;
  int max_ipd = max(ipd);
  std::string pair_feature_status_string = "";
  
  // Loop over loci and find pairs with IPD %in% ipd and on the same 
  // seqname and strand.
  for (int i = 0; i < (nr - 1); i++) {
    j = i + 1;
    while (j <= (nr - 1)) {
      if (seqnames[i] == seqnames[j] and strand[i] == strand[j]) {
        int ipd_ = pos[j] - pos[i];
        if (std::find(ipd.begin(), ipd.end(), ipd_) != ipd.end()) {
          // (i, j) are a pair, so add to id_out, i_out and j_out
          std::string ipd_string = Rcpp::toString(ipd_);
          // pair_feature_status_string: 
          // NA; out/out ("0"); in/out or out/in ("1"); in/in ("2")
          // NA occurs when no feature was supplied to betaCor
          if (LogicalVector::is_na(in_feature[i])) {
            pair_feature_status_string = "NA";
          } else {
            pair_feature_status_string = Rcpp::toString(in_feature[i] + 
                                                        in_feature[j]);
          }
          std::string id_key = ipd_string + strand[i] + 
                               pair_feature_status_string;
          // Look-up id_key in id_map to get the value and store in id_out
          int id_val = id_map[id_key];
          id_out.push_back(id_val);
          // Get (i, j) for the pair
          i_out.push_back(methpat_order[i]);
          j_out.push_back(methpat_order[j]);
          j += 1;
        } else if (ipd_ < max_ipd) {
          // (i, j) not a pair but (i, j + 1) might be, so keep looking.
          j += 1;
        } else {
          // (i, j) not a pair and (i, j + 1) can't be.
          // So move onto pairs with (i + 1) as the first co-ordinate.
          break;
        }
      } else {
        // (i, j) not a pair and (i, j + 1) can't be.
        // So move onto pairs with (i + 1) as the first co-ordinate.
        break;
      }
    }
  }
  
  return List::create(_["ID"] = id_out, _["i"] = i_out, _["j"] = j_out);  
}
