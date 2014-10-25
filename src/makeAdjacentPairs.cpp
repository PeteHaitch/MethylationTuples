#include <Rcpp.h>
using namespace Rcpp;

// NOTE: Could require seqnames and strand to be IntegerVector rather than 
// CharacterVector to possibly save memory but with added minor complication of 
// back-converting integers to seqnames and strand

//' Create adjacent pairs of beta-values.
//' 
//' In the description below, \code{x} is a \code{\link{MethPat}} object 
//' containing 1-tuples. The argument descriptions are chosen for clarity, not 
//' because they are necessarily an efficient way to obtain said arguments.
//'
//' @param methpat_order \code{order(x)}.
//' @param seqnames \code{as.character(seqnames(sort(x)))}.
//' @param strand \code{as.character(sort(x))}.
//' @param pos \code{start(sort(x))}.
//' @param feature_status \code{overlapsAny(x, feature)}.
//' @param betas \code{betaVal(x)}
//' @param id_dt A \code{\link[data.table]{data.table}} mapping the
//' \code{seqnames-strand-IPD-feature_status} combination to an integer ID.
//' 
//' @keywords internal
//' 
//' @return A \code{list} with the \code{ID}, \code{sample}, \code{beta1} and 
//' \code{beta2} of each pair.
// [[Rcpp::export(".makeAdjacentPairsCpp")]]
List makeAdjacentPairs(IntegerVector methpat_order,
                       std::vector<std::string> seqnames,
                       std::vector<std::string> strand, 
                       IntegerVector pos,
                       LogicalVector feature_status,
                       NumericMatrix betas,
                       DataFrame id_dt) {
  
  // Initialise vectors to store results.
  std::vector<int> id_out;
  std::vector<double> beta1_out;
  std::vector<double> beta2_out;
  
  // Reserve memory for output vectors.
  // n is an upper bound on the number of pairs created.
  int nr = seqnames.size();
  int nc = betas.ncol();
  int n = nr * nc;
  id_out.reserve(n);
  beta1_out.reserve(n);
  beta2_out.reserve(n);
  
  // Create id_map from id_dt
  int nid = id_dt.nrows();
  std::map<std::string, int> id_map;
  CharacterVector key = id_dt["KEY"];
  IntegerVector val = id_dt["ID"];
  for (int i = 0; i < nid; i++) {
    String tmp_key = key[i];
    id_map[tmp_key] = val[i];
  }
    
  for (int i = 0; i < nr - 1; i++) {
    if (seqnames[i] == seqnames[i + 1] and strand[i] == strand[i + 1]) {
      int ipd_ = pos[i + 1] - pos[i];
      // (i, i + 1) are a pair, so add to id_out, beta1_out and beta2_out
      std::string ipd_string = Rcpp::toString(ipd_);
      // pair_feature_status: out/out (0); in/out or out/in (1); in/in (2)
      std::string pair_feature_status_string = 
        Rcpp::toString(feature_status[i] + feature_status[i + 1]);
      std::string id_key = seqnames[i] + strand[i] + ipd_string + "_" + 
        pair_feature_status_string;
      // Look-up id_key in id_map to get the value and store in id_out
      int id_val = id_map[id_key];
      // Loop over samples and extract beta values for pair
      for (int k = 0; k < nc; k++) {
        id_out.push_back(id_val);
        // -1 because C++ uses 0-based index.
        beta1_out.push_back(betas(methpat_order[i] - 1, k));
        beta2_out.push_back(betas(methpat_order[i + 1] - 1, k));
      }
    }
  }
  
  // Add sample names to output
  IntegerVector sample_names = rep(seq_len(nc), id_out.size() / nc);
  
  return List::create(_["ID"] = id_out, _["sample"] = sample_names, 
  _["beta1"] = beta1_out, _["beta2"] = beta2_out); 
}
