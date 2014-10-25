#include <Rcpp.h>
using namespace Rcpp;

// TODO: Why am I using std::vector<> in input rather than CharacterVector etc.

// NOTE: Could require seqnames and strand to be IntegerVector rather than 
// CharacterVector to possibly save memory but with added minor complication of 
// back-converting integers to seqnames and strand

//' Create all pairs of beta-values with given IPDs.
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
//' @param ipd An integer vector of IPD, e.g., \code{ipd = 1:100}.
//' @param betas \code{betaVal(x)}
//' @param id_dt A \code{\link[data.table]{data.table}} mapping the
//' \code{seqnames-strand-IPD-feature_status} combination to an integer ID.
//' 
//' @keywords internal
//' 
//' @return A \code{list} with the \code{ID}, \code{sample}, \code{beta1} and 
//' \code{beta2} of each pair.
// [[Rcpp::export(".makeAllPairsCpp")]]
List makeAllPairs(IntegerVector methpat_order,
                  std::vector<std::string> seqnames,
                  std::vector<std::string> strand, 
                  IntegerVector pos,
                  LogicalVector feature_status,
                  IntegerVector ipd,
                  NumericMatrix betas,
                  DataFrame id_dt) {
  
  // Initialise vectors to store results.
  std::vector<int> id_out;
  std::vector<double> beta1_out;
  std::vector<double> beta2_out;
  
  // Reserve memory for output vectors.
  // Hard to estimate this given each loci can be part of multiple pairs and 
  // the number of pairs depends on the genome and the ipd vector.
  // n is an initial guess that assumes each loci is involved in 30 pairs 
  // (which is probably an underestimate but better than no estimate).
  int nr = seqnames.size();
  int nc = betas.ncol();
  int n = nr * nc * 30;
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
  
  // Initialise variables used in the for-loop.
  int j = 0;
  int max_ipd = max(ipd);
  
  // Loop over loci and find pairs with IPD %in% ipd and on the same 
  // seqname and strand.
  for (int i = 0; i < (nr - 1); i++) {
    j = i + 1;
    while (j <= (nr - 1)) {
      if (seqnames[i] == seqnames[j] and strand[i] == strand[j]) {
        int ipd_ = pos[j] - pos[i];
        if (std::find(ipd.begin(), ipd.end(), ipd_) != ipd.end()) {
          // (i, j) are a pair, so add to id_out, beta1_out and beta2_out
          std::string ipd_string = Rcpp::toString(ipd_);
          // pair_feature_status: out/out (0); in/out or out/in (1); in/in (2)
          std::string pair_feature_status_string = 
          Rcpp::toString(feature_status[i] + feature_status[j]);
          std::string id_key = ipd_string + strand[i] + 
                               pair_feature_status_string;
          // Look-up id_key in id_map to get the value and store in id_out
          int id_val = id_map[id_key];
          // Loop over samples and extract beta values for pair
          for (int k = 0; k < nc; k++) {
            id_out.push_back(id_val);
            beta1_out.push_back(betas(methpat_order[i] - 1, k));
            beta2_out.push_back(betas(methpat_order[j] - 1, k));
          }
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
  
  // Add sample names to output
  IntegerVector sample_names = rep(seq_len(nc), id_out.size() / nc);
  
  return List::create(_["ID"] = id_out, _["sample"] = sample_names, 
  _["beta1"] = beta1_out, _["beta2"] = beta2_out);  
}
