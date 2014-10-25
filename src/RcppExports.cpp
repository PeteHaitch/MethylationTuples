// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// makeAdjacentPairs
List makeAdjacentPairs(IntegerVector methpat_order, std::vector<std::string> seqnames, std::vector<std::string> strand, IntegerVector pos, LogicalVector feature_status, NumericMatrix betas, DataFrame id_dt);
RcppExport SEXP MethylationTuples_makeAdjacentPairs(SEXP methpat_orderSEXP, SEXP seqnamesSEXP, SEXP strandSEXP, SEXP posSEXP, SEXP feature_statusSEXP, SEXP betasSEXP, SEXP id_dtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type methpat_order(methpat_orderSEXP );
        Rcpp::traits::input_parameter< std::vector<std::string> >::type seqnames(seqnamesSEXP );
        Rcpp::traits::input_parameter< std::vector<std::string> >::type strand(strandSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP );
        Rcpp::traits::input_parameter< LogicalVector >::type feature_status(feature_statusSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type betas(betasSEXP );
        Rcpp::traits::input_parameter< DataFrame >::type id_dt(id_dtSEXP );
        List __result = makeAdjacentPairs(methpat_order, seqnames, strand, pos, feature_status, betas, id_dt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// makeAllPairs
List makeAllPairs(IntegerVector methpat_order, std::vector<std::string> seqnames, std::vector<std::string> strand, IntegerVector pos, LogicalVector feature_status, IntegerVector ipd, NumericMatrix betas, DataFrame id_dt);
RcppExport SEXP MethylationTuples_makeAllPairs(SEXP methpat_orderSEXP, SEXP seqnamesSEXP, SEXP strandSEXP, SEXP posSEXP, SEXP feature_statusSEXP, SEXP ipdSEXP, SEXP betasSEXP, SEXP id_dtSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type methpat_order(methpat_orderSEXP );
        Rcpp::traits::input_parameter< std::vector<std::string> >::type seqnames(seqnamesSEXP );
        Rcpp::traits::input_parameter< std::vector<std::string> >::type strand(strandSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP );
        Rcpp::traits::input_parameter< LogicalVector >::type feature_status(feature_statusSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type ipd(ipdSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type betas(betasSEXP );
        Rcpp::traits::input_parameter< DataFrame >::type id_dt(id_dtSEXP );
        List __result = makeAllPairs(methpat_order, seqnames, strand, pos, feature_status, ipd, betas, id_dt);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
