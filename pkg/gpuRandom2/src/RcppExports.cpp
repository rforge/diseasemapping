// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_gpuFisher_test
SEXP cpp_gpuFisher_test(Rcpp::S4 srR, Rcpp::S4 scR, Rcpp::S4 resultsR, Rcpp::IntegerMatrix streamsR, Rcpp::IntegerVector max_global_size, Rcpp::IntegerVector max_local_size);
RcppExport SEXP _gpuRandom_cpp_gpuFisher_test(SEXP srRSEXP, SEXP scRSEXP, SEXP resultsRSEXP, SEXP streamsRSEXP, SEXP max_global_sizeSEXP, SEXP max_local_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type srR(srRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type scR(scRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type resultsR(resultsRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type streamsR(streamsRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_global_size(max_global_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_local_size(max_local_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gpuFisher_test(srR, scR, resultsR, streamsR, max_global_size, max_local_size));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gpuRn
SEXP cpp_gpuRn(Rcpp::S4 xR, Rcpp::IntegerMatrix streamsR, IntegerVector max_global_size, IntegerVector max_local_size, std::string random_type, std::string precision_type);
RcppExport SEXP _gpuRandom_cpp_gpuRn(SEXP xRSEXP, SEXP streamsRSEXP, SEXP max_global_sizeSEXP, SEXP max_local_sizeSEXP, SEXP random_typeSEXP, SEXP precision_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type streamsR(streamsRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_global_size(max_global_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_local_size(max_local_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type random_type(random_typeSEXP);
    Rcpp::traits::input_parameter< std::string >::type precision_type(precision_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gpuRn(xR, streamsR, max_global_size, max_local_size, random_type, precision_type));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mrg31k3pCreateStreams
Rcpp::IntegerMatrix cpp_mrg31k3pCreateStreams(int numWorkItems);
RcppExport SEXP _gpuRandom_cpp_mrg31k3pCreateStreams(SEXP numWorkItemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numWorkItems(numWorkItemsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mrg31k3pCreateStreams(numWorkItems));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpuRandom_cpp_gpuFisher_test", (DL_FUNC) &_gpuRandom_cpp_gpuFisher_test, 6},
    {"_gpuRandom_cpp_gpuRn", (DL_FUNC) &_gpuRandom_cpp_gpuRn, 6},
    {"_gpuRandom_cpp_mrg31k3pCreateStreams", (DL_FUNC) &_gpuRandom_cpp_mrg31k3pCreateStreams, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpuRandom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
