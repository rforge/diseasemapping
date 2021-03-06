// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// multiplyLower
SEXP multiplyLower(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 B, int Nglobal0, int Nlocal0, int NlocalCache, std::string type);
RcppExport SEXP _gpuMatrix_multiplyLower(SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP Nglobal0SEXP, SEXP Nlocal0SEXP, SEXP NlocalCacheSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type Nglobal0(Nglobal0SEXP);
    Rcpp::traits::input_parameter< int >::type Nlocal0(Nlocal0SEXP);
    Rcpp::traits::input_parameter< int >::type NlocalCache(NlocalCacheSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplyLower(C, A, B, Nglobal0, Nlocal0, NlocalCache, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpuMatrix_multiplyLower", (DL_FUNC) &_gpuMatrix_multiplyLower, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpuMatrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
