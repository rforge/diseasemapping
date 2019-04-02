// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_cholGpu
SEXP cpp_cholGpu(Rcpp::S4 xR, Rcpp::S4 DR, Rcpp::S4 diagWorkingR, Rcpp::S4 diagTimesRowOfAR, int MCglobal, int MClocal, int localStorage, int colGroupwise, int Ncrossprod, int verbose, std::string kernelR);
RcppExport SEXP _geostatsgpu_cpp_cholGpu(SEXP xRSEXP, SEXP DRSEXP, SEXP diagWorkingRSEXP, SEXP diagTimesRowOfARSEXP, SEXP MCglobalSEXP, SEXP MClocalSEXP, SEXP localStorageSEXP, SEXP colGroupwiseSEXP, SEXP NcrossprodSEXP, SEXP verboseSEXP, SEXP kernelRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DR(DRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type diagWorkingR(diagWorkingRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type diagTimesRowOfAR(diagTimesRowOfARSEXP);
    Rcpp::traits::input_parameter< int >::type MCglobal(MCglobalSEXP);
    Rcpp::traits::input_parameter< int >::type MClocal(MClocalSEXP);
    Rcpp::traits::input_parameter< int >::type localStorage(localStorageSEXP);
    Rcpp::traits::input_parameter< int >::type colGroupwise(colGroupwiseSEXP);
    Rcpp::traits::input_parameter< int >::type Ncrossprod(NcrossprodSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernelR(kernelRSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_cholGpu(xR, DR, diagWorkingR, diagTimesRowOfAR, MCglobal, MClocal, localStorage, colGroupwise, Ncrossprod, verbose, kernelR));
    return rcpp_result_gen;
END_RCPP
}
// gpuNlocal
Rcpp::List gpuNlocal(std::string kernel, std::string functionName, int ctx_id);
RcppExport SEXP _geostatsgpu_gpuNlocal(SEXP kernelSEXP, SEXP functionNameSEXP, SEXP ctx_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< std::string >::type functionName(functionNameSEXP);
    Rcpp::traits::input_parameter< int >::type ctx_id(ctx_idSEXP);
    rcpp_result_gen = Rcpp::wrap(gpuNlocal(kernel, functionName, ctx_id));
    return rcpp_result_gen;
END_RCPP
}
// cpp_maternGpuD
SEXP cpp_maternGpuD(Rcpp::S4 varR, Rcpp::S4 coordsR, Rcpp::S4 DofLDLR, Rcpp::NumericVector param, const int form, Rcpp::IntegerVector numWorkItems, Rcpp::IntegerVector numLocalItems);
RcppExport SEXP _geostatsgpu_cpp_maternGpuD(SEXP varRSEXP, SEXP coordsRSEXP, SEXP DofLDLRSEXP, SEXP paramSEXP, SEXP formSEXP, SEXP numWorkItemsSEXP, SEXP numLocalItemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type varR(varRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type coordsR(coordsRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DofLDLR(DofLDLRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const int >::type form(formSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numWorkItems(numWorkItemsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numLocalItems(numLocalItemsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_maternGpuD(varR, coordsR, DofLDLR, param, form, numWorkItems, numLocalItems));
    return rcpp_result_gen;
END_RCPP
}
// cpp_maternGpuF
SEXP cpp_maternGpuF(Rcpp::S4 varR, Rcpp::S4 coordsR, Rcpp::S4 DofLDLR, Rcpp::NumericVector param, const int form, Rcpp::IntegerVector numWorkItems, Rcpp::IntegerVector numLocalItems);
RcppExport SEXP _geostatsgpu_cpp_maternGpuF(SEXP varRSEXP, SEXP coordsRSEXP, SEXP DofLDLRSEXP, SEXP paramSEXP, SEXP formSEXP, SEXP numWorkItemsSEXP, SEXP numLocalItemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type varR(varRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type coordsR(coordsRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type DofLDLR(DofLDLRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const int >::type form(formSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numWorkItems(numWorkItemsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numLocalItems(numLocalItemsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_maternGpuF(varR, coordsR, DofLDLR, param, form, numWorkItems, numLocalItems));
    return rcpp_result_gen;
END_RCPP
}
// cpp_lu
SEXP cpp_lu(Rcpp::S4 xR, Rcpp::S4 dR);
RcppExport SEXP _geostatsgpu_cpp_lu(SEXP xRSEXP, SEXP dRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type dR(dRSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_lu(xR, dR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_geostatsgpu_cpp_cholGpu", (DL_FUNC) &_geostatsgpu_cpp_cholGpu, 11},
    {"_geostatsgpu_gpuNlocal", (DL_FUNC) &_geostatsgpu_gpuNlocal, 3},
    {"_geostatsgpu_cpp_maternGpuD", (DL_FUNC) &_geostatsgpu_cpp_maternGpuD, 7},
    {"_geostatsgpu_cpp_maternGpuF", (DL_FUNC) &_geostatsgpu_cpp_maternGpuF, 7},
    {"_geostatsgpu_cpp_lu", (DL_FUNC) &_geostatsgpu_cpp_lu, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_geostatsgpu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
