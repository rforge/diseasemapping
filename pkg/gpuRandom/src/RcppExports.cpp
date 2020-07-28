// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// backsolveBatchBackend
SEXP backsolveBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 B, Rcpp::IntegerVector Cstartend, Rcpp::IntegerVector Astartend, Rcpp::IntegerVector Bstartend, const int numbatchB, const int diagIsOne, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_backsolveBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CstartendSEXP, SEXP AstartendSEXP, SEXP BstartendSEXP, SEXP numbatchBSEXP, SEXP diagIsOneSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Cstartend(CstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Astartend(AstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Bstartend(BstartendSEXP);
    Rcpp::traits::input_parameter< const int >::type numbatchB(numbatchBSEXP);
    Rcpp::traits::input_parameter< const int >::type diagIsOne(diagIsOneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(backsolveBatchBackend(C, A, B, Cstartend, Astartend, Bstartend, numbatchB, diagIsOne, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}
// backsolveBatchBackend2
SEXP backsolveBatchBackend2(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 B, const int diagIsOne, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_backsolveBatchBackend2(SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP diagIsOneSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< const int >::type diagIsOne(diagIsOneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(backsolveBatchBackend2(C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}
// cholBatchBackend
void cholBatchBackend(Rcpp::S4 A, Rcpp::S4 D, std::vector<int> Nglobal, std::vector<int> Nlocal, std::vector<int> NlocalCache);
RcppExport SEXP _gpuRandom_cholBatchBackend(SEXP ASEXP, SEXP DSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type D(DSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type NlocalCache(NlocalCacheSEXP);
    cholBatchBackend(A, D, Nglobal, Nlocal, NlocalCache);
    return R_NilValue;
END_RCPP
}
// crossprodBatchBackend
SEXP crossprodBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 D, const int invertD, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_crossprodBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP DSEXP, SEXP invertDSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type invertD(invertDSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(crossprodBatchBackend(C, A, D, invertD, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}
// gemmBatchBackend
SEXP gemmBatchBackend(Rcpp::S4 A, Rcpp::S4 B, Rcpp::S4 C, const int rowbatch, const int Acolbatch, const int Bcolbatch, const int need_transpose, Rcpp::IntegerVector Nglobal);
RcppExport SEXP _gpuRandom_gemmBatchBackend(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP rowbatchSEXP, SEXP AcolbatchSEXP, SEXP BcolbatchSEXP, SEXP need_transposeSEXP, SEXP NglobalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type rowbatch(rowbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type Acolbatch(AcolbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type Bcolbatch(BcolbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type need_transpose(need_transposeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    rcpp_result_gen = Rcpp::wrap(gemmBatchBackend(A, B, C, rowbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gpuFisher_test
SEXP cpp_gpuFisher_test(Rcpp::S4 xR, Rcpp::S4 resultsR, int B, Rcpp::S4 streamsR, Rcpp::IntegerVector max_global_size, Rcpp::IntegerVector max_local_size);
RcppExport SEXP _gpuRandom_cpp_gpuFisher_test(SEXP xRSEXP, SEXP resultsRSEXP, SEXP BSEXP, SEXP streamsRSEXP, SEXP max_global_sizeSEXP, SEXP max_local_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type resultsR(resultsRSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type streamsR(streamsRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_global_size(max_global_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_local_size(max_local_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gpuFisher_test(xR, resultsR, B, streamsR, max_global_size, max_local_size));
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
// gpuRnBackend
SEXP gpuRnBackend(Rcpp::S4 x, Rcpp::S4 streams, IntegerVector max_global_size, std::string random_type);
RcppExport SEXP _gpuRandom_gpuRnBackend(SEXP xSEXP, SEXP streamsSEXP, SEXP max_global_sizeSEXP, SEXP random_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type streams(streamsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type max_global_size(max_global_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type random_type(random_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(gpuRnBackend(x, streams, max_global_size, random_type));
    return rcpp_result_gen;
END_RCPP
}
// cpp_gpu_qqnorm
SEXP cpp_gpu_qqnorm(Rcpp::S4 outR, double mu, double sigma, int lowertail, Rcpp::IntegerVector max_global_size, Rcpp::IntegerVector max_local_size);
RcppExport SEXP _gpuRandom_cpp_gpu_qqnorm(SEXP outRSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP lowertailSEXP, SEXP max_global_sizeSEXP, SEXP max_local_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type outR(outRSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type lowertail(lowertailSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_global_size(max_global_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type max_local_size(max_local_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_gpu_qqnorm(outR, mu, sigma, lowertail, max_global_size, max_local_size));
    return rcpp_result_gen;
END_RCPP
}
// maternBatchBackend
void maternBatchBackend(Rcpp::S4 var, Rcpp::S4 coords, Rcpp::S4 param, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal);
RcppExport SEXP _gpuRandom_maternBatchBackend(SEXP varSEXP, SEXP coordsSEXP, SEXP paramSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type var(varSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    maternBatchBackend(var, coords, param, Nglobal, Nlocal);
    return R_NilValue;
END_RCPP
}
// multiplyLowerDiagonalBatchBackend
SEXP multiplyLowerDiagonalBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 D, Rcpp::S4 B, const int diagIsOne, std::string transformD, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_multiplyLowerDiagonalBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP DSEXP, SEXP BSEXP, SEXP diagIsOneSEXP, SEXP transformDSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< const int >::type diagIsOne(diagIsOneSEXP);
    Rcpp::traits::input_parameter< std::string >::type transformD(transformDSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplyLowerDiagonalBatchBackend(C, A, D, B, diagIsOne, transformD, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}
// multiplyDiagonalBatchBackend
SEXP multiplyDiagonalBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 B, const int inverse, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal);
RcppExport SEXP _gpuRandom_multiplyDiagonalBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP inverseSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< const int >::type inverse(inverseSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplyDiagonalBatchBackend(C, A, B, inverse, Nglobal, Nlocal));
    return rcpp_result_gen;
END_RCPP
}
// multiplyLowerBatchBackend
SEXP multiplyLowerBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 B, const int diagIsOne, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_multiplyLowerBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP diagIsOneSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< const int >::type diagIsOne(diagIsOneSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(multiplyLowerBatchBackend(C, A, B, diagIsOne, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpuRandom_backsolveBatchBackend", (DL_FUNC) &_gpuRandom_backsolveBatchBackend, 11},
    {"_gpuRandom_backsolveBatchBackend2", (DL_FUNC) &_gpuRandom_backsolveBatchBackend2, 7},
    {"_gpuRandom_cholBatchBackend", (DL_FUNC) &_gpuRandom_cholBatchBackend, 5},
    {"_gpuRandom_crossprodBatchBackend", (DL_FUNC) &_gpuRandom_crossprodBatchBackend, 7},
    {"_gpuRandom_gemmBatchBackend", (DL_FUNC) &_gpuRandom_gemmBatchBackend, 8},
    {"_gpuRandom_cpp_gpuFisher_test", (DL_FUNC) &_gpuRandom_cpp_gpuFisher_test, 6},
    {"_gpuRandom_cpp_mrg31k3pCreateStreams", (DL_FUNC) &_gpuRandom_cpp_mrg31k3pCreateStreams, 1},
    {"_gpuRandom_gpuRnBackend", (DL_FUNC) &_gpuRandom_gpuRnBackend, 4},
    {"_gpuRandom_cpp_gpu_qqnorm", (DL_FUNC) &_gpuRandom_cpp_gpu_qqnorm, 6},
    {"_gpuRandom_maternBatchBackend", (DL_FUNC) &_gpuRandom_maternBatchBackend, 5},
    {"_gpuRandom_multiplyLowerDiagonalBatchBackend", (DL_FUNC) &_gpuRandom_multiplyLowerDiagonalBatchBackend, 9},
    {"_gpuRandom_multiplyDiagonalBatchBackend", (DL_FUNC) &_gpuRandom_multiplyDiagonalBatchBackend, 6},
    {"_gpuRandom_multiplyLowerBatchBackend", (DL_FUNC) &_gpuRandom_multiplyLowerBatchBackend, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpuRandom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
