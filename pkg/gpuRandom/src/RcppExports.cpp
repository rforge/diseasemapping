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
void cholBatchBackend(Rcpp::S4 A, Rcpp::S4 D, Rcpp::IntegerVector Astartend, Rcpp::IntegerVector Dstartend, const int numbatchD, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, Rcpp::IntegerVector NlocalCache);
RcppExport SEXP _gpuRandom_cholBatchBackend(SEXP ASEXP, SEXP DSEXP, SEXP AstartendSEXP, SEXP DstartendSEXP, SEXP numbatchDSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type D(DSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Astartend(AstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Dstartend(DstartendSEXP);
    Rcpp::traits::input_parameter< const int >::type numbatchD(numbatchDSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type NlocalCache(NlocalCacheSEXP);
    cholBatchBackend(A, D, Astartend, Dstartend, numbatchD, Nglobal, Nlocal, NlocalCache);
    return R_NilValue;
END_RCPP
}
// crossprodBatchBackend
SEXP crossprodBatchBackend(Rcpp::S4 C, Rcpp::S4 A, Rcpp::S4 D, const int invertD, Rcpp::IntegerVector Cstartend, Rcpp::IntegerVector Astartend, Rcpp::IntegerVector Dstartend, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, const int NlocalCache);
RcppExport SEXP _gpuRandom_crossprodBatchBackend(SEXP CSEXP, SEXP ASEXP, SEXP DSEXP, SEXP invertDSEXP, SEXP CstartendSEXP, SEXP AstartendSEXP, SEXP DstartendSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type invertD(invertDSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Cstartend(CstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Astartend(AstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Dstartend(DstartendSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< const int >::type NlocalCache(NlocalCacheSEXP);
    rcpp_result_gen = Rcpp::wrap(crossprodBatchBackend(C, A, D, invertD, Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache));
    return rcpp_result_gen;
END_RCPP
}
// gemmBatchBackend
SEXP gemmBatchBackend(Rcpp::S4 A, Rcpp::S4 B, Rcpp::S4 C, const int Arowbatch, const int Browbatch, const int Acolbatch, const int Bcolbatch, const int need_transpose, Rcpp::IntegerVector Nglobal);
RcppExport SEXP _gpuRandom_gemmBatchBackend(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP ArowbatchSEXP, SEXP BrowbatchSEXP, SEXP AcolbatchSEXP, SEXP BcolbatchSEXP, SEXP need_transposeSEXP, SEXP NglobalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type Arowbatch(ArowbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type Browbatch(BrowbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type Acolbatch(AcolbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type Bcolbatch(BcolbatchSEXP);
    Rcpp::traits::input_parameter< const int >::type need_transpose(need_transposeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    rcpp_result_gen = Rcpp::wrap(gemmBatchBackend(A, B, C, Arowbatch, Browbatch, Acolbatch, Bcolbatch, need_transpose, Nglobal));
    return rcpp_result_gen;
END_RCPP
}
// gemmBatch2backend
SEXP gemmBatch2backend(Rcpp::S4 A, Rcpp::S4 B, Rcpp::S4 C, Rcpp::IntegerVector transposeABC, Rcpp::IntegerVector submatrixA, Rcpp::IntegerVector submatrixB, Rcpp::IntegerVector submatrixC, Rcpp::IntegerVector batches, Rcpp::IntegerVector workgroupSize, Rcpp::IntegerVector NlocalCache, const int verbose);
RcppExport SEXP _gpuRandom_gemmBatch2backend(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP transposeABCSEXP, SEXP submatrixASEXP, SEXP submatrixBSEXP, SEXP submatrixCSEXP, SEXP batchesSEXP, SEXP workgroupSizeSEXP, SEXP NlocalCacheSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type transposeABC(transposeABCSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type submatrixA(submatrixASEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type submatrixB(submatrixBSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type submatrixC(submatrixCSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type batches(batchesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type workgroupSize(workgroupSizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type NlocalCache(NlocalCacheSEXP);
    Rcpp::traits::input_parameter< const int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gemmBatch2backend(A, B, C, transposeABC, submatrixA, submatrixB, submatrixC, batches, workgroupSize, NlocalCache, verbose));
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
// likfitGpu_Backend
void likfitGpu_Backend(Rcpp::S4 coordsGpuR, Rcpp::S4 bigparamsBatchR, Rcpp::S4 yXR, Rcpp::S4 betasR, Rcpp::S4 variancesR, Rcpp::S4 jacobianR, Rcpp::S4 finalLogLikR, int n, int p, int groupsize, int colbatch, int form, Rcpp::IntegerVector workgroupSize, Rcpp::IntegerVector localSize, Rcpp::IntegerVector NlocalCache);
RcppExport SEXP _gpuRandom_likfitGpu_Backend(SEXP coordsGpuRSEXP, SEXP bigparamsBatchRSEXP, SEXP yXRSEXP, SEXP betasRSEXP, SEXP variancesRSEXP, SEXP jacobianRSEXP, SEXP finalLogLikRSEXP, SEXP nSEXP, SEXP pSEXP, SEXP groupsizeSEXP, SEXP colbatchSEXP, SEXP formSEXP, SEXP workgroupSizeSEXP, SEXP localSizeSEXP, SEXP NlocalCacheSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type coordsGpuR(coordsGpuRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type bigparamsBatchR(bigparamsBatchRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type yXR(yXRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type betasR(betasRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type variancesR(variancesRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type jacobianR(jacobianRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type finalLogLikR(finalLogLikRSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type groupsize(groupsizeSEXP);
    Rcpp::traits::input_parameter< int >::type colbatch(colbatchSEXP);
    Rcpp::traits::input_parameter< int >::type form(formSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type workgroupSize(workgroupSizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type localSize(localSizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type NlocalCache(NlocalCacheSEXP);
    likfitGpu_Backend(coordsGpuR, bigparamsBatchR, yXR, betasR, variancesR, jacobianR, finalLogLikR, n, p, groupsize, colbatch, form, workgroupSize, localSize, NlocalCache);
    return R_NilValue;
END_RCPP
}
// logfactsumBackend
SEXP logfactsumBackend(Rcpp::S4 xR, Rcpp::IntegerVector numWorkItems);
RcppExport SEXP _gpuRandom_logfactsumBackend(SEXP xRSEXP, SEXP numWorkItemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numWorkItems(numWorkItemsSEXP);
    rcpp_result_gen = Rcpp::wrap(logfactsumBackend(xR, numWorkItems));
    return rcpp_result_gen;
END_RCPP
}
// rowsumBackend
void rowsumBackend(Rcpp::S4 xR, Rcpp::S4 SumR, std::string type, int log);
RcppExport SEXP _gpuRandom_rowsumBackend(SEXP xRSEXP, SEXP SumRSEXP, SEXP typeSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type SumR(SumRSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type log(logSEXP);
    rowsumBackend(xR, SumR, type, log);
    return R_NilValue;
END_RCPP
}
// matrix_matrix_sumBackend
void matrix_matrix_sumBackend(Rcpp::S4 aR, Rcpp::S4 bR, Rcpp::S4 sumR, Rcpp::IntegerVector numWorkItems);
RcppExport SEXP _gpuRandom_matrix_matrix_sumBackend(SEXP aRSEXP, SEXP bRSEXP, SEXP sumRSEXP, SEXP numWorkItemsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type aR(aRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type bR(bRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sumR(sumRSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numWorkItems(numWorkItemsSEXP);
    matrix_matrix_sumBackend(aR, bR, sumR, numWorkItems);
    return R_NilValue;
END_RCPP
}
// matrix_vector_sumBackend
void matrix_vector_sumBackend(Rcpp::S4 matrixR, Rcpp::S4 vectorR, Rcpp::S4 sumR, const int byrow, Rcpp::IntegerVector numWorkItems);
RcppExport SEXP _gpuRandom_matrix_vector_sumBackend(SEXP matrixRSEXP, SEXP vectorRSEXP, SEXP sumRSEXP, SEXP byrowSEXP, SEXP numWorkItemsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type matrixR(matrixRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type vectorR(vectorRSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type sumR(sumRSEXP);
    Rcpp::traits::input_parameter< const int >::type byrow(byrowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numWorkItems(numWorkItemsSEXP);
    matrix_vector_sumBackend(matrixR, vectorR, sumR, byrow, numWorkItems);
    return R_NilValue;
END_RCPP
}
// fillParamsExtra
void fillParamsExtra(Rcpp::S4 param);
RcppExport SEXP _gpuRandom_fillParamsExtra(SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type param(paramSEXP);
    fillParamsExtra(param);
    return R_NilValue;
END_RCPP
}
// maternBatchBackend
void maternBatchBackend(Rcpp::S4 var, Rcpp::S4 coords, Rcpp::S4 param, Rcpp::IntegerVector Nglobal, Rcpp::IntegerVector Nlocal, int startrow, int numberofrows, int verbose);
RcppExport SEXP _gpuRandom_maternBatchBackend(SEXP varSEXP, SEXP coordsSEXP, SEXP paramSEXP, SEXP NglobalSEXP, SEXP NlocalSEXP, SEXP startrowSEXP, SEXP numberofrowsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type var(varSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nglobal(NglobalSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Nlocal(NlocalSEXP);
    Rcpp::traits::input_parameter< int >::type startrow(startrowSEXP);
    Rcpp::traits::input_parameter< int >::type numberofrows(numberofrowsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    maternBatchBackend(var, coords, param, Nglobal, Nlocal, startrow, numberofrows, verbose);
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
    {"_gpuRandom_cholBatchBackend", (DL_FUNC) &_gpuRandom_cholBatchBackend, 8},
    {"_gpuRandom_crossprodBatchBackend", (DL_FUNC) &_gpuRandom_crossprodBatchBackend, 10},
    {"_gpuRandom_gemmBatchBackend", (DL_FUNC) &_gpuRandom_gemmBatchBackend, 9},
    {"_gpuRandom_gemmBatch2backend", (DL_FUNC) &_gpuRandom_gemmBatch2backend, 11},
    {"_gpuRandom_cpp_gpuFisher_test", (DL_FUNC) &_gpuRandom_cpp_gpuFisher_test, 6},
    {"_gpuRandom_cpp_mrg31k3pCreateStreams", (DL_FUNC) &_gpuRandom_cpp_mrg31k3pCreateStreams, 1},
    {"_gpuRandom_gpuRnBackend", (DL_FUNC) &_gpuRandom_gpuRnBackend, 4},
    {"_gpuRandom_cpp_gpu_qqnorm", (DL_FUNC) &_gpuRandom_cpp_gpu_qqnorm, 6},
    {"_gpuRandom_likfitGpu_Backend", (DL_FUNC) &_gpuRandom_likfitGpu_Backend, 15},
    {"_gpuRandom_logfactsumBackend", (DL_FUNC) &_gpuRandom_logfactsumBackend, 2},
    {"_gpuRandom_rowsumBackend", (DL_FUNC) &_gpuRandom_rowsumBackend, 4},
    {"_gpuRandom_matrix_matrix_sumBackend", (DL_FUNC) &_gpuRandom_matrix_matrix_sumBackend, 4},
    {"_gpuRandom_matrix_vector_sumBackend", (DL_FUNC) &_gpuRandom_matrix_vector_sumBackend, 5},
    {"_gpuRandom_fillParamsExtra", (DL_FUNC) &_gpuRandom_fillParamsExtra, 1},
    {"_gpuRandom_maternBatchBackend", (DL_FUNC) &_gpuRandom_maternBatchBackend, 8},
    {"_gpuRandom_multiplyLowerDiagonalBatchBackend", (DL_FUNC) &_gpuRandom_multiplyLowerDiagonalBatchBackend, 9},
    {"_gpuRandom_multiplyDiagonalBatchBackend", (DL_FUNC) &_gpuRandom_multiplyDiagonalBatchBackend, 6},
    {"_gpuRandom_multiplyLowerBatchBackend", (DL_FUNC) &_gpuRandom_multiplyLowerBatchBackend, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpuRandom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
