
#include "lgmlikFit.hpp"
//#include "gpuRandom.hpp"
//#include<Rmath.h>

using namespace Rcpp;
using namespace viennacl; 
using namespace viennacl::linalg;


//#define M_PI  3.141592653589793238463
/*
 const double PI  =3.141592653589793238463;
 const float  PI_F=3.14159265358979f;*/
/*#define M_LN_SQRT_2PI 0.918938533204672741780329736406  /* log(sqrt(2*pi))
 == log(2*pi)/2 */


template<typename T> 
void likfitGpu_0(viennacl::matrix<T> &Vbatch,
                 viennacl::matrix<T> &coordsGpu, 
                 viennacl::matrix_range<viennacl::matrix<T> > &paramsBatch, 
                 viennacl::matrix<T> &diagMat, 
                 viennacl::vector_base<T> &logD,
                 viennacl::matrix<T> &ab,
                 viennacl::matrix<T> &yX, //y1,y2,y3,X               
                 viennacl::matrix<T> &temp2,
                 viennacl::matrix<T> &aTDa,   // ssqY
                 viennacl::matrix<T> &betas, //a vclmatrix  //given by the user or provided from formula, default=null
                 viennacl::matrix<T> &temp00,
                 viennacl::matrix<T> &temp0,
                 viennacl::matrix<T> &temp1,
                 viennacl::matrix<T> &ssqBeta0,
                 viennacl::matrix<T> &ssqBeta, //  
                 viennacl::matrix<T> &one,          
                 viennacl::matrix<T> &diagP,
                 viennacl::matrix<T> &temp3,
                 viennacl::matrix<T> &nine0,
                 viennacl::matrix<T> &nine,      // ssqX
                 viennacl::matrix<T> &two,
                 viennacl::vector_base<T> &logP,      
                 viennacl::matrix<T> &Qinverse,
                 viennacl::matrix<T> &identity,
                 viennacl::matrix<T> &QPQinverse,
                 viennacl::matrix<T> &betahat,  // p*rowbatch   colbatch                            
                 viennacl::matrix_range<viennacl::matrix<T>> &variances,   // must be rowbatch * colbatch matrix !!!
                 viennacl::vector_base<T> &logD_plusP,
                 viennacl::matrix<T> &form_temp0,
                 viennacl::matrix<T> &form_temp,
                 viennacl::matrix<T> &form_temp1,
                 viennacl::matrix<T> &form_temp2,
                 viennacl::matrix<T> &jacobian,    
                 viennacl::matrix_range<viennacl::matrix<T> > &LogLik,  
                 int n, 
                 int p, 
                 int rowbatch, 
                 int colbatch,   
                 int form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 Rcpp::IntegerVector workgroupSize,
                 Rcpp::IntegerVector localSize,
                 Rcpp::IntegerVector localSizechol,
                 Rcpp::IntegerVector NlocalCache,
                 const int ctx_id){
  
  //if((form == 1 && betas.size1==betas.size2==0) || (form == 3 && betas.size1==betas.size2==0) )  
  
  //  Rcpp::Rcout << "\n" << "need betas" << "\n";
  //get matern matrix
  maternBatchVcl(Vbatch, coordsGpu, paramsBatch, workgroupSize, localSize, ctx_id, 0, rowbatch, 0);
  //#Vbatch=LDL^T, cholesky decomposition
  Rcpp::IntegerVector Astartend =IntegerVector::create(0,n,0,n);
  Rcpp::IntegerVector Dstartend =IntegerVector::create(0,1,0,n);
  cholBatchVcl(Vbatch, diagMat, Astartend, Dstartend, rowbatch, 
               workgroupSize, localSizechol, NlocalCache, ctx_id);
  
  
  //logD <- apply(log(diagMat),1,sum)
  rowsum(diagMat, logD, "row", 1);
  //variances <- vclMatrix(matrix(paramsBatch[,3], nrow=rowbatch, ncol=colbatch, byrow=FALSE), type=type)     # this has to be a vclMatrix not vector
  
  
  //#L(a1,a2,a3, b) = (y1,y2,y3, X)
  Rcpp::IntegerVector Cstartend =IntegerVector::create(0, n,0, colbatch+p);
  //Astartend = IntegerVector::create(0, n,0, n);
  Rcpp::IntegerVector Bstartend =IntegerVector::create(0,n,0,colbatch+p);  
  backsolveBatch(ab, Vbatch, yX, Cstartend, Astartend, Bstartend, 
                 1, 1, workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  //# temp2 = (ab)^T * D^(-1) *ab
  Cstartend=IntegerVector::create(0,colbatch+p,0,colbatch+p);
  Astartend=IntegerVector::create(0,n,0,colbatch+p);
  //Dstartend=IntegerVector::create(0,1,0,n);
  crossprodBatch(temp2, ab, diagMat, 1, Cstartend, Astartend, Dstartend,  
                 workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  /*extract a^TDa from temp2       1x3
   for (j in 1:colbatch){
   for (i in 1:rowbatch){
   aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
   }
   }*/
  
  //extract a^TDa from temp2       1x3
  viennacl::slice s2(0, 1, colbatch);   // {0,1,2}  colbatch=3    ？？？？？？？？？？？
  for (int i = 0; i < rowbatch; i++){// square sub_matrices
    viennacl::slice s1((colbatch+p)*i, 1, colbatch); 
    matrix_slice<viennacl::matrix<T> >  temp2_sub(temp2, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
    viennacl::vector<T> diags = viennacl::row(aTDa, i);
    diags = viennacl::diag(temp2_sub);
  }
  
  /*
   // Extract 5-th row of A, then overwrite with 6-th diagonal:
   viennacl::vector<T> r = viennacl::row(A, 4);
   r = viennacl::row(A, 5);
   viennacl::vector<T> temp2_diag = viennacl::diag(temp2);
   viennacl::range r(0, colbatch); 
   vector_range<vector<T> > aTDa(temp2_diag, r); */
  
  Rcpp::IntegerVector transposeABC = IntegerVector::create(0,0,0);
  Rcpp::IntegerVector batches = IntegerVector::create(rowbatch,1,0,0,1,0); //batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
  Rcpp::IntegerVector workgroupSize_gemm = IntegerVector::create(workgroupSize[0],workgroupSize[1],workgroupSize[2],localSize[0],localSize[1],localSize[2]);
  Rcpp::IntegerVector NlocalCache_gemm = IntegerVector::create(NlocalCache[0],NlocalCache[0]);
  
  
#ifdef UNDEF
  
  if(form == 1 || form == 3){       
    // a^TD^(-1)b * beta = temp00
    Astartend=IntegerVector::create(0,colbatch,(colbatch+p),colbatch,p,(colbatch+p));
    Bstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
    Cstartend=IntegerVector::create(0,colbatch,colbatch,0,colbatch,colbatch);
    gemmBatch2(temp2, betas, temp00, transposeABC,  
               Astartend, Bstartend, Cstartend, batches, 
               workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
    
    
    // extract a^TD^(-1)b * beta from temp00      1x3
    for (int i = 0; i < rowbatch; i++){
      viennacl::slice s1(colbatch*i, 1, colbatch); 
      matrix_slice<viennacl::matrix<T> >  temp00_sub(temp00, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
      viennacl::vector<T> diag2 = viennacl::row(temp0, i);
      diag2 = viennacl::diag(temp00_sub);
    }
    
    
    
    // b * beta = temp1 (n x colbatch), to get ssqBeta
    Astartend=IntegerVector::create(0,n,n,colbatch,p,(colbatch+p));
    //Bstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
    Cstartend=IntegerVector::create(0,n,n,0,colbatch,colbatch);
    //Rcpp::IntegerVector batches = IntegerVector::create(rowbatch,1,0,0,1,0); 
    gemmBatch2( ab, betas, temp1, transposeABC,  
                Astartend, Bstartend, Cstartend, batches, 
                workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
    
    
    
    
    // ssqBeta0 = temp1^T D^(-1) temp1 = beta T*bT D^(-1) b*beta    3x3
    Cstartend=IntegerVector::create(0, colbatch, 0, colbatch);
    Astartend=IntegerVector::create(0, n, 0, colbatch);
    //Dstartend=IntegerVector::create(0, 1, 0, n);
    crossprodBatch( ssqBeta0, temp1, diagMat,  TRUE,  
                    Cstartend, Astartend, Dstartend,  
                    workgroupSize, localSize, NlocalCache[0], ctx_id);
    
    
    
    // extract beta T*bT D^(-1) b*beta from ssqBeta0       1x3
    for (int i = 0; i < rowbatch; i++){
      viennacl::slice s1(colbatch*i, 1, colbatch); 
      matrix_slice<viennacl::matrix<T> >  ssqBeta0_sub(ssqBeta0, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
      viennacl::vector<T> diag3 = viennacl::row(ssqBeta, i);
      diag3 = viennacl::diag(ssqBeta0_sub);
    }
    
    // to get one, see the notes 
    one = aTDa - 2*temp0 + ssqBeta;      // 1x3
  }
  
#endif
  
  // to get hat_sigma^2 when beta is not given by users
  // b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
  Astartend =IntegerVector::create(colbatch, p, colbatch, p);
  Dstartend =IntegerVector::create(0, 1, 0, p);
  cholBatchVcl( temp2, diagP, Astartend, Dstartend,  
                rowbatch, workgroupSize, localSizechol, NlocalCache, ctx_id);
  
  
  
  // temp3 = Q^(-1) * (b^T * D^(-1) *a), backsolve for temp3    2 by 1
  Cstartend =IntegerVector::create(0,p,0,colbatch);
  //Astartend = IntegerVector::create(colbatch, p, colbatch, p);
  Bstartend =IntegerVector::create(colbatch, p, 0, colbatch);  
  backsolveBatch(temp3, temp2, temp2, 
                 Cstartend, Astartend, Bstartend, 
                 rowbatch, 1, workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  // nine0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices      3x3
  Cstartend = IntegerVector::create(0, colbatch, 0, colbatch);
  Astartend = IntegerVector::create(0, p, 0, colbatch);
  //Dstartend =IntegerVector::create(0, 1, 0, p);  
  crossprodBatch(nine0, temp3, diagP, 1,  
                 Cstartend, Astartend, Dstartend,  
                 workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  // extract needed cells from nine0       1x3
  for (int i = 0; i < rowbatch; i++){
    viennacl::slice s1(colbatch*i, 1, colbatch); 
    matrix_slice<viennacl::matrix<T> >  nine0_sub(nine0, s1, s2);  //project(temp2, r, r); //returns a matrix_range as above
    viennacl::vector<T> diag4 = viennacl::row(nine, i);
    diag4 = viennacl::diag(nine0_sub);
  }
  
  two = aTDa - nine;
  
  
  rowsum( diagP, logP, "row",1L);  
  
  
  // calculate betahat = Q^(-T) p^(-1) Q^(-1) * b^T D^(-1)a
  // first: get Q^(-1)    pxp
  Cstartend = IntegerVector::create(0, p, 0, p);
  Astartend = IntegerVector::create(colbatch, p, colbatch, p);
  Bstartend = IntegerVector::create(0, p, 0, p);    
  backsolveBatch(Qinverse, temp2, identity,
                 Cstartend, Astartend, Bstartend, 
                 1,  1,
                 workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  // Q^(-T) P^(-1) Q^(-1) = QPQinverse
  //Cstartend = IntegerVector::create(0, p, 0, p);
  Astartend = IntegerVector::create(0, p, 0, p);
  //Dstartend = IntegerVector::create(0, 1, 0, p);    
  crossprodBatch(QPQinverse, Qinverse, diagP, 1,
                 Cstartend, Astartend, Dstartend,  
                 workgroupSize, localSize, NlocalCache[0], ctx_id);
  
  
  
  // QPQinverse * b^T D^(-1)a = betahat
  //Rcpp::IntegerVector transposeABC = IntegerVector::create(0,0,0);
  Astartend=IntegerVector::create(0,p,p,0,p,p);
  Bstartend=IntegerVector::create(colbatch, p,(colbatch+p), 0, colbatch, (colbatch+p));
  Cstartend=IntegerVector::create(0,p,p,0,colbatch,colbatch);
  batches = IntegerVector::create(rowbatch,1,0,0,0,0); //batches: nRow, nCol, recycleArow, recycleAcol, recycleBrow, col
  gemmBatch2(QPQinverse, temp2, betahat,
             transposeABC,  Astartend, Bstartend, Cstartend,  
             batches, workgroupSize_gemm, NlocalCache_gemm, 0, ctx_id);
  
  
  
  
  
  
  if(form == 1 || form == 4){     
    // get form_temp <- n*log(variances) + logD    , form_temp = n*log(sigma^2)+log |D|
    viennacl::matrix<T> nlog_variances(rowbatch, colbatch);  
    nlog_variances = n * element_log(variances);
    matrix_vector_sum( nlog_variances, logD, form_temp, 1, workgroupSize, ctx_id);  // form_temp
    
  } 
  
  if (form ==5 || form ==6){
    logD_plusP = logD + logP;   
  }
  
  
  T nlog2pi  = (T) n*2*M_LN_SQRT_2PI;   //  M_LN_SQRT_2PI == log(2*pi)/2  
  
  
#ifdef UNDEF
  
  
  if (form==1){ //loglik
    //form_temp + one/variances + n*log(2*PI)+jacobian 
    //form_temp0 = form_temp + element_div(one, variances);
    //matrix_vector_sum(form_temp0, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp0 = form_temp + element_div(one, variances) + jacobian;
    
    matrix_scalar_sum(form_temp0, nlog2pi, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form ==2){ //ml  =n*log(two/n) +logD +jacobian + n*log(2*PI) + n 
    form_temp0 = n * element_log(two/n);
    
    matrix_vector_sum(form_temp0, logD, form_temp, 1L, workgroupSize, ctx_id);   
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    
    form_temp2 = form_temp + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi+n, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==3){ // mlFixSigma
    //n*log(one/n)+logD + n*log(2*PI) + n+jacobian
    form_temp0 = n * element_log(one/n);
    
    matrix_vector_sum(form_temp0, logD, form_temp, 1L, workgroupSize, ctx_id);
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 = form_temp + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi+ n, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==4){ //mlFixBeta
    //form_temp + two/variances + n*log(2*PI)+jacobian 
    //form_temp0 = form_temp + element_div(two, variances);
    //matrix_vector_sum(form_temp0, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp0 = form_temp + element_div(two, variances) + jacobian;
    
    matrix_scalar_sum(form_temp0, nlog2pi, form_temp1, workgroupSize,ctx_id);
    
    LogLik = -0.5*form_temp1;
    
  }else if(form==5){ //reml
    //form_temp = (n-p)*log(variances) + logD + logP
    //minusTwoLogLik=form_temp + two/variances + n*log(2*PI)+jacobian
    form_temp0 = (n-p) * element_log(variances);
    matrix_vector_sum(form_temp0, logD_plusP, form_temp, 1L, workgroupSize, ctx_id);    
    //viennacl::matrix<T> plus_two_div_variances = form_temp + element_div(two, variances);
    //matrix_vector_sum(plus_two_div_variances, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 = form_temp + element_div(two, variances) + jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi, form_temp1, workgroupSize, ctx_id);
    
    
    LogLik = -0.5 * form_temp1;
    
  }else if(form==6){ // remlPro
    // minusTwoLogLik= (n-p)*log(two/(n-p)) + logD + logP + jacobian + n*log(2*PI) + n-p , 
    form_temp0 = (n-p) * element_log(two/(n-p));
    matrix_vector_sum(form_temp0, logD_plusP, form_temp, 1L, workgroupSize, ctx_id);   
    //matrix_vector_sum(form_temp, jacobian, form_temp2, 0L, workgroupSize, ctx_id);
    form_temp2 =form_temp+ jacobian;
    
    matrix_scalar_sum(form_temp2, nlog2pi + n-p, form_temp1, workgroupSize, ctx_id);
    
    LogLik = -0.5*form_temp1;
  } 
#endif 
  
}










template<typename T> 
void likfitGpu_1(viennacl::matrix<T> &coordsGpu, 
                 viennacl::matrix<T> &bigparamsBatch, 
                 viennacl::matrix<T> &yX, //y1,y2,y0,X,   already transformed            
                 viennacl::matrix<T> &betas, //a vclmatrix  //given by the user or provided from formula, default=null                          
                 viennacl::matrix<T> &bigvariances,   // must be rowbatch * colbatch matrix !!!
                 viennacl::matrix<T> &jacobian,    // (groupsize, colbatch)
                 viennacl::matrix<T> &ssqBeta,   // (groupsize, colbatch)
                 viennacl::matrix<T> &nine,     // ssqX     (groupsize, colbatch)
                 viennacl::matrix<T> &aTDa,     // ssqY     (groupsize, colbatch)
                 viennacl::vector_base<T> &logD,    //of (groupsize)
                 viennacl::vector_base<T> &logP,     //(groupsize)    
                 viennacl::matrix<T> &betahat,   //(groupsize*p, colbatch) 
                 viennacl::matrix<T> &finalLogLik,
                 const int n, 
                 const int p, 
                 const int groupsize,
                 const int colbatch,
                 const int form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 Rcpp::IntegerVector workgroupSize,
                 Rcpp::IntegerVector localSize,
                 Rcpp::IntegerVector NlocalCache,
                 const int ctx_id){
  
  
  viennacl::matrix<T> Vbatch(groupsize*n, n);
  viennacl::matrix<T> diagMat(groupsize, n);//viennacl::vector_base<T> logD(groupsize);
  viennacl::matrix<T> ab(groupsize*n, colbatch+p);             
  viennacl::matrix<T> temp2((colbatch+p)*groupsize, (colbatch+p));//viennacl::matrix<T> aTDa(groupsize, colbatch); // ssqY
  viennacl::matrix<T> temp00(groupsize*colbatch, colbatch);
  viennacl::matrix<T> temp0(groupsize, colbatch);
  viennacl::matrix<T> temp1(groupsize*n, colbatch);  // b* beta
  viennacl::matrix<T> ssqBeta0(groupsize*colbatch, colbatch); //viennacl::matrix<T> ssqBeta(groupsize, colbatch);
  viennacl::matrix<T> one(groupsize, colbatch);         
  viennacl::matrix<T> diagP(groupsize, p);
  viennacl::matrix<T> temp3(groupsize*p, colbatch);
  viennacl::matrix<T> nine0(colbatch*groupsize, colbatch); //viennacl::matrix<T> nine(groupsize, colbatch); // ssqX
  viennacl::matrix<T> two(groupsize, colbatch); //viennacl::vector_base<T> logP(groupsize);     
  viennacl::matrix<T> Qinverse(groupsize*p, p);
  viennacl::matrix<T> identity = viennacl::identity_matrix<T>(p);
  viennacl::matrix<T> QPQinverse(groupsize*p, p); //viennacl::matrix<T> betahat(groupsize*p, colbatch);  // p*groupsize   colbatch  
  viennacl::vector_base<T> logD_plusP(groupsize);     
  viennacl::matrix<T> form_temp0(groupsize, colbatch);                          
  viennacl::matrix<T> form_temp(groupsize, colbatch); 
  viennacl::matrix<T> form_temp1(groupsize, colbatch); 
  viennacl::matrix<T> form_temp2(groupsize, colbatch);
  
  Rcpp::IntegerVector localSizechol = localSize;
  localSizechol[1]= workgroupSize[1];
  
  
  const int totalsets = bigparamsBatch.size1();
  
  viennacl::range r(0, bigparamsBatch.size2()); 
  viennacl::range r2(0, colbatch);  
  
  
  
  ///////////////////////////Loop starts !!!//////////////////////////////////////////////////////////////////////////
  for (int i=0; i< totalsets/groupsize; i++ ){
    
    viennacl::range ri(i*groupsize, (i+1)*groupsize);  //{0, 2}    
    viennacl::matrix_range<viennacl::matrix<T>> paramsBatch(bigparamsBatch, ri, r);//{m(0,0),m(0,1);m(1,0),m(1,1)} matrix_slice<matrix<double> > ms(m, s, s);         
    viennacl::matrix_range<viennacl::matrix<T>> variances(bigvariances, ri, r2);   
    viennacl::matrix_range<viennacl::matrix<T>> loglik(finalLogLik, ri, r2);
    
    //42
    likfitGpu_0(Vbatch,
                coordsGpu, 
                paramsBatch, 
                diagMat, logD,
                ab,  yX, //y1,y2,y3,X               
                temp2,
                aTDa,   // ssqY
                betas, //a vclmatrix  //given by the user or provided from formula, default=null
                temp00, temp0, temp1,
                ssqBeta0, ssqBeta, //  
                one,          
                diagP,
                temp3, nine0,
                nine,      // ssqX
                two,
                logP,      
                Qinverse, identity, QPQinverse,
                betahat,  // p*rowbatch   colbatch                            
                variances,   // must be rowbatch * colbatch matrix !!!
                logD_plusP,
                form_temp0, form_temp, form_temp1, form_temp2,
                jacobian,
                loglik,
                n,  p, groupsize, colbatch,   form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                workgroupSize, localSize, localSizechol, NlocalCache, ctx_id);
    
    
  }
}      







//#ifdef UNDEF


//15
template<typename T> 
void likfitGpu_Templated(
    Rcpp::S4 coordsGpuR,
    Rcpp::S4 bigparamsBatchR,
    Rcpp::S4 yXR, //y1,y2,y0,X,   already transformed     
    Rcpp::S4 betasR,
    Rcpp::S4 bigvariancesR,  //already in the correct form
    Rcpp::S4 jacobianR, // already calculated, a vector 
    Rcpp::S4 ssqBetaR,
    Rcpp::S4 ssqXR,
    Rcpp::S4 ssqYR,
    Rcpp::S4 logDR,
    Rcpp::S4 logPR,
    Rcpp::S4 betahatR,
    Rcpp::S4 finalLogLikR,   
    const int n, 
    const int p, 
    const int groupsize,
    const int colbatch,
    const int form,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache) {
  
  const bool BisVCL=1;
  const int ctx_id = INTEGER(finalLogLikR.slot(".context_index"))[0]-1;
  std::shared_ptr<viennacl::matrix<T> > coordsGpu = getVCLptr<T>(coordsGpuR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > bigparamsBatch = getVCLptr<T>(bigparamsBatchR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > yX = getVCLptr<T>(yXR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > betas = getVCLptr<T>(betasR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > bigvariances = getVCLptr<T>(bigvariancesR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > jacobian = getVCLptr<T>(jacobianR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > ssqBeta = getVCLptr<T>(ssqBetaR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > ssqX = getVCLptr<T>(ssqXR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > ssqY = getVCLptr<T>(ssqYR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > logD = getVCLVecptr<T>(logDR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<T> > logP = getVCLVecptr<T>(logPR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > betahat = getVCLptr<T>(betahatR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > finalLogLik = getVCLptr<T>(finalLogLikR.slot("address"), BisVCL, ctx_id);
  
  
  likfitGpu_1<T>(*coordsGpu, 
                 *bigparamsBatch, 
                 *yX, //y1,y2,y0,X,   already transformed            
                 *betas, //a vclmatrix  //given by the user or provided from formula, default=null                          
                 *bigvariances,   // must be rowbatch * colbatch matrix !!!
                 *jacobian,
                 *ssqBeta,
                 *ssqX,
                 *ssqY,
                 *logD,
                 *logP,
                 *betahat,
                 *finalLogLik,
                 n, p, groupsize, colbatch, form, // c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)
                 workgroupSize, localSize, 
                 NlocalCache, ctx_id);
  
}










//[[Rcpp::export]]
void likfitGpu_Backend(
    Rcpp::S4 coordsGpuR,
    Rcpp::S4 bigparamsBatchR, // vclMatrix
    Rcpp::S4 yXR, //y1,y2,y0,X,   already transformed     
    Rcpp::S4 betasR,
    Rcpp::S4 bigvariancesR,  //already in the correct form
    Rcpp::S4 jacobianR, // already calculated, a vector 
    Rcpp::S4 ssqBetaR,
    Rcpp::S4 ssqXR,
    Rcpp::S4 ssqYR,
    Rcpp::S4 logDR,
    Rcpp::S4 logPR,
    Rcpp::S4 betahatR,
    Rcpp::S4 finalLogLikR,   
    const int n, 
    const int p, 
    const int groupsize,
    const int colbatch,
    const int form,
    Rcpp::IntegerVector workgroupSize,
    Rcpp::IntegerVector localSize,
    Rcpp::IntegerVector NlocalCache) {
  
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(finalLogLikR));
  std::string precision_type = (std::string) classVarR;
  
  
  
  if(precision_type == "fvclMatrix") {
    likfitGpu_Templated<float>(coordsGpuR, bigparamsBatchR, yXR,  betasR,  bigvariancesR,  jacobianR,    ssqBetaR, ssqXR,  ssqYR, logDR, logPR, betahatR, 
                               finalLogLikR, n, p, groupsize, colbatch, form, workgroupSize, localSize, 
                               NlocalCache);
  } 
  else if (precision_type == "dvclMatrix") {
    likfitGpu_Templated<double>(coordsGpuR, bigparamsBatchR, yXR,  
                                betasR,  bigvariancesR,
                                jacobianR,ssqBetaR,ssqXR,  ssqYR, 
                                logDR, logPR, betahatR, finalLogLikR, 
                                n, p, groupsize, colbatch, form, 
                                workgroupSize, 
                                localSize, 
                                NlocalCache);
  } else {
    Rcpp::warning("class of var must be fvclMatrix or dvclMatrix");
  }
}


//#endif







