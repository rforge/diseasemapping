#include "gpuRandom.hpp"
//#define DEBUG
//#define NOKERNELS

// C = A^T A or A^T D A or A^T D^(-1) A 
// if onlyDiagC = 1, then only compute diagonals of C
// C is an Nmatrix by Ncol matrix


template <typename T> 
std::string crossprodBatchString(
    const int Nrow, 
    const int Ncol,
//    const int Nmatrix, 
    const int NpadC, // ignored if onlyDiagC, use NpadBetweenMatrices in this case
    const int NpadA,
    const int NpadD, // set to zero to omit D
    const int invertD, // set to 1 for A^T D^(-1) A
    const int onlyDiagC, // set to 1 to only compute diagonals of C
    const int NstartC,  // newly added
    const int NstartA,  // new
    const int NstartD,  // new
    const int NpadBetweenMatricesC,
    const int NpadBetweenMatricesA,
    const int NlocalCacheA, // numbers of rows to cache of A
    Rcpp::IntegerVector Nlocal// cache a Nlocal[0] by Nlocal[1] submatrix of C
  ) { 
    
 /*
  * global groups dimension 1 is column, groups dimension 2 is matrix
  * local items dimension 1 is inner, dimension 2 is row
  * 
  */
 
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  const int NrowStop = std::min(NlocalCacheA, Nrow);
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result += 
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"  
//    "#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
    "#define NpadC " + std::to_string(NpadC) + "\n"    
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadD " + std::to_string(NpadD) + "\n"   
    "#define NstartC " + std::to_string(NstartC) + "\n"   
    "#define NstartA " + std::to_string(NstartA) + "\n" 
    "#define NstartD " + std::to_string(NstartD) + "\n" 
    "#define NpadLocal " + std::to_string(Nlocal[1]) + "\n"    
    "#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
    "#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
    "#define NrowStop " + std::to_string(NrowStop) + "\n"    
    "#define NlocalCacheA "  + std::to_string(NlocalCacheA) + "\n\n";   
    
    result +=
    "__kernel void crossprodBatch(\n"
    "global " + typeString+ " *C,\n"
    "const global "+ typeString+ " *A,\n";
    
    if (NpadD >0 ) {
      result +=  "const global " + typeString + " *D,\n";
    }
  
    result += "int Nmatrix) {\n\n"
   
    "local " + typeString + " Acache[" + 
      std::to_string(NlocalCacheA) + "];\n" 
    "local " + typeString + " Dcache[" + 
      std::to_string(NlocalCacheA) + "];\n" 
    "local " + typeString + " Ccache[" + 
      std::to_string(Nlocal[0] * Nlocal[1]) + "];\n" +

  typeString + " Cout, Ctemp;\n"
  "event_t ev;\n"
  "int AHere, CHere;\n"
  "int Dmatrix, Drow, Dcol, DrowNpadC, Dinner, DinnerA, DinnerAcol;\n"
  "int A0Dcol, A0Drow;// location of elements A[0,Dcol] and A[0,Drow]\n" 
  "const int AHereInc = get_num_groups(1)*NpadBetweenMatricesA;\n"
  "const int CHereInc = get_num_groups(1)*NpadBetweenMatricesC;\n"
  "const int DrowNpadCInc = get_local_size(1)*NpadC;\n"
 
  "const int localIndex = get_local_id(0) * get_local_size(1) + get_local_id(1);\n"
  "const int NlocalTotal = get_local_size(1)*get_local_size(0);\n"
  "const int cacheIndex = get_local_id(1)+NpadLocal*get_local_id(0);\n";

    if(onlyDiagC) {
      result +=    "const int doLocalSum = (localIndex==0);\n"
      "const int DinnerAinc = NlocalTotal*NpadA;\n";
    } else {
      result +=    "const int doLocalSum = (get_local_id(0)==0);\n"
                  "const int DinnerAinc = get_local_size(0)*NpadA;\n";
    }
    
    
  if(NpadD) {
    result += "int DHere;\n"
    "const int DHereInc = get_num_groups(1)*NpadD;\n";
  }
  
      result +=  "\n\n"
  "for(Dmatrix = get_group_id(1),\n"
  "    AHere = Dmatrix * NpadBetweenMatricesA + NstartA,\n"; // made changes here
  if(NpadD) {
    result +=
    "    DHere = Dmatrix * NpadD + NstartD,\n"; // made changes here
  }
  result +=  
    "    CHere = Dmatrix * NpadBetweenMatricesC + NstartC;\n" // made changes here
  "  Dmatrix < Nmatrix;\n"
  "  Dmatrix += get_num_groups(1),\n"
  "    AHere += AHereInc,\n";
  if(NpadD) {
    result += "    DHere += DHereInc,\n";
  }
  result +=  
  "    CHere += CHereInc\n"
  "  ){\n";

  
  if(NpadD && NrowStop) {  
    result +=  "\n"
    "  /* cache first NrowStop elements of D for this matrix*/\n"
    "  ev=async_work_group_copy(\n"
    "    Dcache, &D[DHere],\n"
    "    NrowStop, 0);\n"
    "  wait_group_events (1, &ev);\n"
    "  barrier(CLK_LOCAL_MEM_FENCE);\n";
    
  }
  
  result +=  "\n"
  "  for(Dcol = get_group_id(0),\n"
  "      A0Dcol = AHere + Dcol;\n"
  "    Dcol < Ncol;\n"
  "    Dcol += get_num_groups(0),\n"
  "    A0Dcol += get_num_groups(0)){\n\n";
  
  
if(NrowStop) {
  result +=  "\n"
  "    // cache A[1:NrowStop, Dcol]\n"
  "    ev=async_work_group_strided_copy (\n"
  "      Acache, &A[A0Dcol],\n"
  "      NrowStop, NpadA, 0);\n"
  "    wait_group_events (1, &ev);\n";
  if(onlyDiagC) {
    result +=
    "    // square the cached A\n"
    "    for(Dinner = localIndex;\n"
    "        Dinner < NrowStop;\n"
    "        Dinner+=NlocalTotal\n"
    "    ){\n";   
  result +=
    "      Acache[Dinner] *= Acache[Dinner];\n";
  result +=
  "    }\n\n";
  }
if(NpadD) {
    result +=
    "    // multiply by D\n"
  "    for(Dinner = localIndex;\n"
  "        Dinner < NrowStop;\n"
  "        Dinner+=NlocalTotal\n"
  "    ){\n";   
  if(invertD) {
  result +=
  "      Acache[Dinner] /= Dcache[Dinner];\n";
  } else {
    result +=
      "      Acache[Dinner] *= Dcache[Dinner];\n";
  }
  result +=
    "    }\n\n";
} // NpadD
} // NrowStop



if(onlyDiagC) {
    result += 
      "     Drow = Dcol;\n"
      "     A0Drow = AHere + Drow;\n"
      "     DrowNpadC = CHere + Drow * NpadC;\n\n";
      
}  else {
    result += 
      "   for(Drow = Dcol + get_local_id(1),\n"
      "       DrowNpadC = CHere + Drow * NpadC,\n"
      "       A0Drow = AHere + Drow;\n"
      "     Drow < Ncol;\n"
      "     Drow += get_local_size(1),\n" 
      "       DrowNpadC += DrowNpadCInc,\n"
      "       A0Drow +=  get_local_size(1)\n"
      "   ){\n\n";
}
  
  result += 
    "      Cout=0.0;\n";
  
  result += 
      "      // cached parts\n";
  
if(onlyDiagC) {
  result +=      "      for(Dinner = localIndex,\n";
} else {
  result +=  "      for(Dinner = get_local_id(0),\n";
}
result +=
      "          DinnerA = A0Drow + Dinner*NpadA;\n"
      "        Dinner < NrowStop;\n"
      "        Dinner += NlocalTotal,\n"
      "          DinnerA += DinnerAinc\n"
      "      ){\n";

if(onlyDiagC) {    
    result += 
      "          Cout += Acache[Dinner];\n"; 
} else {
  result += 
    "          Cout += A[DinnerA] * Acache[Dinner];\n";
    //    "          Cout += A[Dmatrix * NpadBetweenMatricesA + Dinner*NpadA + Drow] * Acache[Dinner];\n";
  }

result += "      }// Dinner\n";

  result +=
    "      // un-cached parts\n";

  if(onlyDiagC) {
  result +=
    "      for(Dinner = NrowStop + localIndex,\n";
} else {
  result +=
  "      for(Dinner = NrowStop + get_local_id(0),\n";
}

result +=
  "          DinnerA = A0Drow + Dinner*NpadA,\n"
    "          DinnerAcol = A0Dcol + Dinner*NpadA;\n"
    "        Dinner < Nrow;\n";

if(onlyDiagC){
} else {
  result +=  "        Dinner += get_local_size(0),\n";
}

result +=    "          DinnerA += DinnerAinc,\n"
    "          DinnerAcol += DinnerAinc\n"
    "      ){\n";

  if(NpadD) {
    if(invertD) {
      result += 
        " Cout += A[DinnerA] * A[DinnerAcol] / D[DHere+Dinner];\n";
    } else {
      result += 
//        "      Cout += A[Dmatrix * NpadBetweenMatricesA + Dinner*NpadA + Dcol] * A[Dmatrix * NpadBetweenMatricesA + Dinner*NpadA + Drow] * D[DHere+Dinner];\n"
//  "      Cout += A[AHere + Dinner*NpadA + Dcol] * A[AHere + Dinner*NpadA + Drow] * D[DHere+Dinner];\n"
//  "      Cout += A[A0Drow + Dinner*NpadA] * A[A0Dcol + Dinner*NpadA] * D[DHere+Dinner];\n"
  "      Cout += A[DinnerA] * A[DinnerAcol] * D[DHere+Dinner];\n";
    }
  } else {
    result += 
      "      Cout += A[DinnerA] * A[DinnerAcol];\n";
  }  
    
    result += 
      "      }// Dinner\n";
    result +=       
      "      Ccache[cacheIndex] = Cout;\n"
      "      barrier(CLK_LOCAL_MEM_FENCE);\n";
    
      result +=
        "      if(doLocalSum){\n";
      if(onlyDiagC) {
        result +=
          "        for(Dinner = 1;Dinner < NlocalTotal;Dinner++){\n"
          "          Ccache[cacheIndex] += Ccache[cacheIndex + Dinner];\n"
          "        }\n";
        result += "          C[CHere + Dcol] = Ccache[cacheIndex];\n";
        // DEBUG
//        result += "          C[CHere + Dcol] = CHere + Dcol;\n";
      } else {
      result +=
      "        for(Dinner = 1;Dinner < get_local_size(0);Dinner++){\n"
      "          Ccache[cacheIndex] += Ccache[cacheIndex + Dinner * NpadLocal];\n"
      "        }\n";
        result += "          C[DrowNpadC + Dcol] = Ccache[cacheIndex];\n";
      } 
      
      result +=
        "      }//doLocalSum \n"
"      barrier(CLK_LOCAL_MEM_FENCE);\n";

    
    if(!onlyDiagC) {
      result += 
      "    }// Drow\n";
    }
    result += 
      "  }// Dcol\n";
    result += 
      "}// Dmatrix\n";
    result += 
      "}// function";
  
  return(result);
}

/*
std::vector<int> Nglobal,
std::vector<int> Nlocal,*/





template <typename T> 
void crossprodBatch(
    viennacl::matrix<T> &C,  // must be a batch of square matrices 
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &D,
    const int invertD,
    Rcpp::IntegerVector Cstartend,
    Rcpp::IntegerVector Astartend,
    Rcpp::IntegerVector Dstartend,  
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal,
    const int NlocalCache, 
    const int ctx_id) {
  
  const int Ncol = Astartend[3];
  const int Nmatrix = C.size1()/C.size2();
  const int Nrow = Astartend[1];
  
  const int NstartC = C.internal_size2() * Cstartend[0] + Cstartend[2];
  const int NstartA = A.internal_size2() * Astartend[0] + Astartend[2];
  const int NstartD = D.internal_size2() * Dstartend[0] + Dstartend[2];
  
  
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
//  cl_device_type type_check = ctx.current_device().type();
  
 /* const int NpadC, 
  const int NpadA,
  const int NpadD, // set to zero to omit D
  const int invertD, // set to 1 for A^T D^(-1) A
  const int NpadBetweenMatricesC,
  const int NpadBetweenMatricesA,
  */
  std::string clString =  crossprodBatchString<T>(  
    Nrow, 
    Ncol, // ncol
//    Nmatrix,
    C.internal_size2(), 
    A.internal_size2(), 
    D.internal_size2(),
    invertD, // A^T D^(-1) A
    0, // don't only compute diagonals of C
    NstartC,
    NstartA,
    NstartD,
    C.internal_size2()*C.size2(),//NpadBetweenMatricesC,
    A.internal_size2()*A.size1()/Nmatrix,//NpadBetweenMatricesA,
    NlocalCache,
    Nlocal);
  
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif  
  
  
  
  viennacl::ocl::program & my_prog = ctx.add_program(clString, "my_kernel");
#ifndef NOKERNELS  
  viennacl::ocl::kernel & multiplyKernel = my_prog.get_kernel("crossprodBatch");
  
  multiplyKernel.global_work_size(0, Nglobal[0]);
  multiplyKernel.global_work_size(1, Nglobal[1]);

  multiplyKernel.local_work_size(0, Nlocal[0]);
  multiplyKernel.local_work_size(1, Nlocal[1]);
  
  // diagonals and diagTimesRowOfA
  viennacl::ocl::enqueue(multiplyKernel( C, A, D, Nmatrix));
#endif  
}



template <typename T> 
SEXP crossprodBatchTyped(
    Rcpp::S4 CR,
    Rcpp::S4 AR,
    Rcpp::S4 DR,
    const int invertD,
    Rcpp::IntegerVector Cstartend,
    Rcpp::IntegerVector Astartend,
    Rcpp::IntegerVector Dstartend,  
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal, 
    const int NlocalCache) {
  /*
  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);*/
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  
  
  std::shared_ptr<viennacl::matrix<T> > 
    AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    DG = getVCLptr<T>(DR.slot("address"), BisVCL, ctx_id);
  
  crossprodBatch<T>(*CG, *AG, *DG, invertD, Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache, ctx_id);
  
  return Rcpp::wrap(0L);
  
}


//' Multiply crossproduct matrices
//' 
//' Computes C = t(A) D A
//'
//' @param C output matrices, stacked row-wise
//' @param A rectangular matrices
//' @param D rectangular matrix, columns are diagonals
//' @param invertD set to 1 for C = t(A) D^(-1) A
//' @param Nglobal vector of number of global work items
//' @param Nlocal vector of number of local work items
//' @param NlocalCache elements in local cache
//' @export
// [[Rcpp::export]]
SEXP crossprodBatchBackend(
    Rcpp::S4 C,
    Rcpp::S4 A,
    Rcpp::S4 D,
    const int invertD,
    Rcpp::IntegerVector Cstartend,
    Rcpp::IntegerVector Astartend,
    Rcpp::IntegerVector Dstartend, 
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal, 
    const int NlocalCache) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  
  
  if(precision_type == "fvclMatrix") {
    result = crossprodBatchTyped<float>(C, A, D, invertD,  Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = crossprodBatchTyped<double>(C, A, D, invertD, Cstartend, Astartend, Dstartend, Nglobal, Nlocal, NlocalCache);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
  
}

