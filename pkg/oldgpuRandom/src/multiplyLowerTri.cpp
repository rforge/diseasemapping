#include "gpuRandom.hpp"
//#define DEBUG

// C = A B, A diagonal

template <typename T> 
std::string multiplyDiagonalBatchString(
    const int sameB,
    const int Nrow, 
    const int Ncol,
    const int Nmatrix, 
    const int NpadC, 
    const int NpadA, 
    const int NpadB,
    const int NpadBetweenMatricesC,
    const int NpadBetweenMatricesB,
    const int NlocalCacheA) { // get_local_size(0)
  
// Dmatrix, Drow, Dcol
// C[Drow,Dcol,Dmatrix] = A[Dmatrix,Drow] * B[Drow,Dcol,Dmatrix]
  // work items are Drow, Dcol, Dmatrix
  // local size is x, y, 1
  // work group caches A[Dmatrix, Drow[local0seq]]
  
  std::string typeString = openclTypeString<T>();
  std::string result = "";
  
  if(typeString == "double") {
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"    
    "#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
    "#define NpadC " + std::to_string(NpadC) + "\n"    
    "#define NpadA " + std::to_string(NpadA) + "\n"    
    "#define NpadB " + std::to_string(NpadB) + "\n"    
    "#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
    "#define NpadBetweenMatricesB " + std::to_string(NpadBetweenMatricesB) + "\n"    
    "#define NlocalCacheA "  + std::to_string(NlocalCacheA) + "\n\n"    
    
    "__kernel void multiplyDiagonalBatch(\n"
    "	__global " + typeString+ " *C,\n"
    "	__global "+ typeString+ " *A,\n"
    "	__global "+ typeString+ " *B) {\n\n" +
      
    "int DrowGlobal, DrowLocal, Dcol,Dmatrix, Diter;\n"
    "	local "+ typeString+ " Acache[NlocalCacheA];\n" 
    " int AHere, BHere, CHere;\n" 
    " "+ typeString+ " AforThisWorkitem;\n" 
    "const int doCacheA = (get_local_id(1) == 0);\n"

  "for(Dmatrix = get_global_id(2); Dmatrix < Nmatrix; Dmatrix += get_global_size(2)){\n"

  "  AHere = Dmatrix*NpadA;\n"

  
  
  "for(DrowGlobal = get_global_id(0);DrowGlobal < Nrow; DrowGlobal += get_global_size(0)){\n"

//  " async_work_group_copy(Acache, &A[AHere+Diter], get_local_size(0), 0);\n"
  "if(doCacheA==0){Acache[get_local_id(0)] = A[AHere + get_global_id(0)];}\n"
  "barrier(CLK_LOCAL_MEM_FENCE);\n"
 
  "  CHere = Dmatrix*NpadBetweenMatricesC + DrowGlobal * NpadC;\n";
    
  if(!sameB){
    result +=   "  BHere = Dmatrix*NpadBetweenMatricesB + DrowGlobal * NpadB;\n";
  } else {
    result +=   "  BHere = DrowGlobal * NpadB;\n";
  }
  
  result +=
  " AforThisWorkitem = Acache[get_local_id(0)];\n"
//  " AforThisWorkitem = A[AHere + DrowGlobal];\n"
    
  " for(Dcol = get_global_id(1); Dcol < Ncol; Dcol += get_global_size(1)){\n"
    
  "   C[CHere+Dcol] = B[BHere+Dcol] * AforThisWorkitem;\n"
  
  " }\n"// Dcol
  "}\n"// Drow
  "}\n"//Dmatrix
  "}\n";//kernel
  return(result);
}


// C = A B, A lower triangular

template <typename T> 
std::string multiplyLowerBatchString(
  const int sameB,
  const int diagIsOne,
  const int Nrow, 
  const int Ncol,
  const int Nmatrix, 
  const int NpadC, 
  const int NpadA, 
  const int NpadB,
  const int NpadD, // set to zero to omit D
  const std::string transformD,
  const int NpadBetweenMatricesC,
  const int NpadBetweenMatricesA,
  const int NpadBetweenMatricesB,
  const int NlocalCacheA,  // greater than Nlocal(0), smaller than Nrow
  const int NlocalCacheB,
  const int NcolInCache
  ) {

  std::string typeString = openclTypeString<T>();
  std::string result = "";
  const int rowsToCache = std::floor(NlocalCacheB/Ncol);

  if(typeString == "double") {
  	result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
  }
  result = result + 
"\n#define Nrow " + std::to_string(Nrow) + "\n"    
"#define Ncol " + std::to_string(Ncol) + "\n"    
"#define DinnerStop " + std::to_string(std::min(Nrow, 
                             rowsToCache)) + "\n"
"#define Nmatrix " + std::to_string(Nmatrix) + "\n"    
"#define NpadC " + std::to_string(NpadC) + "\n"    
"#define NpadA " + std::to_string(NpadA) + "\n"    
"#define NpadB " + std::to_string(NpadB) + "\n"
"#define NpadD " + std::to_string(NpadD) + "\n"
"#define NpadBetweenMatricesC " + std::to_string(NpadBetweenMatricesC) + "\n"    
"#define NpadBetweenMatricesA " + std::to_string(NpadBetweenMatricesA) + "\n"    
"#define NpadBetweenMatricesB " + std::to_string(NpadBetweenMatricesB) + "\n"    
"#define NlocalCache " + std::to_string(NlocalCacheB) + "\n"    
"#define NlocalCacheA " + std::to_string(NlocalCacheA) + "\n"
"#define NcolInCache " + std::to_string(NcolInCache) + "\n\n";
  

// global work items are rows, columns, matrices
// local work items anything, anything, 1
result += "__kernel void multiplyLowerBatch(\n"
  "	__global " + typeString+ " *C,\n"
  "	__global "+ typeString+ " *A,\n";

if(NpadD) {
  result += "	__global "+ typeString+ " *D,\n";
}

result +=  "	__global "+ typeString+ " *B) {\n\n" +
  typeString + " Dout;\n"
  "int AHere, BHere=0, CHere, BcacheHere=0;\n"
  "local "+ typeString+ " Acache[NlocalCacheA];\n" 
  "local "+ typeString+ " Bcache[NlocalCache];\n"
  "int Dmatrix, Drow, Dcol, Dinner, DrowNpadC, DrowNpadA, DcolInCache, DinnerCache;\n"
  "const int incInCacheLocal = get_local_size(0)*NcolInCache;\n"
  "const int incInCacheGlobal = get_global_size(0)*NcolInCache;\n"
  "const int doCacheA = (get_local_id(1) == 0);\n";

if(NpadD) {
  result += "local "+ typeString+ " Dcache[NlocalCacheA];\n"
  "int DHere;\n";
}


result +=  "\n\n"
"for(Dmatrix = get_global_id(2); Dmatrix < Nmatrix; Dmatrix += get_global_size(2)){\n"
  "  AHere = Dmatrix*NpadBetweenMatricesA;\n"
  "  CHere = Dmatrix*NpadBetweenMatricesC;\n";
if(NpadD) {
  result += "  DHere = Dmatrix * NpadD;\n";
}

  if(!sameB){
    // need to cache this B 
    result +=  
      "  BHere = Dmatrix*NpadBetweenMatricesB;\n";
  }
  
  result += "\n  // cache first rows of B for Dmatrix\n";
  result +=  
    "  for(Drow = get_local_id(0),BcacheHere = get_local_id(0)*NcolInCache;\n"
    "      Drow < DinnerStop; Drow += get_local_size(0),BcacheHere += incInCacheLocal){\n"
    
    "    DrowNpadA = BHere+Drow*NpadB;\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  
  if(NpadD) {
    result +=  "    if(doCacheA) {\n"
    "      Dcache[get_local_id(0)] = " +
      transformD + "(D[DHere + Drow]);\n"
    "    }\n"
    "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  }
    result +=     
      "    for(Dcol = get_global_id(1), DcolInCache=get_local_id(1);\n"
      "        Dcol < Ncol; Dcol += get_global_size(1), DcolInCache += get_local_size(1) ){\n";
    //     "Bcache[DcolInCache + Drow*NcolInCache] = B[Dcol + Drow * NpadB];\n"
    if(NpadD) {
      result +=     
        "      Bcache[BcacheHere + DcolInCache] = B[Dcol +DrowNpadA] * Dcache[get_local_id(0)];\n";
    } else {
      result +=     
        "      Bcache[BcacheHere + DcolInCache] = B[Dcol +DrowNpadA];\n";
    }
result +=
    "    } // Dcol\n" 
    "    barrier(CLK_LOCAL_MEM_FENCE);\n"
    "  } // Drow\n" 
    "  barrier(CLK_LOCAL_MEM_FENCE);\n";
    
result += 
  
"\n// looped through rows which are all cached\n"
"  for(Drow = get_global_id(0),BcacheHere = get_global_id(0)*NcolInCache;\n"
"        Drow < DinnerStop; Drow+=get_global_size(0),BcacheHere += incInCacheGlobal){\n"
  "    DrowNpadA= AHere + Drow * NpadA;\n"
  "    DrowNpadC= CHere+Drow * NpadC;\n"
  "    for(Dcol = get_global_id(1), DcolInCache=get_local_id(1);\n"
  "        Dcol < Ncol; Dcol += get_global_size(1),DcolInCache += get_local_size(1)){\n";
  
if(diagIsOne) {
  result += "      Dout = Bcache[BcacheHere + DcolInCache];\n";
  result +=   "      for(Dinner = 0,DinnerCache=0; Dinner < Drow; Dinner++,DinnerCache += NcolInCache){\n";
  
} else {
  result +=   "      Dout = 0.0;\n"
    "      for(Dinner = 0,DinnerCache=0; Dinner <= Drow; Dinner++,DinnerCache += NcolInCache){\n";
} 

  result +=  
  "        barrier(CLK_LOCAL_MEM_FENCE);\n"
  "        if(doCacheA) {\n"
  "           Acache[get_local_id(0)] = A[Dinner + DrowNpadA];\n"
  "        }\n"
	"        barrier(CLK_LOCAL_MEM_FENCE);\n"

	"        Dout += Acache[get_local_id(0)] * Bcache[DinnerCache + DcolInCache];\n"

  "      }// Dinner\n";

result +=  "      C[Dcol + DrowNpadC] = Dout;\n"
  "    }// Dcol\n" 
  "  } //Drow\n"
  "  barrier(CLK_LOCAL_MEM_FENCE);\n";
result +=  "\n// rows which are not all cached\n"

  
  //"\nfor(Drow = DinnerStop + get_global_id(0) ; Drow < Nrow; Drow+=get_global_size(0)){\n"
"  for(1 ; Drow < Nrow; Drow+=get_global_size(0),BcacheHere += incInCacheGlobal){\n"
  "    DrowNpadA= AHere + Drow * NpadA;\n"
  "    DrowNpadC= CHere + Drow * NpadC;\n";

if(NpadD) {
    result +=  
      "    barrier(CLK_LOCAL_MEM_FENCE);\n"
      "    if(doCacheA) {\n"
      "      Dcache[get_local_id(0)] = " +
                transformD + "(D[DHere + Drow]);\n"
      "    }\n"
      "    barrier(CLK_LOCAL_MEM_FENCE);\n";
  }

  result += "    for(Dcol = get_global_id(1), DcolInCache=get_local_id(1);\n"
"        Dcol < Ncol; Dcol += get_global_size(1),DcolInCache += get_local_size(1)){\n";

  result += "    // last row of B\n";
  if(diagIsOne) {
    if(NpadD) {
      result += "      Dout = B[BHere + Dcol + Drow * NpadB] * Dcache[get_local_id(0)];\n";
    } else {
      result += "      Dout = B[BHere + Dcol + Drow * NpadB];\n";
    }
  } else {
    result += "      Dout = 0.0;\n";
  }

  result += "    // cached rows of B\n";
  result +=  
  "      for(Dinner = 0,DinnerCache=0; Dinner < DinnerStop; Dinner++,DinnerCache += NcolInCache){\n"
  
  "        barrier(CLK_LOCAL_MEM_FENCE);\n"
  "        if(doCacheA) {\n"
	"          Acache[get_local_id(0)] = A[Dinner + DrowNpadA];\n"
	"        }\n"
	"        barrier(CLK_LOCAL_MEM_FENCE);\n"
	  
  "        Dout += Acache[get_local_id(0)] * Bcache[DinnerCache + DcolInCache];\n"
	"      }\n" // Dinner
  "      barrier(CLK_LOCAL_MEM_FENCE);\n"

  "\n// un-cached rows\n";
  if(diagIsOne) {
//	  result += "for(Dinner = NlocalCache; Dinner < Drow; Dinner++){\n";
     result += "    for(1; Dinner < Drow; Dinner++){\n";
	} else {
//	  result += "for(Dinner = NlocalCache; Dinner <= Drow; Dinner++){\n";
	   result += "    for(1; Dinner <= Drow; Dinner++){\n";
	}

	result += "    barrier(CLK_LOCAL_MEM_FENCE);\n"
  "    if(doCacheA) {\n";

		if(NpadD) {
	  result += 
	    "      Acache[get_local_id(0)] = " + transformD + 
	    "(D[DHere + Dinner]) * A[Dinner + DrowNpadA];\n";
	} else {
	  result += 
	    "      Acache[get_local_id(0)] = A[Dinner + DrowNpadA];\n";
	}
  result +=	"    }//do cacheA\n"
	"    barrier(CLK_LOCAL_MEM_FENCE);\n";

	 result += "// last row of B\n"
	   "     Dout += Acache[get_local_id(0)] * B[BHere+Dcol + Dinner * NpadB];\n";
	 
	result += "      }// Dinner\n" 
	"\n    C[Dcol + DrowNpadC] = Dout;\n";
	result +=	"    }// Dcol\n";
  result +=	"  }//Drow\n";

  result += "}// Dmatrix\n" 
  "}";

  return(result);
}




template <typename T> 
void multiplyLowerDiagonalBatch(
    viennacl::matrix<T> &C,
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &D,
    viennacl::matrix<T> &B,
    const int diagIsOne,
    const std::string transformD,
    std::vector<int> Nglobal,
    std::vector<int> Nlocal,
    const int NlocalCache, 
    const int ctx_id) {
  
  
  const int Nrow = A.size2(), Ncol = B.size2();
  const int Nmatrix = C.size1()/Nrow;
  const int Nrounds = std::ceil( static_cast<T>(Ncol) / static_cast<T>(Nglobal[1]));
  

  
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
//  cl_device_type type_check = ctx.current_device().type();
  
  std::string clString =  multiplyLowerBatchString<T>(  
    Nrow == (int) B.size1(),
    diagIsOne,
    Nrow, 
    Ncol, // ncol
    Nmatrix,
    C.internal_size2(), 
    A.internal_size2(), 
    B.internal_size2(),
    D.internal_size2(),
    transformD,
    C.internal_size2()*Nrow,//NpadBetweenMatricesC,
    A.internal_size2()*Nrow,//NpadBetweenMatricesA,
    B.internal_size2()*Nrow,//NpadBetweenMatricesB,
    Nlocal[0], 
    std::min(Nrow*Ncol,NlocalCache), //NcacheB
    Nlocal[1]*Nrounds);  //NcolsInCache
    
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif  
  


    viennacl::ocl::program & my_prog = ctx.add_program(
    clString, "my_kernel");
  
  viennacl::ocl::kernel & multiplyKernel = my_prog.get_kernel("multiplyLowerBatch");
  
  multiplyKernel.global_work_size(0, Nglobal[0]);
  multiplyKernel.global_work_size(1, Nglobal[1]);
  multiplyKernel.global_work_size(2, Nglobal[2]);
  
  multiplyKernel.local_work_size(0, Nlocal[0]);
  multiplyKernel.local_work_size(1, Nlocal[1]);
  multiplyKernel.local_work_size(2, 1L);//Nlocal[2]);

    // diagonals and diagTimesRowOfA
  viennacl::ocl::enqueue(multiplyKernel(
      C, A, D, B));

    }


template <typename T> 
void multiplyLowerBatch(
	viennacl::matrix<T> &C,
	viennacl::matrix<T> &A,
	viennacl::matrix<T> &B,
	const int diagIsOne,
	std::vector<int> Nglobal,
	std::vector<int> Nlocal, 
	const int NlocalCache, 
	const int ctx_id) {

  
  const int Nrow = A.size2(), Ncol = B.size2();
  
  const int Nmatrix = C.size1()/Nrow;

  	// the context
	viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));

//	cl_device_type type_check = ctx.current_device().type();

  std::string clString =  multiplyLowerBatchString<T>(  
      Nrow == (int) B.size1(),
      diagIsOne,
      Nrow, 
      Ncol,
      Nmatrix,
      (int) C.internal_size2(), 
      (int) A.internal_size2(), 
      (int) B.internal_size2(),
      0L, // NpadD
      "ignored",//transformD
      ((int) C.internal_size2())*Nrow,//NpadBetweenMatricesC,
       ( (int) A.internal_size2())*Nrow,//NpadBetweenMatricesA,
       ((int) B.internal_size2())*Nrow,//NpadBetweenMatricesB,
      Nlocal[0], 
      std::min(Nrow,NlocalCache), Ncol);
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif

 
	viennacl::ocl::program & my_prog = ctx.add_program(
	  clString, "my_kernel");

	viennacl::ocl::kernel & multiplyKernel = my_prog.get_kernel("multiplyLowerBatch");

	multiplyKernel.global_work_size(0, Nglobal[0]);
	multiplyKernel.global_work_size(1, Nglobal[1]);
  multiplyKernel.global_work_size(2, Nglobal[2]);

  multiplyKernel.local_work_size(0, Nlocal[0]);
  multiplyKernel.local_work_size(1, Nlocal[1]);
  multiplyKernel.local_work_size(2, 1L);//Nlocal[2]);
  
		// diagonals and diagTimesRowOfA
		viennacl::ocl::enqueue(multiplyKernel(C, A, B));

}



template <typename T> 
void multiplyDiagonalBatch(
    viennacl::matrix<T> &C,
    viennacl::matrix<T> &A,
    viennacl::matrix<T> &B,
    std::vector<int> Nglobal,   // work items are Drow, Dcol, Dmatrix
    std::vector<int> Nlocal,   // local size is x, y, 1
    const int ctx_id) {
  
  const int Nrow = A.size2();
  const int Nmatrix = ((int) C.size1())/Nrow;
  // the context
  viennacl::ocl::context ctx(viennacl::ocl::get_context(ctx_id));
  
//  cl_device_type type_check = ctx.current_device().type();
  
  std::string clString =  multiplyDiagonalBatchString<T>(  
    ((int) B.size1()) == Nrow,
      Nrow, 
      (int) B.size2(), // Ncol
      Nmatrix, 
      (int) C.internal_size2(), 
      (int) A.internal_size2(), 
      (int) B.internal_size2(),
      ((int) C.internal_size2())*Nrow,//NpadBetweenMatricesC,
       ((int) B.internal_size2())*Nrow,//NpadBetweenMatricesB,
      Nlocal[0]);

  viennacl::ocl::program & my_prog = ctx.add_program(
    clString, "my_kernel");
  
#ifdef DEBUG
  
  Rcpp::Rcout << clString << "\n\n";
  
#endif  
  

  viennacl::ocl::kernel & multiplyKernel = my_prog.get_kernel("multiplyDiagonalBatch");
  
  multiplyKernel.global_work_size(0, Nglobal[0]);
  multiplyKernel.global_work_size(1, Nglobal[1]);
  multiplyKernel.global_work_size(2, Nglobal[2]);
  
  multiplyKernel.local_work_size(0, Nlocal[0]);
  multiplyKernel.local_work_size(1, Nlocal[1]);
  multiplyKernel.local_work_size(2, 1L);//Nlocal[2]);
  
  // diagonals and diagTimesRowOfA
  viennacl::ocl::enqueue(multiplyKernel(C, A, B));
  
}



template <typename T> 
SEXP multiplyLowerDiagonalBatchTyped(
    Rcpp::S4 CR,
    Rcpp::S4 AR,
    Rcpp::S4 DR,
    Rcpp::S4 BR,
    const int diagIsOne,
    std::string transformD,
    Rcpp::IntegerVector NglobalR,
    Rcpp::IntegerVector NlocalR, 
    const int NlocalCache) {

  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  
  
  std::shared_ptr<viennacl::matrix<T> > 
    AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    BG = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    DG = getVCLptr<T>(DR.slot("address"), BisVCL, ctx_id);
  
  multiplyLowerDiagonalBatch<T>(*CG, *AG, *DG, *BG, 
                                diagIsOne, transformD, 
                        Nglobal, Nlocal, NlocalCache, ctx_id);
  
  return Rcpp::wrap(0L);

}

template <typename T> 
SEXP multiplyLowerBatchTyped(
	Rcpp::S4 CR,
	Rcpp::S4 AR,
	Rcpp::S4 BR,
	const int diagIsOne,
	Rcpp::IntegerVector NglobalR,
	Rcpp::IntegerVector NlocalR, 
	const int NlocalCache) {

  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);

	const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
	const bool BisVCL=1;


	
	std::shared_ptr<viennacl::matrix<T> > 
		AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > 
		BG = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
	std::shared_ptr<viennacl::matrix<T> > 
		CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);

	multiplyLowerBatch<T>(*CG, *AG, *BG, diagIsOne, 
                       Nglobal, Nlocal, NlocalCache, ctx_id);	

	return Rcpp::wrap(0L);
}

template <typename T> 
SEXP multiplyDiagonalBatchTyped(
    Rcpp::S4 CR,
    Rcpp::S4 AR,
    Rcpp::S4 BR,
    Rcpp::IntegerVector NglobalR,
    Rcpp::IntegerVector NlocalR) {
  
  std::vector<int> Nglobal = Rcpp::as<std::vector<int> >(NglobalR);
  std::vector<int> Nlocal = Rcpp::as<std::vector<int> >(NlocalR);
  
  const int ctx_id = INTEGER(CR.slot(".context_index"))[0]-1;
  const bool BisVCL=1;
  
  
  std::shared_ptr<viennacl::matrix<T> > 
    AG = getVCLptr<T>(AR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    BG = getVCLptr<T>(BR.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::matrix<T> > 
    CG = getVCLptr<T>(CR.slot("address"), BisVCL, ctx_id);
  
  multiplyDiagonalBatch<T>(*CG, *AG, *BG, Nglobal, Nlocal, ctx_id);	
  
  return Rcpp::wrap(0L);
}
// [[Rcpp::export]]
SEXP multiplyLowerDiagonalBatchBackend(
    Rcpp::S4 C,
    Rcpp::S4 A,
    Rcpp::S4 D,
    Rcpp::S4 B,
    const int diagIsOne,
    std::string transformD,
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal,
    const int NlocalCache) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  
  
  if(precision_type == "fvclMatrix") {
    result = multiplyLowerDiagonalBatchTyped<float>(C, A, D, B, diagIsOne, transformD, Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = multiplyLowerDiagonalBatchTyped<double>(C, A, D, B, diagIsOne, transformD,Nglobal, Nlocal,NlocalCache);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
}

// [[Rcpp::export]]
SEXP multiplyDiagonalBatchBackend(
    Rcpp::S4 C,
    Rcpp::S4 A,
    Rcpp::S4 B,
    Rcpp::IntegerVector Nglobal,
    Rcpp::IntegerVector Nlocal) {
  
  SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  
  if(precision_type == "fvclMatrix") {
    result = multiplyDiagonalBatchTyped<float>(C, A, B, Nglobal, Nlocal);
  } else if (precision_type == "dvclMatrix") {
    result = multiplyDiagonalBatchTyped<double>(C, A, B, Nglobal, Nlocal);
  } else {
    result = Rcpp::wrap(1L);
  }
  return(result);
  
}

//' Multiply lower triangular matrices
//' 
//' Multiplies a lower triangular matrix by a rectangular matrix
//'
//' @param C output matrices, stacked row-wise
//' @param A lower triangular matrices
//' @param B rectangular matrix or matrices
//' @param Nglobal vector of number of global work items
//' @param Nlocal vector of number of local work items
//' @param NlocalCache elements in local cache
//' @export
// [[Rcpp::export]]
SEXP multiplyLowerBatchBackend(
	Rcpp::S4 C,
	Rcpp::S4 A,
	Rcpp::S4 B,
	const int diagIsOne,
	Rcpp::IntegerVector Nglobal,
	Rcpp::IntegerVector Nlocal, 
	const int NlocalCache) {

	SEXP result;
  
  Rcpp::traits::input_parameter< std::string >::type classVarR(RCPP_GET_CLASS(C));
  std::string precision_type = (std::string) classVarR;
  

  if(precision_type == "fvclMatrix") {
    result = multiplyLowerBatchTyped<float>(
      C, A, B, diagIsOne, 
      Nglobal, Nlocal, NlocalCache);
  } else if (precision_type == "dvclMatrix") {
    result = multiplyLowerBatchTyped<double>(
      C, A, B, diagIsOne, 
      Nglobal, Nlocal, NlocalCache);
	} else {
		result = Rcpp::wrap(1L);
	}
	return(result);

}
