#include "gpuRandom.hpp"





std::string colsumRowsumString(
    const int Nrow, 
    const int Ncol,
    const int NpadCol) { 
  
  std::string typeString = "int";
  std::string typeStringSum = "double"; // type of the sum of log factorial

  std::string result = "";
  
  if(typeStringSum == "double") {
    
    result += "\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
    
  }

  
  result +=
    "\n#define Nrow " + std::to_string(Nrow) + "\n"    
    "#define Ncol " + std::to_string(Ncol) + "\n"
    "#define NpadCol " + std::to_string(NpadCol) + "\n";    
  
    result += mrg31k3pString();
  
  result += 
    "\n\n__kernel void colsumRowsum(\n"
    "  __global " + typeString + "* x,"  
    "  __global " + typeString + "* rowSum,"  
    "  __global " + typeString + "* colSum"
    "){\n\n";  
  
  
 
  result += "int Drow, Dcol, Dindex;\n";
  result += typeString + " result;\n";
  
  result + "if(get_global_id(1) { // sum columns\n"
  "  for(Drow = get_global_id(0); Drow < Nrow; Drow++){\n"
  "    result = 0;\n"
  "    for(Dcol = 0, Dindex = Drow*NpadCol; Dcol < Ncol; Dcol++, Dindex++){\n"
  "       result += x[Dindex];\n"    
  "    } // end loop through columns\n"
  "    colSum[Dcol] = result;\n"
  "  } // end loop through rows\n"
  "} else { // sum rows\n"
  "  for(Dcol = get_global_id(0); Dcol < Ncol; Dcol++){\n"
  "    result = 0;\n"
  "    for(Drow = 0, Dindex = Dcol; Drow < Nrow; Drow++, Dindex+= NpadCol){\n"
  "       result += x[Dindex];\n"    
  "    } // end loop through columns\n"
  "    rowSum[Drow] = result;\n"
  "  } // end loop through columns\n"

  "}\n\n";
  result += 
    "}//kernel\n";

  result += 
    "\n\n__kernel void sumLfactorial(\n"
    "  __global " + typeString + "* x,"  
    "  __global " + typeStringSum + "* result"  
    "){\n\n";  

// TO DO: use groups and local memory

  result += "int Drow, Dcol, Dindex;\n";
  result += typeStringSum + " Dresult=0;\n";

  result += 
 "  for(Drow = get_global_id(0); Drow < Nrow; Drow+=get_global_size(0)){\n"
  "    for(Dcol = get_global_id(1), Dindex = Drow*NpadCol;" 
  "        Dcol < Ncol; Dcol+=get_global_size(1), Dindex++){\n"
  "       Dresult += lgamma(1 + x[Dindex]);"
  "    } // end loop through columns\n"
  "  } // end loop through rows\n";

  result += 
  "result[get_global_id(1) + get_global_id(0)*get_global_size(1)] = resultD;\n";

  result += 
    "}//sumLfactorial kernel\n";

  return(result);
}

//' @export
// [[Rcpp::export]]
SEXP colsumRowsumBackend(
    Rcpp::S4  x,
    Rcpp::S4  rowSum,
    Rcpp::S4  colSum,   
    Rcpp::IntegerVector numWorkItems) {
  
  double result;
   
  const bool BisVCL=1;
  const int ctx_id = INTEGER(x.slot(".context_index"))[0]-1;
  
  std::shared_ptr<viennacl::matrix<int> > xVcl = getVCLptr<int>(x.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<int> > rowSumVcl = getVCLVecptr<int>(rowSum.slot("address"), BisVCL, ctx_id);
  std::shared_ptr<viennacl::vector_base<int> > colSumVcl = getVCLVecptr<int>(colSum.slot("address"), BisVCL, ctx_id);

  std::string sumKernelString = colsumRowsumString(
    (*xVcl).size1(), (*xVcl).size2(),
    (*xVcl).internal_size2() );

    // the context
  viennacl::ocl::switch_context(ctx_id);
  viennacl::ocl::program & my_prog = viennacl::ocl::current_context().add_program(sumKernelString, "my_kernel");
  
  
  viennacl::ocl::kernel &sumKernel = my_prog.get_kernel("colsumRowsum");


  sumKernel.global_work_size(0, numWorkItems[0]);
  sumKernel.global_work_size(1, numWorkItems[1]);
  
  sumKernel.local_work_size(0, 1L);
  sumKernel.local_work_size(1, 1L);

  viennacl::ocl::enqueue(sumKernel(*xVcl, *rowSumVcl, *colSumVcl) );

  viennacl::ocl::kernel &sumLfactorialKernel = my_prog.get_kernel("sumLfactorial");
  sumLfactorialKernel.global_work_size(0, numWorkItems[0]);
  sumLfactorialKernel.global_work_size(1, numWorkItems[1]);
  
  sumLfactorialKernel.local_work_size(0, 1L);
  sumLfactorialKernel.local_work_size(1, 1L);

  viennacl::vector_base<double> logFactorial(numWorkItems[0] * numWorkItems[1]);

  viennacl::ocl::enqueue(sumLfactorialKernel(*xVcl, logFactorial) );


  result = viennacl::linalg::sum(logFactorial);

  return(Rcpp::wrap(result));

}




