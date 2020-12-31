#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix ItemSimilSparseMat(NumericMatrix x, int dim, int DAMP);

TEST(rrecsys_deepstate_test,ItemSimilSparseMat_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/ItemSimilSparseMat/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::ofstream dim_stream;
  int dim  = RcppDeepState_int();
  dim_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/ItemSimilSparseMat/inputs/dim");
  dim_stream << dim;
  std::cout << "dim values: "<< dim << std::endl;
  dim_stream.close();
  std::ofstream DAMP_stream;
  int DAMP  = RcppDeepState_int();
  DAMP_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/ItemSimilSparseMat/inputs/DAMP");
  DAMP_stream << DAMP;
  std::cout << "DAMP values: "<< DAMP << std::endl;
  DAMP_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    ItemSimilSparseMat(x,dim,DAMP);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
