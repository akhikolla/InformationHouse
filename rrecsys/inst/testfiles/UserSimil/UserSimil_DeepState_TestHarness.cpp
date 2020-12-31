#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix UserSimil(NumericMatrix x, int damp);

TEST(rrecsys_deepstate_test,UserSimil_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix x  = RcppDeepState_NumericMatrix();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/UserSimil/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::ofstream damp_stream;
  int damp  = RcppDeepState_int();
  damp_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rrecsys/inst/testfiles/UserSimil/inputs/damp");
  damp_stream << damp;
  std::cout << "damp values: "<< damp << std::endl;
  damp_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    UserSimil(x,damp);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
