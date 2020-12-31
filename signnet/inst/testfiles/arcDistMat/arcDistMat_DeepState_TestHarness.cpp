#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix arcDistMat(NumericMatrix X, double r);

TEST(signnet_deepstate_test,arcDistMat_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix X  = RcppDeepState_NumericMatrix();
  qs::c_qsave(X,"/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/arcDistMat/inputs/X.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "X values: "<< X << std::endl;
  std::ofstream r_stream;
  double r  = RcppDeepState_double();
  r_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/arcDistMat/inputs/r");
  r_stream << r;
  std::cout << "r values: "<< r << std::endl;
  r_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    arcDistMat(X,r);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
