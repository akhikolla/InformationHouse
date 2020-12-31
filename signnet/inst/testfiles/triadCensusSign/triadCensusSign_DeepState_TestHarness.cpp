#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

IntegerVector triadCensusSign(NumericMatrix A, int n);

TEST(signnet_deepstate_test,triadCensusSign_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix A  = RcppDeepState_NumericMatrix();
  qs::c_qsave(A,"/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/triadCensusSign/inputs/A.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "A values: "<< A << std::endl;
  std::ofstream n_stream;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/triadCensusSign/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    triadCensusSign(A,n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
