#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List indexToCoord(IntegerVector V, const int N);

TEST(propr_deepstate_test,indexToCoord_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector V  = RcppDeepState_IntegerVector();
  qs::c_qsave(V,"/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/indexToCoord/inputs/V.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "V values: "<< V << std::endl;
  std::ofstream N_stream;
  int N  = RcppDeepState_int();
  N_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/indexToCoord/inputs/N");
  N_stream << N;
  std::cout << "N values: "<< N << std::endl;
  N_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    indexToCoord(V,N);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
