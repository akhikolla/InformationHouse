#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

int setBitNumber(int n);

TEST(mosum_deepstate_test,setBitNumber_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream n_stream;
  int n  = RcppDeepState_int();
  n_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/setBitNumber/inputs/n");
  n_stream << n;
  std::cout << "n values: "<< n << std::endl;
  n_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    setBitNumber(n);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
