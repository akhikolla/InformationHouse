#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

List labRcpp(int nfeats);

TEST(propr_deepstate_test,labRcpp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream nfeats_stream;
  int nfeats  = RcppDeepState_int();
  nfeats_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/propr/inst/testfiles/labRcpp/inputs/nfeats");
  nfeats_stream << nfeats;
  std::cout << "nfeats values: "<< nfeats << std::endl;
  nfeats_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    labRcpp(nfeats);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
