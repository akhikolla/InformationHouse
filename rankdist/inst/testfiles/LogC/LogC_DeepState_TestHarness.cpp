#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double LogC(NumericVector fai);

TEST(rankdist_deepstate_test,LogC_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector fai  = RcppDeepState_NumericVector();
  qs::c_qsave(fai,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/LogC/inputs/fai.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "fai values: "<< fai << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    LogC(fai);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
