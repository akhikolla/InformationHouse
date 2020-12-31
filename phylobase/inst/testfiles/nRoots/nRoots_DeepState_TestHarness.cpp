#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

int nRoots(Rcpp::IntegerVector ances);

TEST(phylobase_deepstate_test,nRoots_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector ances  = RcppDeepState_IntegerVector();
  qs::c_qsave(ances,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/nRoots/inputs/ances.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ances values: "<< ances << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    nRoots(ances);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
