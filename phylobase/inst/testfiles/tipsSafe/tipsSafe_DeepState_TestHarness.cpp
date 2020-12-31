#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::IntegerVector tipsSafe(Rcpp::IntegerVector ances, Rcpp::IntegerVector desc);

TEST(phylobase_deepstate_test,tipsSafe_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector ances  = RcppDeepState_IntegerVector();
  qs::c_qsave(ances,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/tipsSafe/inputs/ances.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ances values: "<< ances << std::endl;
  IntegerVector desc  = RcppDeepState_IntegerVector();
  qs::c_qsave(desc,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/tipsSafe/inputs/desc.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "desc values: "<< desc << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    tipsSafe(ances,desc);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
