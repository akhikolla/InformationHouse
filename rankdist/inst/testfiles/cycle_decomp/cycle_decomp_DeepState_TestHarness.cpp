#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

int cycle_decomp(NumericVector comp);

TEST(rankdist_deepstate_test,cycle_decomp_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector comp  = RcppDeepState_NumericVector();
  qs::c_qsave(comp,"/home/akhila/fuzzer_packages/fuzzedpackages/rankdist/inst/testfiles/cycle_decomp/inputs/comp.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "comp values: "<< comp << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    cycle_decomp(comp);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
