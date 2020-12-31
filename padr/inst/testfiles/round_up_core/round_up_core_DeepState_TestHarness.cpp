#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

IntegerVector round_up_core(IntegerVector a, IntegerVector b);

TEST(padr_deepstate_test,round_up_core_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector a  = RcppDeepState_IntegerVector();
  qs::c_qsave(a,"/home/akhila/fuzzer_packages/fuzzedpackages/padr/inst/testfiles/round_up_core/inputs/a.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "a values: "<< a << std::endl;
  IntegerVector b  = RcppDeepState_IntegerVector();
  qs::c_qsave(b,"/home/akhila/fuzzer_packages/fuzzedpackages/padr/inst/testfiles/round_up_core/inputs/b.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "b values: "<< b << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    round_up_core(a,b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
