#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector bootstrapped_timeSeries(IntegerVector cpts, NumericVector x);

TEST(mosum_deepstate_test,bootstrapped_timeSeries_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector cpts  = RcppDeepState_IntegerVector();
  qs::c_qsave(cpts,"/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/bootstrapped_timeSeries/inputs/cpts.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "cpts values: "<< cpts << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/bootstrapped_timeSeries/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    bootstrapped_timeSeries(cpts,x);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
