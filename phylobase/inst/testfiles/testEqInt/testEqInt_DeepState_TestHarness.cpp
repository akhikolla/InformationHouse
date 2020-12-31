#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

Rcpp::List testEqInt(Rcpp::IntegerVector x, Rcpp::IntegerVector y);

TEST(phylobase_deepstate_test,testEqInt_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  IntegerVector x  = RcppDeepState_IntegerVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/testEqInt/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  IntegerVector y  = RcppDeepState_IntegerVector();
  qs::c_qsave(y,"/home/akhila/fuzzer_packages/fuzzedpackages/phylobase/inst/testfiles/testEqInt/inputs/y.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y values: "<< y << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    testEqInt(x,y);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
