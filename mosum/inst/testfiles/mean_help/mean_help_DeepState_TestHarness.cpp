#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double mean_help(NumericVector x, int l, int r);

TEST(mosum_deepstate_test,mean_help_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/mean_help/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  std::ofstream l_stream;
  int l  = RcppDeepState_int();
  l_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/mean_help/inputs/l");
  l_stream << l;
  std::cout << "l values: "<< l << std::endl;
  l_stream.close();
  std::ofstream r_stream;
  int r  = RcppDeepState_int();
  r_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/mosum/inst/testfiles/mean_help/inputs/r");
  r_stream << r;
  std::cout << "r values: "<< r << std::endl;
  r_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    mean_help(x,l,r);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
