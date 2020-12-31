#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericMatrix simulateBS(NumericVector param, int ndays);

TEST(pinbasic_deepstate_test,simulateBS_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector param  = RcppDeepState_NumericVector();
  qs::c_qsave(param,"/home/akhila/fuzzer_packages/fuzzedpackages/pinbasic/inst/testfiles/simulateBS/inputs/param.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "param values: "<< param << std::endl;
  std::ofstream ndays_stream;
  int ndays  = RcppDeepState_int();
  ndays_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/pinbasic/inst/testfiles/simulateBS/inputs/ndays");
  ndays_stream << ndays;
  std::cout << "ndays values: "<< ndays << std::endl;
  ndays_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    simulateBS(param,ndays);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
