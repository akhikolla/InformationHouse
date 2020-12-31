#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double arcDist(NumericVector x, NumericVector y, double r);

TEST(signnet_deepstate_test,arcDist_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector x  = RcppDeepState_NumericVector();
  qs::c_qsave(x,"/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/arcDist/inputs/x.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "x values: "<< x << std::endl;
  NumericVector y  = RcppDeepState_NumericVector();
  qs::c_qsave(y,"/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/arcDist/inputs/y.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "y values: "<< y << std::endl;
  std::ofstream r_stream;
  double r  = RcppDeepState_double();
  r_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/signnet/inst/testfiles/arcDist/inputs/r");
  r_stream << r;
  std::cout << "r values: "<< r << std::endl;
  r_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    arcDist(x,y,r);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
