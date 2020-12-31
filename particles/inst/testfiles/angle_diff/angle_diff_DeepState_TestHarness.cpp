#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector angle_diff(NumericMatrix a, NumericMatrix b);

TEST(particles_deepstate_test,angle_diff_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericMatrix a  = RcppDeepState_NumericMatrix();
  qs::c_qsave(a,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/angle_diff/inputs/a.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "a values: "<< a << std::endl;
  NumericMatrix b  = RcppDeepState_NumericMatrix();
  qs::c_qsave(b,"/home/akhila/fuzzer_packages/fuzzedpackages/particles/inst/testfiles/angle_diff/inputs/b.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "b values: "<< b << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    angle_diff(a,b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
