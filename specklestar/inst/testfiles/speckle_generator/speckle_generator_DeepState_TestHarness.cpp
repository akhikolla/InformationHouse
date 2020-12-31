// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// speckle_generator_DeepState_TestHarness_generation.cpp and speckle_generator_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

NumericVector speckle_generator(double rho, double theta, double dm, double seeing, double speckle_sigma, double wind);

TEST(specklestar_deepstate_test,speckle_generator_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector rho(1);
  rho[0]  = RcppDeepState_double();
  qs::c_qsave(rho,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/rho.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "rho values: "<< rho << std::endl;
  NumericVector theta(1);
  theta[0]  = RcppDeepState_double();
  qs::c_qsave(theta,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/theta.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "theta values: "<< theta << std::endl;
  NumericVector dm(1);
  dm[0]  = RcppDeepState_double();
  qs::c_qsave(dm,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/dm.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "dm values: "<< dm << std::endl;
  NumericVector seeing(1);
  seeing[0]  = RcppDeepState_double();
  qs::c_qsave(seeing,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/seeing.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "seeing values: "<< seeing << std::endl;
  NumericVector speckle_sigma(1);
  speckle_sigma[0]  = RcppDeepState_double();
  qs::c_qsave(speckle_sigma,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/speckle_sigma.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "speckle_sigma values: "<< speckle_sigma << std::endl;
  NumericVector wind(1);
  wind[0]  = RcppDeepState_double();
  qs::c_qsave(wind,"/home/akhila/fuzzer_packages/fuzzedpackages/specklestar/inst/testfiles/speckle_generator/inputs/wind.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "wind values: "<< wind << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    speckle_generator(rho[0],theta[0],dm[0],seeing[0],speckle_sigma[0],wind[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}