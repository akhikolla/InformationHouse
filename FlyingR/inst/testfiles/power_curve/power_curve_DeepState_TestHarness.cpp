// AUTOMATICALLY GENERATED BY RCPPDEEPSTATE PLEASE DO NOT EDIT BY HAND, INSTEAD EDIT
// power_curve_DeepState_TestHarness_generation.cpp and power_curve_DeepState_TestHarness_checks.cpp

#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double power_curve(double bm, double ws, double wa, double tas, double g, double airDensity, double ipf, double bdc, double ppc);

TEST(FlyingR_deepstate_test,power_curve_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  NumericVector bm(1);
  bm[0]  = RcppDeepState_double();
  qs::c_qsave(bm,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/bm.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "bm values: "<< bm << std::endl;
  NumericVector ws(1);
  ws[0]  = RcppDeepState_double();
  qs::c_qsave(ws,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/ws.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ws values: "<< ws << std::endl;
  NumericVector wa(1);
  wa[0]  = RcppDeepState_double();
  qs::c_qsave(wa,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/wa.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "wa values: "<< wa << std::endl;
  NumericVector tas(1);
  tas[0]  = RcppDeepState_double();
  qs::c_qsave(tas,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/tas.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "tas values: "<< tas << std::endl;
  NumericVector g(1);
  g[0]  = RcppDeepState_double();
  qs::c_qsave(g,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/g.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "g values: "<< g << std::endl;
  NumericVector airDensity(1);
  airDensity[0]  = RcppDeepState_double();
  qs::c_qsave(airDensity,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/airDensity.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "airDensity values: "<< airDensity << std::endl;
  NumericVector ipf(1);
  ipf[0]  = RcppDeepState_double();
  qs::c_qsave(ipf,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/ipf.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ipf values: "<< ipf << std::endl;
  NumericVector bdc(1);
  bdc[0]  = RcppDeepState_double();
  qs::c_qsave(bdc,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/bdc.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "bdc values: "<< bdc << std::endl;
  NumericVector ppc(1);
  ppc[0]  = RcppDeepState_double();
  qs::c_qsave(ppc,"/home/akhila/fuzzer_packages/fuzzedpackages/FlyingR/inst/testfiles/power_curve/inputs/ppc.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "ppc values: "<< ppc << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    power_curve(bm[0],ws[0],wa[0],tas[0],g[0],airDensity[0],ipf[0],bdc[0],ppc[0]);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}