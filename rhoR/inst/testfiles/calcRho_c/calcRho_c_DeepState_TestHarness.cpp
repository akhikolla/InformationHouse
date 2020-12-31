#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double calcRho_c(double x, double OcSBaserate, int testSetLength, double testSetBaserateInflation, int OcSLength, int replicates, double ScSKappaThreshold, double ScSKappaMin, double ScSPrecisionMin, double ScSPrecisionMax, NumericMatrix KPs);

TEST(rhoR_deepstate_test,calcRho_c_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream x_stream;
  double x  = RcppDeepState_double();
  x_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/x");
  x_stream << x;
  std::cout << "x values: "<< x << std::endl;
  x_stream.close();
  std::ofstream OcSBaserate_stream;
  double OcSBaserate  = RcppDeepState_double();
  OcSBaserate_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/OcSBaserate");
  OcSBaserate_stream << OcSBaserate;
  std::cout << "OcSBaserate values: "<< OcSBaserate << std::endl;
  OcSBaserate_stream.close();
  std::ofstream testSetLength_stream;
  int testSetLength  = RcppDeepState_int();
  testSetLength_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/testSetLength");
  testSetLength_stream << testSetLength;
  std::cout << "testSetLength values: "<< testSetLength << std::endl;
  testSetLength_stream.close();
  std::ofstream testSetBaserateInflation_stream;
  double testSetBaserateInflation  = RcppDeepState_double();
  testSetBaserateInflation_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/testSetBaserateInflation");
  testSetBaserateInflation_stream << testSetBaserateInflation;
  std::cout << "testSetBaserateInflation values: "<< testSetBaserateInflation << std::endl;
  testSetBaserateInflation_stream.close();
  std::ofstream OcSLength_stream;
  int OcSLength  = RcppDeepState_int();
  OcSLength_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/OcSLength");
  OcSLength_stream << OcSLength;
  std::cout << "OcSLength values: "<< OcSLength << std::endl;
  OcSLength_stream.close();
  std::ofstream replicates_stream;
  int replicates  = RcppDeepState_int();
  replicates_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/replicates");
  replicates_stream << replicates;
  std::cout << "replicates values: "<< replicates << std::endl;
  replicates_stream.close();
  std::ofstream ScSKappaThreshold_stream;
  double ScSKappaThreshold  = RcppDeepState_double();
  ScSKappaThreshold_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/ScSKappaThreshold");
  ScSKappaThreshold_stream << ScSKappaThreshold;
  std::cout << "ScSKappaThreshold values: "<< ScSKappaThreshold << std::endl;
  ScSKappaThreshold_stream.close();
  std::ofstream ScSKappaMin_stream;
  double ScSKappaMin  = RcppDeepState_double();
  ScSKappaMin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/ScSKappaMin");
  ScSKappaMin_stream << ScSKappaMin;
  std::cout << "ScSKappaMin values: "<< ScSKappaMin << std::endl;
  ScSKappaMin_stream.close();
  std::ofstream ScSPrecisionMin_stream;
  double ScSPrecisionMin  = RcppDeepState_double();
  ScSPrecisionMin_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/ScSPrecisionMin");
  ScSPrecisionMin_stream << ScSPrecisionMin;
  std::cout << "ScSPrecisionMin values: "<< ScSPrecisionMin << std::endl;
  ScSPrecisionMin_stream.close();
  std::ofstream ScSPrecisionMax_stream;
  double ScSPrecisionMax  = RcppDeepState_double();
  ScSPrecisionMax_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/ScSPrecisionMax");
  ScSPrecisionMax_stream << ScSPrecisionMax;
  std::cout << "ScSPrecisionMax values: "<< ScSPrecisionMax << std::endl;
  ScSPrecisionMax_stream.close();
  NumericMatrix KPs  = RcppDeepState_NumericMatrix();
  qs::c_qsave(KPs,"/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/calcRho_c/inputs/KPs.qs",
		"high", "zstd", 1, 15, true, 1);
  std::cout << "KPs values: "<< KPs << std::endl;
  std::cout << "input ends" << std::endl;
  try{
    calcRho_c(x,OcSBaserate,testSetLength,testSetBaserateInflation,OcSLength,replicates,ScSKappaThreshold,ScSKappaMin,ScSPrecisionMin,ScSPrecisionMax,KPs);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
