#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

double recall(double kappa, double BR, double P);

TEST(rhoR_deepstate_test,recall_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream kappa_stream;
  double kappa  = RcppDeepState_double();
  kappa_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/recall/inputs/kappa");
  kappa_stream << kappa;
  std::cout << "kappa values: "<< kappa << std::endl;
  kappa_stream.close();
  std::ofstream BR_stream;
  double BR  = RcppDeepState_double();
  BR_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/recall/inputs/BR");
  BR_stream << BR;
  std::cout << "BR values: "<< BR << std::endl;
  BR_stream.close();
  std::ofstream P_stream;
  double P  = RcppDeepState_double();
  P_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/recall/inputs/P");
  P_stream << P;
  std::cout << "P values: "<< P << std::endl;
  P_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    recall(kappa,BR,P);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
