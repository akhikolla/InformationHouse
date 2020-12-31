#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <qs.h>
#include <DeepState.hpp>

bool check_BRK_combo(double BR, double P, double K);

TEST(rhoR_deepstate_test,check_BRK_combo_test){
  RInside R;
  std::cout << "input starts" << std::endl;
  std::ofstream BR_stream;
  double BR  = RcppDeepState_double();
  BR_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/check_BRK_combo/inputs/BR");
  BR_stream << BR;
  std::cout << "BR values: "<< BR << std::endl;
  BR_stream.close();
  std::ofstream P_stream;
  double P  = RcppDeepState_double();
  P_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/check_BRK_combo/inputs/P");
  P_stream << P;
  std::cout << "P values: "<< P << std::endl;
  P_stream.close();
  std::ofstream K_stream;
  double K  = RcppDeepState_double();
  K_stream.open("/home/akhila/fuzzer_packages/fuzzedpackages/rhoR/inst/testfiles/check_BRK_combo/inputs/K");
  K_stream << K;
  std::cout << "K values: "<< K << std::endl;
  K_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    check_BRK_combo(BR,P,K);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
